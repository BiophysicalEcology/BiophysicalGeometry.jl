module BiophysicalGeometryWhiskExt

# Compiles models to wasm with Whisk.jl, behind one entry point with a fixed interface: settings in, and into
# buffers the page owns, the body's triangles to draw, its areas and masses, and its shadow. See `compile_wasm`.

using BiophysicalGeometry, Unitful, Whisk
const BG = BiophysicalGeometry

# An array over memory owned by someone else: wasm memory that the page allocated and reads back.
struct PtrArray{T,N} <: AbstractArray{T,N}
    ptr::Ptr{T}
    dims::NTuple{N,Int}
end
Base.size(a::PtrArray) = a.dims
Base.IndexStyle(::Type{<:PtrArray}) = IndexLinear()
Base.getindex(a::PtrArray, i::Int) = unsafe_load(a.ptr, i)
Base.setindex!(a::PtrArray, v, i::Int) = unsafe_store!(a.ptr, v, i)

# Where each number goes, counting from 1: the triangles written, the totals, the shadow's area and extent, and
# then each part's exposed area and mass.
const NUMBERS = (; triangle_count = 1, outer_area = 2, skin_area = 3, hidden_area = 4, volume = 5,
                   shadow_area = 6, shadow_extent = 7, parts = 11)
const MAX_PARTS = 64
const TRIANGLE_FLOATS = 10   # three corners, in metres, then the part's number from 0

# The entry point: the models, and their settings with defaults, are in its type, so that it compiles as a plain
# function and each model's branch is concrete.
struct Entry{Models,Defaults} end

read_settings(::NamedTuple{K,NTuple{N,Float64}}, address) where {K,N} =
    NamedTuple{K}(ntuple(i -> unsafe_load(Ptr{Float64}(address), i), Val(N)))

# Build model `index` (from 1) from the settings at `settings`, and write it into the buffers at `triangles`
# (`capacity` triangles), `numbers` and `shadow` (`n × n`), with the sun in the direction `(sx, sy, sz)`. Returns its
# number of parts, or 0 for no such model.
function (::Entry{Models,Defaults})(index::Int, settings::Int, triangles::Int, capacity::Int, numbers::Int,
                                     shadow::Int, n::Int, sx::Float64, sy::Float64, sz::Float64) where {Models,Defaults}
    out = (PtrArray(Ptr{Float32}(triangles), (TRIANGLE_FLOATS, capacity)),
           PtrArray(Ptr{Float64}(numbers), (NUMBERS.parts - 1 + 2 * MAX_PARTS,)),
           PtrArray(Ptr{Bool}(shadow), (n, n)), (sx, sy, sz))
    counts = map(Models, Defaults, ntuple(identity, Val(length(Models)))) do model, defaults, k
        k == index ? write_body!(model(read_settings(defaults, settings)), out...) : 0
    end
    max(counts...)
end

# Drawn from tiles a third as fine as the meshes the shadow is rasterised from: smooth enough, and few enough
# triangles to redraw as the body is turned.
coarse(tile::BG.Tile) = BG.Tile(tile.point, max(2, cld(tile.nu, 3)), max(2, cld(tile.nv, 3)))

function write_body!(body::CompositeBody, triangles, numbers, shadow, sun)
    nparts = length(body.parts)
    numbers[NUMBERS.triangle_count] = 0.0
    capacity = size(triangles, 2)
    part_tiles = map(tiles -> map(coarse, tiles), BG._part_tiles(body))
    map(part_tiles, ntuple(identity, Val(nparts))) do tiles, k
        BG.foreach_triangle(tiles) do p1, p2, p3
            t = Int(numbers[NUMBERS.triangle_count]) + 1
            if t <= capacity
                for (c, p) in enumerate((p1, p2, p3)), i in 1:3
                    triangles[3 * (c - 1) + i, t] = p[i]
                end
                triangles[TRIANGLE_FLOATS, t] = k - 1
                numbers[NUMBERS.triangle_count] = t
            end
        end
    end
    covered = values(BG.covered_areas(body.parts, body.joins))
    map(values(body.parts), covered, ntuple(identity, Val(nparts))) do part, hidden, k
        numbers[NUMBERS.parts + 2 * (k - 1)] = ustrip(u"m^2", total_area(part) - hidden)
        numbers[NUMBERS.parts + 2 * (k - 1) + 1] = ustrip(u"kg", BG.mass(shape(part)))
    end
    numbers[NUMBERS.outer_area] = ustrip(u"m^2", total_area(body))
    numbers[NUMBERS.skin_area] = ustrip(u"m^2", skin_area(body))
    numbers[NUMBERS.hidden_area] = ustrip(u"m^2", +(covered...))
    numbers[NUMBERS.volume] = ustrip(u"m^3", flesh_volume(body))
    area, g = BG._silhouette_rasterized!(shadow, body, sun)
    numbers[NUMBERS.shadow_area] = ustrip(u"m^2", area)
    for (i, x) in enumerate((g.x0, g.x1, g.y0, g.y1))
        numbers[NUMBERS.shadow_extent + i - 1] = x
    end
    return nparts
end

function BiophysicalGeometry.compile_wasm(models::NamedTuple)
    callables = map(first, values(models))
    defaults = map(m -> map(Float64, last(m)), values(models))
    for (name, f, d) in zip(keys(models), callables, defaults)
        Base.issingletontype(typeof(f)) ||
            throw(ArgumentError("model `$name` must be a function or callable without fields, to compile; got $(typeof(f))"))
        body = f(d)
        body isa CompositeBody || throw(ArgumentError("model `$name` must make a CompositeBody; got $(typeof(body))"))
        length(body.parts) <= MAX_PARTS || throw(ArgumentError("model `$name` has more than $MAX_PARTS parts"))
    end
    spec = (;
        models = [(; name = string(name), settings = collect(string.(keys(d))), defaults = collect(values(d)),
                     parts = collect(string.(keys(f(d).parts))))
                  for (name, f, d) in zip(keys(models), callables, defaults)],
        triangle_floats = TRIANGLE_FLOATS, max_parts = MAX_PARTS, numbers = NUMBERS)
    file = tempname() * ".wasm"
    Whisk.compile_to_wasm(Entry{callables,defaults}(),
                          Tuple{Int,Int,Int,Int,Int,Int,Int,Float64,Float64,Float64};
                          out_path = file, name = "biophysical_run", optimize = true)
    wasm = read(file)
    rm(file; force = true)
    return (; wasm, spec)
end
BiophysicalGeometry.compile_wasm(model, defaults::NamedTuple) = compile_wasm((; model = (model, defaults)))

function BiophysicalGeometry.compile_wasm(path::AbstractString, models::NamedTuple; name::AbstractString = "model")
    compiled = compile_wasm(models)
    mkpath(path)
    write(joinpath(path, "$name.wasm"), compiled.wasm)
    open(io -> json(io, compiled.spec), joinpath(path, "$name.json"), "w")
    cp(joinpath(@__DIR__, "biophysical.mjs"), joinpath(path, "biophysical.mjs"); force = true)
    return path
end
BiophysicalGeometry.compile_wasm(path::AbstractString, model, defaults::NamedTuple; name::AbstractString = "model") =
    compile_wasm(path, (; model = (model, defaults)); name)

# Just enough JSON for the spec.
json(io, x::AbstractString) = (print(io, '"'); escape_string(io, x, '"'); print(io, '"'))
json(io, x::Symbol) = json(io, string(x))
json(io, x::Real) = print(io, isfinite(x) ? x : "null")
json(io, x::Union{AbstractVector,Tuple}) =
    (print(io, '['); for (i, v) in enumerate(x); i > 1 && print(io, ','); json(io, v); end; print(io, ']'))
function json(io, x::NamedTuple)
    print(io, '{')
    for (i, (k, v)) in enumerate(pairs(x))
        i > 1 && print(io, ',')
        json(io, k); print(io, ':'); json(io, v)
    end
    print(io, '}')
end

end
