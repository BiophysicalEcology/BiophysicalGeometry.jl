# The builder page's entry point, compiled to wasm by `build.jl`. It builds one of the animals from the page's
# settings with the package's machine-level builder, and writes what the page shows into buffers in wasm memory
# that the page owns: the triangles to draw, the numbers for the tables, and the shadow.

using BiophysicalGeometry, Unitful
const BG = BiophysicalGeometry
const AB = BG.AnimalBuilder

# An array over memory owned by someone else: here, wasm memory the page allocated and reads back.
struct PtrArray{T,N} <: AbstractArray{T,N}
    ptr::Ptr{T}
    dims::NTuple{N,Int}
end
Base.size(a::PtrArray) = a.dims
Base.IndexStyle(::Type{<:PtrArray}) = IndexLinear()
Base.getindex(a::PtrArray, i::Int) = unsafe_load(a.ptr, i)
Base.setindex!(a::PtrArray, v, i::Int) = unsafe_store!(a.ptr, v, i)

const NSETTINGS = length(AB.SETTING_NAMES)
read_settings(address) = NamedTuple{AB.SETTING_NAMES}(ntuple(i -> unsafe_load(Ptr{Float64}(address), i), Val(NSETTINGS)))

# The numbers, in order: the triangles written, the totals, the shadow and its extent, then each part's exposed
# area and mass.
const TRIANGLE_COUNT, OUTER_AREA, SKIN_AREA, HIDDEN_AREA, VOLUME, MEEH, SHADOW_AREA, SHADOW_EXTENT, PARTS = 1, 2, 3, 4, 5, 6, 7, 8, 12
const TRIANGLE_FLOATS = 10   # three corners, then the part's number

"""
    animal_builder(animal, settings, triangles, capacity, numbers, shadow, n, zenith, azimuth) -> Int

Build animal number `animal` of `AnimalBuilder.ANIMALS` from the `NSETTINGS` `Float64`s at address `settings`.
Write up to `capacity` of its posed triangles at `triangles` (`Float32`s, `TRIANGLE_FLOATS` each, in metres), its
numbers at `numbers` (`Float64`s), and its shadow toward the sun at `zenith` and `azimuth` (radians) into the
`n × n` `Bool`s at `shadow`. Returns the number of parts.
"""
function animal_builder(animal::Int, settings::Int, triangles::Int, capacity::Int, numbers::Int, shadow::Int,
                        n::Int, zenith::Float64, azimuth::Float64)
    s = read_settings(settings)
    out = (PtrArray(Ptr{Float32}(triangles), (TRIANGLE_FLOATS, capacity)), PtrArray(Ptr{Float64}(numbers), (PARTS + 2 * 32,)),
           PtrArray(Ptr{Bool}(shadow), (n, n)), (sin(zenith) * cos(azimuth), sin(zenith) * sin(azimuth), cos(zenith)))
    animal == 1 && return output!(AB.build(AB.Dog(), s), s, out...)
    animal == 2 && return output!(AB.build(AB.Mouse(), s), s, out...)
    animal == 3 && return output!(AB.build(AB.Elephant(), s), s, out...)
    animal == 4 && return output!(AB.build(AB.Human(), s), s, out...)
    animal == 5 && return output!(AB.build(AB.Kangaroo(), s), s, out...)
    animal == 6 && return output!(AB.build(AB.Tyrannosaur(), s), s, out...)
    animal == 7 && return output!(AB.build(AB.Giraffe(), s), s, out...)
    animal == 8 && return output!(AB.build(AB.Cow(), s), s, out...)
    animal == 9 && return output!(AB.build(AB.Bird(), s), s, out...)
    animal == 10 && return output!(AB.build(AB.Seal(), s), s, out...)
    return 0
end

coarse(tile::BG.Tile) = BG.Tile(tile.point, max(2, cld(tile.nu, 3)), max(2, cld(tile.nv, 3)))

function output!(body, s, triangles, numbers, shadow, sun)
    nparts = length(body.parts)
    numbers[TRIANGLE_COUNT] = 0.0
    capacity = size(triangles, 2)
    # Drawn from tiles a third as fine as the meshes the shadow is rasterised from: enough to look smooth, and few
    # enough triangles for the page to redraw as the animal is turned.
    part_tiles = map(tiles -> map(coarse, tiles), BG._part_tiles(body))
    map(part_tiles, ntuple(identity, Val(nparts))) do tiles, k
        BG.foreach_triangle(tiles) do p1, p2, p3
            t = Int(numbers[TRIANGLE_COUNT]) + 1
            if t <= capacity
                for (c, p) in enumerate((p1, p2, p3)), i in 1:3
                    triangles[3 * (c - 1) + i, t] = p[i]
                end
                triangles[TRIANGLE_FLOATS, t] = k - 1
                numbers[TRIANGLE_COUNT] = t
            end
        end
    end
    covered = values(BG.covered_areas(body.parts, body.joins))
    map(values(body.parts), covered, ntuple(identity, Val(nparts))) do part, hidden, k
        numbers[PARTS + 2 * (k - 1)] = ustrip(u"m^2", total_area(part) - hidden)
        numbers[PARTS + 2 * (k - 1) + 1] = ustrip(u"kg", BG.mass(shape(part)))
    end
    total = ustrip(u"m^2", total_area(body))
    numbers[OUTER_AREA] = total
    numbers[SKIN_AREA] = ustrip(u"m^2", skin_area(body))
    numbers[HIDDEN_AREA] = ustrip(u"m^2", +(covered...))
    numbers[VOLUME] = ustrip(u"m^3", flesh_volume(body))
    numbers[MEEH] = total / s.mass^(2 / 3)
    area, g = BG._silhouette_rasterized!(shadow, body, sun)
    numbers[SHADOW_AREA] = ustrip(u"m^2", area)
    for (i, x) in enumerate((g.x0, g.x1, g.y0, g.y1))
        numbers[SHADOW_EXTENT + i - 1] = x
    end
    return nparts
end
