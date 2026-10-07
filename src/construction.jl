# Shapes are built from any sufficient set of keywords. Mass, volume, dimensions
# and ratios are all products of powers of each other, so in log space they form
# a small linear system, solved here by least squares.
#
# There are two ways in. The keyword constructors users call, `Cylinder(; ...)`,
# check the keywords with `check_shape` and throw errors that say what is wrong.
# The machine-level ones, `Cylinder(Unchecked(); ...)`, only solve: they have no
# checks and nothing to throw, so they run on GPUs and in wasm, and bad input
# gives bad numbers.
#
# Keywords are never looked up by name. Each keyword of a shape has a descriptor
# saying which unknown it fixes; the given keywords are merged onto a NamedTuple
# of `nothing` with the same names, so they line up with the descriptors by
# position, and everything after that is dispatch on the descriptor and value types.

struct MassKeyword end
struct DensityKeyword end
struct VolumeKeyword end
struct DimensionKeyword
    index::Int
end
# ratio = dimension[numerator] / (scale · dimension[denominator])
struct RatioKeyword
    numerator::Int
    denominator::Int
    scale::Float64
end
RatioKeyword(r::RatioKeyword) = r
RatioKeyword(r::Tuple) = RatioKeyword(r...)

"""
    ShapeSpec{shape}(powers, log_constant, ratios)

`powers` is a `NamedTuple` of dimensions and their powers in
`volume = exp(log_constant) · ∏ dimension^power`. `ratios` is a `NamedTuple`
of `(i, j, scale)`, each meaning `ratio = dimᵢ / (scale · dimⱼ)`.
"""
struct ShapeSpec{S,P<:NamedTuple,R<:NamedTuple}
    powers::P
    log_constant::Float64
    ratios::R
end
function ShapeSpec{S}(powers::P, log_constant, ratios::NamedTuple) where {S,P<:NamedTuple}
    r = _map_values(RatioKeyword, ratios)
    ShapeSpec{S,P,typeof(r)}(powers, log_constant, r)
end

# The spec of the half of a shape, with the same dimensions, has half the volume.
# Base's NamedTuple `map` doesn't specialise on the function, so map the values.
_map_values(f, nt::NamedTuple{K}) where {K} = NamedTuple{K}(map(f, values(nt)))

half_spec(spec::ShapeSpec{S}) where {S} = ShapeSpec{Half{S}}(spec.powers, spec.log_constant - log(2), spec.ratios)

_ndims(::ShapeSpec{S,<:NamedTuple{<:Any,<:NTuple{N,Any}}}) where {S,N} = N
_dimension_keywords(::NamedTuple{D,<:NTuple{N,Any}}) where {D,N} = NamedTuple{D}(ntuple(DimensionKeyword, Val(N)))

shape_keywords(spec::ShapeSpec) =
    merge((; mass = MassKeyword(), density = DensityKeyword(), volume = VolumeKeyword()),
          _dimension_keywords(spec.powers), spec.ratios)

# Unknowns: the log of each dimension, then of the volume, mass and density.
_n(spec) = Val(_ndims(spec) + 3)
_e(::Val{N}, i) where {N} = SVector(ntuple(j -> Float64(j == i), Val(N)))
_index(spec, ::VolumeKeyword) = _ndims(spec) + 1
_index(spec, ::MassKeyword) = _ndims(spec) + 2
_index(spec, ::DensityKeyword) = _ndims(spec) + 3
_index(spec, k::DimensionKeyword) = k.index

_si(x, unit) = Float64(ustrip(unit, x))

# A keyword with no matching descriptor gives a NamedTuple with extra names.
_known(::NamedTuple{K}, given::NamedTuple{K}, ::ShapeSpec) where {K} = given
_known(keywords, given, ::ShapeSpec{S}) where {S} =
    throw(ArgumentError("$S takes any sufficient set of the keywords $(keys(keywords)), got $(keys(given))"))

# Each row is (coefficients, right-hand side). A keyword not given adds nothing.
_row(spec::ShapeSpec, ::Any, ::Nothing) = (zero(_e(_n(spec), 1)), 0.0)
_row(spec::ShapeSpec, k::MassKeyword, x::Unitful.Mass) = (_e(_n(spec), _index(spec, k)), log(_si(x, u"kg")))
_row(spec::ShapeSpec, k::DensityKeyword, x::Unitful.Density) = (_e(_n(spec), _index(spec, k)), log(_si(x, u"kg/m^3")))
_row(spec::ShapeSpec, k::VolumeKeyword, x::Unitful.Volume) = (_e(_n(spec), _index(spec, k)), log(_si(x, u"m^3")))
_row(spec::ShapeSpec, k::DimensionKeyword, x::Unitful.Length) = (_e(_n(spec), k.index), log(_si(x, u"m")))
_row(spec::ShapeSpec, k::RatioKeyword, x::Real) =
    (_e(_n(spec), k.numerator) - _e(_n(spec), k.denominator), log(Float64(x) * k.scale))
_row(::ShapeSpec{S}, k, x) where {S} = throw(ArgumentError("$S keyword got $(repr(x)), which needs " *
    (k isa RatioKeyword ? "to be a plain number" : "units of $(_unit(k))")))

_unit(::MassKeyword) = u"kg"
_unit(::DensityKeyword) = u"kg/m^3"
_unit(::VolumeKeyword) = u"m^3"
_unit(::DimensionKeyword) = u"m"
_unit(::RatioKeyword) = NoUnits

_has_units(::MassKeyword, ::Unitful.Mass) = true
_has_units(::DensityKeyword, ::Unitful.Density) = true
_has_units(::VolumeKeyword, ::Unitful.Volume) = true
_has_units(::DimensionKeyword, ::Unitful.Length) = true
_has_units(::RatioKeyword, ::Real) = true
_has_units(k, x) = false

# The system the keywords make: the volume and density equations, then a row per keyword.
function _system(spec::ShapeSpec, keywords, given)
    n = _n(spec)
    volume = _e(n, _index(spec, VolumeKeyword())) -
        sum(ntuple(i -> values(spec.powers)[i] * _e(n, i), Val(_ndims(spec))); init = zero(_e(n, 1)))
    density = _e(n, _index(spec, MassKeyword())) - _e(n, _index(spec, DensityKeyword())) -
        _e(n, _index(spec, VolumeKeyword()))
    rows = ((volume, spec.log_constant), (density, 0.0),
            map((k, g) -> _row(spec, k, g), values(keywords), values(given))...)
    (transpose(hcat(map(first, rows)...)), SVector(map(last, rows)))
end

_given(spec, keywords, kw) = _known(keywords, merge(_map_values(_ -> nothing, keywords), kw), spec)

"""
    check_shape(spec, kw)

Throw an `ArgumentError` that says what is wrong if the keywords `kw` don't fix a
shape of `spec`: a keyword it doesn't take, one without the right units or that
isn't positive, too few to fix the shape, or ones that contradict each other.
"""
function check_shape(spec::ShapeSpec{S}, kw::NamedTuple) where {S}
    keywords = shape_keywords(spec)
    given = _given(spec, keywords, kw)
    for (name, x) in pairs(kw)
        _has_units(keywords[name], x) || throw(ArgumentError(
            "$S `$name` must be $(_unit(keywords[name]) == NoUnits ? "a plain number" : "a quantity like $(1.0 * _unit(keywords[name]))"), got $(repr(x))"))
        isfinite(ustrip(x)) && ustrip(x) > 0 ||
            throw(ArgumentError("$S `$name` must be positive and finite, got $(repr(x))"))
    end
    A, b = _system(spec, keywords, given)
    # The rows are integer, so the normal matrix is singular exactly when something is free.
    M = transpose(A) * A
    others = filter(k -> !(k in keys(kw)), keys(keywords))
    abs(det(M)) > 0.5 || throw(ArgumentError(
        "$S needs more keywords: $(join(keys(kw), ", ")) don't fix it. Add some of $(join(others, ", "))"))
    maximum(abs, A * (M \ (transpose(A) * b)) - b) < 1e-8 ||
        throw(ArgumentError("$S keywords contradict each other: $(join(keys(kw), ", ")) can't all hold"))
    nothing
end

_solved(spec, k::MassKeyword, logs) = exp(logs[_index(spec, k)]) * u"kg"
_solved(spec, k::DensityKeyword, logs) = exp(logs[_index(spec, k)]) * u"kg/m^3"
_solved(spec, k::VolumeKeyword, logs) = exp(logs[_index(spec, k)]) * u"m^3"
_solved(spec, k::DimensionKeyword, logs) = exp(logs[k.index]) * u"m"
_solved(spec, k::RatioKeyword, logs) = exp(logs[k.numerator] - logs[k.denominator]) / k.scale

_pick(::Nothing, solved) = solved
_pick(given, solved) = given

"""
    _resolve_shape(spec, kw) -> NamedTuple

Solve for every keyword of a shape (mass, density, volume, dimensions and ratios)
from the keywords `kw`. Given values come back unchanged; solved ones are in SI.
Nothing is checked: see `check_shape`.
"""
function _resolve_shape(spec::ShapeSpec, kw::NamedTuple)
    keywords = shape_keywords(spec)
    given = _given(spec, keywords, kw)
    A, b = _system(spec, keywords, given)
    _resolved(spec, keywords, given, _solve_unchecked(transpose(A) * A, transpose(A) * b))
end

_resolved(spec, keywords::NamedTuple{K}, given, logs) where {K} =
    NamedTuple{K}(map((k, g) -> _pick(g, _solved(spec, k, logs)), values(keywords), values(given)))

# Solve M x = r by LU, without the singularity check that StaticArrays' `\` throws
# from. A singular M gives non-finite x.
function _solve_unchecked(M, r)
    F = lu(M; check = false)
    L, U = F.L, F.U
    x = r[F.p]
    for i in 2:length(x), j in 1:i - 1
        x = Base.setindex(x, x[i] - L[i, j] * x[j], i)
    end
    for i in length(x):-1:1
        s = x[i]
        for j in i + 1:length(x)
            s -= U[i, j] * x[j]
        end
        x = Base.setindex(x, s / U[i, i], i)
    end
    x
end

"""
    Unchecked()

Builds a shape, layer or body at machine level: `Cylinder(Unchecked(); mass, density,
axis_ratio_b)` solves the keywords like `Cylinder(; ...)` but checks nothing and
throws nothing, so it can run on a GPU or in wasm. Bad input gives bad numbers.
"""
struct Unchecked end
