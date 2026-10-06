# Shapes are built from any sufficient set of keywords. Mass, volume, dimensions
# and ratios are all products of powers of each other, so in log space they form
# a small linear system, solved here by least squares.

"""
    ShapeSpec(dimensions, powers, log_constant, ratios)

`volume = exp(log_constant) · ∏ dimensionᵢ^powersᵢ`, and each entry of `ratios`
is `name => (i, j, scale)` meaning `name = dimᵢ / (scale · dimⱼ)`.
"""
struct ShapeSpec{Dims,Ratios,P,R}
    powers::P
    log_constant::Float64
    ratios::R
end
ShapeSpec(dims, powers, log_constant, ratios) =
    ShapeSpec{dims,map(first, ratios),typeof(powers),typeof(map(last, ratios))}(
        powers, log_constant, map(last, ratios))

_si(x, unit) = Float64(ustrip(unit, x))

# Unknowns: the log of each dimension, then of the volume, mass and density.
_nunknowns(::ShapeSpec{D}) where {D} = Val(length(D) + 3)
_e(::Val{N}, i) where {N} = SVector(ntuple(j -> Float64(j == i), Val(N)))

_unit(::ShapeSpec, ::Val{:mass}) = u"kg"
_unit(::ShapeSpec, ::Val{:density}) = u"kg/m^3"
_unit(::ShapeSpec, ::Val{:volume}) = u"m^3"
_unit(::ShapeSpec{D,R}, ::Val{k}) where {D,R,k} = k in D ? u"m" : k in R ? NoUnits : nothing

function _check_keyword(name, spec::ShapeSpec{D,R}, key::Val{k}, x) where {D,R,k}
    unit = _unit(spec, key)
    unit === nothing && throw(ArgumentError("$name has no keyword `$k`; use any sufficient set of: " *
                                            join((:mass, :density, :volume, D..., R...), ", ")))
    dimension(x) == dimension(unit) || throw(ArgumentError(
        "`$k` must be $(unit == NoUnits ? "a plain number" : "a quantity like $unit"), got $(repr(x))"))
    isfinite(ustrip(x)) && ustrip(x) > 0 ||
        throw(ArgumentError("`$k` must be positive and finite, got $(repr(x))"))
    nothing
end

_row(spec::ShapeSpec{D}, ::Val{:volume}, x) where {D} = (_e(_nunknowns(spec), length(D) + 1), log(_si(x, u"m^3")))
_row(spec::ShapeSpec{D}, ::Val{:mass}, x) where {D} = (_e(_nunknowns(spec), length(D) + 2), log(_si(x, u"kg")))
_row(spec::ShapeSpec{D}, ::Val{:density}, x) where {D} = (_e(_nunknowns(spec), length(D) + 3), log(_si(x, u"kg/m^3")))
function _row(spec::ShapeSpec{D,R}, ::Val{k}, x) where {D,R,k}
    n = _nunknowns(spec)
    i = findfirst(==(k), D)
    i === nothing || return (_e(n, i), log(_si(x, u"m")))
    (a, b, scale) = spec.ratios[findfirst(==(k), R)]
    (_e(n, a) - _e(n, b), log(x * scale))
end

_pick(kw::NamedTuple{K}, ::Val{k}, solved) where {K,k} = k in K ? kw[k] : solved

"""
    _resolve_shape(name, spec, kw) -> NamedTuple{(:mass, :density, ratios...)}

Solve for a shape's mass, density and axis ratios from the keywords `kw`. Given
values come back unchanged; solved ones are in kg and kg/m³.
"""
function _resolve_shape(name, spec::ShapeSpec{D,R}, kw::NamedTuple{K}) where {D,R,K}
    keys = map(Val, K)
    map((key, x) -> _check_keyword(name, spec, key, x), keys, values(kw))

    n = _nunknowns(spec)
    iV, iM, iρ = length(D) + 1, length(D) + 2, length(D) + 3
    volume = _e(n, iV) - sum(ntuple(i -> spec.powers[i] * _e(n, i), Val(length(D))); init = zero(_e(n, 1)))
    rows = ((volume, spec.log_constant), (_e(n, iM) - _e(n, iρ) - _e(n, iV), 0.0),
            map((key, x) -> _row(spec, key, x), keys, values(kw))...)
    A = transpose(hcat(map(first, rows)...))
    b = SVector(map(last, rows))

    # The rows are integer, so the normal matrix is singular exactly when something is free.
    M = transpose(A) * A
    abs(det(M)) > 0.5 || throw(ArgumentError(
        "$name is under-determined: given $(join(K, ", ")); give more of: " *
        join(filter(k -> !(k in K), (:mass, :density, :volume, D..., R...)), ", ")))
    x = M \ (transpose(A) * b)
    maximum(abs, A * x - b) < 1e-8 || throw(ArgumentError(
        "$name inputs are inconsistent: $(join(K, ", ")) don't fit together"))

    ratios = map(map(Val, R), spec.ratios) do key, (i, j, scale)
        _pick(kw, key, exp(x[i] - x[j]) / scale)
    end
    NamedTuple{(:mass, :density, R...)}((_pick(kw, Val(:mass), exp(x[iM]) * u"kg"),
                                         _pick(kw, Val(:density), exp(x[iρ]) * u"kg/m^3"), ratios...))
end

# Marker for the shapes' only inner constructor, so they can't be built positionally.
struct Resolved end
const RESOLVED = Resolved()
