# Keyword construction for shapes.
#
# Every shape is built from keywords — `mass`, `density`, `volume`, its own
# dimensions and its own axis ratios — and whatever isn't given is solved for.
# Each shape struct stores only mass, density and its dimensionless ratios; the
# solve here gets any sufficient set of keywords down to those.
#
# All the relations are products of powers:
#
#     mass   = density · volume
#     volume = k · ∏ dimensionᵢ^pᵢ          (e.g. π·radius²·length)
#     ratio  = dimensionᵢ / (s · dimensionⱼ)  (e.g. length / (2·radius))
#
# so in log space they are linear, and one small linear solve covers every
# combination of inputs: given values fix some unknowns, the relations tie the
# rest together, and the rank tells us whether the inputs are enough. Each shape
# only declares its dimensions, its volume monomial and its ratios.

"""
    _ShapeSpec(dimensions, powers, log_constant, ratios)

How a shape's dimensions relate to its volume and ratios. `dimensions` are the
keyword names, `volume = exp(log_constant) · ∏ dimensionᵢ^powersᵢ`, and each
entry of `ratios` is `name => (i, j, scale)` meaning `name = dimᵢ / (scale · dimⱼ)`.
"""
struct _ShapeSpec{D,P,R}
    dimensions::D
    powers::P
    log_constant::Float64
    ratios::R
end

# SI units the solve works in. Inputs are converted to these; plain numbers are
# taken to be in them already.
const _MASS_UNIT = u"kg"
const _LENGTH_UNIT = u"m"
const _VOLUME_UNIT = u"m^3"
const _DENSITY_UNIT = u"kg/m^3"

_to_si(x::Unitful.AbstractQuantity, unit) = Float64(ustrip(unit, x))
_to_si(x::Real, _) = Float64(x)

function _check_positive(name, x)
    x isa Union{Real,Unitful.AbstractQuantity} ||
        throw(ArgumentError("`$name` must be a number, got $(repr(x))"))
    v = ustrip(x)
    isfinite(v) && v > 0 || throw(ArgumentError("`$name` must be positive and finite, got $x"))
    return x
end

"""
    _resolve_shape(name, spec, kw) -> (; mass, density, ratios)

Solve for a shape's mass, density and axis ratios (in `spec.ratios` order) from
the keywords `kw` (a NamedTuple of the inputs actually given). Inputs that are
given come back unchanged; solved ones are in SI units if any input was
unitful, plain numbers otherwise.
"""
function _resolve_shape(name::AbstractString, spec::_ShapeSpec, kw::NamedTuple)
    dims = spec.dimensions
    ratio_names = map(first, spec.ratios)
    known = (:mass, :density, :volume, dims..., ratio_names...)
    for k in keys(kw)
        k in known || throw(ArgumentError(
            "$name has no keyword `$k`; use any sufficient set of: $(join(known, ", "))"))
        _check_positive(k, kw[k])
    end

    # Fast path: the stored fields themselves. No solve, so the given values
    # (and their types — units, AD duals) pass straight through.
    if haskey(kw, :mass) && haskey(kw, :density) && all(r -> haskey(kw, r), ratio_names) &&
            !haskey(kw, :volume) && !any(d -> haskey(kw, d), dims)
        return (; mass = kw.mass, density = kw.density, ratios = map(r -> kw[r], ratio_names))
    end

    # Unknowns: log of each dimension, then log volume, log mass, log density.
    n = length(dims)
    iV, iM, iρ = n + 1, n + 2, n + 3
    rows = Vector{Float64}[]
    rhs = Float64[]
    row() = zeros(n + 3)
    # volume = k · ∏ dimᵢ^pᵢ
    r = row(); r[iV] = 1; for i in 1:n; r[i] = -spec.powers[i]; end
    push!(rows, r); push!(rhs, spec.log_constant)
    # mass = density · volume
    r = row(); r[iM] = 1; r[iρ] = -1; r[iV] = -1
    push!(rows, r); push!(rhs, 0.0)
    for (i, d) in enumerate(dims)
        haskey(kw, d) || continue
        r = row(); r[i] = 1; push!(rows, r); push!(rhs, log(_to_si(kw[d], _LENGTH_UNIT)))
    end
    for (rname, (i, j, s)) in spec.ratios
        haskey(kw, rname) || continue
        r = row(); r[i] = 1; r[j] = -1; push!(rows, r); push!(rhs, log(Float64(kw[rname]) * s))
    end
    for (k, i, unit) in ((:volume, iV, _VOLUME_UNIT), (:mass, iM, _MASS_UNIT), (:density, iρ, _DENSITY_UNIT))
        haskey(kw, k) || continue
        r = row(); r[i] = 1; push!(rows, r); push!(rhs, log(_to_si(kw[k], unit)))
    end

    A = reduce(vcat, permutedims.(rows))
    F = svd(A)
    tol = maximum(F.S) * 1e-10
    rank = count(>(tol), F.S)
    given = isempty(kw) ? "nothing" : join(keys(kw), ", ")
    rank < n + 3 && throw(ArgumentError(
        "$name is under-determined: given $given; give $(n + 3 - rank) more of: " *
        join(filter(k -> !haskey(kw, k), known), ", ")))
    x = A \ rhs
    resid = A * x - rhs
    maximum(abs, resid) < 1e-8 || throw(ArgumentError(
        "$name inputs are inconsistent: $given don't fit together"))

    unitful = any(v -> v isa Unitful.AbstractQuantity, values(kw))
    wrap(v, unit) = unitful ? v * unit : v
    mass = haskey(kw, :mass) ? kw.mass : wrap(exp(x[iM]), _MASS_UNIT)
    density = haskey(kw, :density) ? kw.density : wrap(exp(x[iρ]), _DENSITY_UNIT)
    ratios = map(spec.ratios) do (rname, (i, j, s))
        haskey(kw, rname) ? kw[rname] : exp(x[i] - x[j]) / s
    end
    return (; mass, density, ratios)
end

# Marker for the shape structs' only inner constructor: shapes are built by
# keyword, never positionally (positional arguments carry no meaning).
struct _Resolved end
const _RESOLVED = _Resolved()
