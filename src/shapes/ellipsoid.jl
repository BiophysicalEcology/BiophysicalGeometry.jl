"""
    Ellipsoid(mass, density, b, c, pole_a_truncation=0.0) <: AbstractShape

A spheroid (b = c) with polar axis `a = axis_ratio_b * b` along `+x` — prolate
for `axis_ratio_b > 1`, oblate for `axis_ratio_b < 1`. With
`pole_a_truncation > 0`, the `+x` end is sliced flat: the cut plane sits at
`x = (1 - pole_a_truncation) * a`. `pole_a_truncation = 0` is a full ellipsoid;
`pole_a_truncation = 1` cuts through the centre (half ellipsoid).

The flat face exposed by the cut has semi-axes
`(b * sqrt(1 - (1 - pole_a_truncation)^2), c * sqrt(1 - (1 - pole_a_truncation)^2))`.

Surface area is reported for the full (untruncated) ellipsoid — this slightly
over-counts (by the removed spherical cap) for thin truncations.
"""
struct Ellipsoid{M,D,B,C,T} <: AbstractEllipsoidal
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
    pole_a_truncation::T
end
Ellipsoid(mass, density, b, c) = Ellipsoid(mass, density, b, c, 0.0)

# x-position of the truncated pole_a (as a fraction of a). 1.0 = full ellipsoid.
_pole_a_x_ratio(s::Ellipsoid) = 1 - s.pole_a_truncation
# Radial scale at the truncated pole (so cut disc has y/z extents = scale * b/c).
_pole_a_radial_scale(s::Ellipsoid) = sqrt(max(0.0, 1 - _pole_a_x_ratio(s)^2))

# Exact spheroid surface area (b = c, the equatorial radius; a the polar axis).
# The formula uses sqrt / asin / log on dimensionless eccentricity, so this is
# the one place we step out of Unitful — ustrip once, compute, re-wrap.
# Prolate (a > c): eccentricity about the major axis a, standard asin form.
# Oblate (a < c): eccentricity about the equatorial major axis c, log form.
function _spheroid_area(a, b, c)
    am = ustrip(u"m", a); bm = ustrip(u"m", b); cm = ustrip(u"m", c)
    area = if abs(am - cm) < max(am, cm) * 1e-9
        4 * π * bm^2                       # sphere limit: a == c ⇒ e → 0
    elseif am > cm
        e = sqrt(am^2 - cm^2) / am
        2 * π * bm^2 + 2 * π * (am * bm / e) * asin(e)
    else
        e = sqrt(cm^2 - am^2) / cm
        2 * π * bm^2 + π * (am^2 / e) * log((1 + e) / (1 - e))
    end
    area * u"m^2"
end

# Canonicalise to m^3 to avoid weird unit ratios from mass/density (e.g. g/(kg/m^3))
# that propagate into cbrt() and trip up Enzyme's typeunstablerules path.
_ellipsoid_volume(shape::Ellipsoid) = uconvert(u"m^3", body_volume(shape))
_ellipsoid_fat_volume(shape::Ellipsoid, fat_layer::FatLayer) =
    uconvert(u"m^3", fat_volume(shape, fat_layer))

# b-semi-minor from an enclosed volume. Identical to the equivalent-sphere radius
# of the volume scaled by the axis ratio, so the cube-root lives once (geometry.jl).
_ellipsoid_b(shape::Ellipsoid, volume) = _sphere_radius(volume / shape.axis_ratio_b)

function _skin_level(shape::Ellipsoid, volume)
    b = _ellipsoid_b(shape, volume)          # c = b; a = axis_ratio_b · b
    a = shape.axis_ratio_b * b
    (; dims = (; a_semi_major_skin = a, b_semi_minor_skin = b, c_semi_minor_skin = b),
       area = _spheroid_area(a, b, b))
end
function _fibrous_level(shape::Ellipsoid, skin, thickness)
    a = skin.a_semi_major_skin + thickness
    b = skin.b_semi_minor_skin + thickness
    c = skin.c_semi_minor_skin + thickness
    (; dims = (; a_semi_major_fibrous = a, b_semi_minor_fibrous = b, c_semi_minor_fibrous = c),
       area = _spheroid_area(a, b, c))
end

# Smooth Heaviside on a length: ≈ 1 for fat ≫ ε, ≈ 0 for fat ≪ -ε, smooth
# everywhere. Used to blend the "no fat" and "with fat" geometry branches.
@inline _smooth_step_meters(fat) = 0.5 * (1 + fat / sqrt(fat*fat + (1e-9u"m")^2))

# Ellipsoid fat is a *uniform* shell: each skin semi-axis is the flesh semi-axis
# plus one fat thickness (so the skin is not a scaled flesh ellipsoid). Skin dims
# are therefore coupled to the fat solve, so the fat paths override the
# generic orchestrator rather than compute skin from volume independently.
#
# The "not enough fat" fallback (full-volume axes, fat = 0) and the flesh + fat
# axes are smooth-blended rather than switched with an if/else, which would be
# a value discontinuity at raw_fat = 0 — bad for reverse-mode AD. A single
# smooth Heaviside drives both the axes blend and the effective fat value, so
# they stay consistent at every raw_fat.
function _ellipsoid_fat_skin(shape::Ellipsoid, fat_layer::FatLayer)
    volume = _ellipsoid_volume(shape)
    fat_v = _ellipsoid_fat_volume(shape, fat_layer)
    flesh_v = volume - fat_v
    b_flesh = _ellipsoid_b(shape, flesh_v)
    a_flesh = shape.axis_ratio_b * b_flesh
    raw_fat = prolate_fat_layer(flesh_v, fat_v, shape.axis_ratio_b, b_flesh)
    b_full = _ellipsoid_b(shape, volume)
    a_full = shape.axis_ratio_b * b_full
    w = _smooth_step_meters(raw_fat)
    fat = w * raw_fat                                   # smooth-clamped to ≥ 0
    a = w * (a_flesh + raw_fat) + (1 - w) * a_full
    b = w * (b_flesh + raw_fat) + (1 - w) * b_full
    dims = (; a_semi_major_skin = a, b_semi_minor_skin = b, c_semi_minor_skin = b)
    (; volume, dims, fat, area = _spheroid_area(a, b, b))
end

function geometry(shape::Ellipsoid, fat_layer::FatLayer)
    s = _ellipsoid_fat_skin(shape, fat_layer)
    Geometry(s.volume, merge(s.dims, (; fat = s.fat)), SurfaceAreas(; total = s.area))
end
function geometry(shape::Ellipsoid, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    s = _ellipsoid_fat_skin(shape, fat_layer)
    fibrous = _fibrous_level(shape, s.dims, fibrous_layer.thickness)
    Geometry(s.volume, merge(s.dims, fibrous.dims, (; fat = s.fat)),
             SurfaceAreas(; total = fibrous.area, skin = s.area,
                          convection = _convective_area(fibrous_layer, s.area)))
end

# Fat thickness calculation
#
# The Newton solve is unitless. `prolate_fat_layer` is the dimensional boundary:
# it ustrips volumes and radii once, solves, and re-wraps the answer.

function prolate_fat_layer(flesh_volume, fat_volume, axis_ratio_b, semi_minor_flesh)
    fat_m = _prolate_fat_layer_m(
        ustrip(u"m^3", flesh_volume),
        ustrip(u"m^3", fat_volume),
        axis_ratio_b,
        ustrip(u"m", semi_minor_flesh),
    )
    return max(0.0, fat_m) * u"m"
end

function _prolate_fat_layer_m(flesh_volume, fat_volume, axis_ratio_b, semi_minor_flesh)
    # Find uniform fat thickness X such that the outer prolate spheroid
    # (semi-axes a+X, b+X, b+X, a = axis_ratio_b*b) has volume V_total.
    # Solves f(X) = (a+X)(b+X)² = (3/4π)*V_total via Newton's method.
    # f is strictly monotone (f'(X) > 0 for all X > -b), so there is exactly
    # one non-negative root. Fixed 10 iterations reaches machine precision
    # without branching, which keeps Enzyme AD well-behaved. Replaces the
    # Cardano formula which silently returns the wrong root when the
    # discriminant is negative (casus irreducibilis).
    b      = semi_minor_flesh
    a      = axis_ratio_b * b
    target = (3.0 / (4.0 * π)) * (flesh_volume + fat_volume)
    X      = 0.0
    for _ in 1:10
        bX = b + X
        aX = a + X
        X -= (aX * bX^2 - target) / (bX^2 + 2 * aX * bX)
    end
    return X
end

# Surface area

function surface_area(shape::Ellipsoid, body::AbstractBody)
    _spheroid_area(body.geometry.length.a_semi_major_skin,
                   body.geometry.length.b_semi_minor_skin,
                   body.geometry.length.c_semi_minor_skin)
end

# Silhouette area

# For an ellipsoid (a, b, c) the silhouette projected along direction d
# is an ellipse of area π·sqrt(b²c²·d_x² + a²c²·d_y² + a²b²·d_z²). With
# the sun in the equatorial plane at angle θ from the long (a) axis,
# d = (cos θ, sin θ, 0) gives π·c·sqrt(b²·cos²θ + a²·sin²θ).
function silhouette(shape::Ellipsoid, a, b, c, θ)
    π * c * sqrt(b^2 * cos(θ)^2 + a^2 * sin(θ)^2)
end
function silhouette(sh::Ellipsoid, ::AbstractInsulationLayer, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = π * d.a * d.b, parallel = π * d.b * d.c)
end
function silhouette(sh::Ellipsoid, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette(sh, d.a, d.b, d.c, θ)
end

# Radius accessors — shared by every ellipsoidal shape (`Ellipsoid`,
# `HalfEllipsoid`); all store the same `b_semi_minor_skin` /
# `b_semi_minor_fibrous` fields, so the dispatch lives once on the family type.

_skin_radius(::AbstractEllipsoidal, length) = length.b_semi_minor_skin
_fibrous_radius(::AbstractEllipsoidal, length) = length.b_semi_minor_fibrous

# Composition

attachment_surfaces(::Ellipsoid) = (PoleA, PoleB, Equator)

# Outer (insulation-aware) semi-axes (a, b, c). Insulation-dispatched.
outer_dims(sh::Ellipsoid, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Ellipsoid, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_skin,
     b = body.geometry.length.b_semi_minor_skin,
     c = body.geometry.length.c_semi_minor_skin)
outer_dims(::Ellipsoid, ::FibrousLayer, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_fibrous,
     b = body.geometry.length.b_semi_minor_fibrous,
     c = body.geometry.length.c_semi_minor_fibrous)

# Skin-level semi-axes — used for flesh-anchored attachment positions.
function _ellipsoid_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.a_semi_major_skin, gl.b_semi_minor_skin, gl.c_semi_minor_skin)
end

# Notional area for pole attachments. For a truncated pole_a this is the
# actual flat-disc area exposed by the cut; for an untruncated pole it's
# the cross-sectional disc through the pole (used as a loose upper bound).
function surface_area(sh::Ellipsoid, body::AbstractBody, ::PoleA)
    d = outer_dims(sh, body)
    s = _pole_a_radial_scale(sh)
    s == 0 ? π * d.b * d.c : π * (s * d.b) * (s * d.c)
end
function surface_area(sh::Ellipsoid, body::AbstractBody, ::PoleB)
    d = outer_dims(sh, body)
    π * d.b * d.c
end

# Loose bound for equator joins: full ellipsoid surface area.
surface_area(::Ellipsoid, body::AbstractBody, ::Equator) =
    body.geometry.area.total

# For untruncated pole_a, coordinate fields are irrelevant (point pole).
# For truncated pole_a, PoleA may be bare (centre) or PoleA(r, φ) on the disc.
validate_range(::Ellipsoid, ::AbstractBody, ::PoleA) = nothing
validate_range(::Ellipsoid, ::AbstractBody, ::PoleB) = nothing
validate_range(::Ellipsoid, ::AbstractBody, ::Equator) = nothing

function surface_point(sh::Ellipsoid, body::AbstractBody, ::PoleA)
    a, _, _ = _ellipsoid_skin(body)
    (a * _pole_a_x_ratio(sh), zero(a), zero(a))
end
function surface_point(::Ellipsoid, body::AbstractBody, ::PoleB)
    a, _, _ = _ellipsoid_skin(body)
    (-a, zero(a), zero(a))
end
function surface_point(::Ellipsoid, body::AbstractBody, loc::Equator)
    _, b, c = _ellipsoid_skin(body)
    (zero(b), b * cos(loc.φ), c * sin(loc.φ))
end

surface_normal(::Ellipsoid, ::AbstractBody, ::PoleA) = ( 1.0, 0.0, 0.0)
surface_normal(::Ellipsoid, ::AbstractBody, ::PoleB) = (-1.0, 0.0, 0.0)
surface_normal(::Ellipsoid, ::AbstractBody, loc::Equator) =
    (0.0, cos(loc.φ), sin(loc.φ))

function surface_centroid(sh::Ellipsoid, body::AbstractBody, ::PoleA)
    a, _, _ = _ellipsoid_skin(body); (a * _pole_a_x_ratio(sh), zero(a), zero(a))
end
function surface_centroid(::Ellipsoid, body::AbstractBody, ::PoleB)
    a, _, _ = _ellipsoid_skin(body); (-a, zero(a), zero(a))
end
function surface_centroid(::Ellipsoid, body::AbstractBody, ::Equator)
    _, b, _ = _ellipsoid_skin(body); (zero(b), b, zero(b))  # arbitrary point on equator
end
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::PoleA) = ( 1.0, 0.0, 0.0)
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::PoleB) = (-1.0, 0.0, 0.0)
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::Equator) = (0.0, 1.0, 0.0)
