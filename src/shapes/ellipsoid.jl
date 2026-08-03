"""
    Ellipsoid <: AbstractShape

An ellipsoidal organism shape.
"""
struct Ellipsoid{M,D,B,C} <: AbstractShape
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
end

# Canonicalise to m^3 to avoid weird unit ratios from mass/density (e.g. g/(kg/m^3))
# that propagate into cbrt() and trip up Enzyme's typeunstablerules path on the
# fat-layer branch where the if/else introduces Union-typed locals.
_ellipsoid_volume(shape::Ellipsoid) = uconvert(u"m^3", body_volume(shape))
_ellipsoid_fat_volume(shape::Ellipsoid, fat_layer::FatLayer) =
    uconvert(u"m^3", fat_volume(shape, fat_layer))

function _ellipsoid_skin_axes(shape::Ellipsoid, volume)
    b = cbrt((3 / 4) * volume / (π * shape.axis_ratio_b))
    c = b
    a = b * shape.axis_ratio_b
    return (a, b, c)
end

# Smooth Heaviside on a length: ≈ 1 for fat ≫ ε, ≈ 0 for fat ≪ -ε, smooth
# everywhere. Used to blend the "no fat" and "with fat" geometry branches.
@inline _smooth_step_meters(fat) = 0.5 * (1 + fat / sqrt(fat*fat + (1e-9u"m")^2))

# "a" is the axisymmetric (polar) axis, "c" (=b) the equatorial radius.
# Prolate (a>c): eccentricity of the meridian ellipse about its major axis a.
# Oblate (a<c): eccentricity about the equatorial major axis c instead.
_ellipsoid_eccentricity(a, c) = a >= c ? sqrt(a^2 - c^2) / a : sqrt(c^2 - a^2) / c

# Pass the (a, b, c, e) tuple into the unitless surface_area method.
_ellipsoid_surface_args((a, b, c)::Tuple) = (
    ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c), _ellipsoid_eccentricity(a, c),
)

function geometry(shape::Ellipsoid, ::Naked)
    volume = _ellipsoid_volume(shape)
    (a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin) = _ellipsoid_skin_axes(shape, volume)
    total = surface_area(shape, _ellipsoid_surface_args((a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin))...)
    return Geometry(volume, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin), SurfaceAreas(; total))
end
function geometry(shape::Ellipsoid, fibrous_layer::FibrousLayer)
    volume = _ellipsoid_volume(shape)
    skin_axes = _ellipsoid_skin_axes(shape, volume)
    (a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin) = skin_axes
    a_semi_major_fibrous = a_semi_major_skin + fibrous_layer.thickness
    b_semi_minor_fibrous = b_semi_minor_skin + fibrous_layer.thickness
    c_semi_minor_fibrous = c_semi_minor_skin + fibrous_layer.thickness
    fibrous_axes = (a_semi_major_fibrous, b_semi_minor_fibrous, c_semi_minor_fibrous)
    areas = fibrous_areas(shape, fibrous_layer,
        _ellipsoid_surface_args(skin_axes),
        _ellipsoid_surface_args(fibrous_axes))
    return Geometry(volume,
        (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
           a_semi_major_fibrous, b_semi_minor_fibrous, c_semi_minor_fibrous),
        areas)
end
function geometry(shape::Ellipsoid, fat_layer::FatLayer)
    volume = _ellipsoid_volume(shape)
    fat_v = _ellipsoid_fat_volume(shape, fat_layer)
    flesh_v = volume - fat_v
    (a_flesh, b_flesh, c_flesh) = _ellipsoid_skin_axes(shape, flesh_v)
    raw_fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_v),
        ustrip(u"m^3", fat_v),
        shape.axis_ratio_b,
        ustrip(u"m", b_flesh))
    # Smooth-blend the "not enough fat" fallback (full-volume axes, fat=0)
    # and the normal flesh+fat axes. The previous if/else was a value
    # discontinuity at raw_fat = 0 — bad for reverse-mode AD. A single
    # smooth Heaviside drives both the axes blend and the effective_fat
    # value, so they stay consistent at every raw_fat.
    (a_full, b_full, c_full) = _ellipsoid_skin_axes(shape, volume)
    _w = _smooth_step_meters(raw_fat)
    fat = _w * raw_fat                                  # smooth-clamped to ≥ 0
    a_semi_major_skin = _w * (a_flesh + raw_fat) + (1 - _w) * (a_full + zero(raw_fat))
    b_semi_minor_skin = _w * (b_flesh + raw_fat) + (1 - _w) * (b_full + zero(raw_fat))
    c_semi_minor_skin = _w * (c_flesh + raw_fat) + (1 - _w) * (c_full + zero(raw_fat))
    skin_axes = (a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin)
    total = surface_area(shape, _ellipsoid_surface_args(skin_axes)...)
    return Geometry(volume, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Ellipsoid, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = _ellipsoid_volume(shape)
    fat_v = _ellipsoid_fat_volume(shape, fat_layer)
    flesh_v = volume - fat_v
    (a_flesh, b_flesh, c_flesh) = _ellipsoid_skin_axes(shape, flesh_v)
    raw_fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_v),
        ustrip(u"m^3", fat_v),
        shape.axis_ratio_b,
        ustrip(u"m", b_flesh))
    (a_full, b_full, c_full) = _ellipsoid_skin_axes(shape, volume)
    _w = _smooth_step_meters(raw_fat)
    fat = _w * raw_fat
    a_semi_major_skin = _w * (a_flesh + raw_fat) + (1 - _w) * (a_full + zero(raw_fat))
    b_semi_minor_skin = _w * (b_flesh + raw_fat) + (1 - _w) * (b_full + zero(raw_fat))
    c_semi_minor_skin = _w * (c_flesh + raw_fat) + (1 - _w) * (c_full + zero(raw_fat))
    skin_axes = (a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin)
    a_semi_major_fibrous = a_semi_major_skin + fibrous_layer.thickness
    b_semi_minor_fibrous = b_semi_minor_skin + fibrous_layer.thickness
    c_semi_minor_fibrous = c_semi_minor_skin + fibrous_layer.thickness
    fibrous_axes = (a_semi_major_fibrous, b_semi_minor_fibrous, c_semi_minor_fibrous)
    areas = fibrous_areas(shape, fibrous_layer,
        _ellipsoid_surface_args(skin_axes),
        _ellipsoid_surface_args(fibrous_axes))
    return Geometry(volume,
        (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
           a_semi_major_fibrous, b_semi_minor_fibrous, c_semi_minor_fibrous, fat),
        areas)
end

# Fat thickness calculation

function prolate_fat_layer(
    flesh_volume,
    fat_volume,
    axis_ratio_b,
    semi_minor_flesh
)
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
    return max(0.0, X) * u"m"
end

# Surface area

function surface_area(shape::Ellipsoid, body::AbstractBody)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    return surface_area(shape, a, b, c)
end
function surface_area(shape::Ellipsoid, a, b, c)
    p =  1.6075
    return(4 * π * (((a ^ p * b ^ p + a ^ p * c ^ p + b ^ p * c ^ p)) / 3) ^ (1 / p))
end
# Exact spheroid surface area (b=c, the equatorial radius; a the polar axis).
# Prolate (a>c): standard prolate-spheroid formula. Oblate (a<c): the
# analogous formula (log replaces asin).
function surface_area(shape::Ellipsoid, a, b, c, e)
    return a >= c ?
        (2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)) * u"m^2" :
        (2 * π * b ^ 2 + π * (a ^ 2 / e) * log((1 + e) / (1 - e))) * u"m^2"
end

# Silhouette area

# Project the 3D ellipsoid onto the plane perpendicular to the view direction
# by setting up a 2D quadratic form (a2·x² + 2·twohh·xy + b2·y² = 1),
# diagonalising via a rotation by θ2, and reading off the semi-axes. The view
# direction is parameterised by two angles: the zenith θ passed in, plus a
# body-orientation angle that the original Julia port hardcoded as the bare
# numeric `90`. Julia treats `cos(90)` / `sin(90)` as 90 *radians* (≈ 5156°),
# almost certainly a port bug — the value was meant in degrees, so we use
# `cosd(90)` / `sind(90)` here.
#
# KNOWN ISSUES with this implementation (kept as-is pending clarification of the
# original convention):
#
#   1. With the body-orientation angle hardcoded, the function does not produce
#      π·a·b at θ = 0 (the footprint of a horizontal ellipsoid viewed straight
#      down), which is what most callers will expect from
#      `silhouette_area(body, ZenithAngleVarying(), 0)`. The orientation
#      convention here is something other than "long axis horizontal, θ from
#      vertical".
#   2. The body-orientation angle is not exposed in the API. To support arbitrary
#      body tilts the function would need to accept a second angle.
#
# A correct standard formula for the simple "horizontal prolate ellipsoid,
# θ = zenith from vertical" case is
#     π · a · b · c · √((sin θ/a)² + (cos θ/c)²)
# with endpoints π·a·b at θ = 0 and π·b·c at θ = π/2.
function silhouette_area(shape::Ellipsoid, a, b, c, θ)
    a2 = cosd(90) ^ 2 * (cos(θ) ^ 2 / a ^ 2 + sin(θ) ^ 2 / b ^ 2) + sind(90) ^ 2 / c ^ 2
    twohh = 2 * cosd(90) * sind(90) * cos(θ) * (1 / b ^ 2 - 1 / a ^ 2)
    b2 = sin(θ) ^ 2 / a ^ 2 + cos(θ) ^ 2 / b ^ 2
    θ2 = 0.5 * atan(twohh, a2 - b2)
    sps = sin(θ2)
    cps = cos(θ2)
    a3 = cps * (a2 * cps + twohh * sps) + b2 * sps * sps
    b3 = sps * (a2 * sps - twohh * cps) + b2 * cps * cps
    semax1 = 1 / sqrt(a3)
    semax2 = 1 / sqrt(b3)
    π * semax1 * semax2
end

silhouette_area(shape::Ellipsoid, ins::AbstractInsulationLayer, body::AbstractBody, θ) =
    silhouette_area(shape, _ellipsoid_outer_axes(ins, body)..., θ)
function silhouette_area(shape::Ellipsoid, ins::AbstractInsulationLayer, body::AbstractBody)
    (a, b, c) = _ellipsoid_outer_axes(ins, body)
    return (; normal=π * a * b, parallel=π * b * c)
end

_ellipsoid_outer_axes(::Union{Naked,FatLayer}, body) = (
    body.geometry.length.a_semi_major_skin,
    body.geometry.length.b_semi_minor_skin,
    body.geometry.length.c_semi_minor_skin,
)
_ellipsoid_outer_axes(::Union{FibrousLayer,CompositeInsulation}, body) = (
    body.geometry.length.a_semi_major_fibrous,
    body.geometry.length.b_semi_minor_fibrous,
    body.geometry.length.c_semi_minor_fibrous,
)

# Radius accessors

_skin_radius(::Ellipsoid, length) = length.b_semi_minor_skin
_fibrous_radius(::Ellipsoid, length) = length.b_semi_minor_fibrous
