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

_ellipsoid_eccentricity(a, c) = sqrt(a^2 - c^2) / a

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
    a_semi_major_fur = a_semi_major_skin + fibrous_layer.thickness
    b_semi_minor_fur = b_semi_minor_skin + fibrous_layer.thickness
    c_semi_minor_fur = c_semi_minor_skin + fibrous_layer.thickness
    fur_axes = (a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur)
    areas = fibrous_areas(shape, fibrous_layer,
        _ellipsoid_surface_args(skin_axes),
        _ellipsoid_surface_args(fur_axes))
    return Geometry(volume,
        (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
           a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur),
        areas)
end
function geometry(shape::Ellipsoid, fat_layer::FatLayer)
    volume = _ellipsoid_volume(shape)
    fat_v = _ellipsoid_fat_volume(shape, fat_layer)
    flesh_v = volume - fat_v
    (a_flesh, b_flesh, c_flesh) = _ellipsoid_skin_axes(shape, flesh_v)
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_v),
        ustrip(u"m^3", fat_v),
        shape.axis_ratio_b,
        ustrip(u"m", b_flesh))
    skin_axes = if fat <= 0.0u"m"
        _ellipsoid_skin_axes(shape, volume)
    else
        (a_flesh + fat, b_flesh + fat, c_flesh + fat)
    end
    (a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin) = skin_axes
    total = surface_area(shape, _ellipsoid_surface_args(skin_axes)...)
    return Geometry(volume, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Ellipsoid, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = _ellipsoid_volume(shape)
    fat_v = _ellipsoid_fat_volume(shape, fat_layer)
    flesh_v = volume - fat_v
    (a_flesh, b_flesh, c_flesh) = _ellipsoid_skin_axes(shape, flesh_v)
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_v),
        ustrip(u"m^3", fat_v),
        shape.axis_ratio_b,
        ustrip(u"m", b_flesh))
    skin_axes = if fat <= 0.0u"m"
        _ellipsoid_skin_axes(shape, volume)
    else
        (a_flesh + fat, b_flesh + fat, c_flesh + fat)
    end
    (a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin) = skin_axes
    a_semi_major_fur = a_semi_major_skin + fibrous_layer.thickness
    b_semi_minor_fur = b_semi_minor_skin + fibrous_layer.thickness
    c_semi_minor_fur = c_semi_minor_skin + fibrous_layer.thickness
    fur_axes = (a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur)
    areas = fibrous_areas(shape, fibrous_layer,
        _ellipsoid_surface_args(skin_axes),
        _ellipsoid_surface_args(fur_axes))
    return Geometry(volume,
        (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
           a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur, fat),
        areas)
end

# Fat thickness calculation

function prolate_fat_layer(
    flesh_volume,
    fat_volume,
    axis_ratio_b,
    semi_minor_flesh
)
    # Flesh is approximated as a prolate spheroid:
    # Volume = 4/3 π * A * B * C
    # C = B, A = axis_ratio_b * B
    # → B = ((3 * volume) / (4 * axis_ratio_b * π))^(1//3)
    # Fat thickness X is root of cubic: A X^3 + B X^2 + C X + D = 0
    A = 1.0
    B = axis_ratio_b * semi_minor_flesh + 2 * semi_minor_flesh
    C = 2 * axis_ratio_b * semi_minor_flesh^2 + semi_minor_flesh^2
    D = axis_ratio_b * semi_minor_flesh^3 - (( (fat_volume + flesh_volume) * 3.0 ) / (4.0 * π))

    # Components of cubic formula

    T1a = (-B)^3 / (27 * A^3)
    T1b = (B * C) / (6 * A^2)
    T1c = D / (2 * A)
    T1  = T1a + T1b - T1c

    T2a = T1^2
    T2b = ( (C / (3*A)) - (B^2) / (9 * A^2) )^3

    # Prevent sqrt of negative number
    T2 = (T2a + T2b >= 0) ? sqrt(T2a + T2b) : 0.0

    T3 = B / (3*A)

    # Cube roots with sign handling

    function signed_cuberoot(x)
        # cbrt'(0) is +∞, which is mathematically correct but propagates as
        # NaN through Enzyme when chained with very-small (Cardano formula
        # discriminant) values. The cubic root of an effectively-zero
        # discriminant is itself effectively zero with effectively-zero
        # contribution to the cubic's solution, so clamp the singular point
        # to zero to keep the derivative finite.
        abs(x) < 1e-12 && return zero(x)
        x < 0 ? -((-x)^(1//3)) : x^(1//3)
    end

    root1 = signed_cuberoot(T1 + T2)
    root2 = signed_cuberoot(T1 - T2)

    fat = (root1 + root2 - T3) * u"m"

    # If negative, not enough fat to cover spheroid
    return max(0.0u"m", fat)
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
function surface_area(shape::Ellipsoid, a, b, c, e)
    (2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)) * u"m^2"
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
    body.geometry.length.a_semi_major_fur,
    body.geometry.length.b_semi_minor_fur,
    body.geometry.length.c_semi_minor_fur,
)

# Radius accessors

_skin_radius(::Ellipsoid, length) = length.b_semi_minor_skin
_fur_radius(::Ellipsoid, length) = length.b_semi_minor_fur
