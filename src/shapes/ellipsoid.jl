"""
    Ellipsoid(mass, density, b, c, pole_a_truncation=0.0) <: AbstractShape

A prolate ellipsoid (b = c) with semi-major axis along `+x`. With
`pole_a_truncation > 0`, the `+x` end is sliced flat: the cut plane sits at
`x = (1 - pole_a_truncation) * a`. `pole_a_truncation = 0` is a full ellipsoid;
`pole_a_truncation = 1` cuts through the centre (half ellipsoid).

The flat face exposed by the cut has semi-axes
`(b * sqrt(1 - (1 - pole_a_truncation)^2), c * sqrt(1 - (1 - pole_a_truncation)^2))`.

Surface area is reported for the full (untruncated) ellipsoid — this slightly
over-counts (by the removed spherical cap) for thin truncations.
"""
mutable struct Ellipsoid{M,D,B,C,T} <: AbstractShape
    mass::M
    density::D
    b::B
    c::C
    pole_a_truncation::T
end
Ellipsoid(mass, density, b, c) = Ellipsoid(mass, density, b, c, 0.0)

# x-position of the truncated pole_a (as a fraction of a). 1.0 = full ellipsoid.
_pole_a_x_ratio(s::Ellipsoid) = 1 - s.pole_a_truncation
# Radial scale at the truncated pole (so cut disc has y/z extents = scale * b/c).
_pole_a_radial_scale(s::Ellipsoid) = sqrt(max(0.0, 1 - _pole_a_x_ratio(s)^2))

# Prolate-spheroid surface area. The formula uses sqrt / asin on
# dimensionless eccentricity, so this is the one place we step out of
# Unitful — ustrip once, compute, re-wrap. `a, b, c` are lengths; result is
# an area.
function _prolate_area(a, b, c)
    am = ustrip(u"m", a); bm = ustrip(u"m", b); cm = ustrip(u"m", c)
    e = sqrt(am^2 - cm^2) / am
    (2 * π * bm^2 + 2 * π * (am * bm / e) * asin(e)) * u"m^2"
end

function geometry(shape::Ellipsoid, ::Naked)
    volume = shape.mass / shape.density
    b_semi_minor_skin = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    total = _prolate_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin)
    characteristic_dimension = volume^(1 / 3)
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin),
                    SurfaceAreas(; total))
end
function geometry(shape::Ellipsoid, fur::Fur)
    volume = shape.mass / shape.density
    b_semi_minor_skin = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    total = _prolate_area(a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur)
    skin = _prolate_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
                       a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur),
                    SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Ellipsoid, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (((3 / 4) * flesh_volume) / (π * shape.b)) ^ (1 / 3)
    c_flesh = b_flesh # assuming c = b
    a_flesh = shape.b * b_flesh
    fat = prolate_fat_layer(flesh_volume, fat_volume, shape.b, b_flesh)
    if fat <= 0.0u"m"
        b_semi_minor_skin = (((3 / 4) * volume) / (π * shape.b)) ^ (1 / 3)
        c_semi_minor_skin = b_semi_minor_skin
        a_semi_major_skin = shape.b * b_semi_minor_skin
    else
        a_semi_major_skin = a_flesh + fat
        b_semi_minor_skin = b_flesh + fat
        c_semi_minor_skin = c_flesh + fat
    end
    total = _prolate_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin)
    characteristic_dimension = volume^(1 / 3)
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, fat),
                    SurfaceAreas(; total))
end
function geometry(shape::Ellipsoid, fur::Fur, fat::Fat)
    # TODO reduce duplication here
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (((3 / 4) * flesh_volume) / (π * shape.b)) ^ (1 / 3)
    c_flesh = b_flesh # assuming c = b
    a_flesh = shape.b * b_flesh
    fat = prolate_fat_layer(flesh_volume, fat_volume, shape.b, b_flesh)
    if fat <= 0.0u"m"
        b_semi_minor_skin = (((3 / 4) * volume) / (π * shape.b)) ^ (1 / 3)
        c_semi_minor_skin = b_semi_minor_skin
        a_semi_major_skin = shape.b * b_semi_minor_skin
    else
        a_semi_major_skin = a_flesh + fat
        b_semi_minor_skin = b_flesh + fat
        c_semi_minor_skin = c_flesh + fat
    end
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    total = _prolate_area(a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur)
    skin = _prolate_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
                       a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur, fat),
                    SurfaceAreas(; total, skin, convection))
end

# Fat thickness calculation
#
# The cubic solve uses fractional powers, so the numeric core (`_prolate_fat_layer_m`)
# is unitless. `prolate_fat_layer` is the dimensional boundary: it ustrips
# volumes and radii once, calls the numeric core, and re-wraps the answer.

function prolate_fat_layer(flesh_volume, fat_volume, shape_b, semi_minor_flesh)
    fat_m = _prolate_fat_layer_m(
        ustrip(u"m^3", flesh_volume),
        ustrip(u"m^3", fat_volume),
        shape_b,
        ustrip(u"m", semi_minor_flesh),
    )
    return max(0.0u"m", fat_m * u"m")
end

function _prolate_fat_layer_m(flesh_volume, fat_volume, shape_b, semi_minor_flesh)
    # Flesh is approximated as a prolate spheroid:
    # Volume = 4/3 π * A * B * C   with C = B, A = shape_b * B
    # Fat thickness X is root of cubic: A X^3 + B X^2 + C X + D = 0
    A = 1.0
    B = shape_b * semi_minor_flesh + 2 * semi_minor_flesh
    C = 2 * shape_b * semi_minor_flesh^2 + semi_minor_flesh^2
    D = shape_b * semi_minor_flesh^3 - (((fat_volume + flesh_volume) * 3.0) / (4.0 * π))

    T1a = (-B)^3 / (27 * A^3)
    T1b = (B * C) / (6 * A^2)
    T1c = D / (2 * A)
    T1 = T1a + T1b - T1c

    T2a = T1^2
    T2b = ((C / (3*A)) - (B^2) / (9 * A^2))^3

    # Prevent sqrt of negative number
    T2 = (T2a + T2b >= 0) ? sqrt(T2a + T2b) : 0.0

    T3 = B / (3*A)

    signed_cuberoot(x) = x < 0 ? -((-x)^(1/3)) : x^(1/3)
    root1 = signed_cuberoot(T1 + T2)
    root2 = signed_cuberoot(T1 - T2)

    root1 + root2 - T3
end

# Surface area

function surface_area(shape::Ellipsoid, body::AbstractBody)
    _prolate_area(body.geometry.length.a_semi_major_skin,
                  body.geometry.length.b_semi_minor_skin,
                  body.geometry.length.c_semi_minor_skin)
end

# Silhouette area 

# For an ellipsoid (a, b, c) the silhouette projected along direction d
# is an ellipse of area π·sqrt(b²c²·d_x² + a²c²·d_y² + a²b²·d_z²). With
# the sun in the equatorial plane at angle θ from the long (a) axis,
# d = (cos θ, sin θ, 0) gives π·c·sqrt(b²·cos²θ + a²·sin²θ).
function silhouette_area(shape::Ellipsoid, a, b, c, θ)
    π * c * sqrt(b^2 * cos(θ)^2 + a^2 * sin(θ)^2)
end
function silhouette_area(sh::Ellipsoid, ::AbstractInsulation, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = π * d.a * d.b, parallel = π * d.b * d.c)
end
function silhouette_area(sh::Ellipsoid, ::AbstractInsulation, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette_area(sh, d.a, d.b, d.c, θ)
end

# Radius

skin_radius(shape::Ellipsoid, insulation::AbstractInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor_skin

# naked
insulation_radius(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
flesh_radius(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor_skin

# fur
insulation_radius(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor_fur
flesh_radius(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor_skin

# fat
insulation_radius(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor_skin
flesh_radius(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat

# fur and fat
insulation_radius(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor_fur
flesh_radius(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::Ellipsoid) = (PoleA, PoleB, Equator)

# Outer (insulation-aware) semi-axes (a, b, c). Insulation-dispatched.
outer_dims(sh::Ellipsoid, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Ellipsoid, ::Union{Naked,Fat}, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_skin,
     b = body.geometry.length.b_semi_minor_skin,
     c = body.geometry.length.c_semi_minor_skin)
outer_dims(::Ellipsoid, ::Fur, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_fur,
     b = body.geometry.length.b_semi_minor_fur,
     c = body.geometry.length.c_semi_minor_fur)

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
