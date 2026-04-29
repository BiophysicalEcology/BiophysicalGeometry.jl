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

function geometry(shape::Ellipsoid, ::Naked)
    volume = shape.mass / shape.density
    b_semi_minor_skin = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    e = ((ustrip(u"m", a_semi_major_skin) ^ 2 - ustrip(u"m", c_semi_minor_skin) ^ 2) ^ (1 / 2)) / ustrip(u"m", a_semi_major_skin)
    total = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", b_semi_minor_skin), e)
    characteristic_dimension = volume^(1 / 3) # b_semi_minor_skin * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin), SurfaceAreas(; total))
end
function geometry(shape::Ellipsoid, fur::Fur)
    volume = shape.mass / shape.density
    b_semi_minor_skin = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    e = ((ustrip(u"m", a_semi_major_skin) ^ 2 - ustrip(u"m", c_semi_minor_skin) ^ 2) ^ (1 / 2)) / ustrip(u"m", a_semi_major_skin)
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    e_fur = ((a_semi_major_fur ^ 2 - c_semi_minor_fur ^ 2) ^ (1 / 2)) / a_semi_major_fur
    total = surface_area(shape, ustrip(u"m", a_semi_major_fur),  ustrip(u"m", b_semi_minor_fur),  ustrip(u"m", c_semi_minor_fur),  e_fur)
    skin = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", c_semi_minor_skin), e)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    characteristic_dimension = volume^(1 / 3) + fur.thickness # b_semi_minor_fur * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Ellipsoid, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (((3 / 4) * flesh_volume) / (π * shape.b)) ^ (1 / 3)
    c_flesh = b_flesh # assuming c = b
    a_flesh = shape.b * b_flesh
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_volume), 
        ustrip(u"m^3", fat_volume), 
        shape.b, 
        ustrip(u"m", b_flesh))
    if fat <= 0.0u"m"
        b_semi_minor_skin = (((3 / 4) * volume) / (π * shape.b)) ^ (1 / 3)
        c_semi_minor_skin = b_semi_minor_skin # assuming c = b
        a_semi_major_skin = shape.b * b_semi_minor_skin
    else
        a_semi_major_skin = a_flesh + fat
        b_semi_minor_skin = b_flesh + fat
        c_semi_minor_skin = c_flesh + fat
    end
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    total = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", c_semi_minor_skin), e)
    characteristic_dimension = volume^(1 / 3) # b_semi_minor_skin * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, fat), SurfaceAreas(; total))
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
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_volume), 
        ustrip(u"m^3", fat_volume), 
        shape.b, 
        ustrip(u"m", b_flesh))
    if fat <= 0.0u"m"
        b_semi_minor_skin = (((3 / 4) * volume) / (π * shape.b)) ^ (1 / 3)
        c_semi_minor_skin = b_semi_minor_skin # assuming c = b
        a_semi_major_skin = shape.b * b_semi_minor_skin
    else
        a_semi_major_skin = a_flesh + fat
        b_semi_minor_skin = b_flesh + fat
        c_semi_minor_skin = c_flesh + fat
    end
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    e = ((a_semi_major_skin ^ 2 - c_semi_minor_skin ^ 2) ^ (1 / 2)) / a_semi_major_skin
    e_fur = ((a_semi_major_fur ^ 2 - c_semi_minor_fur ^ 2) ^ (1 / 2)) / a_semi_major_fur
    total = surface_area(shape, ustrip(u"m", a_semi_major_fur),  ustrip(u"m", b_semi_minor_fur),  ustrip(u"m", c_semi_minor_fur),  e_fur)
    skin = surface_area(shape, ustrip(u"m", a_semi_major_skin), ustrip(u"m", b_semi_minor_skin), ustrip(u"m", c_semi_minor_skin), e)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    characteristic_dimension = volume^(1 / 3) + fur.thickness # b_semi_minor_fur * 2
    return Geometry(volume, characteristic_dimension, (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur, fat), SurfaceAreas(; total, skin, convection))
end

# Fat thickness calculation

function prolate_fat_layer(
    flesh_volume,
    fat_volume,
    shape_b,
    semi_minor_flesh
)
    # Flesh is approximated as a prolate spheroid:
    # Volume = 4/3 π * A * B * C
    # C = B, A = shape_b * B
    # → B = ((3 * volume) / (4 * shape_b * π))^(1/3)
    # Fat thickness X is root of cubic: A X^3 + B X^2 + C X + D = 0
    A = 1.0
    B = shape_b * semi_minor_flesh + 2 * semi_minor_flesh
    C = 2 * shape_b * semi_minor_flesh^2 + semi_minor_flesh^2
    D = shape_b * semi_minor_flesh^3 - (( (fat_volume + flesh_volume) * 3.0 ) / (4.0 * π))

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
        x < 0 ? -((-x)^(1/3)) : x^(1/3)
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
    #e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    #2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
    p =  1.6075
    return(4 * π * (((a ^ p * b ^ p + a ^ p * c ^ p + b ^ p * c ^ p)) / 3) ^ (1 / p))
end
function surface_area(shape::Ellipsoid, a, b, c, e)
    (2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)) * u"m^2"
end

# Silhouette area 

# For an ellipsoid (a, b, c) the silhouette projected along direction d
# is an ellipse of area π·sqrt(b²c²·d_x² + a²c²·d_y² + a²b²·d_z²). With
# the sun in the equatorial plane at angle θ from the long (a) axis,
# d = (cos θ, sin θ, 0) gives π·c·sqrt(b²·cos²θ + a²·sin²θ).
function silhouette_area(shape::Ellipsoid, a, b, c, θ)
    π * c * sqrt(b^2 * cos(θ)^2 + a^2 * sin(θ)^2)
end
function silhouette_area(shape::Ellipsoid, insulation::Union{Naked,Fat}, body::AbstractBody)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    normal = π * a * b
    parallel = π * b * c
    return (; normal, parallel)
end
function silhouette_area(shape::Ellipsoid, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    normal = π * a * b
    parallel = π * b * c
    return (; normal, parallel)
end
function silhouette_area(shape::Ellipsoid, ::Naked, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    return silhouette_area(shape, a, b, c, θ)
end
function silhouette_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    return silhouette_area(shape, a, b, c, θ)
end
function silhouette_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    return silhouette_area(shape, a, b, c, θ)
end
function silhouette_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    return silhouette_area(shape, a, b, c, θ)
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

attachment_surfaces(::Ellipsoid) = (:pole_a, :pole_b, :equator)

# Outer semi-axes (a, b, c). Returns the fur-level if furred, else skin.
# Used by surface_area; attachment positions use _ellipsoid_skin instead.
function _ellipsoid_outer(body::AbstractBody)
    gl = body.geometry.length
    if haskey(gl, :a_semi_major_fur)
        (gl.a_semi_major_fur, gl.b_semi_minor_fur, gl.c_semi_minor_fur)
    else
        (gl.a_semi_major_skin, gl.b_semi_minor_skin, gl.c_semi_minor_skin)
    end
end

# Skin-level semi-axes — used for flesh-anchored attachment positions.
function _ellipsoid_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.a_semi_major_skin, gl.b_semi_minor_skin, gl.c_semi_minor_skin)
end

# Notional area for pole attachments. For a truncated pole_a this is the
# actual flat-disc area exposed by the cut; for an untruncated pole it's
# the cross-sectional disc through the pole (used as a loose upper bound).
function surface_area(sh::Ellipsoid, body::AbstractBody, ::Val{:pole_a})
    _, b, c = _ellipsoid_outer(body)
    s = _pole_a_radial_scale(sh)
    s == 0 ? π * b * c : π * (s * b) * (s * c)
end
function surface_area(::Ellipsoid, body::AbstractBody, ::Val{:pole_b})
    _, b, c = _ellipsoid_outer(body)
    π * b * c
end

# Loose bound for equator joins: full ellipsoid surface area.
surface_area(::Ellipsoid, body::AbstractBody, ::Val{:equator}) =
    body.geometry.area.total

function validate_position(sh::Ellipsoid, ::AbstractBody, ::Val{:pole_a}, pos)
    if sh.pole_a_truncation == 0
        isempty(keys(pos)) || error(":pole_a position must be empty (;); got $(keys(pos))")
    else
        # Truncated pole_a is a flat disc — accept (;) (centre) or (r, φ).
        isempty(keys(pos)) || issetequal(keys(pos), (:r, :φ)) ||
            error(":pole_a (truncated) needs (;) or (r, φ); got $(keys(pos))")
    end
end
validate_position(::Ellipsoid, body::AbstractBody, ::Val{:pole_b}, pos) =
    validate_position(shape(body), body, Val(:pole_a), pos)
function validate_position(::Ellipsoid, ::AbstractBody, ::Val{:equator}, pos)
    issetequal(keys(pos), (:φ,)) ||
        error(":equator position needs (φ,); got $(keys(pos))")
end

function surface_point(sh::Ellipsoid, body::AbstractBody, ::Val{:pole_a}, _)
    a, _, _ = _ellipsoid_skin(body)
    (a * _pole_a_x_ratio(sh), zero(a), zero(a))
end
function surface_point(::Ellipsoid, body::AbstractBody, ::Val{:pole_b}, _)
    a, _, _ = _ellipsoid_skin(body)
    (-a, zero(a), zero(a))
end
function surface_point(::Ellipsoid, body::AbstractBody, ::Val{:equator}, pos)
    _, b, c = _ellipsoid_skin(body)
    (zero(b), b * cos(pos.φ), c * sin(pos.φ))
end

surface_normal(::Ellipsoid, ::AbstractBody, ::Val{:pole_a}, _) = ( 1.0, 0.0, 0.0)
surface_normal(::Ellipsoid, ::AbstractBody, ::Val{:pole_b}, _) = (-1.0, 0.0, 0.0)
surface_normal(::Ellipsoid, ::AbstractBody, ::Val{:equator}, pos) =
    (0.0, cos(pos.φ), sin(pos.φ))

function surface_centroid(sh::Ellipsoid, body::AbstractBody, ::Val{:pole_a})
    a, _, _ = _ellipsoid_skin(body); (a * _pole_a_x_ratio(sh), zero(a), zero(a))
end
function surface_centroid(::Ellipsoid, body::AbstractBody, ::Val{:pole_b})
    a, _, _ = _ellipsoid_skin(body); (-a, zero(a), zero(a))
end
function surface_centroid(::Ellipsoid, body::AbstractBody, ::Val{:equator})
    _, b, _ = _ellipsoid_skin(body); (zero(b), b, zero(b))  # arbitrary point on equator
end
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::Val{:pole_a}) = ( 1.0, 0.0, 0.0)
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::Val{:pole_b}) = (-1.0, 0.0, 0.0)
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::Val{:equator}) = (0.0, 1.0, 0.0)
