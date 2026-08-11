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
mutable struct Ellipsoid{M,D,B,C,T} <: AbstractEllipsoidal
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

# Prolate-spheroid surface area. The formula uses sqrt / asin on
# dimensionless eccentricity, so this is the one place we step out of
# Unitful — ustrip once, compute, re-wrap. `a, b, c` are lengths; result is
# an area.
function _prolate_area(a, b, c)
    am = ustrip(u"m", a); bm = ustrip(u"m", b); cm = ustrip(u"m", c)
    area = if abs(am - cm) < am * 1e-9
        4 * π * bm^2                       # sphere limit: a == c ⇒ e → 0, asin(e)/e → 1
    else
        e = sqrt(am^2 - cm^2) / am
        2 * π * bm^2 + 2 * π * (am * bm / e) * asin(e)
    end
    area * u"m^2"
end

# b-semi-minor from an enclosed volume. Identical to the equivalent-sphere radius
# of the volume scaled by the axis ratio, so the cube-root lives once (geometry.jl).
_ellipsoid_b(shape::Ellipsoid, volume) = _sphere_radius(volume / shape.axis_ratio_b)

function _skin_level(shape::Ellipsoid, volume)
    b = _ellipsoid_b(shape, volume)          # c = b; a = axis_ratio_b · b
    a = shape.axis_ratio_b * b
    (; dims = (; a_semi_major_skin = a, b_semi_minor_skin = b, c_semi_minor_skin = b),
       area = _prolate_area(a, b, b))
end
function _fibrous_level(shape::Ellipsoid, skin, thickness)
    a = skin.a_semi_major_skin + thickness
    b = skin.b_semi_minor_skin + thickness
    c = skin.c_semi_minor_skin + thickness
    (; dims = (; a_semi_major_fibrous = a, b_semi_minor_fibrous = b, c_semi_minor_fibrous = c),
       area = _prolate_area(a, b, c))
end
# Ellipsoid fat is a *uniform* shell: each skin semi-axis is the flesh semi-axis
# plus one fat thickness (so the skin is not a scaled flesh ellipsoid). Skin dims
# are therefore coupled to the cubic fat solve, so the fat paths override the
# generic orchestrator rather than compute skin from volume independently.
function _ellipsoid_fat_skin(shape::Ellipsoid, fat::FatLayer)
    volume = _body_volume(shape)
    flesh_volume = _flesh_volume(shape, fat)
    b_flesh = _ellipsoid_b(shape, flesh_volume)
    fat_thickness = prolate_fat_layer(flesh_volume, volume - flesh_volume, shape.axis_ratio_b, b_flesh)
    if fat_thickness <= 0.0u"m"
        b = _ellipsoid_b(shape, volume)
        a = shape.axis_ratio_b * b
    else
        b = b_flesh + fat_thickness
        a = shape.axis_ratio_b * b_flesh + fat_thickness
    end
    dims = (; a_semi_major_skin = a, b_semi_minor_skin = b, c_semi_minor_skin = b)
    (; volume, dims, fat = fat_thickness, area = _prolate_area(a, b, b))
end

function geometry(shape::Ellipsoid, fat::FatLayer)
    s = _ellipsoid_fat_skin(shape, fat)
    Geometry(s.volume, _characteristic_length(s.volume),
             merge(s.dims, (; fat = s.fat)), SurfaceAreas(; total = s.area))
end
function geometry(shape::Ellipsoid, fur::FibrousLayer, fat::FatLayer)
    s = _ellipsoid_fat_skin(shape, fat)
    fibrous = _fibrous_level(shape, s.dims, fur.thickness)
    Geometry(s.volume, _characteristic_length(s.volume) + fur.thickness,
             merge(s.dims, fibrous.dims, (; fat = s.fat)),
             SurfaceAreas(; total = fibrous.area, skin = s.area,
                          convection = _convective_area(fur, s.area)))
end

# FatLayer thickness calculation
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
    # FatLayer thickness X is root of cubic: A X^3 + B X^2 + C X + D = 0
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

# Radius — shared by every ellipsoidal shape (`Ellipsoid`, `HalfEllipsoid`);
# all store the same `b_semi_minor_skin` / `b_semi_minor_fibrous` / `fat`
# fields, so the dispatch lives once on the family type.

skin_radius(::AbstractEllipsoidal, ::AbstractInsulationLayer, body) = body.geometry.length.b_semi_minor_skin

# naked
insulation_radius(::AbstractEllipsoidal, ::Naked, body) = body.geometry.length.b_semi_minor_skin
flesh_radius(::AbstractEllipsoidal, ::Naked, body) = body.geometry.length.b_semi_minor_skin

# fur
insulation_radius(::AbstractEllipsoidal, ::FibrousLayer, body) = body.geometry.length.b_semi_minor_fibrous
flesh_radius(::AbstractEllipsoidal, ::FibrousLayer, body) = body.geometry.length.b_semi_minor_skin

# fat
insulation_radius(::AbstractEllipsoidal, ::FatLayer, body) = body.geometry.length.b_semi_minor_skin
flesh_radius(::AbstractEllipsoidal, ::FatLayer, body) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat

# fur and fat
insulation_radius(::AbstractEllipsoidal, ::CompositeInsulation, body) = body.geometry.length.b_semi_minor_fibrous
flesh_radius(::AbstractEllipsoidal, ::CompositeInsulation, body) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat

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
