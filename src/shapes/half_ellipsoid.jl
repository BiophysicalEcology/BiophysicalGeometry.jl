"""
    HalfEllipsoid <: AbstractShape

A half-ellipsoid organism shape — a prolate ellipsoid (b = c) cut by the
plane perpendicular to the minor axis through the centre. Dimensions are
sized so that two half-ellipsoids of equal mass joined at their flat faces
reconstruct the full `Ellipsoid(2*mass, density, b, c)`.

In local coordinates the long axis lies along `+x`, and the flat face is
the elliptical disc at `z = 0` with the dome in `z ≥ 0`.

The flat-face area is reported at *skin* level even when fur is present,
so two halves with mixed insulation can join cleanly with `FullCover`.
"""
mutable struct HalfEllipsoid{M,D,B,C} <: AbstractShape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::HalfEllipsoid, ::Naked)
    volume = shape.mass / shape.density
    b_semi_minor_skin = (3 * volume / (2π * shape.b))^(1/3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    a = ustrip(u"m", a_semi_major_skin)
    b = ustrip(u"m", b_semi_minor_skin)
    c = ustrip(u"m", c_semi_minor_skin)
    half_dome = _half_ellipsoid_dome(a, b, c)
    flat = π * b * c
    total = (half_dome + flat) * u"m^2"
    characteristic_dimension = volume^(1/3)
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin),
                    SurfaceAreas(; total))
end
function geometry(shape::HalfEllipsoid, fur::Fur)
    volume = shape.mass / shape.density
    b_semi_minor_skin = (3 * volume / (2π * shape.b))^(1/3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    as = ustrip(u"m", a_semi_major_skin)
    bs = ustrip(u"m", b_semi_minor_skin)
    cs = ustrip(u"m", c_semi_minor_skin)
    af = ustrip(u"m", a_semi_major_fur)
    bf = ustrip(u"m", b_semi_minor_fur)
    cf = ustrip(u"m", c_semi_minor_fur)
    flat = π * bs * cs   # at skin level for FullCover correctness
    total = (_half_ellipsoid_dome(af, bf, cf) + flat) * u"m^2"
    skin  = (_half_ellipsoid_dome(as, bs, cs) + flat) * u"m^2"
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1/3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
                       a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur),
                    SurfaceAreas(; total, skin, convection))
end
function geometry(shape::HalfEllipsoid, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (3 * flesh_volume / (2π * shape.b))^(1/3)
    fat_thickness = (((3 * volume / (2π * shape.b))^(1/3)) - b_flesh)
    if fat_thickness <= 0.0u"m"
        b_semi_minor_skin = (3 * volume / (2π * shape.b))^(1/3)
    else
        b_semi_minor_skin = b_flesh + fat_thickness
    end
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    a = ustrip(u"m", a_semi_major_skin)
    b = ustrip(u"m", b_semi_minor_skin)
    c = ustrip(u"m", c_semi_minor_skin)
    total = (_half_ellipsoid_dome(a, b, c) + π * b * c) * u"m^2"
    characteristic_dimension = volume^(1/3)
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin, fat=fat_thickness),
                    SurfaceAreas(; total))
end
function geometry(shape::HalfEllipsoid, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (3 * flesh_volume / (2π * shape.b))^(1/3)
    fat_thickness = ((3 * volume / (2π * shape.b))^(1/3)) - b_flesh
    if fat_thickness <= 0.0u"m"
        b_semi_minor_skin = (3 * volume / (2π * shape.b))^(1/3)
    else
        b_semi_minor_skin = b_flesh + fat_thickness
    end
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    a_semi_major_fur = a_semi_major_skin + fur.thickness
    b_semi_minor_fur = b_semi_minor_skin + fur.thickness
    c_semi_minor_fur = c_semi_minor_skin + fur.thickness
    as = ustrip(u"m", a_semi_major_skin)
    bs = ustrip(u"m", b_semi_minor_skin)
    cs = ustrip(u"m", c_semi_minor_skin)
    af = ustrip(u"m", a_semi_major_fur)
    bf = ustrip(u"m", b_semi_minor_fur)
    cf = ustrip(u"m", c_semi_minor_fur)
    flat = π * bs * cs
    total = (_half_ellipsoid_dome(af, bf, cf) + flat) * u"m^2"
    skin  = (_half_ellipsoid_dome(as, bs, cs) + flat) * u"m^2"
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1/3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
                       a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur, fat=fat_thickness),
                    SurfaceAreas(; total, skin, convection))
end

# Dome area for one half of a prolate spheroid (a > c, b = c assumed).
# Inputs are dimensionless metres; output is dimensionless m².
function _half_ellipsoid_dome(a, b, c)
    if abs(a - c) < a * 1e-9
        return 2 * π * b^2  # half-sphere
    end
    e = sqrt(a^2 - c^2) / a
    return π * b^2 + π * (a * b / e) * asin(e)
end

# Silhouette area: half of full ellipsoid's, so two halves reproduce the full.
function silhouette_area(::HalfEllipsoid, a, b, c, θ)
    a2 = cos(90)^2 * (cos(θ)^2 / a^2 + sin(θ)^2 / b^2) + sin(90)^2 / c^2
    twohh = 2 * cos(90) * sin(90) * cos(θ) * (1 / b^2 - 1 / a^2)
    b2 = sin(θ)^2 / a^2 + cos(θ)^2 / b^2
    θ2 = 0.5 * atan(twohh, a2 - b2)
    sps = sin(θ2); cps = cos(θ2)
    a3 = cps * (a2 * cps + twohh * sps) + b2 * sps * sps
    b3 = sps * (a2 * sps - twohh * cps) + b2 * cps * cps
    semax1 = 1 / sqrt(a3)
    semax2 = 1 / sqrt(b3)
    π * semax1 * semax2 / 2
end

function silhouette_area(::HalfEllipsoid, ::Union{Naked,Fat}, body::AbstractBody)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    c = body.geometry.length.c_semi_minor_skin
    (; normal = π * a * b / 2, parallel = π * b * c / 2)
end
function silhouette_area(::HalfEllipsoid, ::Union{Fur,CompositeInsulation}, body::AbstractBody)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    (; normal = π * a * b / 2, parallel = π * b * c / 2)
end
function silhouette_area(sh::HalfEllipsoid, ::Naked, body::AbstractBody, θ)
    silhouette_area(sh,
        body.geometry.length.a_semi_major_skin,
        body.geometry.length.b_semi_minor_skin,
        body.geometry.length.c_semi_minor_skin, θ)
end
function silhouette_area(sh::HalfEllipsoid, ::Fur, body::AbstractBody, θ)
    silhouette_area(sh,
        body.geometry.length.a_semi_major_fur,
        body.geometry.length.b_semi_minor_fur,
        body.geometry.length.c_semi_minor_fur, θ)
end
function silhouette_area(sh::HalfEllipsoid, ::Fat, body::AbstractBody, θ)
    silhouette_area(sh,
        body.geometry.length.a_semi_major_skin,
        body.geometry.length.b_semi_minor_skin,
        body.geometry.length.c_semi_minor_skin, θ)
end
function silhouette_area(sh::HalfEllipsoid, ::CompositeInsulation, body::AbstractBody, θ)
    silhouette_area(sh,
        body.geometry.length.a_semi_major_fur,
        body.geometry.length.b_semi_minor_fur,
        body.geometry.length.c_semi_minor_fur, θ)
end

# Radii — same dispatch as Ellipsoid.
skin_radius(::HalfEllipsoid, ::AbstractInsulation, body) = body.geometry.length.b_semi_minor_skin

insulation_radius(::HalfEllipsoid, ::Naked, body) = body.geometry.length.b_semi_minor_skin
flesh_radius(::HalfEllipsoid, ::Naked, body) = body.geometry.length.b_semi_minor_skin

insulation_radius(::HalfEllipsoid, ::Fur, body) = body.geometry.length.b_semi_minor_fur
flesh_radius(::HalfEllipsoid, ::Fur, body) = body.geometry.length.b_semi_minor_skin

insulation_radius(::HalfEllipsoid, ::Fat, body) = body.geometry.length.b_semi_minor_skin
flesh_radius(::HalfEllipsoid, ::Fat, body) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat

insulation_radius(::HalfEllipsoid, ::CompositeInsulation, body) = body.geometry.length.b_semi_minor_fur
flesh_radius(::HalfEllipsoid, ::CompositeInsulation, body) = body.geometry.length.b_semi_minor_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::HalfEllipsoid) = (:dome, :flat)

# Outer (insulation-aware) semi-axes of the dome — drives surface_area.
function _halfellip_outer(body::AbstractBody)
    gl = body.geometry.length
    if haskey(gl, :a_semi_major_fur)
        (gl.a_semi_major_fur, gl.b_semi_minor_fur, gl.c_semi_minor_fur)
    else
        (gl.a_semi_major_skin, gl.b_semi_minor_skin, gl.c_semi_minor_skin)
    end
end

# Skin-level semi-axes — drives flesh-anchored attachment positions.
function _halfellip_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.a_semi_major_skin, gl.b_semi_minor_skin, gl.c_semi_minor_skin)
end

# Loose bound for dome attachments.
surface_area(::HalfEllipsoid, body::AbstractBody, ::Val{:dome}) =
    body.geometry.area.total

# Flat at skin level so two halves with mixed insulation join cleanly.
function surface_area(::HalfEllipsoid, body::AbstractBody, ::Val{:flat})
    π * body.geometry.length.b_semi_minor_skin * body.geometry.length.c_semi_minor_skin
end

function validate_position(::HalfEllipsoid, ::AbstractBody, ::Val{:dome}, pos)
    issetequal(keys(pos), (:α, :β)) ||
        error(":dome needs (α, β); got $(keys(pos))")
    0 ≤ pos.α ≤ π || error(":dome α out of range [0, π]: $(pos.α)")
    0 ≤ pos.β ≤ π || error(":dome β out of range [0, π]: $(pos.β)")
end
function validate_position(::HalfEllipsoid, body::AbstractBody, ::Val{:flat}, pos)
    issetequal(keys(pos), (:x, :y)) || error(":flat needs (x, y); got $(keys(pos))")
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    abs(pos.x) ≤ a || error(":flat x out of range ±$a: $(pos.x)")
    abs(pos.y) ≤ b || error(":flat y out of range ±$b: $(pos.y)")
    rel = (pos.x / a)^2 + (pos.y / b)^2
    rel ≤ 1 + 1e-9 || error(":flat (x,y) outside ellipse boundary")
end

function surface_point(::HalfEllipsoid, body::AbstractBody, ::Val{:dome}, pos)
    a, b, c = _halfellip_skin(body)
    (a * cos(pos.α), b * sin(pos.α) * cos(pos.β), c * sin(pos.α) * sin(pos.β))
end
function surface_point(::HalfEllipsoid, body::AbstractBody, ::Val{:flat}, pos)
    (pos.x, pos.y, zero(pos.x))
end

function surface_normal(::HalfEllipsoid, body::AbstractBody, ::Val{:dome}, pos)
    a, b, c = _halfellip_skin(body)
    nx = cos(pos.α) / ustrip(u"m", a)
    ny = sin(pos.α) * cos(pos.β) / ustrip(u"m", b)
    nz = sin(pos.α) * sin(pos.β) / ustrip(u"m", c)
    n = sqrt(nx^2 + ny^2 + nz^2)
    (nx / n, ny / n, nz / n)
end
surface_normal(::HalfEllipsoid, ::AbstractBody, ::Val{:flat}, _) = (0.0, 0.0, -1.0)

# Centroids — flesh level.
function surface_centroid(::HalfEllipsoid, body::AbstractBody, ::Val{:dome})
    a, _, c = _halfellip_skin(body)
    (zero(a), zero(a), c)  # top of the dome
end
function surface_centroid(::HalfEllipsoid, body::AbstractBody, ::Val{:flat})
    a = body.geometry.length.a_semi_major_skin
    (zero(a), zero(a), zero(a))
end
surface_centroid_normal(::HalfEllipsoid, ::AbstractBody, ::Val{:dome}) = (0.0, 0.0, 1.0)
surface_centroid_normal(::HalfEllipsoid, ::AbstractBody, ::Val{:flat}) = (0.0, 0.0, -1.0)
