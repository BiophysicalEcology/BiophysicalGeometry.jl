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

# Dome area for one half of a prolate spheroid (b = c assumed). The
# formula uses sqrt / asin on dimensionless eccentricity, so this is the
# one place we step out of Unitful — ustrip once, compute, re-wrap.
function _half_dome_area(a, b, c)
    am = ustrip(u"m", a); bm = ustrip(u"m", b); cm = ustrip(u"m", c)
    dome = if abs(am - cm) < am * 1e-9
        2 * π * bm^2  # half-sphere
    else
        e = sqrt(am^2 - cm^2) / am
        π * bm^2 + π * (am * bm / e) * asin(e)
    end
    dome * u"m^2"
end

function geometry(shape::HalfEllipsoid, ::Naked)
    volume = shape.mass / shape.density
    b_semi_minor_skin = (3 * volume / (2π * shape.b))^(1/3)
    c_semi_minor_skin = b_semi_minor_skin
    a_semi_major_skin = b_semi_minor_skin * shape.b
    dome = _half_dome_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin)
    flat = π * b_semi_minor_skin * c_semi_minor_skin
    total = dome + flat
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
    flat = π * b_semi_minor_skin * c_semi_minor_skin  # skin level for FullCover
    total = _half_dome_area(a_semi_major_fur,  b_semi_minor_fur,  c_semi_minor_fur)  + flat
    skin  = _half_dome_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin) + flat
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
    total = _half_dome_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin) +
            π * b_semi_minor_skin * c_semi_minor_skin
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
    flat  = π * b_semi_minor_skin * c_semi_minor_skin
    total = _half_dome_area(a_semi_major_fur,  b_semi_minor_fur,  c_semi_minor_fur)  + flat
    skin  = _half_dome_area(a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin) + flat
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1/3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; a_semi_major_skin, b_semi_minor_skin, c_semi_minor_skin,
                       a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur, fat=fat_thickness),
                    SurfaceAreas(; total, skin, convection))
end

# Silhouette area: half of full ellipsoid's, so two halves reproduce the
# full. For an ellipsoid (a, b, c) the silhouette projected along direction
# d is an ellipse of area π·sqrt(b²c²·d_x² + a²c²·d_y² + a²b²·d_z²). With
# the sun in the equatorial plane at angle θ from the long (a) axis,
# d = (cos θ, sin θ, 0) gives π·c·sqrt(b²·cos²θ + a²·sin²θ).
function silhouette_area(::HalfEllipsoid, a, b, c, θ)
    π * c * sqrt(b^2 * cos(θ)^2 + a^2 * sin(θ)^2) / 2
end

function silhouette_area(sh::HalfEllipsoid, ::AbstractInsulation, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = π * d.a * d.b / 2, parallel = π * d.b * d.c / 2)
end
function silhouette_area(sh::HalfEllipsoid, ::AbstractInsulation, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette_area(sh, d.a, d.b, d.c, θ)
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

attachment_surfaces(::HalfEllipsoid) = (Dome, Flat)

# Outer (insulation-aware) semi-axes of the dome. Insulation-dispatched.
outer_dims(sh::HalfEllipsoid, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::HalfEllipsoid, ::Union{Naked,Fat}, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_skin,
     b = body.geometry.length.b_semi_minor_skin,
     c = body.geometry.length.c_semi_minor_skin)
outer_dims(::HalfEllipsoid, ::Fur, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_fur,
     b = body.geometry.length.b_semi_minor_fur,
     c = body.geometry.length.c_semi_minor_fur)

# Skin-level semi-axes — drives flesh-anchored attachment positions.
function _halfellip_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.a_semi_major_skin, gl.b_semi_minor_skin, gl.c_semi_minor_skin)
end

# Loose bound for dome attachments.
surface_area(::HalfEllipsoid, body::AbstractBody, ::Dome) =
    body.geometry.area.total

# Flat at skin level so two halves with mixed insulation join cleanly.
function surface_area(::HalfEllipsoid, body::AbstractBody, ::Flat)
    π * body.geometry.length.b_semi_minor_skin * body.geometry.length.c_semi_minor_skin
end

function validate_range(::HalfEllipsoid, ::AbstractBody, loc::Dome)
    0 ≤ loc.α ≤ π || error("Dome α out of range [0, π]: $(loc.α)")
    0 ≤ loc.β ≤ π || error("Dome β out of range [0, π]: $(loc.β)")
end
# For HalfEllipsoid, Flat uses (u=x, v=y).
function validate_range(::HalfEllipsoid, body::AbstractBody, loc::Flat)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    abs(loc.u) ≤ a || error("Flat x out of range ±$a: $(loc.u)")
    abs(loc.v) ≤ b || error("Flat y out of range ±$b: $(loc.v)")
    rel = (loc.u / a)^2 + (loc.v / b)^2
    rel ≤ 1 + 1e-9 || error("Flat (x,y) outside ellipse boundary")
end

function surface_point(::HalfEllipsoid, body::AbstractBody, loc::Dome)
    a, b, c = _halfellip_skin(body)
    (a * cos(loc.α), b * sin(loc.α) * cos(loc.β), c * sin(loc.α) * sin(loc.β))
end
function surface_point(::HalfEllipsoid, body::AbstractBody, loc::Flat)
    (loc.u, loc.v, zero(loc.u))
end

function surface_normal(::HalfEllipsoid, body::AbstractBody, loc::Dome)
    a, b, c = _halfellip_skin(body)
    # Gradient of (x/a)² + (y/b)² + (z/c)² = 1 is (2x/a², 2y/b², 2z/c²).
    # Scale each component by a common length (b·c, a·c, a·b) so units
    # cancel: multiply top and bottom by a·b·c to get a unitless vector.
    nx = cos(loc.α)            * (b * c)
    ny = sin(loc.α) * cos(loc.β) * (a * c)
    nz = sin(loc.α) * sin(loc.β) * (a * b)
    n = sqrt(nx^2 + ny^2 + nz^2)
    (nx / n, ny / n, nz / n)
end
surface_normal(::HalfEllipsoid, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)

# Centroids — flesh level.
function surface_centroid(::HalfEllipsoid, body::AbstractBody, ::Dome)
    a, _, c = _halfellip_skin(body)
    (zero(a), zero(a), c)  # top of the dome
end
function surface_centroid(::HalfEllipsoid, body::AbstractBody, ::Flat)
    a = body.geometry.length.a_semi_major_skin
    (zero(a), zero(a), zero(a))
end
surface_centroid_normal(::HalfEllipsoid, ::AbstractBody, ::Dome) = (0.0, 0.0, 1.0)
surface_centroid_normal(::HalfEllipsoid, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)
