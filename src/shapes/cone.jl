"""
    Cone(mass, density, b, top_ratio=0.0) <: AbstractShape

A right circular cone or truncated cone (frustum). `b` is the length to
base-diameter ratio. `top_ratio` is the ratio of apex (top) radius to base
radius — `0` makes a sharp cone, values in `(0, 1)` make a frustum.

Local frame: axis along `+z`, base disc at `z = 0` with radius `radius_skin`,
top disc at `z = length_skin` with radius `top_ratio * radius_skin`.
Insulation expands radii and length as for `Cylinder`; attachment positions
stay at flesh level.
"""
mutable struct Cone{M,D,B,T} <: AbstractCylindrical
    mass::M
    density::D
    axis_ratio_b::B
    top_ratio::T
end
Cone(mass, density, b) = Cone(mass, density, b, 0.0)

# Volume of a frustum = (π/3) · L · (R² + R·r + r²) with r = top_ratio·R.
# With L = 2·b·R: V = (2π/3) · b · R³ · (1 + t + t²).
_cone_volume_factor(t) = 1 + t + t^2

function _cone_radius(volume, b, t)
    (3 * volume / (2π * b * _cone_volume_factor(t)))^(1/3)
end

function geometry(shape::Cone, ::Naked)
    volume = shape.mass / shape.density
    t = shape.top_ratio
    radius_skin = _cone_radius(volume, shape.axis_ratio_b, t)
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1/3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin),
                    SurfaceAreas(; total))
end
function geometry(shape::Cone, fur::FibrousLayer)
    volume = shape.mass / shape.density
    t = shape.top_ratio
    radius_skin = _cone_radius(volume, shape.axis_ratio_b, t)
    radius_fibrous = radius_skin + fur.thickness
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    length_fibrous = length_skin + 2 * fur.thickness
    total = surface_area(shape, radius_fibrous, length_fibrous)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1/3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fibrous, length_skin, length_fibrous),
                    SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Cone, fat::FatLayer)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    t = shape.top_ratio
    radius_skin = _cone_radius(volume, shape.axis_ratio_b, t)
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    radius_flesh = _cone_radius(flesh_volume, shape.axis_ratio_b, t)
    fat_thickness = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1/3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin, fat=fat_thickness),
                    SurfaceAreas(; total))
end
function geometry(shape::Cone, fur::FibrousLayer, fat::FatLayer)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    t = shape.top_ratio
    radius_skin = _cone_radius(volume, shape.axis_ratio_b, t)
    radius_fibrous = radius_skin + fur.thickness
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    length_fibrous = length_skin + 2 * fur.thickness
    radius_flesh = _cone_radius(flesh_volume, shape.axis_ratio_b, t)
    fat_thickness = radius_skin - radius_flesh
    total = surface_area(shape, radius_fibrous, length_fibrous)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1/3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fibrous, length_skin, length_fibrous, fat=fat_thickness),
                    SurfaceAreas(; total, skin, convection))
end

# Aggregate surface area for a frustum: base disc + top disc + slant surface.
# Slant length s = sqrt((R - r)² + L²) where r = t·R.
function surface_area(shape::Cone, R, L)
    t = shape.top_ratio
    r = t * R
    s = sqrt((R - r)^2 + L^2)
    π * R^2 + π * r^2 + π * (R + r) * s
end

# Silhouette: same shape as the cylinder pattern, but the lateral profile is a
# triangle (r·L) rather than a rectangle (2·r·L). `outer_dims` picks skin- vs
# fibrous-level (r, L) so one wrapper per arity covers all insulation kinds.
function silhouette_area(::Cone, r, L, θ)
    r * L * sin(θ) + π * r^2 * cos(θ)
end
function silhouette_area(sh::Cone, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette_area(sh, d.r, d.L, θ)
end
function silhouette_area(sh::Cone, ::AbstractInsulationLayer, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = d.r * d.L, parallel = π * d.r^2)
end

# Radii come from the shared `AbstractCylindrical` dispatch in cylinder.jl.

# Composition
#
# `EndA` is the base disc (z=0), `EndB` is the top disc (z=length_skin,
# radius = top_ratio * radius_skin), `Lateral` is the slant surface.

attachment_surfaces(::Cone) = (EndA, EndB, Lateral)

# Outer (insulation-aware) dimensions.
outer_dims(sh::Cone, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Cone, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (r = body.geometry.length.radius_skin, L = body.geometry.length.length_skin)
outer_dims(::Cone, ::FibrousLayer, body::AbstractBody) =
    (r = body.geometry.length.radius_fibrous, L = body.geometry.length.length_fibrous)

function surface_area(sh::Cone, body::AbstractBody, ::EndA)
    d = outer_dims(sh, body); π * d.r^2
end
function surface_area(sh::Cone, body::AbstractBody, ::EndB)
    d = outer_dims(sh, body); π * (sh.top_ratio * d.r)^2
end
function surface_area(sh::Cone, body::AbstractBody, ::Lateral)
    d = outer_dims(sh, body)
    r = sh.top_ratio * d.r
    s = sqrt((d.r - r)^2 + d.L^2)
    π * (d.r + r) * s
end

function validate_range(::Cone, body::AbstractBody, loc::EndA)
    R = body.geometry.length.radius_skin
    loc.r ≥ zero(loc.r) && loc.r ≤ R || error("EndA r out of range [0, $R]: $(loc.r)")
end
function validate_range(shape::Cone, body::AbstractBody, loc::EndB)
    Rt = shape.top_ratio * body.geometry.length.radius_skin
    Rt > zero(Rt) || error("EndB has zero radius (top_ratio=0); use a Disc(0) only")
    loc.r ≥ zero(loc.r) && loc.r ≤ Rt || error("EndB r out of range [0, $Rt]: $(loc.r)")
end
function validate_range(::Cone, body::AbstractBody, loc::Lateral)
    L = body.geometry.length.length_skin
    loc.z ≥ zero(loc.z) && loc.z ≤ L || error("Lateral z out of range [0, $L]: $(loc.z)")
end

# Attachment positions at flesh (skin) level.
surface_point(::Cone, body::AbstractBody, loc::EndA) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), zero(loc.r))
surface_point(::Cone, body::AbstractBody, loc::EndB) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), body.geometry.length.length_skin)
function surface_point(shape::Cone, body::AbstractBody, loc::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    t = shape.top_ratio
    rz = R * (1 - (1 - t) * loc.z / L)
    (rz * cos(loc.φ), rz * sin(loc.φ), loc.z)
end

surface_normal(::Cone, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_normal(::Cone, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
function surface_normal(shape::Cone, body::AbstractBody, loc::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    # Slant surface tangent direction has Δr = -(R - r) over Δz = L → outward
    # normal is (radial)·(L/s) + (axial)·((R-r)/s). All quantities are
    # lengths; each ratio is unitless — no ustrip needed.
    Δr = R * (1 - shape.top_ratio)
    s = sqrt(Δr^2 + L^2)
    (L/s * cos(loc.φ), L/s * sin(loc.φ), Δr/s)
end

# Centroids (FullCover).
function surface_centroid(::Cone, body::AbstractBody, ::EndA)
    R = body.geometry.length.radius_skin; (zero(R), zero(R), zero(R))
end
function surface_centroid(::Cone, body::AbstractBody, ::EndB)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), zero(R), L)
end
function surface_centroid(shape::Cone, body::AbstractBody, ::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    rz = R * (1 + shape.top_ratio) / 2  # midpoint radius
    (rz, zero(R), L/2)
end
surface_centroid_normal(::Cone, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_centroid_normal(::Cone, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
function surface_centroid_normal(shape::Cone, body::AbstractBody, ::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    Δr = R * (1 - shape.top_ratio)
    s = sqrt(Δr^2 + L^2)
    (L/s, 0.0, Δr/s)
end
