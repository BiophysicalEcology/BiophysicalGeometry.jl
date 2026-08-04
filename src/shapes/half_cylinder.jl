"""
    HalfCylinder <: AbstractShape

A half-cylinder organism shape — a cylinder cut along a plane through its
axis. Dimensions are sized so that two half-cylinders of equal mass joined
at their flat faces reconstruct the full `Cylinder(2*mass, density, b)`.

In local coordinates the axis lies along `+z`, the curved surface occupies
`y ≥ 0` half-space, and the flat face is at `y = 0`.

The flat-face area is reported at *skin* level even when fur is present,
so that two halves (mixed insulation) join cleanly with `FullCover`.
"""
mutable struct HalfCylinder{M,D,B} <: AbstractCylindrical
    mass::M
    density::D
    axis_ratio_b::B
end

function geometry(shape::HalfCylinder, ::Naked)
    volume = shape.mass / shape.density
    radius_skin = (volume / (π * shape.axis_ratio_b))^(1 / 3)
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1 / 3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin),
                    SurfaceAreas(; total))
end
function geometry(shape::HalfCylinder, fur::FibrousLayer)
    volume = shape.mass / shape.density
    radius_skin = (volume / (π * shape.axis_ratio_b))^(1 / 3)
    radius_fibrous = radius_skin + fur.thickness
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    length_fibrous = length_skin + 2 * fur.thickness
    flat_area = 2 * radius_skin * length_skin
    total = π * radius_fibrous * length_fibrous + π * radius_fibrous^2 + flat_area
    skin = π * radius_skin * length_skin + π * radius_skin^2 + flat_area
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fibrous, length_skin, length_fibrous),
                    SurfaceAreas(; total, skin, convection))
end
function geometry(shape::HalfCylinder, fat::FatLayer)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (π * shape.axis_ratio_b))^(1 / 3)
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    radius_flesh = (flesh_volume / (π * shape.axis_ratio_b))^(1 / 3)
    fat_thickness = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1 / 3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin, fat=fat_thickness),
                    SurfaceAreas(; total))
end
function geometry(shape::HalfCylinder, fur::FibrousLayer, fat::FatLayer)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (π * shape.axis_ratio_b))^(1 / 3)
    radius_fibrous = radius_skin + fur.thickness
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    length_fibrous = length_skin + 2 * fur.thickness
    radius_flesh = (flesh_volume / (π * shape.axis_ratio_b))^(1 / 3)
    fat_thickness = radius_skin - radius_flesh
    flat_area = 2 * radius_skin * length_skin
    total = π * radius_fibrous * length_fibrous + π * radius_fibrous^2 + flat_area
    skin = π * radius_skin * length_skin + π * radius_skin^2 + flat_area
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fibrous, length_skin, length_fibrous, fat=fat_thickness),
                    SurfaceAreas(; total, skin, convection))
end

# Aggregate surface area (naked).
surface_area(::HalfCylinder, r, L) = π * r * L + π * r^2 + 2 * r * L

# Silhouette area: half of the equivalent full cylinder's, so two halves
# joined back together reproduce the full cylinder's silhouettes.
function silhouette_area(::HalfCylinder, r, L, θ)
    r * L * sin(θ) + (π * r^2 / 2) * cos(θ)
end
function silhouette_area(::HalfCylinder, ::Naked, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    silhouette_area(shape(body), r, L, θ)
end
function silhouette_area(::HalfCylinder, ::FatLayer, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    silhouette_area(shape(body), r, L, θ)
end
function silhouette_area(::HalfCylinder, ::FibrousLayer, body::AbstractBody, θ)
    r = body.geometry.length.radius_fibrous
    L = body.geometry.length.length_fibrous
    silhouette_area(shape(body), r, L, θ)
end
function silhouette_area(::HalfCylinder, ::CompositeInsulation, body::AbstractBody, θ)
    r = body.geometry.length.radius_fibrous
    L = body.geometry.length.length_fibrous
    silhouette_area(shape(body), r, L, θ)
end
function silhouette_area(::HalfCylinder, ::Union{Naked,FatLayer}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (; normal = r * L, parallel = π * r^2 / 2)
end
function silhouette_area(::HalfCylinder, ::Union{FibrousLayer,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fibrous
    L = body.geometry.length.length_fibrous
    (; normal = r * L, parallel = π * r^2 / 2)
end

# Radii — same dispatch as Cylinder.
skin_radius(::HalfCylinder, ::AbstractInsulationLayer, body) = body.geometry.length.radius_skin

insulation_radius(::HalfCylinder, ::Naked, body) = body.geometry.length.radius_skin
flesh_radius(::HalfCylinder, ::Naked, body) = body.geometry.length.radius_skin

insulation_radius(::HalfCylinder, ::FibrousLayer, body) = body.geometry.length.radius_fibrous
flesh_radius(::HalfCylinder, ::FibrousLayer, body) = body.geometry.length.radius_skin

insulation_radius(::HalfCylinder, ::FatLayer, body) = body.geometry.length.radius_skin
flesh_radius(::HalfCylinder, ::FatLayer, body) = body.geometry.length.radius_skin - body.geometry.length.fat

insulation_radius(::HalfCylinder, ::CompositeInsulation, body) = body.geometry.length.radius_fibrous
flesh_radius(::HalfCylinder, ::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::HalfCylinder) = (EndA, EndB, Lateral, Flat)

# Outer (insulation-aware) dimensions. `r` matches insulation_radius(body).
outer_dims(sh::HalfCylinder, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::HalfCylinder, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (r = body.geometry.length.radius_skin, L = body.geometry.length.length_skin)
outer_dims(::HalfCylinder, ::FibrousLayer, body::AbstractBody) =
    (r = body.geometry.length.radius_fibrous, L = body.geometry.length.length_fibrous)

# FibrousLayer axial half-thickness (0 for naked); fur extends past the flesh by this on each end.
_halfcyl_fur_pad(body) =
    (outer_dims(shape(body), body).L - body.geometry.length.length_skin) / 2

# Per-surface areas. Outer (insulation-aware) for the curved/end faces;
# skin level for the flat face so two halves join cleanly under FullCover.
function surface_area(::HalfCylinder, body::AbstractBody, ::EndA)
    R = insulation_radius(body)
    π * R^2 / 2
end
surface_area(sh::HalfCylinder, body::AbstractBody, ::EndB) =
    surface_area(sh, body, EndA())
function surface_area(sh::HalfCylinder, body::AbstractBody, ::Lateral)
    d = outer_dims(sh, body)
    π * d.r * d.L
end
function surface_area(::HalfCylinder, body::AbstractBody, ::Flat)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    2 * r * L
end

# Attachment positions are anchored at flesh (skin) level so joins are
# flesh-to-flesh: local z=0 is the flesh end at -z, z=length_skin is the
# flesh end at +z. Lateral attachment radius = radius_skin. The fur axial
# overhang past z=0 / z=length_skin is part of the outer mesh and counted
# in surface_area, but bears no attachments.

function validate_range(::HalfCylinder, body::AbstractBody, loc::EndA)
    R = body.geometry.length.radius_skin
    loc.r ≥ zero(loc.r) && loc.r ≤ R || error("EndA r out of range [0, $R]: $(loc.r)")
    0 ≤ loc.φ ≤ π || error("EndA φ out of range [0, π]: $(loc.φ)")
end
validate_range(sh::HalfCylinder, body::AbstractBody, loc::EndB) =
    validate_range(sh, body, EndA(loc.r, loc.φ))
function validate_range(::HalfCylinder, body::AbstractBody, loc::Lateral)
    L = body.geometry.length.length_skin
    loc.z ≥ zero(loc.z) && loc.z ≤ L || error("Lateral z out of range [0, $L]: $(loc.z)")
    0 ≤ loc.φ ≤ π || error("Lateral φ out of range [0, π]: $(loc.φ)")
end
# For HalfCylinder, Flat uses (u=z, v=x).
function validate_range(::HalfCylinder, body::AbstractBody, loc::Flat)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    loc.u ≥ zero(loc.u) && loc.u ≤ L || error("Flat z out of range [0, $L]: $(loc.u)")
    abs(loc.v) ≤ r || error("Flat x out of range ±$r: $(loc.v)")
end

surface_point(::HalfCylinder, body::AbstractBody, loc::EndA) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), zero(loc.r))
surface_point(::HalfCylinder, body::AbstractBody, loc::EndB) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), body.geometry.length.length_skin)
function surface_point(::HalfCylinder, body::AbstractBody, loc::Lateral)
    R = body.geometry.length.radius_skin
    (R * cos(loc.φ), R * sin(loc.φ), loc.z)
end
function surface_point(::HalfCylinder, body::AbstractBody, loc::Flat)
    (loc.v, zero(loc.v), loc.u)
end

surface_normal(::HalfCylinder, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_normal(::HalfCylinder, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_normal(::HalfCylinder, ::AbstractBody, loc::Lateral) =
    (cos(loc.φ), sin(loc.φ), 0.0)
surface_normal(::HalfCylinder, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)

# Centroids — also at flesh level.
function surface_centroid(::HalfCylinder, body::AbstractBody, ::EndA)
    R = body.geometry.length.radius_skin; (zero(R), R/2, zero(R))
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::EndB)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), R/2, L)
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), R, L/2)
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Flat)
    L = body.geometry.length.length_skin
    (zero(L), zero(L), L/2)
end
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Lateral) = (0.0, 1.0, 0.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)
