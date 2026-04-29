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
mutable struct HalfCylinder{M,D,B} <: AbstractShape
    mass::M
    density::D
    b::B
end

function geometry(shape::HalfCylinder, ::Naked)
    volume = shape.mass / shape.density
    radius_skin = (volume / (π * shape.b))^(1 / 3)
    length_skin = 2 * shape.b * radius_skin
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1 / 3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin),
                    SurfaceAreas(; total))
end
function geometry(shape::HalfCylinder, fur::Fur)
    volume = shape.mass / shape.density
    radius_skin = (volume / (π * shape.b))^(1 / 3)
    radius_fur = radius_skin + fur.thickness
    length_skin = 2 * shape.b * radius_skin
    length_fur = length_skin + 2 * fur.thickness
    flat_area = 2 * radius_skin * length_skin
    total = π * radius_fur * length_fur + π * radius_fur^2 + flat_area
    skin  = π * radius_skin * length_skin + π * radius_skin^2 + flat_area
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fur, length_skin, length_fur),
                    SurfaceAreas(; total, skin, convection))
end
function geometry(shape::HalfCylinder, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (π * shape.b))^(1 / 3)
    length_skin = 2 * shape.b * radius_skin
    radius_flesh = (flesh_volume / (π * shape.b))^(1 / 3)
    fat_thickness = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1 / 3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin, fat=fat_thickness),
                    SurfaceAreas(; total))
end
function geometry(shape::HalfCylinder, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (π * shape.b))^(1 / 3)
    radius_fur = radius_skin + fur.thickness
    length_skin = 2 * shape.b * radius_skin
    length_fur = length_skin + 2 * fur.thickness
    radius_flesh = (flesh_volume / (π * shape.b))^(1 / 3)
    fat_thickness = radius_skin - radius_flesh
    flat_area = 2 * radius_skin * length_skin
    total = π * radius_fur * length_fur + π * radius_fur^2 + flat_area
    skin  = π * radius_skin * length_skin + π * radius_skin^2 + flat_area
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fur, length_skin, length_fur, fat=fat_thickness),
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
function silhouette_area(::HalfCylinder, ::Fat, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    silhouette_area(shape(body), r, L, θ)
end
function silhouette_area(::HalfCylinder, ::Fur, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    L = body.geometry.length.length_fur
    silhouette_area(shape(body), r, L, θ)
end
function silhouette_area(::HalfCylinder, ::CompositeInsulation, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    L = body.geometry.length.length_fur
    silhouette_area(shape(body), r, L, θ)
end
function silhouette_area(::HalfCylinder, ::Union{Naked,Fat}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (; normal = r * L, parallel = π * r^2 / 2)
end
function silhouette_area(::HalfCylinder, ::Union{Fur,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fur
    L = body.geometry.length.length_fur
    (; normal = r * L, parallel = π * r^2 / 2)
end

# Radii — same dispatch as Cylinder.
skin_radius(::HalfCylinder, ::AbstractInsulation, body) = body.geometry.length.radius_skin

insulation_radius(::HalfCylinder, ::Naked, body) = body.geometry.length.radius_skin
flesh_radius(::HalfCylinder, ::Naked, body) = body.geometry.length.radius_skin

insulation_radius(::HalfCylinder, ::Fur, body) = body.geometry.length.radius_fur
flesh_radius(::HalfCylinder, ::Fur, body) = body.geometry.length.radius_skin

insulation_radius(::HalfCylinder, ::Fat, body) = body.geometry.length.radius_skin
flesh_radius(::HalfCylinder, ::Fat, body) = body.geometry.length.radius_skin - body.geometry.length.fat

insulation_radius(::HalfCylinder, ::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(::HalfCylinder, ::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::HalfCylinder) = (:end_a, :end_b, :lateral, :flat)

_halfcyl_outer_length(body) =
    haskey(body.geometry.length, :length_fur) ? body.geometry.length.length_fur : body.geometry.length.length_skin

# Fur axial half-thickness (0 for naked); fur extends past the flesh by this on each end.
_halfcyl_fur_pad(body) = haskey(body.geometry.length, :length_fur) ?
    (body.geometry.length.length_fur - body.geometry.length.length_skin) / 2 :
    zero(body.geometry.length.length_skin)

# Per-surface areas. Outer (insulation-aware) for the curved/end faces;
# skin level for the flat face so two halves join cleanly under FullCover.
function surface_area(::HalfCylinder, body::AbstractBody, ::Val{:end_a})
    R = insulation_radius(body)
    π * R^2 / 2
end
surface_area(sh::HalfCylinder, body::AbstractBody, ::Val{:end_b}) =
    surface_area(sh, body, Val(:end_a))
function surface_area(::HalfCylinder, body::AbstractBody, ::Val{:lateral})
    R = insulation_radius(body)
    π * R * _halfcyl_outer_length(body)
end
function surface_area(::HalfCylinder, body::AbstractBody, ::Val{:flat})
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    2 * r * L
end

# Attachment positions are anchored at flesh (skin) level so joins are
# flesh-to-flesh: local z=0 is the flesh end at -z, z=length_skin is the
# flesh end at +z. Lateral attachment radius = radius_skin. The fur axial
# overhang past z=0 / z=length_skin is part of the outer mesh and counted
# in surface_area, but bears no attachments.

function validate_position(::HalfCylinder, body::AbstractBody, ::Val{:end_a}, pos)
    issetequal(keys(pos), (:r, :φ)) || error(":end_a needs (r, φ); got $(keys(pos))")
    R = body.geometry.length.radius_skin
    pos.r ≥ zero(pos.r) && pos.r ≤ R || error(":end_a r out of range [0, $R]: $(pos.r)")
    0 ≤ pos.φ ≤ π || error(":end_a φ out of range [0, π]: $(pos.φ)")
end
validate_position(::HalfCylinder, body::AbstractBody, ::Val{:end_b}, pos) =
    validate_position(shape(body), body, Val(:end_a), pos)
function validate_position(::HalfCylinder, body::AbstractBody, ::Val{:lateral}, pos)
    issetequal(keys(pos), (:z, :φ)) || error(":lateral needs (z, φ); got $(keys(pos))")
    L = body.geometry.length.length_skin
    pos.z ≥ zero(pos.z) && pos.z ≤ L || error(":lateral z out of range [0, $L]: $(pos.z)")
    0 ≤ pos.φ ≤ π || error(":lateral φ out of range [0, π]: $(pos.φ)")
end
function validate_position(::HalfCylinder, body::AbstractBody, ::Val{:flat}, pos)
    issetequal(keys(pos), (:z, :x)) || error(":flat needs (z, x); got $(keys(pos))")
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    pos.z ≥ zero(pos.z) && pos.z ≤ L || error(":flat z out of range [0, $L]: $(pos.z)")
    abs(pos.x) ≤ r || error(":flat x out of range ±$r: $(pos.x)")
end

surface_point(::HalfCylinder, body::AbstractBody, ::Val{:end_a}, pos) =
    (pos.r * cos(pos.φ), pos.r * sin(pos.φ), zero(pos.r))
surface_point(::HalfCylinder, body::AbstractBody, ::Val{:end_b}, pos) =
    (pos.r * cos(pos.φ), pos.r * sin(pos.φ), body.geometry.length.length_skin)
function surface_point(::HalfCylinder, body::AbstractBody, ::Val{:lateral}, pos)
    R = body.geometry.length.radius_skin
    (R * cos(pos.φ), R * sin(pos.φ), pos.z)
end
function surface_point(::HalfCylinder, body::AbstractBody, ::Val{:flat}, pos)
    (pos.x, zero(pos.x), pos.z)
end

surface_normal(::HalfCylinder, ::AbstractBody, ::Val{:end_a}, _) = (0.0, 0.0, -1.0)
surface_normal(::HalfCylinder, ::AbstractBody, ::Val{:end_b}, _) = (0.0, 0.0,  1.0)
surface_normal(::HalfCylinder, ::AbstractBody, ::Val{:lateral}, pos) =
    (cos(pos.φ), sin(pos.φ), 0.0)
surface_normal(::HalfCylinder, ::AbstractBody, ::Val{:flat}, _) = (0.0, -1.0, 0.0)

# Centroids — also at flesh level.
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Val{:end_a})
    R = body.geometry.length.radius_skin; (zero(R), R/2, zero(R))
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Val{:end_b})
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), R/2, L)
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Val{:lateral})
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), R, L/2)
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Val{:flat})
    L = body.geometry.length.length_skin
    (zero(L), zero(L), L/2)
end
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Val{:end_a}) = (0.0, 0.0, -1.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Val{:end_b}) = (0.0, 0.0,  1.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Val{:lateral}) = (0.0, 1.0, 0.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Val{:flat}) = (0.0, -1.0, 0.0)
