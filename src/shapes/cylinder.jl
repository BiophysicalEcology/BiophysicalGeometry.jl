"""
    Cylinder <: AbstractShape

A cylindrical organism shape.
"""
mutable struct Cylinder{M,D,B} <: AbstractCylindrical
    mass::M
    density::D
    axis_ratio_b::B
end

# Radial dimension from an enclosed volume; used for both skin and flesh radii.
_cylinder_radius(shape::Cylinder, volume) = (volume / (shape.axis_ratio_b * π * 2))^(1 / 3)

function _skin_level(shape::Cylinder, volume)
    radius_skin = _cylinder_radius(shape, volume)
    length_skin = shape.axis_ratio_b * radius_skin * 2
    (; dims = (; radius_skin, length_skin), area = surface_area(shape, radius_skin, length_skin))
end
function _fibrous_level(shape::Cylinder, skin, thickness)
    radius_fibrous = skin.radius_skin + thickness
    length_fibrous = skin.length_skin + thickness * 2
    (; dims = (; radius_fibrous, length_fibrous), area = surface_area(shape, radius_fibrous, length_fibrous))
end
_fat_thickness(shape::Cylinder, skin, flesh_volume, fat_volume) =
    skin.radius_skin - _cylinder_radius(shape, flesh_volume)

# Surface area

# function surface_area(shape::Cylinder, body::AbstractBody)
#     r = body.geometry.length.radius
#     l = body.geometry.length.length
#     surface_area(shape, r, l)
# end
surface_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2

# Silhouette area. `outer_dims` selects skin- vs fibrous-level (r, L) by
# insulation, so a single insulation-dispatched wrapper per arity covers
# all four insulation kinds.
silhouette_area(shape::Cylinder, r, l, θ) = 2 * r * l * sin(θ) + π * r^2 * cos(θ)
function silhouette_area(sh::Cylinder, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette_area(sh, d.r, d.L, θ)
end
function silhouette_area(sh::Cylinder, ::AbstractInsulationLayer, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = 2 * d.r * d.L, parallel = π * d.r^2)
end

# Radius — shared by every cylindrical shape (`Cylinder`, `HalfCylinder`,
# `Cone`); all store the same `radius_skin` / `radius_fibrous` / `fat`
# fields, so the dispatch lives once on the family type.

skin_radius(::AbstractCylindrical, ::AbstractInsulationLayer, body) = body.geometry.length.radius_skin

# naked
insulation_radius(::AbstractCylindrical, ::Naked, body) = body.geometry.length.radius_skin
flesh_radius(::AbstractCylindrical, ::Naked, body) = body.geometry.length.radius_skin

# fur
insulation_radius(::AbstractCylindrical, ::FibrousLayer, body) = body.geometry.length.radius_fibrous
flesh_radius(::AbstractCylindrical, ::FibrousLayer, body) = body.geometry.length.radius_skin

# fat
insulation_radius(::AbstractCylindrical, ::FatLayer, body) = body.geometry.length.radius_skin
flesh_radius(::AbstractCylindrical, ::FatLayer, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# fur and fat
insulation_radius(::AbstractCylindrical, ::CompositeInsulation, body) = body.geometry.length.radius_fibrous
flesh_radius(::AbstractCylindrical, ::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::Cylinder) = (EndA, EndB, Lateral)

# Outer (insulation-aware) dimensions. Insulation-dispatched; no runtime
# field-lookup. `r` matches insulation_radius(body); `L` is the axial extent
# of the outer surface (skin + fur padding at each end).
outer_dims(sh::Cylinder, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Cylinder, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (r = body.geometry.length.radius_skin, L = body.geometry.length.length_skin)
outer_dims(::Cylinder, ::FibrousLayer, body::AbstractBody) =
    (r = body.geometry.length.radius_fibrous, L = body.geometry.length.length_fibrous)

# Surface AREAS report the actual outer (insulation-aware) area. This drives
# composition's patch-fits-surface validation and FullCover lookup.
surface_area(::Cylinder, body::AbstractBody, ::EndA) =
    π * insulation_radius(body)^2
surface_area(::Cylinder, body::AbstractBody, ::EndB) =
    π * insulation_radius(body)^2
function surface_area(sh::Cylinder, body::AbstractBody, ::Lateral)
    d = outer_dims(sh, body)
    2 * π * d.r * d.L
end

# Attachment POSITIONS are anchored at flesh (skin) level so joins meet
# flesh-to-flesh. Local frame: z=0 is the flesh end at -z, z=length_skin is
# the flesh end at +z, and the lateral cylindrical attachment surface has
# radius = radius_skin. The fur axial overhang past z=0 and z=length_skin
# is part of the outer mesh (and counted in surface_area) but never bears
# attachments.

function validate_range(::Cylinder, body::AbstractBody, loc::EndA)
    R = skin_radius(body)
    loc.r ≥ zero(loc.r) && loc.r ≤ R ||
        error("EndA r out of range [0, $R]: got $(loc.r)")
end
validate_range(sh::Cylinder, body::AbstractBody, loc::EndB) =
    validate_range(sh, body, EndA(loc.r, loc.φ))

function validate_range(::Cylinder, body::AbstractBody, loc::Lateral)
    L = body.geometry.length.length_skin
    loc.z ≥ zero(loc.z) && loc.z ≤ L ||
        error("Lateral z out of range [0, $L]: got $(loc.z)")
end

surface_point(::Cylinder, body::AbstractBody, loc::EndA) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), zero(loc.r))
surface_point(::Cylinder, body::AbstractBody, loc::EndB) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), body.geometry.length.length_skin)
function surface_point(::Cylinder, body::AbstractBody, loc::Lateral)
    R = skin_radius(body)
    (R * cos(loc.φ), R * sin(loc.φ), loc.z)
end

surface_normal(::Cylinder, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_normal(::Cylinder, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_normal(::Cylinder, ::AbstractBody, loc::Lateral) =
    (cos(loc.φ), sin(loc.φ), 0.0)

# Centroids (used by FullCover) — also at flesh level.
function surface_centroid(::Cylinder, body::AbstractBody, ::EndA)
    R = skin_radius(body); (zero(R), zero(R), zero(R))
end
function surface_centroid(::Cylinder, body::AbstractBody, ::EndB)
    (zero(body.geometry.length.length_skin), zero(body.geometry.length.length_skin),
     body.geometry.length.length_skin)
end
function surface_centroid(::Cylinder, body::AbstractBody, ::Lateral)
    R = skin_radius(body); L = body.geometry.length.length_skin
    (R, zero(R), L/2)
end
surface_centroid_normal(::Cylinder, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_centroid_normal(::Cylinder, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_centroid_normal(::Cylinder, ::AbstractBody, ::Lateral) = (1.0, 0.0, 0.0)
