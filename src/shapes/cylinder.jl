"""
    Cylinder <: AbstractShape

A cylindrical organism shape.
"""
mutable struct Cylinder{M,D,B} <: AbstractShape
    mass::M
    density::D
    b::B
end

function geometry(shape::Cylinder, ::Naked)
    volume = shape.mass / shape.density
    radius_skin = (volume / (shape.b * π * 2))^(1 / 3)
    length_skin = shape.b * radius_skin * 2
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1 / 3) # radius_skin * 2
    return Geometry(volume, characteristic_dimension, (; length_skin, radius_skin), SurfaceAreas(; total))
end
function geometry(shape::Cylinder, fur::Fur)
    volume = shape.mass / shape.density
    radius_skin = (volume / (shape.b * π * 2))^(1 / 3)
    radius_fur = radius_skin + fur.thickness
    length_skin = shape.b * radius_skin * 2
    length_fur = shape.b * radius_skin * 2 + fur.thickness * 2
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness #radius_fur * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, length_skin, length_fur), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Cylinder, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (shape.b * π * 2))^(1 / 3)
    length_skin = shape.b * radius_skin * 2
    radius_flesh = (flesh_volume / (shape.b * π * 2))^(1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1 / 3) # radius_skin * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, length_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Cylinder, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (shape.b * π * 2))^(1 / 3)
    radius_fur = radius_skin + fur.thickness
    length_skin = shape.b * radius_skin * 2
    length_fur = shape.b * radius_skin * 2 + fur.thickness * 2
    radius_flesh = (flesh_volume / (shape.b * π * 2))^(1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness # radius_fur * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, length_skin, length_fur, fat), SurfaceAreas(; total, skin, convection))
end

# Surface area

# function surface_area(shape::Cylinder, body::AbstractBody)
#     r = body.geometry.length.radius
#     l = body.geometry.length.length
#     surface_area(shape, r, l)
# end
surface_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2

# Silhouette area

function silhouette_area(shape::Cylinder, ::Naked, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    l = body.geometry.length.length_skin
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::Fat, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    l = body.geometry.length.length_skin
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::Fur, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::Union{Naked,Fat}, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    l = body.geometry.length.length_skin
    return silhouette_area(shape, r, l, θ)
end

function silhouette_area(shape::Cylinder, insulation::CompositeInsulation, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    return silhouette_area(shape, r, l, θ)
end
silhouette_area(shape::Cylinder, r, l, θ) = 2 * r * l * sin(θ) + π * r^2 * cos(θ)
function silhouette_area(shape::Cylinder, insulation::Union{Naked,Fat}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    l = body.geometry.length.length_skin
    normal = 2 * r * l
    parallel = π *r^2
    return (; normal, parallel)
end
function silhouette_area(shape::Cylinder, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    normal = 2 * r * l
    parallel = π *r^2
    return (; normal, parallel)
end

# Radius

skin_radius(shape::Cylinder, insulation::AbstractInsulation, body) = body.geometry.length.radius_skin

# naked
insulation_radius(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius_skin
flesh_radius(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius_skin

# fur
insulation_radius(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_fur
flesh_radius(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_skin

# fat
insulation_radius(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius_skin
flesh_radius(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# fur and fat
insulation_radius(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::Cylinder) = (:end_a, :end_b, :lateral)

# Outer length: includes fur if present.
_cylinder_outer_length(body::AbstractBody) =
    haskey(body.geometry.length, :length_fur) ? body.geometry.length.length_fur : body.geometry.length.length_skin

# Surface AREAS report the actual outer (insulation-aware) area. This drives
# composition's patch-fits-surface validation and FullCover lookup.
surface_area(::Cylinder, body::AbstractBody, ::Val{:end_a}) =
    π * insulation_radius(body)^2
surface_area(::Cylinder, body::AbstractBody, ::Val{:end_b}) =
    π * insulation_radius(body)^2
surface_area(::Cylinder, body::AbstractBody, ::Val{:lateral}) =
    2 * π * insulation_radius(body) * _cylinder_outer_length(body)

# Attachment POSITIONS are anchored at flesh (skin) level so joins meet
# flesh-to-flesh. Local frame: z=0 is the flesh end at -z, z=length_skin is
# the flesh end at +z, and the lateral cylindrical attachment surface has
# radius = radius_skin. The fur axial overhang past z=0 and z=length_skin
# is part of the outer mesh (and counted in surface_area) but never bears
# attachments.

function validate_position(::Cylinder, body::AbstractBody, ::Val{:end_a}, pos)
    issetequal(keys(pos), (:r, :φ)) || error(":end_a position needs (r, φ); got $(keys(pos))")
    R = skin_radius(body)
    pos.r ≥ zero(pos.r) && pos.r ≤ R ||
        error(":end_a r out of range [0, $R]: got $(pos.r)")
end
validate_position(::Cylinder, body::AbstractBody, ::Val{:end_b}, pos) =
    validate_position(shape(body), body, Val(:end_a), pos)

function validate_position(::Cylinder, body::AbstractBody, ::Val{:lateral}, pos)
    issetequal(keys(pos), (:z, :φ)) || error(":lateral position needs (z, φ); got $(keys(pos))")
    L = body.geometry.length.length_skin
    pos.z ≥ zero(pos.z) && pos.z ≤ L ||
        error(":lateral z out of range [0, $L]: got $(pos.z)")
end

surface_point(::Cylinder, body::AbstractBody, ::Val{:end_a}, pos) =
    (pos.r * cos(pos.φ), pos.r * sin(pos.φ), zero(pos.r))
surface_point(::Cylinder, body::AbstractBody, ::Val{:end_b}, pos) =
    (pos.r * cos(pos.φ), pos.r * sin(pos.φ), body.geometry.length.length_skin)
function surface_point(::Cylinder, body::AbstractBody, ::Val{:lateral}, pos)
    R = skin_radius(body)
    (R * cos(pos.φ), R * sin(pos.φ), pos.z)
end

surface_normal(::Cylinder, ::AbstractBody, ::Val{:end_a}, _) = (0.0, 0.0, -1.0)
surface_normal(::Cylinder, ::AbstractBody, ::Val{:end_b}, _) = (0.0, 0.0,  1.0)
surface_normal(::Cylinder, ::AbstractBody, ::Val{:lateral}, pos) =
    (cos(pos.φ), sin(pos.φ), 0.0)

# Centroids (used by FullCover) — also at flesh level.
function surface_centroid(::Cylinder, body::AbstractBody, ::Val{:end_a})
    R = skin_radius(body); (zero(R), zero(R), zero(R))
end
function surface_centroid(::Cylinder, body::AbstractBody, ::Val{:end_b})
    (zero(body.geometry.length.length_skin), zero(body.geometry.length.length_skin),
     body.geometry.length.length_skin)
end
function surface_centroid(::Cylinder, body::AbstractBody, ::Val{:lateral})
    R = skin_radius(body); L = body.geometry.length.length_skin
    (R, zero(R), L/2)
end
surface_centroid_normal(::Cylinder, ::AbstractBody, ::Val{:end_a}) = (0.0, 0.0, -1.0)
surface_centroid_normal(::Cylinder, ::AbstractBody, ::Val{:end_b}) = (0.0, 0.0,  1.0)
surface_centroid_normal(::Cylinder, ::AbstractBody, ::Val{:lateral}) = (1.0, 0.0, 0.0)
