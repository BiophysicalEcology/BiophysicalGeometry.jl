"""
    Sphere <: AbstractShape

A spherical organism shape.
"""
mutable struct Sphere{M,D} <: AbstractShape
    mass::M
    density::D
end

function geometry(shape::Sphere, ::Naked)
    volume = shape.mass / shape.density
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    total = surface_area(shape, radius_skin)
    characteristic_dimension = volume^(1 / 3) #radius_skin * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin), SurfaceAreas(; total))
end
function geometry(shape::Sphere, fur::Fur)
    volume = shape.mass / shape.density
    radius_skin = ((3 / 4)* volume / π) ^ (1 / 3)
    radius_fur = radius_skin + fur.thickness
    total = surface_area(shape, radius_fur)
    skin = surface_area(shape, radius_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness # radius_fur * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Sphere, fat::Fat)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    radius_flesh = ((3 / 4) * flesh_volume / π) ^ (1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin)
    characteristic_dimension = volume^(1 / 3) #radius_skin * 2 #
    return Geometry(volume, characteristic_dimension, (; radius_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Sphere, fur::Fur, fat::Fat)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    radius_fur = radius_skin + fur.thickness
    radius_flesh = ((3 / 4) * flesh_volume / π) ^ (1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_fur)
    skin = surface_area(shape, radius_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness #radius_fur * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, fat), SurfaceAreas(; total, skin, convection))
end

# Surface area

function surface_area(shape::Sphere, body::AbstractBody)
    r = body.geometry.length_skin / 2
    return surface_area(shape, r)
end
function surface_area(shape::Sphere, r)
    4 * π * r ^ 2
end

# Silhouette area

silhouette_area(shape::Sphere, r) = π * r ^ 2
function silhouette_area(shape::Sphere, insulation::Union{Naked,Fat}, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    return silhouette_area(shape, r)
end
function silhouette_area(shape::Sphere, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    return silhouette_area(shape, r)
end
function silhouette_area(shape::Sphere, insulation::Union{Naked,Fat}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    area = silhouette_area(shape, r)
    normal = area
    parallel = area
    return (; normal, parallel)
end
function silhouette_area(shape::Sphere, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fur
    area = silhouette_area(shape, r)
    normal = area
    parallel = area
    return (; normal, parallel)
end

# Radius

skin_radius(shape::Sphere, insulation::AbstractInsulation, body) = body.geometry.length.radius_skin

# naked
insulation_radius(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius_skin
flesh_radius(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius_skin

# fur
insulation_radius(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_fur
flesh_radius(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_skin

# fat
insulation_radius(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius_skin
flesh_radius(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# fur and fat
insulation_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::Sphere) = (:radial,)

surface_area(::Sphere, body::AbstractBody, ::Val{:radial}) =
    4 * π * insulation_radius(body)^2

function validate_position(::Sphere, body::AbstractBody, ::Val{:radial}, pos)
    issetequal(keys(pos), (:θ, :φ)) ||
        error(":radial position needs (θ, φ); got $(keys(pos))")
end

function surface_point(::Sphere, body::AbstractBody, ::Val{:radial}, pos)
    R = skin_radius(body)
    (R * sin(pos.θ) * cos(pos.φ), R * sin(pos.θ) * sin(pos.φ), R * cos(pos.θ))
end
surface_normal(::Sphere, ::AbstractBody, ::Val{:radial}, pos) =
    (sin(pos.θ) * cos(pos.φ), sin(pos.θ) * sin(pos.φ), cos(pos.θ))

function surface_centroid(::Sphere, body::AbstractBody, ::Val{:radial})
    R = skin_radius(body); (R, zero(R), zero(R))  # arbitrary point
end
surface_centroid_normal(::Sphere, ::AbstractBody, ::Val{:radial}) = (1.0, 0.0, 0.0)
