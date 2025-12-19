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
    return Geometry(volume, characteristic_dimension, (; radius_skin), (; total))
end

function geometry(shape::Sphere, fur::Fur)
    volume = shape.mass / shape.density
    radius_skin = ((3 / 4)* volume / π) ^ (1 / 3)
    radius_fur = radius_skin + fur.thickness
    total = surface_area(shape, radius_fur)
    skin = surface_area(shape, radius_skin)
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness # radius_fur * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur), (; total, skin, convection))
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
    return Geometry(volume, characteristic_dimension, (; radius_skin, fat), (; total))
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
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness #radius_fur * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, fat), (; total, skin, convection))
end

# area functions

function surface_area(shape::Sphere, body::AbstractBody)
    r = body.geometry.length_skin / 2
    return surface_area(shape, r)
end
function surface_area(shape::Sphere, r)
    4 * π * r ^ 2
end

# silhouette area functions

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

# area and radii functions

# naked
total_area(shape::Sphere, insulation::Naked, body) = body.geometry.area.total
skin_area(shape::Sphere, insulation::Naked, body) = body.geometry.area.total
evaporation_area(shape::Sphere, insulation::Naked, body) = body.geometry.area.total

skin_radius(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius_skin
insulation_radius(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius_skin
flesh_radius(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius_skin

# fur
total_area(shape::Sphere, insulation::Fur, body) = body.geometry.area.total
skin_area(shape::Sphere, insulation::Fur, body) = body.geometry.area.skin
evaporation_area(shape::Sphere, insulation::Fur, body) = body.geometry.area.convection

skin_radius(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_skin
insulation_radius(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_fur
flesh_radius(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_skin

# fat
total_area(shape::Sphere, insulation::Fat, body) = body.geometry.area.total
skin_area(shape::Sphere, insulation::Fat, body) = body.geometry.area.total
evaporation_area(shape::Sphere, insulation::Fat, body) = body.geometry.area.total

skin_radius(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius_skin
insulation_radius(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius_skin
flesh_radius(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# fur and fat
total_area(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.area.total
skin_area(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.area.skin
evaporation_area(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.area.convection

skin_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin
insulation_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat