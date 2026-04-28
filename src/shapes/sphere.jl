"""
    Sphere <: AbstractShape

A spherical organism shape.
"""
struct Sphere{M,D} <: AbstractShape
    mass::M
    density::D
end

function geometry(shape::Sphere, ::Naked)
    volume = shape.mass / shape.density
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    total = surface_area(shape, radius_skin)
    return Geometry(volume, (; radius_skin), SurfaceAreas(; total))
end
function geometry(shape::Sphere, fibrous_layer::FibrousLayer)
    volume = shape.mass / shape.density
    radius_skin = ((3 / 4)* volume / π) ^ (1 / 3)
    radius_fur = radius_skin + fibrous_layer.thickness
    total = surface_area(shape, radius_fur)
    skin = surface_area(shape, radius_skin)
    area_hair = insulation_area(fibrous_layer.fibre_diameter, fibrous_layer.fibre_density, skin)
    convection = skin - area_hair
    return Geometry(volume, (; radius_skin, radius_fur), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Sphere, fat_layer::FatLayer)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat_layer.fraction
    fat_volume = fat_mass / fat_layer.density
    flesh_volume = volume - fat_volume
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    radius_flesh = ((3 / 4) * flesh_volume / π) ^ (1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin)
    return Geometry(volume, (; radius_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Sphere, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat_layer.fraction
    fat_volume = fat_mass / fat_layer.density
    flesh_volume = volume - fat_volume
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    radius_fur = radius_skin + fibrous_layer.thickness
    radius_flesh = ((3 / 4) * flesh_volume / π) ^ (1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_fur)
    skin = surface_area(shape, radius_skin)
    area_hair = insulation_area(fibrous_layer.fibre_diameter, fibrous_layer.fibre_density, skin)
    convection = skin - area_hair
    return Geometry(volume, (; radius_skin, radius_fur, fat), SurfaceAreas(; total, skin, convection))
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
function silhouette_area(shape::Sphere, insulation::Union{Naked,FatLayer}, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    return silhouette_area(shape, r)
end
function silhouette_area(shape::Sphere, insulation::Union{FibrousLayer,CompositeInsulation}, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    return silhouette_area(shape, r)
end
function silhouette_area(shape::Sphere, insulation::Union{Naked,FatLayer}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    area = silhouette_area(shape, r)
    normal = area
    parallel = area
    return (; normal, parallel)
end
function silhouette_area(shape::Sphere, insulation::Union{FibrousLayer,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fur
    area = silhouette_area(shape, r)
    normal = area
    parallel = area
    return (; normal, parallel)
end

# Radius

skin_radius(shape::Sphere, insulation::AbstractInsulationLayer, body) = body.geometry.length.radius_skin

# naked
insulation_radius(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius_skin
flesh_radius(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius_skin

# fibrous layer
insulation_radius(shape::Sphere, insulation::FibrousLayer, body) = body.geometry.length.radius_fur
flesh_radius(shape::Sphere, insulation::FibrousLayer, body) = body.geometry.length.radius_skin

# fat layer
insulation_radius(shape::Sphere, insulation::FatLayer, body) = body.geometry.length.radius_skin
flesh_radius(shape::Sphere, insulation::FatLayer, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# fibrous + fat layer
insulation_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat
