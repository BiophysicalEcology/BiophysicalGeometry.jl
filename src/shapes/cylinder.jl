"""
    Cylinder <: AbstractShape

A cylindrical organism shape.
"""
struct Cylinder{M,D,B} <: AbstractShape
    mass::M
    density::D
    axis_ratio_b::B
end

function geometry(shape::Cylinder, ::Naked)
    volume = shape.mass / shape.density
    radius_skin = (volume / (shape.axis_ratio_b * π * 2))^(1 / 3)
    length_skin = shape.axis_ratio_b * radius_skin * 2
    total = surface_area(shape, radius_skin, length_skin)
    return Geometry(volume, (; length_skin, radius_skin), SurfaceAreas(; total))
end
function geometry(shape::Cylinder, fibrous_layer::FibrousLayer)
    volume = shape.mass / shape.density
    radius_skin = (volume / (shape.axis_ratio_b * π * 2))^(1 / 3)
    radius_fur = radius_skin + fibrous_layer.thickness
    length_skin = shape.axis_ratio_b * radius_skin * 2
    length_fur = shape.axis_ratio_b * radius_skin * 2 + fibrous_layer.thickness * 2
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fibrous_layer.fibre_diameter, fibrous_layer.fibre_density, skin)
    convection = skin - area_hair
    return Geometry(volume, (; radius_skin, radius_fur, length_skin, length_fur), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Cylinder, fat_layer::FatLayer)
    fat_mass = shape.mass * fat_layer.fraction
    fat_volume = fat_mass / fat_layer.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (shape.axis_ratio_b * π * 2))^(1 / 3)
    length_skin = shape.axis_ratio_b * radius_skin * 2
    radius_flesh = (flesh_volume / (shape.axis_ratio_b * π * 2))^(1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin, length_skin)
    return Geometry(volume, (; radius_skin, length_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Cylinder, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    fat_mass = shape.mass * fat_layer.fraction
    fat_volume = fat_mass / fat_layer.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (shape.axis_ratio_b * π * 2))^(1 / 3)
    radius_fur = radius_skin + fibrous_layer.thickness
    length_skin = shape.axis_ratio_b * radius_skin * 2
    length_fur = shape.axis_ratio_b * radius_skin * 2 + fibrous_layer.thickness * 2
    radius_flesh = (flesh_volume / (shape.axis_ratio_b * π * 2))^(1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fibrous_layer.fibre_diameter, fibrous_layer.fibre_density, skin)
    convection = skin - area_hair
    return Geometry(volume, (; radius_skin, radius_fur, length_skin, length_fur, fat), SurfaceAreas(; total, skin, convection))
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
function silhouette_area(shape::Cylinder, insulation::FatLayer, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    l = body.geometry.length.length_skin
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::FibrousLayer, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::Union{Naked,FatLayer}, body::AbstractBody, θ)
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
function silhouette_area(shape::Cylinder, insulation::Union{Naked,FatLayer}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    l = body.geometry.length.length_skin
    normal = 2 * r * l
    parallel = π *r^2
    return (; normal, parallel)
end
function silhouette_area(shape::Cylinder, insulation::Union{FibrousLayer,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    normal = 2 * r * l
    parallel = π *r^2
    return (; normal, parallel)
end

# Radius

skin_radius(shape::Cylinder, insulation::AbstractInsulationLayer, body) = body.geometry.length.radius_skin

# naked
insulation_radius(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius_skin
flesh_radius(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius_skin

# fur
insulation_radius(shape::Cylinder, insulation::FibrousLayer, body) = body.geometry.length.radius_fur
flesh_radius(shape::Cylinder, insulation::FibrousLayer, body) = body.geometry.length.radius_skin

# fat
insulation_radius(shape::Cylinder, insulation::FatLayer, body) = body.geometry.length.radius_skin
flesh_radius(shape::Cylinder, insulation::FatLayer, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# fur and fat
insulation_radius(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat
