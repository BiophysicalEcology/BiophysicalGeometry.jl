"""
    Plate <: AbstractShape

A flat plate-shaped organism shape.
"""
mutable struct Plate{M,D,B,C} <: AbstractShape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Plate, ::Naked)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c 
    total = surface_area(shape, length_skin, width_skin, height_skin)
    characteristic_dimension = volume^(1 / 3) # width_skin * 2
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin), (; total))
end

function geometry(shape::Plate, fur::Fur)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c 
    length_fur = length_skin + fur.thickness * 2
    width_fur = width_skin + fur.thickness * 2
    height_fur = height_skin + fur.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length_skin, width_skin, height_skin)
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    characteristic_dimension = volume^(1 / 3) + fur.thickness # width_fur * 2
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), (; total, skin, convection))
end

function geometry(shape::Plate, fat::Fat)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c 
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    length_flesh = (flesh_volume * shape.b * shape.c)^(1 / 3)
    fat = (length_skin - length_flesh) / 2
    total = surface_area(shape, length_skin, width_skin, height_skin)
    characteristic_dimension = volume^(1 / 3) # width_skin * 2 
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin, fat), (; total))
end

function geometry(shape::Plate, fur::Fur, fat::Fat)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    length_flesh = (flesh_volume * shape.b * shape.c)^(1 / 3)
    fat = (length_skin - length_flesh) / 2
    length_fur = length_skin + fur.thickness * 2
    width_fur = width_skin + fur.thickness * 2
    height_fur = height_skin + fur.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length_skin, width_skin, height_skin)
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness #width_fur * 2 
    return Geometry(volume, characteristic_dimension, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), (; total, skin, convection))
end

# surface area functions

function surface_area(shape::Plate, body)
    length = body.geometry.length.length_skin
    width = body.geometry.length.width_skin
    height = body.geometry.length.height_skin
    surface_area(shape, length, width, height)
end
surface_area(shape::Plate, length, width, height) = length * width * 2 + length * height * 2 + width * height * 2

# silhouette area functions

function silhouette_area(shape::Plate, insulation::Union{Naked,Fat}, body::AbstractBody)
    length = body.geometry.length.length_skin
    width = body.geometry.length.width_skin
    height = body.geometry.length.height_skin
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end

function silhouette_area(shape::Plate, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    length = body.geometry.length.length_fur
    width = body.geometry.length.width_fur
    height = body.geometry.length.height_fur
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end

# area and radii functions

# naked
total_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total
skin_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total
evaporation_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total

skin_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2
insulation_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2

# fur
total_area(shape::Plate, insulation::Fur, body) = body.geometry.area.total
skin_area(shape::Plate, insulation::Fur, body) = body.geometry.area.skin
evaporation_area(shape::Plate, insulation::Fur, body) = body.geometry.area.convection

skin_radius(shape::Plate, insulation::Fur, body) = body.geometry.length.width_skin / 2
insulation_radius(shape::Plate, insulation::Fur, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::Fur, body) = body.geometry.length.width_skin / 2

# fat
total_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total
skin_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total
evaporation_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total

skin_radius(shape::Plate, insulation::Fat, body) = body.geometry.length.width_skin / 2
insulation_radius(shape::Plate, insulation::Fat, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::Fat, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat

# fur plus fat
total_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.total
skin_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.skin
evaporation_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.convection

skin_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_skin / 2
insulation_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat