"""
    Plate <: AbstractShape

A flat plate-shaped organism shape.
"""
struct Plate{M,D,B,C} <: AbstractShape
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
end

function geometry(shape::Plate, ::Naked)
    volume = shape.mass / shape.density
    length_skin = cbrt(volume * shape.axis_ratio_b * shape.axis_ratio_c)
    width_skin = length_skin / shape.axis_ratio_b
    height_skin = length_skin / shape.axis_ratio_c 
    total = surface_area(shape, length_skin, width_skin, height_skin)
    return Geometry(volume, (; length_skin, width_skin, height_skin), SurfaceAreas(; total))
end
function geometry(shape::Plate, fibrous_layer::FibrousLayer)
    volume = shape.mass / shape.density
    length_skin = cbrt(volume * shape.axis_ratio_b * shape.axis_ratio_c)
    width_skin = length_skin / shape.axis_ratio_b
    height_skin = length_skin / shape.axis_ratio_c
    length_fur = length_skin + fibrous_layer.thickness * 2
    width_fur = width_skin + fibrous_layer.thickness * 2
    height_fur = height_skin + fibrous_layer.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length_skin, width_skin, height_skin)
    area_hair = insulation_area(fibrous_layer.fibre_diameter, fibrous_layer.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    return Geometry(volume, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Plate, fat_layer::FatLayer)
    volume = shape.mass / shape.density
    length_skin = cbrt(volume * shape.axis_ratio_b * shape.axis_ratio_c)
    width_skin = length_skin / shape.axis_ratio_b
    height_skin = length_skin / shape.axis_ratio_c
    fat_mass = shape.mass * fat_layer.fraction
    fat_volume = fat_mass / fat_layer.density
    flesh_volume = volume - fat_volume
    width_flesh = cbrt(flesh_volume * shape.axis_ratio_b * shape.axis_ratio_c) / shape.axis_ratio_b
    fat = (width_skin - width_flesh) / 2
    total = surface_area(shape, length_skin, width_skin, height_skin)
    return Geometry(volume, (; length_skin, width_skin, height_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Plate, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = shape.mass / shape.density
    length_skin = cbrt(volume * shape.axis_ratio_b * shape.axis_ratio_c)
    width_skin = length_skin / shape.axis_ratio_b
    height_skin = length_skin / shape.axis_ratio_c
    fat_mass = shape.mass * fat_layer.fraction
    fat_volume = fat_mass / fat_layer.density
    flesh_volume = volume - fat_volume
    width_flesh = cbrt(flesh_volume * shape.axis_ratio_b * shape.axis_ratio_c) / shape.axis_ratio_b
    fat = (width_skin - width_flesh) / 2
    length_fur = length_skin + fibrous_layer.thickness * 2
    width_fur = width_skin + fibrous_layer.thickness * 2
    height_fur = height_skin + fibrous_layer.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length_skin, width_skin, height_skin)
    area_hair = insulation_area(fibrous_layer.fibre_diameter, fibrous_layer.fibre_density, skin)
    convection = skin - area_hair
    return Geometry(volume, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), SurfaceAreas(; total, skin, convection))
end

# Surface area

function surface_area(shape::Plate, body)
    length = body.geometry.length.length_skin
    width = body.geometry.length.width_skin
    height = body.geometry.length.height_skin
    surface_area(shape, length, width, height)
end
surface_area(shape::Plate, length, width, height) = length * width * 2 + length * height * 2 + width * height * 2

# Silhouette area

function silhouette_area(shape::Plate, insulation::Union{Naked,FatLayer}, body::AbstractBody)
    length = body.geometry.length.length_skin
    width = body.geometry.length.width_skin
    height = body.geometry.length.height_skin
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end
function silhouette_area(shape::Plate, insulation::Union{FibrousLayer,CompositeInsulation}, body::AbstractBody)
    length = body.geometry.length.length_fur
    width = body.geometry.length.width_fur
    height = body.geometry.length.height_fur
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end

# characteristic dimension functions

shortest_outer_dim(::Plate, ::Union{Naked,FatLayer}, geom) =
    min(geom.length.length_skin, geom.length.width_skin, geom.length.height_skin)
shortest_outer_dim(::Plate, ::Union{FibrousLayer,CompositeInsulation}, geom) =
    min(geom.length.length_fur, geom.length.width_fur, geom.length.height_fur)

# radii functions

skin_radius(shape::Plate, insulation::AbstractInsulationLayer, body) = body.geometry.length.width_skin / 2

# naked
insulation_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2

# fur
insulation_radius(shape::Plate, insulation::FibrousLayer, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::FibrousLayer, body) = body.geometry.length.width_skin / 2

# fat
insulation_radius(shape::Plate, insulation::FatLayer, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::FatLayer, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat

# fur plus fat
insulation_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat
