"""
    Plate <: AbstractShape

A flat plate-shaped organism shape.
"""
struct Plate{M,D,B,C} <: AbstractShape
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
    return Geometry(volume, (; length_skin, width_skin, height_skin), SurfaceAreas(; total))
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
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    return Geometry(volume, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Plate, fat::Fat)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c 
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    width_flesh = (flesh_volume * shape.b * shape.c)^(1 / 3) / shape.b
    fat = (width_skin - width_flesh) / 2
    total = surface_area(shape, length_skin, width_skin, height_skin)
    return Geometry(volume, (; length_skin, width_skin, height_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Plate, fur::Fur, fat::Fat)
    volume = shape.mass / shape.density
    length_skin = (volume * shape.b * shape.c)^(1 / 3)
    width_skin = length_skin / shape.b
    height_skin = length_skin / shape.c
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    width_flesh = (flesh_volume * shape.b * shape.c)^(1 / 3) / shape.b
    fat = (width_skin - width_flesh) / 2
    length_fur = length_skin + fur.thickness * 2
    width_fur = width_skin + fur.thickness * 2
    height_fur = height_skin + fur.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length_skin, width_skin, height_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
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

# characteristic dimension functions

shortest_outer_dim(::Plate, ::Union{Naked,Fat}, geom) =
    min(geom.length.length_skin, geom.length.width_skin, geom.length.height_skin)
shortest_outer_dim(::Plate, ::Union{Fur,CompositeInsulation}, geom) =
    min(geom.length.length_fur, geom.length.width_fur, geom.length.height_fur)

# radii functions

skin_radius(shape::Plate, insulation::AbstractInsulation, body) = body.geometry.length.width_skin / 2

# naked
insulation_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::Naked, body) = body.geometry.length.width_skin / 2

# fur
insulation_radius(shape::Plate, insulation::Fur, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::Fur, body) = body.geometry.length.width_skin / 2

# fat
insulation_radius(shape::Plate, insulation::Fat, body) = body.geometry.length.width_skin / 2
flesh_radius(shape::Plate, insulation::Fat, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat

# fur plus fat
insulation_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_fur / 2
flesh_radius(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_skin / 2 - body.geometry.length.fat
