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

function _plate_skin_dims(shape::Plate, volume)
    length_skin = cbrt(volume * shape.axis_ratio_b * shape.axis_ratio_c)
    width_skin = length_skin / shape.axis_ratio_b
    height_skin = length_skin / shape.axis_ratio_c
    return (length_skin, width_skin, height_skin)
end

function _plate_fur_dims(skin_dims, thickness)
    return ntuple(i -> skin_dims[i] + thickness * 2, 3)
end

function geometry(shape::Plate, ::Naked)
    volume = body_volume(shape)
    (length_skin, width_skin, height_skin) = _plate_skin_dims(shape, volume)
    total = surface_area(shape, length_skin, width_skin, height_skin)
    return Geometry(volume, (; length_skin, width_skin, height_skin), SurfaceAreas(; total))
end
function geometry(shape::Plate, fibrous_layer::FibrousLayer)
    volume = body_volume(shape)
    skin_dims = _plate_skin_dims(shape, volume)
    (length_skin, width_skin, height_skin) = skin_dims
    fur_dims = _plate_fur_dims(skin_dims, fibrous_layer.thickness)
    (length_fur, width_fur, height_fur) = fur_dims
    areas = fibrous_areas(shape, fibrous_layer, skin_dims, fur_dims)
    fat = 0.0u"m"
    return Geometry(volume, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), areas)
end
function geometry(shape::Plate, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_v = volume - fat_volume(shape, fat_layer)
    (length_skin, width_skin, height_skin) = _plate_skin_dims(shape, volume)
    width_flesh = cbrt(flesh_v * shape.axis_ratio_b * shape.axis_ratio_c) / shape.axis_ratio_b
    fat = (width_skin - width_flesh) / 2
    total = surface_area(shape, length_skin, width_skin, height_skin)
    return Geometry(volume, (; length_skin, width_skin, height_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Plate, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_v = volume - fat_volume(shape, fat_layer)
    skin_dims = _plate_skin_dims(shape, volume)
    (length_skin, width_skin, height_skin) = skin_dims
    width_flesh = cbrt(flesh_v * shape.axis_ratio_b * shape.axis_ratio_c) / shape.axis_ratio_b
    fat = (width_skin - width_flesh) / 2
    fur_dims = _plate_fur_dims(skin_dims, fibrous_layer.thickness)
    (length_fur, width_fur, height_fur) = fur_dims
    areas = fibrous_areas(shape, fibrous_layer, skin_dims, fur_dims)
    return Geometry(volume, (; length_skin, width_skin, height_skin, length_fur, width_fur, height_fur, fat), areas)
end

# Surface area

function surface_area(shape::Plate, body::AbstractBody)
    length = body.geometry.length.length_skin
    width = body.geometry.length.width_skin
    height = body.geometry.length.height_skin
    surface_area(shape, length, width, height)
end
surface_area(shape::Plate, length, width, height) =
    length * width * 2 + length * height * 2 + width * height * 2

# Silhouette area

function silhouette_area(shape::Plate, ins::AbstractInsulationLayer, body::AbstractBody)
    (length, width, height) = _plate_outer_dims(ins, body)
    sides = (length * width, length * height, height * width)
    return (; normal=max(sides...), parallel=min(sides...))
end

_plate_outer_dims(::Union{Naked,FatLayer}, body) = _plate_skin_outer_lengths(body.geometry.length)
_plate_outer_dims(::Union{FibrousLayer,CompositeInsulation}, body) = _plate_fur_outer_lengths(body.geometry.length)

_plate_skin_outer_lengths(length) = (length.length_skin, length.width_skin, length.height_skin)
_plate_fur_outer_lengths(length) = (length.length_fur, length.width_fur, length.height_fur)

# Characteristic dimension

shortest_outer_dim(::Plate, ::Union{Naked,FatLayer}, geom) =
    min(_plate_skin_outer_lengths(geom.length)...)
shortest_outer_dim(::Plate, ::Union{FibrousLayer,CompositeInsulation}, geom) =
    min(_plate_fur_outer_lengths(geom.length)...)

# Radius accessors — Plate uses width/2 as the equivalent radius

_skin_radius(::Plate, length) = length.width_skin / 2
_fur_radius(::Plate, length) = length.width_fur / 2
