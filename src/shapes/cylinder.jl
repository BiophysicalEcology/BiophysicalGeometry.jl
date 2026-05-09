"""
    Cylinder <: AbstractShape

A cylindrical organism shape.
"""
struct Cylinder{M,D,B} <: AbstractShape
    mass::M
    density::D
    axis_ratio_b::B
end

_cylinder_radius(shape::Cylinder, volume) = cbrt(volume / (shape.axis_ratio_b * π * 2))

function geometry(shape::Cylinder, ::Naked)
    volume = body_volume(shape)
    radius_skin = _cylinder_radius(shape, volume)
    length_skin = shape.axis_ratio_b * radius_skin * 2
    total = surface_area(shape, radius_skin, length_skin)
    return Geometry(volume, (; length_skin, radius_skin), SurfaceAreas(; total))
end
function geometry(shape::Cylinder, fibrous_layer::FibrousLayer)
    volume = body_volume(shape)
    radius_skin = _cylinder_radius(shape, volume)
    radius_fur = radius_skin + fibrous_layer.thickness
    length_skin = shape.axis_ratio_b * radius_skin * 2
    length_fur = length_skin + fibrous_layer.thickness * 2
    areas = fibrous_areas(shape, fibrous_layer, (radius_skin, length_skin), (radius_fur, length_fur))
    return Geometry(volume, (; radius_skin, radius_fur, length_skin, length_fur), areas)
end
function geometry(shape::Cylinder, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_volume = volume - fat_volume(shape, fat_layer)
    radius_skin = _cylinder_radius(shape, volume)
    length_skin = shape.axis_ratio_b * radius_skin * 2
    radius_flesh = _cylinder_radius(shape, flesh_volume)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin, length_skin)
    return Geometry(volume, (; radius_skin, length_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Cylinder, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_volume = volume - fat_volume(shape, fat_layer)
    radius_skin = _cylinder_radius(shape, volume)
    radius_fur = radius_skin + fibrous_layer.thickness
    length_skin = shape.axis_ratio_b * radius_skin * 2
    length_fur = length_skin + fibrous_layer.thickness * 2
    radius_flesh = _cylinder_radius(shape, flesh_volume)
    fat = radius_skin - radius_flesh
    areas = fibrous_areas(shape, fibrous_layer, (radius_skin, length_skin), (radius_fur, length_fur))
    return Geometry(volume, (; radius_skin, radius_fur, length_skin, length_fur, fat), areas)
end

# Surface area

surface_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2

# Silhouette area

silhouette_area(shape::Cylinder, r, l, θ) = 2 * r * l * sin(θ) + π * r^2 * cos(θ)

silhouette_area(shape::Cylinder, ins::AbstractInsulationLayer, body::AbstractBody, θ) =
    silhouette_area(shape, _cylinder_outer_dims(ins, body)..., θ)
function silhouette_area(shape::Cylinder, ins::AbstractInsulationLayer, body::AbstractBody)
    (r, l) = _cylinder_outer_dims(ins, body)
    return (; normal=2 * r * l, parallel=π * r^2)
end

_cylinder_outer_dims(::Union{Naked,FatLayer}, body) =
    (body.geometry.length.radius_skin, body.geometry.length.length_skin)
_cylinder_outer_dims(::Union{FibrousLayer,CompositeInsulation}, body) =
    (body.geometry.length.radius_fur, body.geometry.length.length_fur)

# Radius accessors

_skin_radius(::Cylinder, length) = length.radius_skin
_fur_radius(::Cylinder, length) = length.radius_fur
