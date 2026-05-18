"""
    Sphere <: AbstractShape

A spherical organism shape.
"""
struct Sphere{M,D} <: AbstractShape
    mass::M
    density::D
end

_sphere_radius(volume) = cbrt((3 / 4) * volume / π)

function geometry(shape::Sphere, ::Naked)
    volume = body_volume(shape)
    radius_skin = _sphere_radius(volume)
    total = surface_area(shape, radius_skin)
    return Geometry(volume, (; radius_skin), SurfaceAreas(; total))
end
function geometry(shape::Sphere, fibrous_layer::FibrousLayer)
    volume = body_volume(shape)
    radius_skin = _sphere_radius(volume)
    radius_fibrous = radius_skin + fibrous_layer.thickness
    areas = fibrous_areas(shape, fibrous_layer, (radius_skin,), (radius_fibrous,))
    return Geometry(volume, (; radius_skin, radius_fibrous), areas)
end
function geometry(shape::Sphere, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_volume = volume - fat_volume(shape, fat_layer)
    radius_skin = _sphere_radius(volume)
    radius_flesh = _sphere_radius(flesh_volume)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin)
    return Geometry(volume, (; radius_skin, fat), SurfaceAreas(; total))
end
function geometry(shape::Sphere, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_volume = volume - fat_volume(shape, fat_layer)
    radius_skin = _sphere_radius(volume)
    radius_fibrous = radius_skin + fibrous_layer.thickness
    radius_flesh = _sphere_radius(flesh_volume)
    fat = radius_skin - radius_flesh
    areas = fibrous_areas(shape, fibrous_layer, (radius_skin,), (radius_fibrous,))
    return Geometry(volume, (; radius_skin, radius_fibrous, fat), areas)
end

# Surface area

function surface_area(shape::Sphere, body::AbstractBody)
    r = body.geometry.length_skin / 2
    return surface_area(shape, r)
end
surface_area(shape::Sphere, r) = 4 * π * r ^ 2

# Silhouette area

silhouette_area(shape::Sphere, r) = π * r ^ 2

silhouette_area(shape::Sphere, ins::AbstractInsulationLayer, body::AbstractBody, θ) =
    silhouette_area(shape, _sphere_outer_radius(ins, body))
function silhouette_area(shape::Sphere, ins::AbstractInsulationLayer, body::AbstractBody)
    area = silhouette_area(shape, _sphere_outer_radius(ins, body))
    return (; normal=area, parallel=area)
end

_sphere_outer_radius(::Union{Naked,FatLayer}, body) = body.geometry.length.radius_skin
_sphere_outer_radius(::Union{FibrousLayer,CompositeInsulation}, body) = body.geometry.length.radius_fibrous

# Radius accessors

_skin_radius(::Sphere, length) = length.radius_skin
_fibrous_radius(::Sphere, length) = length.radius_fibrous
