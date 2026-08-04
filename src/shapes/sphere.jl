"""
    Sphere <: AbstractShape

A spherical organism shape.
"""
mutable struct Sphere{M,D} <: AbstractSpherical
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
function geometry(shape::Sphere, fur::FibrousLayer)
    volume = shape.mass / shape.density
    radius_skin = ((3 / 4)* volume / π) ^ (1 / 3)
    radius_fibrous = radius_skin + fur.thickness
    total = surface_area(shape, radius_fibrous)
    skin = surface_area(shape, radius_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness # radius_fibrous * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fibrous), SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Sphere, fat::FatLayer)
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
function geometry(shape::Sphere, fur::FibrousLayer, fat::FatLayer)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    radius_fibrous = radius_skin + fur.thickness
    radius_flesh = ((3 / 4) * flesh_volume / π) ^ (1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_fibrous)
    skin = surface_area(shape, radius_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1 / 3) + fur.thickness #radius_fibrous * 2
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fibrous, fat), SurfaceAreas(; total, skin, convection))
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
    r = body.geometry.length.radius_fibrous
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
    r = body.geometry.length.radius_fibrous
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

# fur
insulation_radius(shape::Sphere, insulation::FibrousLayer, body) = body.geometry.length.radius_fibrous
flesh_radius(shape::Sphere, insulation::FibrousLayer, body) = body.geometry.length.radius_skin

# fat
insulation_radius(shape::Sphere, insulation::FatLayer, body) = body.geometry.length.radius_skin
flesh_radius(shape::Sphere, insulation::FatLayer, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# fur and fat
insulation_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_fibrous
flesh_radius(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::Sphere) = (Radial,)

surface_area(::Sphere, body::AbstractBody, ::Radial) =
    4 * π * insulation_radius(body)^2

validate_range(::Sphere, ::AbstractBody, ::Radial) = nothing

function surface_point(::Sphere, body::AbstractBody, loc::Radial)
    R = skin_radius(body)
    (R * sin(loc.θ) * cos(loc.φ), R * sin(loc.θ) * sin(loc.φ), R * cos(loc.θ))
end
surface_normal(::Sphere, ::AbstractBody, loc::Radial) =
    (sin(loc.θ) * cos(loc.φ), sin(loc.θ) * sin(loc.φ), cos(loc.θ))

function surface_centroid(::Sphere, body::AbstractBody, ::Radial)
    R = skin_radius(body); (R, zero(R), zero(R))  # arbitrary point
end
surface_centroid_normal(::Sphere, ::AbstractBody, ::Radial) = (1.0, 0.0, 0.0)
