"""
    Sphere <: AbstractShape

A spherical organism shape.
"""
mutable struct Sphere{M,D} <: AbstractSpherical
    mass::M
    density::D
end

function _skin_level(shape::Sphere, volume)
    radius_skin = _sphere_radius(volume)
    (; dims = (; radius_skin), area = surface_area(shape, radius_skin))
end
function _fibrous_level(shape::Sphere, skin, thickness)
    radius_fibrous = skin.radius_skin + thickness
    (; dims = (; radius_fibrous), area = surface_area(shape, radius_fibrous))
end
_fat_thickness(shape::Sphere, skin, flesh_volume, fat_volume) =
    skin.radius_skin - _sphere_radius(flesh_volume)

# Surface area

function surface_area(shape::Sphere, body::AbstractBody)
    r = body.geometry.length_skin / 2
    return surface_area(shape, r)
end
function surface_area(shape::Sphere, r)
    4 * π * r ^ 2
end

# Silhouette area

silhouette(shape::Sphere, r) = π * r ^ 2
function silhouette(shape::Sphere, insulation::Union{Naked,FatLayer}, body::AbstractBody, θ)
    r = body.geometry.length.radius_skin
    return silhouette(shape, r)
end
function silhouette(shape::Sphere, insulation::Union{FibrousLayer,CompositeInsulation}, body::AbstractBody, θ)
    r = body.geometry.length.radius_fibrous
    return silhouette(shape, r)
end
function silhouette(shape::Sphere, insulation::Union{Naked,FatLayer}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    area = silhouette(shape, r)
    normal = area
    parallel = area
    return (; normal, parallel)
end
function silhouette(shape::Sphere, insulation::Union{FibrousLayer,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fibrous
    area = silhouette(shape, r)
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
