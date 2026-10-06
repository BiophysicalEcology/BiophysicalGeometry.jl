"""
    Sphere <: AbstractShape

A spherical organism shape.
"""
struct Sphere{M,D} <: AbstractSpherical
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

surface_area(shape::Sphere, r) = 4 * π * r ^ 2

# Silhouette area

silhouette(shape::Sphere, r) = π * r ^ 2

silhouette(shape::Sphere, ins::AbstractInsulationLayer, body::AbstractBody, θ) =
    silhouette(shape, _sphere_outer_radius(ins, body))
function silhouette(shape::Sphere, ins::AbstractInsulationLayer, body::AbstractBody)
    area = silhouette(shape, _sphere_outer_radius(ins, body))
    return (; normal=area, parallel=area)
end

_sphere_outer_radius(::Union{Naked,FatLayer}, body) = body.geometry.length.radius_skin
_sphere_outer_radius(::Union{FibrousLayer,CompositeInsulation}, body) = body.geometry.length.radius_fibrous

# Radius accessors

_skin_radius(::Sphere, length) = length.radius_skin
_fibrous_radius(::Sphere, length) = length.radius_fibrous

# Composition

attachment_surfaces(::Sphere) = (Radial,)

# Outer (insulation-aware) radius, matching insulation_radius(body).
outer_dims(::Sphere, body::AbstractBody) = (r = insulation_radius(body),)

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
