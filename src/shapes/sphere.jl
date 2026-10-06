"""
    Sphere(; mass, density, volume, radius) <: AbstractShape

A spherical organism shape centred on the origin. Give any two of `mass`,
`density`, `volume` / `radius` and the rest is solved for. `radius` is at skin level.
"""
struct Sphere{M,D} <: AbstractSpherical
    mass::M
    density::D
    Sphere(::Resolved, mass::M, density::D) where {M,D} = new{M,D}(mass, density)
end

# volume = (4π/3)·radius³
const SPHERE_SPEC = ShapeSpec((:radius,), (3,), log(4π / 3), ())

function Sphere(; kw...)
    s = _resolve_shape("Sphere", SPHERE_SPEC, NamedTuple(kw))
    Sphere(RESOLVED, s.mass, s.density)
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
outer_dims(::Sphere, body::AbstractBody) = (radius = insulation_radius(body),)

surface_area(::Sphere, body::AbstractBody, ::Radial) =
    4 * π * insulation_radius(body)^2

validate_range(::Sphere, ::AbstractBody, ::Radial) = nothing

function surface_point(::Sphere, body::AbstractBody, loc::Radial)
    R = skin_radius(body)
    (R * sin(loc.polar) * cos(loc.azimuth), R * sin(loc.polar) * sin(loc.azimuth), R * cos(loc.polar))
end
surface_normal(::Sphere, ::AbstractBody, loc::Radial) =
    (sin(loc.polar) * cos(loc.azimuth), sin(loc.polar) * sin(loc.azimuth), cos(loc.polar))

function surface_centroid(::Sphere, body::AbstractBody, ::Radial)
    R = skin_radius(body); (R, zero(R), zero(R))  # arbitrary point
end
surface_centroid_normal(::Sphere, ::AbstractBody, ::Radial) = (1.0, 0.0, 0.0)
