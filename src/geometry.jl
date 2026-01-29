# TODO: remove AbstractGeometryModel
abstract type AbstractGeometryModel end

"""
    AbstractShape

Abstract supertype for the shape of the organism being modelled.
"""
abstract type AbstractShape <: AbstractGeometryModel end

"""
    AbstractInsulation

Abstract supertype for the insulation of the organism being modelled.
"""
abstract type AbstractInsulation <: AbstractGeometryModel end

"""
    Naked <: AbstractInsulation

    Naked()

Insulation trait for an organism without fur.
"""
struct Naked <: AbstractInsulation end

"""
    Fur <: AbstractInsulation
    
    Fur(thickness)

Insulation trait for an organism with fur.
"""
struct Fur{T,D,R} <: AbstractInsulation
    thickness::T
    fibre_diameter::D
    fibre_density::R
end

# TODO
# """
#     Feathers <: AbstractInsulation
    
#     Feathers(thickness)

# Insulation trait for an organism with feathers.
# """
# struct Feathers{T,D,R} <: AbstractInsulation
#     thickness::T
#     fibre_diameter::D
#     fibre_density::R
# end

"""
    Fat <: AbstractInsulation
    
    Fat(fraction, density)

Insulation trait for an organism with fat.
"""
struct Fat{F,D} <: AbstractInsulation
    fraction::F
    density::D
end

"""
    CompositeInsulation <: AbstractInsulation

    CompositeInsulation(layers)

A composite of insulation layers (e.g., fur, fat) for an organism.
"""
struct CompositeInsulation{T<:Tuple} <: AbstractInsulation
    layers::T
end
CompositeInsulation(i::AbstractInsulation) = CompositeInsulation((i,))
CompositeInsulation(is::AbstractInsulation...) = CompositeInsulation((is...,))

geometry(shape, ins::CompositeInsulation) = geometry(shape, ins.layers...)

abstract type AbstractGeometryPars end

"""
    SurfaceAreas

    SurfaceAreas(; total, skin=total, convection=total, ventral=nothing)

Surface areas of an organism for heat exchange calculations.

# Keywords

- `total`: Total outer surface area (including insulation if present)
- `skin`: Skin surface area (under insulation)
- `convection`: Area available for convection (skin minus hair coverage)
- `ventral`: Ventral surface area (for ground contact, optional)
"""
@kwdef struct SurfaceAreas{T,S,C,V}
    total::T
    skin::S = total
    convection::C = total
    ventral::V = nothing
end

abstract type SolarOrientation <: AbstractGeometryPars end

struct NormalToSun <: SolarOrientation end
struct ParallelToSun <: SolarOrientation end
struct Intermediate <: SolarOrientation end

"""
    Geometry

    Geometry(volume, characteristic_dimension, length, area)

The geometry of an organism.
"""
struct Geometry{V,C,L,A<:SurfaceAreas} <: AbstractGeometryPars
    volume::V
    characteristic_dimension::C
    length::L
    area::A
end

"""
    AbstractBody

Abstract supertype for organism bodies.
"""
abstract type AbstractBody <: AbstractGeometryPars end

shape(body::AbstractBody) = body.shape
insulation(body::AbstractBody) = body.insulation
geometry(body::AbstractBody) = body.geometry
surface_area(body::AbstractBody) = surface_area(shape(body), body)

"""
    Body <: AbstractBody

    Body(shape::AbstractShape, insulation::AbstractInsulation)
    Body(shape::AbstractShape, insulation::AbstractInsulation, geometry::AbstractGeometryPars)

Physical dimensions of a body or body part that may or may note be insulated.
"""
struct Body{S<:AbstractShape, I<:AbstractInsulation, G<:AbstractGeometryPars} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end
Body(shape::AbstractShape, insulation::AbstractInsulation) =
    Body(shape, insulation, geometry(shape, insulation))

# Surface areas

total_area(body::AbstractBody) = total_area(shape(body), insulation(body), body)
skin_area(body::AbstractBody) = skin_area(shape(body), insulation(body), body)
evaporation_area(body::AbstractBody) = evaporation_area(shape(body), insulation(body), body)

# Fallbacks - mostly these are the same for all shapes
total_area(shape::AbstractShape, insulation::AbstractInsulation, body::AbstractBody) = body.geometry.area.total
skin_area(shape::AbstractShape, insulation::AbstractInsulation, body::AbstractBody) = body.geometry.area.total
evaporation_area(shape::AbstractShape, insulation::AbstractInsulation, body::AbstractBody) = body.geometry.area.total

# CompositeInsulation uses the outer layer
total_area(shape::AbstractShape, ins::CompositeInsulation, body::AbstractBody) =
    total_area(shape, outer_insulation(ins), body)
skin_area(shape::AbstractShape, ins::CompositeInsulation, body::AbstractBody) =
    skin_area(shape, outer_insulation(ins), body)
evaporation_area(shape::AbstractShape, ins::CompositeInsulation, body::AbstractBody) =
    evaporation_area(shape, outer_insulation(ins), body)

# Silhouette area

"""
    silhouette_area(body::AbstractBody, θ)
    silhouette_area(body::AbstractBody, orientation::SolarOrientation)

Calculates the silhouette (projected) area of an object given a 
solar zenith angle `θ` or a fixed [`SolarOrientation`](@ref) such as 
[`NormalToSun`](@ref), [`ParallelToSun`](@ref), or [`Intermediate`](@ref).
"""
silhouette_area(body::AbstractBody, θ) = silhouette_area(shape(body), insulation(body), body, θ)
silhouette_area(body::AbstractBody) = silhouette_area(shape(body), insulation(body), body)
# Orientation-specific implementations
silhouette_area(body::AbstractBody, ::NormalToSun) = silhouette_area(body).normal
silhouette_area(body::AbstractBody, ::ParallelToSun) = silhouette_area(body).parallel
silhouette_area(body::AbstractBody, ::Intermediate) =
    (silhouette_area(body).normal + silhouette_area(body).parallel) * 0.5

# Insulation area

# TODO make this 'insulation_area'
function hair_area(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end

# Volume

flesh_volume(body::AbstractBody) = flesh_volume(insulation(body), body)
function flesh_volume(ins::Union{Fat, CompositeInsulation}, body)
    fat = inner_insulation(body.insulation)
    if body.geometry.length.fat > zero(body.geometry.length.fat)
        body.geometry.volume - body.shape.mass * fat.fraction / fat.density
    else
        body.geometry.volume
    end
end
flesh_volume(ins::Fur, body) = body.geometry.volume
flesh_volume(ins::Naked, body) = body.geometry.volume

# Radius

skin_radius(body::AbstractBody) = skin_radius(shape(body), insulation(body), body)
insulation_radius(body::AbstractBody) = insulation_radius(shape(body), insulation(body), body)
flesh_radius(body::AbstractBody) = flesh_radius(shape(body), insulation(body), body)


# Helpers for handling CompositeInsulation

# for composite insulation cases (fat and fur/feathers)
outer_insulation(ins::AbstractInsulation) = ins
function outer_insulation(ins::CompositeInsulation)
    # find fur layer if present
    fur_layer = findlast(i -> i isa Fur, ins.layers)
    if fur_layer !== nothing
        ins.layers[fur_layer]
    else
        # otherwise the last layer
        ins.layers[end]
    end
end

inner_insulation(ins::AbstractInsulation) = ins
function inner_insulation(ins::CompositeInsulation)
    # find fur layer if present
    fat_layer = findfirst(i -> i isa Fat, ins.layers)
    if fat_layer !== nothing
        ins.layers[fat_layer]
    else
        # otherwise the last layer
        ins.layers[end]
    end
end
