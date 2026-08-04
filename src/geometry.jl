# TODO: remove AbstractGeometryModel
abstract type AbstractGeometryModel end

"""
    AbstractShape

Abstract supertype for the shape of the organism being modelled.
"""
abstract type AbstractShape <: AbstractGeometryModel end

# Physics-relevant family intermediates between `AbstractShape` and the
# concrete shapes. Thermal consumers (HeatExchange) dispatch on these —
# one method per family covers every concrete shape in it, and a new
# concrete shape joins its family with no new physics methods.
# `TriMesh` deliberately has no family: a thermal solve on a mesh errors
# with a missing-method dispatch rather than silently using the wrong
# correlation.

"""
    AbstractCylindrical <: AbstractShape

Family of axial, circular-cross-section shapes (`Cylinder`, `HalfCylinder`,
`Cone`) that share cylindrical thermal correlations.
"""
abstract type AbstractCylindrical <: AbstractShape end

"""
    AbstractSpherical <: AbstractShape

Family of spherical shapes (`Sphere`) sharing spherical thermal correlations.
"""
abstract type AbstractSpherical <: AbstractShape end

"""
    AbstractEllipsoidal <: AbstractShape

Family of ellipsoidal shapes (`Ellipsoid`, `HalfEllipsoid`) sharing
ellipsoidal thermal correlations.
"""
abstract type AbstractEllipsoidal <: AbstractShape end

"""
    AbstractSlab <: AbstractShape

Family of flat-plate shapes (`Plate`) sharing slab thermal correlations.
"""
abstract type AbstractSlab <: AbstractShape end

"""
    AbstractInsulationLayer

Abstract supertype for the insulation of the organism being modelled.
"""
abstract type AbstractInsulationLayer <: AbstractGeometryModel end

"""
    AbstractPorousLayer <: AbstractInsulationLayer

Family of fibrous insulations (fur, feathers, hair) sharing radiative /
conductive fibre-bed physics.
"""
abstract type AbstractPorousLayer <: AbstractInsulationLayer end

"""
    AbstractSolidLayer <: AbstractInsulationLayer

Family of solid insulating layers (fat, muscle, skin) sharing
lumped-conductive-shell physics.
"""
abstract type AbstractSolidLayer <: AbstractInsulationLayer end

"""
    Naked <: AbstractInsulationLayer

    Naked()

Insulation trait for an organism without fur.
"""
struct Naked <: AbstractInsulationLayer end

"""
    FibrousLayer <: AbstractPorousLayer

    FibrousLayer(thickness, fibre_diameter, fibre_density)

A porous insulation layer characterised by fibre geometry (fur, feathers, hair, wool).
"""
struct FibrousLayer{T,D,R} <: AbstractPorousLayer
    thickness::T
    fibre_diameter::D
    fibre_density::R
end

# TODO
# """
#     Feathers <: AbstractInsulationLayer

#     Feathers(thickness)

# Insulation trait for an organism with feathers.
# """
# struct Feathers{T,D,R} <: AbstractPorousLayer
#     thickness::T
#     fibre_diameter::D
#     fibre_density::R
# end

"""
    FatLayer <: AbstractSolidLayer

    FatLayer(fraction, density)

A solid insulation layer of subcutaneous fat.
"""
struct FatLayer{F,D} <: AbstractSolidLayer
    fraction::F
    density::D
end

"""
    CompositeInsulation <: AbstractInsulationLayer

    CompositeInsulation(layers)

A composite of insulation layers (e.g., fur, fat) for an organism.
"""
struct CompositeInsulation{T<:Tuple} <: AbstractInsulationLayer
    layers::T
end
CompositeInsulation(i::AbstractInsulationLayer) = CompositeInsulation((i,))
CompositeInsulation(is::AbstractInsulationLayer...) = CompositeInsulation((is...,))

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
struct ZenithAngleVarying <: SolarOrientation end

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

    Body(shape::AbstractShape, insulation::AbstractInsulationLayer)
    Body(shape::AbstractShape, insulation::AbstractInsulationLayer, geometry::AbstractGeometryPars)

Physical dimensions of a body or body part that may or may note be insulated.
"""
struct Body{S<:AbstractShape, I<:AbstractInsulationLayer, G<:AbstractGeometryPars} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end

Body(shape::AbstractShape, insulation::AbstractInsulationLayer) =
    Body(shape, insulation, geometry(shape, insulation))

# Surface areas

total_area(body::AbstractBody) = total_area(shape(body), insulation(body), body)
skin_area(body::AbstractBody) = skin_area(shape(body), insulation(body), body)
evaporation_area(body::AbstractBody) = evaporation_area(shape(body), insulation(body), body)

# Fallbacks - mostly these are the same for all shapes
total_area(shape::AbstractShape, insulation::AbstractInsulationLayer, body::AbstractBody) = body.geometry.area.total
skin_area(shape::AbstractShape, insulation::AbstractInsulationLayer, body::AbstractBody) = body.geometry.area.skin
evaporation_area(shape::AbstractShape, insulation::AbstractInsulationLayer, body::AbstractBody) = body.geometry.area.convection

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

# Generic 3-arg fallback: zenith angle ignored for fixed orientations
silhouette_area(body::AbstractBody, o::SolarOrientation, zenith_angle) = silhouette_area(body, o)

# ZenithAngleVarying: compute from zenith angle using shape-specific 4-arg dispatch;
# falls back to Intermediate() for shapes that don't implement silhouette_area(shape, ins, body, θ)
function silhouette_area(body::AbstractBody, ::ZenithAngleVarying, zenith_angle)
    sh = shape(body)
    ins = insulation(body)
    θ = uconvert(u"rad", zenith_angle)
    if hasmethod(silhouette_area, (typeof(sh), typeof(ins), typeof(body), typeof(θ)))
        return silhouette_area(sh, ins, body, θ)
    end
    return silhouette_area(body, Intermediate())
end

# Insulation area

function insulation_area(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end

# TODO make this 'insulation_area'
function hair_area(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end

# Volume

flesh_volume(body::AbstractBody) = flesh_volume(insulation(body), body)
function flesh_volume(ins::Union{FatLayer, CompositeInsulation}, body)
    fat = inner_insulation(body.insulation)
    if body.geometry.length.fat > zero(body.geometry.length.fat)
        body.geometry.volume - body.shape.mass * fat.fraction / fat.density
    else
        body.geometry.volume
    end
end
flesh_volume(ins::FibrousLayer, body) = body.geometry.volume
flesh_volume(ins::Naked, body) = body.geometry.volume

# Radius

skin_radius(body::AbstractBody) = skin_radius(shape(body), insulation(body), body)
insulation_radius(body::AbstractBody) = insulation_radius(shape(body), insulation(body), body)
flesh_radius(body::AbstractBody) = flesh_radius(shape(body), insulation(body), body)

# Helpers for handling CompositeInsulation

# for composite insulation cases (fat and fur/feathers)
outer_insulation(ins::AbstractInsulationLayer) = ins
function outer_insulation(ins::CompositeInsulation)
    # find fur layer if present
    fur_layer = findlast(i -> i isa FibrousLayer, ins.layers)
    if fur_layer !== nothing
        ins.layers[fur_layer]
    else
        # otherwise the last layer
        ins.layers[end]
    end
end

inner_insulation(ins::AbstractInsulationLayer) = ins
function inner_insulation(ins::CompositeInsulation)
    # find fur layer if present
    fat_layer = findfirst(i -> i isa FatLayer, ins.layers)
    if fat_layer !== nothing
        ins.layers[fat_layer]
    else
        # otherwise the last layer
        ins.layers[end]
    end
end
