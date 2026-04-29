# TODO: remove AbstractGeometryModel
abstract type AbstractGeometryModel end

"""
    AbstractShape

Abstract supertype for the shape of the organism being modelled.
"""
abstract type AbstractShape <: AbstractGeometryModel end

"""
    AbstractInsulationLayer

Abstract supertype for all insulation layers of an organism being modelled.
"""
abstract type AbstractInsulationLayer end

"""
    AbstractSolidLayer <: AbstractInsulationLayer

Abstract supertype for solid (non-porous) insulation layers such as subcutaneous fat,
chitin (arthropod cuticle), or keratin (scales, scutes). Users may define additional
solid layers by subtyping this.
"""
abstract type AbstractSolidLayer <: AbstractInsulationLayer end

"""
    AbstractPorousLayer <: AbstractInsulationLayer

Abstract supertype for porous (fibre/air-filled) insulation layers such as fur, feathers,
or hair. The porous structure is characterised by fibre geometry rather than bulk
fraction. Users may define additional porous layers by subtyping this.
"""
abstract type AbstractPorousLayer <: AbstractInsulationLayer end

"""
    Naked <: AbstractInsulationLayer

    Naked()

Insulation trait for an organism without any insulation layer.
"""
struct Naked <: AbstractInsulationLayer end

"""
    FibrousLayer <: AbstractPorousLayer

    FibrousLayer(thickness, fibre_diameter, fibre_density)

A porous insulation layer characterised by fibre geometry (fur, feathers, wool, hair, etc.).
"""
struct FibrousLayer{T,D,R} <: AbstractPorousLayer
    thickness::T
    fibre_diameter::D
    fibre_density::R
end

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

A composite of insulation layers (e.g., fibrous layer over fat layer) for an organism.
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

"""
    SolarOrientation

Abstract supertype for solar orientation traits used to select how silhouette area is computed.

Concrete subtypes: [`NormalToSun`](@ref), [`ParallelToSun`](@ref), [`Intermediate`](@ref), [`ZenithAngleVarying`](@ref).
"""
abstract type SolarOrientation <: AbstractGeometryPars end

"""
    NormalToSun <: SolarOrientation

Orientation trait: body axis perpendicular to the sun, maximising silhouette area.
"""
struct NormalToSun <: SolarOrientation end

"""
    ParallelToSun <: SolarOrientation

Orientation trait: body axis parallel to the sun, minimising silhouette area.
"""
struct ParallelToSun <: SolarOrientation end

"""
    Intermediate <: SolarOrientation

Orientation trait: silhouette area is the average of [`NormalToSun`](@ref) and [`ParallelToSun`](@ref).
"""
struct Intermediate <: SolarOrientation end

"""
    ZenithAngleVarying <: SolarOrientation

Orientation trait: silhouette area computed from the solar zenith angle via shape-specific dispatch.
Falls back to [`Intermediate`](@ref) for shapes that do not implement zenith-angle silhouette area.
"""
struct ZenithAngleVarying <: SolarOrientation end

"""
    Geometry

    Geometry(volume, length, area)

The geometry of an organism.
"""
struct Geometry{V,L,A<:SurfaceAreas} <: AbstractGeometryPars
    volume::V
    length::L
    area::A
end

"""
    AbstractBody

Abstract supertype for organism bodies.
"""
abstract type AbstractBody <: AbstractGeometryPars end

"""
    Body <: AbstractBody

    Body(shape::AbstractShape, insulation::AbstractInsulationLayer)
    Body(shape::AbstractShape, insulation::AbstractInsulationLayer, geometry::AbstractGeometryPars)

Physical dimensions of a body or body part that may or may not be insulated.
"""
struct Body{S<:AbstractShape, I<:AbstractInsulationLayer, G<:AbstractGeometryPars} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end

Body(shape::AbstractShape, insulation::AbstractInsulationLayer) =
    Body(shape, insulation, geometry(shape, insulation))

"""
    shape(body::AbstractBody) -> AbstractShape

Return the shape of `body`.
"""
shape(body::AbstractBody) = body.shape

"""
    insulation(body::AbstractBody) -> AbstractInsulationLayer

Return the insulation of `body`.
"""
insulation(body::AbstractBody) = body.insulation

"""
    geometry(body::AbstractBody) -> AbstractGeometryPars

Return the geometry of `body`.
"""
geometry(body::AbstractBody) = body.geometry

"""
    surface_area(body::AbstractBody)

Return the outer surface area of `body`.
"""
surface_area(body::AbstractBody) = surface_area(shape(body), body)

# Surface areas

"""
    total_area(body::AbstractBody)

Return the total outer surface area of `body` (including insulation if present).
"""
total_area(body::AbstractBody) = total_area(shape(body), insulation(body), body)

"""
    skin_area(body::AbstractBody)

Return the skin surface area of `body` (beneath any insulation).
"""
skin_area(body::AbstractBody) = skin_area(shape(body), insulation(body), body)

"""
    evaporation_area(body::AbstractBody)

Return the area available for evaporative water loss from `body`.
"""
evaporation_area(body::AbstractBody) = evaporation_area(shape(body), insulation(body), body)

# Fallbacks — mostly the same for all shapes
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

# Insulation area

"""
    insulation_area(fibre_diameter, fibre_density, skin)

Return the total cross-sectional area of insulation fibres (fur/feather) covering `skin` area,
given `fibre_diameter` and `fibre_density` (fibres per unit area).
"""
function insulation_area(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end

# Volume

"""
    flesh_volume(body::AbstractBody)

Return the volume of the flesh (non-fat) component of `body`.
"""
flesh_volume(body::AbstractBody) = flesh_volume(insulation(body), body)
function flesh_volume(ins::Union{FatLayer, CompositeInsulation}, body)
    fat_layer = inner_insulation(body.insulation)
    if body.geometry.length.fat > zero(body.geometry.length.fat)
        body.geometry.volume - body.shape.mass * fat_layer.fraction / fat_layer.density
    else
        body.geometry.volume
    end
end
flesh_volume(ins::FibrousLayer, body) = body.geometry.volume
flesh_volume(ins::Naked, body) = body.geometry.volume

# Radius

"""
    skin_radius(body::AbstractBody)

Return the radius at the skin surface of `body`.
"""
skin_radius(body::AbstractBody) = skin_radius(shape(body), insulation(body), body)

"""
    insulation_radius(body::AbstractBody)

Return the outer radius of the insulation layer of `body`.
"""
insulation_radius(body::AbstractBody) = insulation_radius(shape(body), insulation(body), body)

"""
    flesh_radius(body::AbstractBody)

Return the radius of the flesh (inner) cylinder of `body`.
"""
flesh_radius(body::AbstractBody) = flesh_radius(shape(body), insulation(body), body)

# Helpers for handling CompositeInsulation

"""
    outer_insulation(ins::AbstractInsulationLayer) -> AbstractInsulationLayer

Return the outermost insulation layer. For [`CompositeInsulation`](@ref), returns the
porous layer (e.g. [`FibrousLayer`](@ref)) if one is present, otherwise the last layer.
For other types, returns `ins` itself.
"""
outer_insulation(ins::AbstractInsulationLayer) = ins
function outer_insulation(ins::CompositeInsulation)
    idx = findlast(i -> i isa AbstractPorousLayer, ins.layers)
    if idx !== nothing
        ins.layers[idx]
    else
        ins.layers[end]
    end
end

"""
    inner_insulation(ins::AbstractInsulationLayer) -> AbstractInsulationLayer

Return the innermost insulation layer. For [`CompositeInsulation`](@ref), returns the
solid layer (e.g. [`FatLayer`](@ref)) if one is present, otherwise the last layer.
For other types, returns `ins` itself.
"""
inner_insulation(ins::AbstractInsulationLayer) = ins
function inner_insulation(ins::CompositeInsulation)
    idx = findfirst(i -> i isa AbstractSolidLayer, ins.layers)
    if idx !== nothing
        ins.layers[idx]
    else
        ins.layers[end]
    end
end

# Silhouette area

"""
    silhouette_area(body::AbstractBody, θ)

Return the silhouette (projected) area of `body` at solar zenith angle `θ`.
"""
silhouette_area(body::AbstractBody, θ) = silhouette_area(shape(body), insulation(body), body, θ)

"""
    silhouette_area(body::AbstractBody) -> NamedTuple

Return a named tuple `(normal=..., parallel=...)` of silhouette areas for the two
bounding orientations of `body` relative to the sun.
"""
silhouette_area(body::AbstractBody) = silhouette_area(shape(body), insulation(body), body)

"""
    silhouette_area(body::AbstractBody, orientation::SolarOrientation)
    silhouette_area(body::AbstractBody, orientation::SolarOrientation, zenith_angle)

Return the silhouette (projected) area of `body` for a given [`SolarOrientation`](@ref)
or a zenith angle paired with an orientation trait.

- [`NormalToSun`](@ref): body perpendicular to sun (maximum silhouette area)
- [`ParallelToSun`](@ref): body parallel to sun (minimum silhouette area)
- [`Intermediate`](@ref): mean of normal and parallel areas
- [`ZenithAngleVarying`](@ref): computed from `zenith_angle` via shape-specific dispatch;
  falls back to [`Intermediate`](@ref) if the shape does not implement zenith-angle silhouette area
"""
silhouette_area(body::AbstractBody, ::NormalToSun) = silhouette_area(body).normal
silhouette_area(body::AbstractBody, ::ParallelToSun) = silhouette_area(body).parallel
silhouette_area(body::AbstractBody, ::Intermediate) =
    (silhouette_area(body).normal + silhouette_area(body).parallel) * 0.5

silhouette_area(body::AbstractBody, o::SolarOrientation, ::Any) = silhouette_area(body, o)

function silhouette_area(body::AbstractBody, ::ZenithAngleVarying, zenith_angle)
    sh  = shape(body)
    ins = insulation(body)
    θ   = uconvert(u"rad", zenith_angle)
    if hasmethod(silhouette_area, (typeof(sh), typeof(ins), typeof(body), typeof(θ)))
        return silhouette_area(sh, ins, body, θ)
    end
    return silhouette_area(body, Intermediate())
end
