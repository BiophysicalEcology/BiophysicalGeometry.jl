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
# Shapes with no family (the animal shapes) get a missing-method dispatch
# on a thermal solve rather than silently using the wrong correlation.

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
    Half(parent::AbstractShape) <: AbstractShape

A shape cut through the centre and mirror-joined at a flat face: the dorsal or
ventral half of `parent`. `parent` is the full shape of *double* mass, so a
`Half` reuses all of the parent's radius/flesh/fat/fur math (single source of
truth); only surface area, the flat join face, and the mesh are halved/added.
Family membership comes from the type parameter (`Half{<:AbstractCylindrical}`);
where a family supertype won't match a wrapper (HeatExchange), the `Half`
forwards to its parent. `HalfCylinder`, `HalfEllipsoid` and `HalfSphere` are
constructors returning the matching `Half`.
"""
struct Half{S<:AbstractShape} <: AbstractShape
    parent::S
end

"""
    mass(shape) -> mass

Physical mass of a shape. A `Half` wraps a *double*-mass parent, so its own mass
is half the parent's — never read `shape.mass` directly, as wrappers have no such
field.
"""
mass(s::AbstractShape) = s.mass
mass(h::Half) = mass(h.parent) / 2

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

No insulation.
"""
struct Naked <: AbstractInsulationLayer end

"""
    FibrousLayer <: AbstractPorousLayer

    FibrousLayer(thickness, fibre_diameter, fibre_density)

A layer of fibres outside the skin (fur, feathers, hair, clothing).

- `thickness`: depth of the layer (length)
- `fibre_diameter`: diameter of a fibre (length)
- `fibre_density`: fibres per area of skin (1/area)
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

A layer of subcutaneous fat inside the skin.

- `fraction`: fat mass as a fraction of body mass
- `density`: density of fat
"""
struct FatLayer{F,D} <: AbstractSolidLayer
    fraction::F
    density::D
end

"""
    CompositeInsulation <: AbstractInsulationLayer

    CompositeInsulation(fibrous, fat)

A [`FibrousLayer`](@ref) and a [`FatLayer`](@ref) together, in that order.
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

Supertype for the orientation of the long axis of a body to the sun, used by [`silhouette`](@ref).
"""
abstract type SolarOrientation <: AbstractGeometryPars end

"""
    NormalToSun()

Long axis at right angles to the sun: the largest silhouette.
"""
struct NormalToSun <: SolarOrientation end

"""
    ParallelToSun()

Long axis pointing at the sun: the smallest silhouette.
"""
struct ParallelToSun <: SolarOrientation end

"""
    Intermediate()

The mean of [`NormalToSun`](@ref) and [`ParallelToSun`](@ref).
"""
struct Intermediate <: SolarOrientation end

"""
    ZenithAngleVarying()

Silhouette computed from the zenith angle of the sun, where the shape allows.
"""
struct ZenithAngleVarying <: SolarOrientation end

"""
    Geometry

    Geometry(volume, characteristic_dimension, length, area)

The computed geometry of a [`Body`](@ref).

- `volume`: mass over density
- `characteristic_dimension`: cube root of the volume, plus the depth of any fibrous layer
- `length`: `NamedTuple` of dimensions, with names that depend on the shape and layers
- `area`: [`SurfaceAreas`](@ref)
"""
struct Geometry{V,C,L,A<:SurfaceAreas} <: AbstractGeometryPars
    volume::V
    characteristic_dimension::C
    length::L
    area::A
end

# ── Generic geometry construction ─────────────────────────────────────────
#
# The standard parametric shapes (cylinder, cone, sphere, ellipsoid, plate)
# all build their `Geometry` the same way: enclosed volume from mass/density,
# skin dimensions from that volume, then optional fat (a shell inside the skin,
# derived from the fat mass fraction) and fur (a shell outside the skin, from a
# thickness). Only three things vary per shape, so each provides those as
# primitives and the four insulation combinations are assembled here once:
#
#   _skin_level(shape, volume)             -> (; dims, area)   skin NamedTuple + skin area
#   _fibrous_level(shape, skin_dims, t)    -> (; dims, area)   fur NamedTuple + outer area
#   _fat_thickness(shape, skin_dims, flesh_volume, fat_volume) -> Length
#
# `Half` and the animal shapes are not in `_StandardShape`; they keep their
# own `geometry` methods.

const _StandardShape = Union{AbstractCylindrical, AbstractSpherical, AbstractEllipsoidal, AbstractSlab}

_body_volume(shape::AbstractShape) = shape.mass / shape.density
_flesh_volume(shape::AbstractShape, fat::FatLayer) =
    _body_volume(shape) - shape.mass * fat.fraction / fat.density
_characteristic_length(volume) = volume^(1 / 3)
_convective_area(fur::FibrousLayer, skin_area) =
    skin_area - insulation_area(fur.fibre_diameter, fur.fibre_density, skin_area)

# Equivalent-sphere radius enclosing `volume`; the ellipsoid reuses it on the
# volume scaled by its axis ratio, so the cube-root formula lives in one place.
_sphere_radius(volume) = ((3 / 4) * volume / π)^(1 / 3)

function geometry(shape::_StandardShape, ::Naked)
    volume = _body_volume(shape)
    skin = _skin_level(shape, volume)
    Geometry(volume, _characteristic_length(volume), skin.dims, SurfaceAreas(; total = skin.area))
end
function geometry(shape::_StandardShape, fur::FibrousLayer)
    volume = _body_volume(shape)
    skin = _skin_level(shape, volume)
    fibrous = _fibrous_level(shape, skin.dims, fur.thickness)
    Geometry(volume, _characteristic_length(volume) + fur.thickness,
             merge(skin.dims, fibrous.dims),
             SurfaceAreas(; total = fibrous.area, skin = skin.area,
                          convection = _convective_area(fur, skin.area)))
end
function geometry(shape::_StandardShape, fat::FatLayer)
    volume = _body_volume(shape)
    flesh_volume = _flesh_volume(shape, fat)
    skin = _skin_level(shape, volume)
    fat = _fat_thickness(shape, skin.dims, flesh_volume, volume - flesh_volume)
    Geometry(volume, _characteristic_length(volume),
             merge(skin.dims, (; fat)), SurfaceAreas(; total = skin.area))
end
function geometry(shape::_StandardShape, fur::FibrousLayer, fat::FatLayer)
    volume = _body_volume(shape)
    flesh_volume = _flesh_volume(shape, fat)
    skin = _skin_level(shape, volume)
    fibrous = _fibrous_level(shape, skin.dims, fur.thickness)
    fat = _fat_thickness(shape, skin.dims, flesh_volume, volume - flesh_volume)
    Geometry(volume, _characteristic_length(volume) + fur.thickness,
             merge(skin.dims, fibrous.dims, (; fat)),
             SurfaceAreas(; total = fibrous.area, skin = skin.area,
                          convection = _convective_area(fur, skin.area)))
end

"""
    AbstractBody

Abstract supertype for organism bodies.
"""
abstract type AbstractBody <: AbstractGeometryPars end

"""
    shape(body)

The shape of a body; of the root part for a [`CompositeBody`](@ref).
"""
shape(body::AbstractBody) = body.shape

"""
    insulation(body)

The insulation of a body; of the root part for a [`CompositeBody`](@ref).
"""
insulation(body::AbstractBody) = body.insulation

"""
    geometry(body)
    geometry(shape, insulation)

The [`Geometry`](@ref) of a body, or of a shape with the given insulation.
"""
geometry(body::AbstractBody) = body.geometry
surface_area(body::AbstractBody) = surface_area(shape(body), body)

"""
    Body <: AbstractBody

    Body(shape::AbstractShape, insulation::AbstractInsulationLayer)
    Body(shape::AbstractShape, insulation::AbstractInsulationLayer, geometry::AbstractGeometryPars)

A shape with its insulation and the geometry computed from them.
"""
struct Body{S<:AbstractShape, I<:AbstractInsulationLayer, G<:AbstractGeometryPars} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end

Body(shape::AbstractShape, insulation::AbstractInsulationLayer) =
    Body(shape, insulation, geometry(shape, insulation))

# Surface areas

"""
    total_area(body)

Outer surface area, over any fibrous layer. For a [`CompositeBody`](@ref), summed over parts less joined patches.
"""
total_area(body::AbstractBody) = total_area(shape(body), insulation(body), body)

"""
    skin_area(body)

Surface area of the skin. For a [`CompositeBody`](@ref), summed over parts less joined patches.
"""
skin_area(body::AbstractBody) = skin_area(shape(body), insulation(body), body)

"""
    evaporation_area(body)

Skin area not covered by the bases of fibres. For a [`CompositeBody`](@ref), summed over parts less joined patches.
"""
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
    silhouette(body::AbstractBody, θ)
    silhouette(body::AbstractBody, orientation::SolarOrientation)

Calculates the silhouette (projected) area of an object given a
solar zenith angle `θ` or a fixed [`SolarOrientation`](@ref) such as
[`NormalToSun`](@ref), [`ParallelToSun`](@ref), or [`Intermediate`](@ref).
"""
silhouette(body::AbstractBody, θ) = silhouette(shape(body), insulation(body), body, θ)
silhouette(body::AbstractBody) = silhouette(shape(body), insulation(body), body)

# Orientation-specific implementations
silhouette(body::AbstractBody, ::NormalToSun) = silhouette(body).normal
silhouette(body::AbstractBody, ::ParallelToSun) = silhouette(body).parallel
silhouette(body::AbstractBody, ::Intermediate) =
    (silhouette(body).normal + silhouette(body).parallel) * 0.5

# Generic 3-arg fallback: zenith angle ignored for fixed orientations
silhouette(body::AbstractBody, o::SolarOrientation, zenith_angle) = silhouette(body, o)

# ZenithAngleVarying: compute from zenith angle using shape-specific 4-arg dispatch;
# falls back to Intermediate() for shapes that don't implement silhouette(shape, ins, body, θ)
function silhouette(body::AbstractBody, ::ZenithAngleVarying, zenith_angle)
    sh = shape(body)
    ins = insulation(body)
    θ = uconvert(u"rad", zenith_angle)
    if hasmethod(silhouette, (typeof(sh), typeof(ins), typeof(body), typeof(θ)))
        return silhouette(sh, ins, body, θ)
    end
    return silhouette(body, Intermediate())
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

"""
    flesh_volume(body)

Volume inside the fat layer. For a [`CompositeBody`](@ref), summed over parts.
"""
flesh_volume(body::AbstractBody) = flesh_volume(insulation(body), body)
function flesh_volume(ins::Union{FatLayer, CompositeInsulation}, body)
    fat = inner_insulation(body.insulation)
    if body.geometry.length.fat > zero(body.geometry.length.fat)
        body.geometry.volume - mass(body.shape) * fat.fraction / fat.density
    else
        body.geometry.volume
    end
end
flesh_volume(ins::FibrousLayer, body) = body.geometry.volume
flesh_volume(ins::Naked, body) = body.geometry.volume

# Radius

"""
    skin_radius(body)

Radius at the skin: cylinder or sphere radius, short semi-axis of an ellipsoid, half-width of a plate.
"""
skin_radius(body::AbstractBody) = skin_radius(shape(body), insulation(body), body)

"""
    insulation_radius(body)

Radius at the outside of the fibrous layer; the [`skin_radius`](@ref) if there is none.
"""
insulation_radius(body::AbstractBody) = insulation_radius(shape(body), insulation(body), body)

"""
    flesh_radius(body)

Radius at the inside of the fat layer; the [`skin_radius`](@ref) if there is none.
"""
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
