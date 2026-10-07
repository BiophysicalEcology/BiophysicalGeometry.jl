"""
    AbstractShape

Abstract supertype for the shape of the organism being modelled.
"""
abstract type AbstractShape end

# Physics-relevant family intermediates between `AbstractShape` and the
# concrete shapes. Thermal consumers (HeatExchange) dispatch on these —
# one method per family covers every concrete shape in it, and a new
# concrete shape joins its family with no new physics methods.

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

Abstract supertype for all insulation layers of an organism being modelled.
"""
abstract type AbstractInsulationLayer end

"""
    AbstractSolidLayer <: AbstractInsulationLayer

Abstract supertype for solid (non-porous) insulation layers such as subcutaneous fat,
chitin (arthropod cuticle), or keratin (scales, scutes), sharing lumped-conductive-shell
physics. Users may define additional solid layers by subtyping this.
"""
abstract type AbstractSolidLayer <: AbstractInsulationLayer end

"""
    AbstractPorousLayer <: AbstractInsulationLayer

Abstract supertype for porous (fibre/air-filled) insulation layers such as fur, feathers,
or hair, sharing radiative / conductive fibre-bed physics. The porous structure is
characterised by fibre geometry rather than bulk fraction. Users may define additional
porous layers by subtyping this.
"""
abstract type AbstractPorousLayer <: AbstractInsulationLayer end

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

A [`FibrousLayer`](@ref) and a [`FatLayer`](@ref) together, in either order.
"""
struct CompositeInsulation{T<:Tuple} <: AbstractInsulationLayer
    layers::T
end
CompositeInsulation(i::AbstractInsulationLayer) = CompositeInsulation((i,))
CompositeInsulation(is::AbstractInsulationLayer...) = CompositeInsulation((is...,))

# Shapes define geometry for (fibrous, fat) in that order; accept the layers of
# a composite in either order by putting the porous (outer) layers first.
geometry(shape, ins::CompositeInsulation) =
    geometry(shape, filter(l -> l isa AbstractPorousLayer, ins.layers)...,
             filter(l -> !(l isa AbstractPorousLayer), ins.layers)...)

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

The computed geometry of a [`Body`](@ref).

- `volume`: mass over density
- `length`: `NamedTuple` of dimensions, with names that depend on the shape and layers
- `area`: [`SurfaceAreas`](@ref)
"""
struct Geometry{V,L,A<:SurfaceAreas} <: AbstractGeometryPars
    volume::V
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
# `Half` is not in `StandardShape`; it keeps its own `geometry` methods.

const StandardShape = Union{AbstractCylindrical, AbstractSpherical, AbstractEllipsoidal, AbstractSlab}

_flesh_volume(shape::AbstractShape, fat_layer::FatLayer) =
    body_volume(shape) - fat_volume(shape, fat_layer)
_convective_area(fibrous_layer::FibrousLayer, skin_area) =
    skin_area - insulation_area(fibrous_layer.fibre_diameter, fibrous_layer.fibre_density, skin_area)

# Equivalent-sphere radius enclosing `volume`; the ellipsoid reuses it on the
# volume scaled by its axis ratio, so the cube-root formula lives in one place.
_sphere_radius(volume) = cbrt((3 / 4) * volume / π)

function geometry(shape::StandardShape, ::Naked)
    volume = body_volume(shape)
    skin = _skin_level(shape, volume)
    Geometry(volume, skin.dims, SurfaceAreas(; total = skin.area))
end
function geometry(shape::StandardShape, fibrous_layer::FibrousLayer)
    volume = body_volume(shape)
    skin = _skin_level(shape, volume)
    fibrous = _fibrous_level(shape, skin.dims, fibrous_layer.thickness)
    Geometry(volume, merge(skin.dims, fibrous.dims),
             SurfaceAreas(; total = fibrous.area, skin = skin.area,
                          convection = _convective_area(fibrous_layer, skin.area)))
end
function geometry(shape::StandardShape, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_volume = _flesh_volume(shape, fat_layer)
    skin = _skin_level(shape, volume)
    fat = _fat_thickness(shape, skin.dims, flesh_volume, volume - flesh_volume)
    Geometry(volume, merge(skin.dims, (; fat)), SurfaceAreas(; total = skin.area))
end
function geometry(shape::StandardShape, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    volume = body_volume(shape)
    flesh_volume = _flesh_volume(shape, fat_layer)
    skin = _skin_level(shape, volume)
    fibrous = _fibrous_level(shape, skin.dims, fibrous_layer.thickness)
    fat = _fat_thickness(shape, skin.dims, flesh_volume, volume - flesh_volume)
    Geometry(volume, merge(skin.dims, fibrous.dims, (; fat)),
             SurfaceAreas(; total = fibrous.area, skin = skin.area,
                          convection = _convective_area(fibrous_layer, skin.area)))
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

Return the outer surface area of `body` — the same as [`total_area`](@ref).
"""
surface_area(body::AbstractBody) = total_area(body)

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

# Silhouette area

"""
    silhouette(body::AbstractBody, θ)
    silhouette(body::AbstractBody) -> NamedTuple
    silhouette(body::AbstractBody, orientation::SolarOrientation)
    silhouette(body::AbstractBody, orientation::SolarOrientation, zenith_angle)

Return the silhouette (projected) area of `body` at solar zenith angle `θ`, or for a
fixed [`SolarOrientation`](@ref). With no second argument, returns a named tuple
`(normal=..., parallel=...)` of the two bounding orientations.

- [`NormalToSun`](@ref): body perpendicular to sun (maximum silhouette area)
- [`ParallelToSun`](@ref): body parallel to sun (minimum silhouette area)
- [`Intermediate`](@ref): mean of normal and parallel areas
- [`ZenithAngleVarying`](@ref): computed from `zenith_angle` via shape-specific dispatch;
  falls back to [`Intermediate`](@ref) if the shape does not implement zenith-angle silhouette area
"""
silhouette(body::AbstractBody, θ) = silhouette(shape(body), insulation(body), body, θ)
silhouette(body::AbstractBody) = silhouette(shape(body), insulation(body), body)

silhouette(body::AbstractBody, ::NormalToSun) = silhouette(body).normal
silhouette(body::AbstractBody, ::ParallelToSun) = silhouette(body).parallel
silhouette(body::AbstractBody, ::Intermediate) =
    (silhouette(body).normal + silhouette(body).parallel) * 0.5

# Generic 3-arg fallback: zenith angle ignored for fixed orientations
silhouette(body::AbstractBody, o::SolarOrientation, ::Any) = silhouette(body, o)

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

"""
    insulation_area(fibre_diameter, fibre_density, skin)

Return the total cross-sectional area of insulation fibres (fur, feathers, etc.) covering `skin` area,
given `fibre_diameter` and `fibre_density` (fibres per unit area).
"""
function insulation_area(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end

# Volume

"""
    body_volume(shape::AbstractShape)

Return the body volume `mass / density` for `shape`, in m³ when unitful — so
mixing units (grams with kg/m³) doesn't leave odd units in every length.
"""
body_volume(shape::AbstractShape) = _in_m³(shape.mass / shape.density)
_in_m³(v::Unitful.Volume) = uconvert(u"m^3", v)
_in_m³(v) = v

"""
    fat_volume(shape::AbstractShape, fat_layer::FatLayer)

Return the fat volume implied by a [`FatLayer`](@ref): `mass * fraction / density`.
"""
fat_volume(shape::AbstractShape, fat_layer::FatLayer) =
    shape.mass * fat_layer.fraction / fat_layer.density

"""
    flesh_volume(body::AbstractBody)

Return the volume of the flesh (non-fat) component of `body`.
"""
flesh_volume(body::AbstractBody) = flesh_volume(insulation(body), body)
function flesh_volume(ins::Union{FatLayer, CompositeInsulation}, body)
    fat_layer = inner_insulation(body.insulation)
    if body.geometry.length.fat > zero(body.geometry.length.fat)
        body.geometry.volume - mass(body.shape) * fat_layer.fraction / fat_layer.density
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

# Generic radius dispatch. Each concrete shape only needs to define
# `_skin_radius(shape, length)` and (where insulation is supported)
# `_fibrous_radius(shape, length)` accessors over its `body.geometry.length` NamedTuple.

skin_radius(s::AbstractShape, ::AbstractInsulationLayer, b::AbstractBody) =
    _skin_radius(s, b.geometry.length)

insulation_radius(s::AbstractShape, ::Union{Naked,FatLayer}, b::AbstractBody) =
    _skin_radius(s, b.geometry.length)
insulation_radius(s::AbstractShape, ::Union{FibrousLayer,CompositeInsulation}, b::AbstractBody) =
    _fibrous_radius(s, b.geometry.length)

flesh_radius(s::AbstractShape, ::Union{Naked,FibrousLayer}, b::AbstractBody) =
    _skin_radius(s, b.geometry.length)
flesh_radius(s::AbstractShape, ::Union{FatLayer,CompositeInsulation}, b::AbstractBody) =
    _skin_radius(s, b.geometry.length) - b.geometry.length.fat

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
