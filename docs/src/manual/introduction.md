# Introduction

BiophysicalGeometry.jl computes the areas, lengths and volumes of organisms from geometry, for biophysical
calculations. A heat or water budget needs the surface area that exchanges heat with the air, the silhouette area
that intercepts sunlight, the thickness of fur and fat that heat is conducted through, and the volume that stores
heat. This package gets them from simple shapes, alone or joined into a body with many parts.

```@setup intro
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## Geometry or allometry

There are two ways to get the surface area of an animal of a given mass. An *allometric* equation is fitted to
measurements of many animals, and is the domain of
[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl). A *geometric* calculation assumes
a shape, and is the domain of this package:

```@example intro
using BiophysicalGeometry, Unitful
import BiologicalScaling

mass = 70.0u"kg"
allometric_area = BiologicalScaling.surface_area(BiologicalScaling.EutherianMammal(), mass)
geometric_area = total_area(Body(Cylinder(; mass, density = 1000.0u"kg/m^3", axis_ratio_b = 4.0), Naked()))
allometric_area, geometric_area
```

Allometry needs only a mass, but gives one number for an average animal. Geometry needs a shape, and in return
gives every area and length consistently, for any posture, coat or body part. Both packages export `surface_area`
and `skin_area`, so load one of them with `import`. The tutorial
[Geometry and allometry](../tutorials/scaling.md) compares the two over a range of sizes.

## Shapes, layers and parts

A [`Body`](@ref) is a shape and its layers of insulation. The shape type decides how the geometry is computed, see
[Shapes](shapes.md), and the layers add fat inside the skin and fur or feathers outside it, see
[Layers](layers.md). Several bodies joined at named surfaces make a [`CompositeBody`](@ref), see
[Bodies as graphs](graphs.md).

```@example intro
fur = FibrousLayer(20.0u"mm", 30.0u"μm", 3000u"cm^-2") # hide
fat = FatLayer(0.2, 901.0u"kg/m^3") # hide
layers = CompositeInsulation(fur, fat) # hide
shape_gallery("Sphere" => Body(Sphere(; mass = 2.0u"kg", density = 1000.0u"kg/m^3"), layers), # hide
    "Cylinder" => Body(Cylinder(; mass = 2.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 2.0), layers), # hide
    "Ellipsoid" => Body(Ellipsoid(; mass = 2.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 2.0, axis_ratio_c = 2.0), layers)) # hide
```

## Origins

The package grew out of the `GEOM` routines of
[NicheMapR](https://github.com/mrke/NicheMapR) (Kearney and Porter 2020, Kearney et al. 2021), which return the
dimensions of one shape from a mass, a density and a shape code. The same calculations are here, with types in
place of codes:

| NicheMapR `GEOM_ENDO` | BiophysicalGeometry.jl |
|:--|:--|
| `SHAPE` = 1, 2, 3, 4 | [`Cylinder`](@ref), [`Sphere`](@ref), [`Plate`](@ref), [`Ellipsoid`](@ref) |
| `AMASS`, `ANDENS` | mass and density of the shape |
| `SHAPE_B`, `SHAPE_C` | `axis_ratio_b`, `axis_ratio_c` |
| `ZFUR`, `DHARA`, `RHOARA` | [`FibrousLayer`](@ref) |
| `FATPCT`, `FATDEN`, `SUBQFAT` | [`FatLayer`](@ref) |
| `ORIENT`, `ZEN` | [`silhouette`](@ref) with an angle or a [`SolarOrientation`](@ref) |
| `PCOND`, used by `HomoTherm` as the fraction of area joined to other parts | [`Join`](@ref), with the place and size of the joined patch |

NicheMapR's human model, `HomoTherm` (Kearney et al. 2026), calls `GEOM_ENDO` once per body part, and joins the parts by giving each a
fraction of its area that is in contact with the others. What is new here is that parts have positions and
orientations, so the joined areas, the shading of one part by another and the views between parts are computed
from where the parts are. The tutorial
[A human: comparison with NicheMapR](../tutorials/human.md) builds the same human both ways.

## What the geometry is for

Two packages of the BiophysicalEcology ecosystem are built on this one.

[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) solves the heat budget of a body. It chooses
its heat-transfer equations from the shape family ([`AbstractCylindrical`](@ref), [`AbstractSpherical`](@ref),
[`AbstractEllipsoidal`](@ref), [`AbstractSlab`](@ref)), and takes from the body:

| From this package | Used for |
|:--|:--|
| [`total_area`](@ref), [`skin_area`](@ref), [`evaporation_area`](@ref) | convection, radiation and evaporation |
| [`flesh_radius`](@ref), [`skin_radius`](@ref), [`insulation_radius`](@ref), [`flesh_volume`](@ref) | conduction from the core through fat and fur |
| [`silhouette`](@ref), [`silhouette_factors`](@ref) | absorbed sunlight, and exchange with sky, ground and other parts |
| [`join_area`](@ref), [`internal_distance`](@ref) | conduction of heat between joined parts |

The commonest body of several parts is a shape split into a dorsal and a ventral half. The back and the belly of an
animal differ in their coat and in what they face: sun and sky above, ground below. NicheMapR handles this by
solving a body that is all back in the environment above, and one that is all belly in the environment below, and
averaging the two. Here the two halves are parts of one body, each with its own layers, areas and views, exchanging
heat across the join.

[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl) adds thermoregulation. It
works on the parts of a [`CompositeBody`](@ref), computes each part's view of the sun, sky and ground once per
posture with [`Beam`](@ref), [`Sky`](@ref), [`Ground`](@ref) and [`Horizon`](@ref).

## Units

All inputs and outputs are [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities, in any units of
the right dimension. Ratios, fractions and angles in radians are plain numbers:

```@example intro
small = Body(Sphere(; mass = 20.0u"g", density = 1.0u"g/cm^3"), Naked())
uconvert(u"cm^2", total_area(small)), uconvert(u"mm", skin_radius(small))
```

## BiophysicalGeometry.jl in the BiophysicalEcology ecosystem

BiophysicalGeometry.jl can be used on its own, and is part of the
[BiophysicalEcology](https://github.com/BiophysicalEcology) packages for mechanistic niche modelling, to be brought
together in [NicheMapper.jl](https://github.com/BiophysicalEcology/NicheMapper.jl) (in development):

| Package | Link with BiophysicalGeometry.jl |
|:--|:--|
| [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) | Uses the areas, radii, volumes and joins of a body in its heat budget |
| [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl) | Uses multipart bodies, their views of sun, sky and ground, and changes of posture |
| [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl) | Provides limb aspect ratios and body-part proportions as shape parameters, and allometric areas to compare with |
| [SolarRadiation.jl](https://github.com/BiophysicalEcology/SolarRadiation.jl) | Provides the direction and strength of the sunlight that a silhouette intercepts |
| [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) | Provides the horizon angles used by [`Horizon`](@ref) |
