# Shapes

A shape is a type, such as [`Cylinder`](@ref) or [`Sphere`](@ref), holding a mass, a density and its proportions.
The functions of the package are the same for every shape, and Julia picks the calculation from the type of the
shape it is given. This is *multiple dispatch*: to use a different shape, change the shape and nothing else.

```@setup shapes
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

```@example shapes
using BiophysicalGeometry, Unitful

mass, density = 2.0u"kg", 1000.0u"kg/m^3"
shapes = (Sphere(mass, density), Cylinder(mass, density, 2.0), Ellipsoid(mass, density, 2.0, 1.0),
          Plate(mass, density, 2.0, 4.0), Cone(mass, density, 2.0, 0.4))
[total_area(Body(shape, Naked())) for shape in shapes]
```

```@example shapes
shape_gallery(("$(nameof(typeof(shape)))" => Body(shape, Naked()) for shape in shapes)...; ncols = 5) # hide
```

## From mass to dimensions

The volume of a body is its mass over its density. The shape and its ratios then give the dimensions, and the
dimensions give the areas. [`Body`](@ref) does this once, and keeps the result in its `geometry`:

```@example shapes
body = Body(Cylinder(mass, density, 2.0), Naked())
body.geometry.volume, body.geometry.length, body.geometry.area.total
```

## The shapes

Each shape has its own dimensions, named in `body.geometry.length`, and its own position relative to the axes,
which matters when shapes are joined, see [Surfaces](surfaces.md).

::: tabs

== Sphere

`Sphere(mass, density)`. Centred on the origin.

```@example shapes
Body(Sphere(mass, density), Naked()).geometry.length
```

== Cylinder

`Cylinder(mass, density, axis_ratio_b)`, where `axis_ratio_b` is the length over the diameter. The axis is ``z``,
from ``z = 0`` to the length.

```@example shapes
Body(Cylinder(mass, density, 2.0), Naked()).geometry.length
```

== Ellipsoid

`Ellipsoid(mass, density, axis_ratio_b, axis_ratio_c)`, where `axis_ratio_b` is the long semi-axis over the short
one. The long axis is ``x``, centred on the origin. The ellipsoid is prolate, with two equal short axes, and
`axis_ratio_c` is not yet used.

```@example shapes
Body(Ellipsoid(mass, density, 2.0, 1.0), Naked()).geometry.length
```

== Plate

`Plate(mass, density, axis_ratio_b, axis_ratio_c)`, where the width is the length over `axis_ratio_b` and the
height is the length over `axis_ratio_c`. Length, width and height are along ``x``, ``y`` and ``z``, centred on
the origin.

```@example shapes
Body(Plate(mass, density, 2.0, 4.0), Naked()).geometry.length
```

== Cone

`Cone(mass, density, axis_ratio_b, top_ratio)`, where `axis_ratio_b` is the length over the diameter of the base,
and `top_ratio` the radius of the top over that of the base, see
[Frustums, truncation and halves](truncation.md). The axis is ``z``, with the base at ``z = 0``.

```@example shapes
Body(Cone(mass, density, 2.0, 0.4), Naked()).geometry.length
```

:::

The effect of the ratio on a cylinder of constant mass, from a disc to a rod:

```@example shapes
ratios = exp10.(range(-1, 1.5; length = 60))
areas = [total_area(Body(Cylinder(mass, density, ratio), Naked())) for ratio in ratios]
fig, ax = figure_axis("Length / diameter", "Surface area (cm²)"; xscale = log10)
lines!(ax, ratios, ustrip.(u"cm^2", areas); linewidth = 2, label = "Cylinder")
hlines!(ax, ustrip(u"cm^2", total_area(Body(Sphere(mass, density), Naked()))); linestyle = :dash, color = :black,
        label = "Sphere")
axislegend(ax; position = :ct)
fig
```

The sphere has the least area for its volume, which is why animals curl up in the cold.

## Shape families

Shapes belong to families, which other packages use to choose their equations: the heat lost from a cylinder and
from a cone is computed in the same way by
[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl).

| Family | Shapes |
|:--|:--|
| [`AbstractCylindrical`](@ref) | [`Cylinder`](@ref), [`Cone`](@ref) |
| [`AbstractSpherical`](@ref) | [`Sphere`](@ref) |
| [`AbstractEllipsoidal`](@ref) | [`Ellipsoid`](@ref) |
| [`AbstractSlab`](@ref) | [`Plate`](@ref) |

```@example shapes
Cone(mass, density, 2.0, 0.4) isa AbstractCylindrical
```

A [`Half`](@ref) shape wraps one of these, and belongs with the shape it wraps.

## The same functions for every shape

| Function | Returns |
|:--|:--|
| [`total_area`](@ref) | outer surface area, of the fur if there is any |
| [`skin_area`](@ref) | area of the skin |
| [`evaporation_area`](@ref) | skin area not covered by the bases of fibres |
| [`skin_radius`](@ref), [`flesh_radius`](@ref), [`insulation_radius`](@ref) | radius at each layer, see [Layers](layers.md) |
| [`flesh_volume`](@ref) | volume inside the fat |
| [`silhouette`](@ref) | area projected towards the sun, see [Silhouettes](silhouettes.md) |
| [`shape`](@ref), [`insulation`](@ref), [`geometry`](@ref) | the parts of a body |
| [`mass`](@ref) | mass of a shape |

The radius of a shape that is not round is its half-width: the short semi-axis of an ellipsoid, and half the width
of a plate.
