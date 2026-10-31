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
shapes = (Sphere(; mass, density), Cylinder(; mass, density, axis_ratio_b = 2.0),
          Ellipsoid(; mass, density, axis_ratio_b = 2.0, axis_ratio_c = 3.0),
          Plate(; mass, density, axis_ratio_b = 2.0, axis_ratio_c = 4.0),
          Cone(; mass, density, axis_ratio_b = 2.0, top_ratio = 0.4),
          TriangularPlate(; mass, density, axis_ratio_b = 1.5, axis_ratio_c = 8.0))
[total_area(Body(shape, Naked())) for shape in shapes]
```

```@example shapes
shape_gallery(("$(nameof(typeof(shape)))" => Body(shape, Naked()) for shape in shapes)...; ncols = 6) # hide
```

## Keywords

A shape is made from keywords, never from positional arguments: `mass`, `density` and `volume`, the shape's own
dimensions (such as `length` and `radius`) and its ratios (such as `axis_ratio_b`). Give any set that fixes the
shape and the rest is worked out — weigh an animal, or measure it:

```@example shapes
weighed = Cylinder(; mass, density, axis_ratio_b = 2.0)
l = Body(weighed, Naked()).geometry.length
measured = Cylinder(; length = l.length_skin, radius = l.radius_skin, density)
BiophysicalGeometry.mass(measured), measured.axis_ratio_b
```

Too little, or too much that doesn't agree, is an error that says what to give:

```@example shapes
try Cylinder(; mass) catch e; e end
```

The dimensions given are at the skin, inside any fur.

## From mass to dimensions

The volume of a body is its mass over its density. The shape and its ratios then give the dimensions, and the
dimensions give the areas. [`Body`](@ref) does this once, and keeps the result in its `geometry`:

```@example shapes
body = Body(Cylinder(; mass, density, axis_ratio_b = 2.0), Naked())
body.geometry.volume, body.geometry.length, body.geometry.area.total
```

## The shapes

Each shape has its own dimensions, named in `body.geometry.length`, and its own position relative to the axes,
which matters when shapes are joined, see [Surfaces](surfaces.md). Every shape lies along ``x``, with its width
along ``y`` and its height along ``z``; turn it with a pose.

::: tabs

== Sphere

`Sphere(; mass, density, volume, radius)`. Centred on the origin.

```@example shapes
Body(Sphere(; mass, density), Naked()).geometry.length
```

== Cylinder

`Cylinder(; mass, density, volume, length, radius, axis_ratio_b)`, where `axis_ratio_b` is the length over the
diameter. The axis is ``x``, from ``x = 0`` to the length.

```@example shapes
Body(Cylinder(; mass, density, axis_ratio_b = 2.0), Naked()).geometry.length
```

== Ellipsoid

`Ellipsoid(; mass, density, volume, length, width, height, axis_ratio_b, axis_ratio_c)`, a general ellipsoid with
three independent axes. As for a plate, `axis_ratio_b` is the length over the width and `axis_ratio_c` the length
over the height; give them equal for a spheroid, and below 1 for an oblate one. Length, width and height are the
full extents along ``x``, ``y`` and ``z``, centred on the origin. The area is exact, from Carlson's symmetric
elliptic integral.

```@example shapes
Body(Ellipsoid(; mass, density, axis_ratio_b = 2.0, axis_ratio_c = 3.0), Naked()).geometry.length
```

== Plate

`Plate(; mass, density, volume, length, width, height, axis_ratio_b, axis_ratio_c)`, where the width is the length
over `axis_ratio_b` and the height is the length over `axis_ratio_c`. Length, width and height are along ``x``,
``y`` and ``z``, centred on the origin.

```@example shapes
Body(Plate(; mass, density, axis_ratio_b = 2.0, axis_ratio_c = 4.0), Naked()).geometry.length
```

== Triangular plate

`TriangularPlate(; mass, density, volume, length, width, height, axis_ratio_b, axis_ratio_c)`: a plate cut along its
diagonal, for wings and ears. The right angle is at the origin, with the legs `length` along ``x`` and `width`
along ``y``, and the thickness `height` along ``z``, centred on ``z = 0``. Fur offsets every face outward.

```@example shapes
Body(TriangularPlate(; mass, density, axis_ratio_b = 1.5, axis_ratio_c = 8.0), Naked()).geometry.length
```

== Cone

`Cone(; mass, density, volume, length, radius, axis_ratio_b, top_ratio = 0.0)`, where `axis_ratio_b` is the length
over the diameter of the base, and `top_ratio` the radius of the top over that of the base, see
[Frustums, truncation and halves](truncation.md). The axis is ``x``, with the base at ``x = 0``.

```@example shapes
Body(Cone(; mass, density, axis_ratio_b = 2.0, top_ratio = 0.4), Naked()).geometry.length
```

:::

The effect of the ratio on a cylinder of constant mass, from a disc to a rod:

```@example shapes
ratios = exp10.(range(-1, 1.5; length = 60))
areas = [total_area(Body(Cylinder(; mass, density, axis_ratio_b), Naked())) for axis_ratio_b in ratios]
fig, ax = figure_axis("Length / diameter", "Surface area (cm²)"; xscale = log10)
lines!(ax, ratios, ustrip.(u"cm^2", areas); linewidth = 2, label = "Cylinder")
hlines!(ax, ustrip(u"cm^2", total_area(Body(Sphere(; mass, density), Naked()))); linestyle = :dash, color = :black,
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
| [`AbstractSlab`](@ref) | [`Plate`](@ref), [`TriangularPlate`](@ref) |

```@example shapes
Cone(; mass, density, axis_ratio_b = 2.0, top_ratio = 0.4) isa AbstractCylindrical
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

The radius of a shape that is not round is its half-width: half the width of an ellipsoid or a plate, and the
inradius of a triangular plate.
