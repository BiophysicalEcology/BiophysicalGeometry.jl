# Get started

BiophysicalGeometry.jl builds the body of an organism from simple shapes and returns its areas, lengths and volumes
as [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities.

```julia
using Pkg
Pkg.add(url = "https://github.com/BiophysicalEcology/BiophysicalGeometry.jl")
```

```@setup get_started
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## A body

A body is a shape with its insulation. A shape is sized from a mass, a density and, for most shapes, a ratio of
length to width:

```@example get_started
using BiophysicalGeometry, Unitful

shape = Cylinder(2.0u"kg", 1000.0u"kg/m^3", 3.0)
body = Body(shape, Naked())
total_area(body)
```

The dimensions and volume are in `body.geometry`:

```@example get_started
body.geometry.length, body.geometry.volume
```

## Other shapes

The same functions work for every shape. Here are four bodies of the same mass:

::: tabs

== Sphere

```@example get_started
total_area(Body(Sphere(2.0u"kg", 1000.0u"kg/m^3"), Naked()))
```

== Cylinder

```@example get_started
total_area(Body(Cylinder(2.0u"kg", 1000.0u"kg/m^3", 3.0), Naked()))
```

== Ellipsoid

```@example get_started
total_area(Body(Ellipsoid(2.0u"kg", 1000.0u"kg/m^3", 3.0, 1.0), Naked()))
```

== Plate

```@example get_started
total_area(Body(Plate(2.0u"kg", 1000.0u"kg/m^3", 3.0, 6.0), Naked()))
```

:::

See [Shapes](manual/shapes.md).

## Layers

Fur or feathers sit outside the skin, and fat inside it:

```@example get_started
fur = FibrousLayer(20.0u"mm", 30.0u"μm", 3000u"cm^-2")   # depth, fibre diameter, fibres per area
fat = FatLayer(0.2, 901.0u"kg/m^3")                      # fraction of body mass, density
furry = Body(shape, CompositeInsulation(fur, fat))
flesh_radius(furry), skin_radius(furry), insulation_radius(furry)
```

The total area is now that of the outside of the fur, and the skin area is unchanged:

```@example get_started
total_area(furry), skin_area(furry)
```

See [Layers](manual/layers.md).

## Plotting

With a Makie backend loaded, [`plot_body`](@ref) draws a body with a quarter cut away to show its layers:

```@example get_started
using CairoMakie

plot_body(furry)
```

Makie also exports `Sphere`, `Top` and `Bottom`. When both packages are loaded, import them from this one:
`import BiophysicalGeometry: Sphere, Top, Bottom`.

## More than one part

The back of an animal faces the sun and sky and its belly faces the ground, often with a different coat. The most
common body of more than one part is therefore a shape split into a dorsal and a ventral half, joined over their
flat faces into a [`CompositeBody`](@ref):

```@example get_started
dorsal = Body(HalfEllipsoid(1.0u"kg", 1000.0u"kg/m^3", 3.0, 1.0), FibrousLayer(20.0u"mm", 30.0u"μm", 3000u"cm^-2"))
ventral = Body(HalfEllipsoid(1.0u"kg", 1000.0u"kg/m^3", 3.0, 1.0), FibrousLayer(5.0u"mm", 30.0u"μm", 3000u"cm^-2"))
animal = CompositeBody(;
    parts = (; dorsal, ventral),
    joins = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),),
)
total_area(animal)
```

```@example get_started
composite_views(animal; views = (:oblique, :side, :front)) # hide
```

Each half has its own areas, and its own view of the sky and the ground:

```@example get_started
views = silhouette_factors(animal, Sky(0.5))
views.dorsal.sky, views.ventral.sky
```

Heads, limbs and tails are joined on in the same way. Try it on the interactive page [Build an animal](builder.md),
and see [Bodies as graphs](manual/graphs.md) and the [tutorials](tutorials/first_body.md) for worked examples up to
a cow and a human.
