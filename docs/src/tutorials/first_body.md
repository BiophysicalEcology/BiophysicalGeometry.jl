# A first body

This tutorial builds a small mammal as a single shape, sees how its shape changes its surface area, and gives it
fat and fur.

```@setup first
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## A sphere

The simplest animal is a sphere. All it needs is a mass and a density, close to that of water:

```@example first
using BiophysicalGeometry, Unitful

mass = 500.0u"g"
density = 1000.0u"kg/m^3"
sphere = Body(Sphere(; mass, density), Naked())
uconvert(u"cm", skin_radius(sphere)), uconvert(u"cm^2", total_area(sphere))
```

`Naked()` says that it has no insulation.

## Stretching it

A cylinder with a length three times its diameter is more like a rat. Only the shape changes:

```@example first
cylinder = Body(Cylinder(; mass, density, axis_ratio_b = 3.0), Naked())
uconvert(u"cm^2", total_area(cylinder))
```

It has the same volume and more surface, so it would lose heat faster. Its dimensions:

```@example first
map(x -> uconvert(u"cm", x), cylinder.geometry.length)
```

An ellipsoid has rounded ends:

```@example first
ellipsoid = Body(Ellipsoid(; mass, density, axis_ratio_b = 3.0, axis_ratio_c = 3.0), Naked())
shape_gallery("Sphere" => sphere, "Cylinder" => cylinder, "Ellipsoid" => ellipsoid) # hide
```

```@example first
markdown_table(["Shape", "Surface area", "Relative to sphere"],
               [(name, uconvert(u"cm^2", total_area(body)), uconvert(NoUnits, total_area(body) / total_area(sphere)))
                for (name, body) in ("Sphere" => sphere, "Cylinder" => cylinder, "Ellipsoid" => ellipsoid)])
```

## Fur

Fur is a [`FibrousLayer`](@ref): a depth, the diameter of a hair, and the number of hairs per area of skin.

```@example first
fur = FibrousLayer(8.0u"mm", 20.0u"μm", 8000u"cm^-2")
furry = Body(Cylinder(; mass, density, axis_ratio_b = 3.0), fur)
uconvert(u"cm", skin_radius(furry)), uconvert(u"cm", insulation_radius(furry))
```

The animal is now larger on the outside, where it meets the air and the sun, and the same at the skin:

```@example first
uconvert(u"cm^2", total_area(furry)), uconvert(u"cm^2", skin_area(furry))
```

## Fat

Fat is a [`FatLayer`](@ref): a fraction of the body mass, and a density. It goes under the skin, so the flesh
shrinks and the outside stays the same:

```@example first
fat = FatLayer(0.15, 901.0u"kg/m^3")
fat_and_fur = Body(Cylinder(; mass, density, axis_ratio_b = 3.0), CompositeInsulation(fur, fat))
map(x -> uconvert(u"cm", x), (flesh_radius(fat_and_fur), skin_radius(fat_and_fur), insulation_radius(fat_and_fur)))
```

```@example first
layer_diagram("Naked" => cylinder, "Fur and fat" => fat_and_fur) # hide
```

```@example first
plot_body(fat_and_fur)
```

## In the sun

The silhouette area is what the sun sees. It is largest side-on and smallest end-on:

```@example first
map(x -> uconvert(u"cm^2", x), silhouette(fat_and_fur))
```

These areas, radii and volumes are what a heat budget needs, see
[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl). The [next tutorial](two_parts.md)
gives the animal a back and a belly.
