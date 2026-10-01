# Layers

A body is a set of shells, one inside the other: flesh at the centre, then fat, the skin, and fur or feathers. The
skin is fixed by the mass and density of the shape. Fat is placed inside it and fibres outside it.

```@setup layers
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

```@example layers
using BiophysicalGeometry, Unitful

shape = Cylinder(2.0u"kg", 1000.0u"kg/m^3", 2.0)
fur = FibrousLayer(20.0u"mm", 30.0u"μm", 3000u"cm^-2")   # thickness, fibre diameter, fibres per area of skin
fat = FatLayer(0.2, 901.0u"kg/m^3")                      # fraction of body mass, density
layer_diagram("Naked()" => Body(shape, Naked()), # hide
    "CompositeInsulation(fur, fat)" => Body(shape, CompositeInsulation(fur, fat))) # hide
```

## The layer types

::: tabs

== Naked

No insulation. All three radii are the same.

```@example layers
body = Body(shape, Naked())
flesh_radius(body), skin_radius(body), insulation_radius(body)
```

== FibrousLayer

`FibrousLayer(thickness, fibre_diameter, fibre_density)` is fur, feathers, hair or clothing, outside the skin. Its
three values are described under [Fibres](#Fibres).

```@example layers
body = Body(shape, fur)
flesh_radius(body), skin_radius(body), insulation_radius(body)
```

== FatLayer

`FatLayer(fraction, density)` is subcutaneous fat: a fraction of the body mass, with its own density, in a layer
of even thickness under the skin.

```@example layers
body = Body(shape, fat)
flesh_radius(body), skin_radius(body), insulation_radius(body)
```

== Both

`CompositeInsulation(fur, fat)`, with the fibres first.

```@example layers
body = Body(shape, CompositeInsulation(fur, fat))
flesh_radius(body), skin_radius(body), insulation_radius(body)
```

:::

## What each layer changes

| | [`FibrousLayer`](@ref) | [`FatLayer`](@ref) |
|:--|:--|:--|
| Where | outside the skin | inside the skin |
| Set by | a thickness | a fraction of body mass |
| [`total_area`](@ref) | grows to the outside of the fibres | unchanged |
| [`skin_area`](@ref) | unchanged | unchanged |
| [`evaporation_area`](@ref) | skin area less the bases of the fibres | unchanged |
| [`insulation_radius`](@ref) | skin radius plus the thickness | unchanged |
| [`flesh_radius`](@ref) | unchanged | skin radius less the fat thickness |
| `geometry.length` | gains the outer dimensions | gains `fat`, the thickness |

Fibres also cover the ends of a cylinder, so its outer length grows by twice the thickness:

```@example layers
Body(shape, CompositeInsulation(fur, fat)).geometry.length
```

Fatter and furrier bodies of the same mass:

```@example layers
shape_gallery(("$(Int(100f))% fat, $(d) mm fur" => # hide
    Body(shape, CompositeInsulation(FibrousLayer(d * u"mm", 30.0u"μm", 3000u"cm^-2"), FatLayer(f, 901.0u"kg/m^3"))) # hide
    for (f, d) in ((0.05, 5), (0.2, 20), (0.4, 40)))...) # hide
```

## Fibres

A [`FibrousLayer`](@ref) has three values:

| Value | Meaning | In the example |
|:--|:--|:--|
| `thickness` | depth of the coat, from the skin to its outer surface | 20 mm |
| `fibre_diameter` | diameter of one hair, or one barb of a feather | 30 μm |
| `fibre_density` | number of fibres per area of skin | 3000 cm⁻² |

```@example layers
plot_insulation_properties(fur)
```

Only the thickness changes the dimensions of a body. The diameter and density give the fraction of the skin
covered by the bases of fibres, ``f = \pi (d / 2)^2 N``, which is taken off the [`evaporation_area`](@ref):

```@example layers
covered = π * (fur.fibre_diameter / 2)^2 * fur.fibre_density
body = Body(shape, fur)
uconvert(NoUnits, covered), 1 - evaporation_area(body) / skin_area(body)
```

All three are used again by [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) for the heat
that passes through the coat. Heat is conducted through the still air between the fibres and along the fibres
themselves, so the conductivity of the coat lies between that of air and that of keratin, according to the
fraction of the layer that is fibre (Conley and Porter 1986). For this HeatExchange.jl adds the length of the
fibres, their conductivity and their reflectance. Fibres longer than the coat is deep lie at an angle, which puts
more fibre into the layer, by the ratio of length to thickness. The diameter and density also set how far
radiation penetrates the coat.

## Shells and heat flow

Heat made in the flesh is conducted out through each shell in turn, and the resistance of a shell depends on its
inner and outer radius. [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) reads those radii
from the body:

| Shell | Inner radius | Outer radius |
|:--|:--|:--|
| Flesh, where heat is produced | centre | [`flesh_radius`](@ref) |
| Fat | [`flesh_radius`](@ref) | [`skin_radius`](@ref) |
| Fibres | [`skin_radius`](@ref) | [`insulation_radius`](@ref) |

The shells of each shape family, cut across the long axis and along it:

```@example layers
layers = CompositeInsulation(fur, fat) # hide
layer_sections("Cylinder" => Body(Cylinder(2.0u"kg", 1000.0u"kg/m^3", 2.0), layers), # hide
    "Sphere" => Body(Sphere(2.0u"kg", 1000.0u"kg/m^3"), layers), # hide
    "Ellipsoid" => Body(Ellipsoid(2.0u"kg", 1000.0u"kg/m^3", 2.0, 1.0), layers), # hide
    "Cone" => Body(Cone(2.0u"kg", 1000.0u"kg/m^3", 2.0, 0.4), layers)) # hide
```

Fat is an even layer under the skin. On a cylinder or cone it lines the sides and not the ends. Fibres cover the
whole surface, ends included.

A body has at most one layer of fat and one of fibres. A second layer of flesh, or clothing over fur, cannot yet be
represented.
