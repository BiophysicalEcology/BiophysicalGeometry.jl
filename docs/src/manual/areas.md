# Areas and volumes

The area and volume functions work on a single [`Body`](@ref) and on a [`CompositeBody`](@ref). For a body of
several parts they sum over the parts, and take off the patches hidden by the joins.

```@setup areas
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## One part

```@example areas
using BiophysicalGeometry, Unitful

density = 1000.0u"kg/m^3"
fur = FibrousLayer(10.0u"mm", 30.0u"μm", 3000u"cm^-2")
fat = FatLayer(0.1, 901.0u"kg/m^3")
torso = Body(Cylinder(20.0u"kg", density, 3.0), CompositeInsulation(fur, fat))
(; total = total_area(torso), skin = skin_area(torso), evaporation = evaporation_area(torso),
   flesh_volume = flesh_volume(torso), volume = torso.geometry.volume)
```

| Function | Area or volume |
|:--|:--|
| [`total_area`](@ref) | the outside of the body, over the fur |
| [`skin_area`](@ref) | the skin |
| [`evaporation_area`](@ref) | the skin between the fibres |
| [`surface_area`](@ref)`(shape, body, surface)` | one named surface, such as `EndA()` |
| [`flesh_volume`](@ref) | the volume inside the fat |
| `body.geometry.volume` | the whole volume inside the skin |

## Several parts

Each join hides a patch on both of its parts:

```@example areas
head = Body(Ellipsoid(2.0u"kg", density, 1.5, 1.0), fur)
patch = Disc(4.0u"cm")
join = Join(torso = Attachment(EndB(0.0u"m", 0.0), patch), head = Attachment(PoleB(), patch))
body = CompositeBody(; parts = (; torso, head), joins = (join,))

total_area(body), total_area(torso) + total_area(head) - 2 * join_area(join, body)
```

Volumes add, as nothing is lost at a join:

```@example areas
flesh_volume(body), flesh_volume(torso) + flesh_volume(head)
```

The patch is taken at its full size from the total area, the skin area and the evaporation area alike.

## Halves make a whole

Two halves joined over the whole of their flat faces have the area and volume of the shape they were cut from:

```@example areas
dorsal = Body(HalfCylinder(10.0u"kg", density, 3.0), Naked())
ventral = Body(HalfCylinder(10.0u"kg", density, 3.0), Naked())
halves = CompositeBody(; parts = (; dorsal, ventral), joins = (
    Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),))
whole = Body(Cylinder(20.0u"kg", density, 3.0), Naked())

total_area(halves), total_area(whole)
```

The halves need not match. Thick fur above and thin fur below:

```@example areas
furry_back = Body(HalfCylinder(10.0u"kg", density, 3.0), FibrousLayer(30.0u"mm", 30.0u"μm", 3000u"cm^-2"))
thin_belly = Body(HalfCylinder(10.0u"kg", density, 3.0), FibrousLayer(5.0u"mm", 30.0u"μm", 3000u"cm^-2"))
two_coats = CompositeBody(; parts = (; furry_back, thin_belly), joins = (
    Join(furry_back = Attachment(Flat(), FullCover()), thin_belly = Attachment(Flat(), FullCover())),))
fig = Figure(size = (520, 260)) # hide
draw_parts!(body_axis(fig[1, 1]; decorations = false, azimuth = -π / 2, elevation = π / 2), two_coats) # hide
draw_parts!(body_axis(fig[1, 2]; decorations = false), two_coats) # hide
fig # hide
```

## Part by part

The parts of a composite are ordinary bodies, so any function can be mapped over them:

```@example areas
map(total_area, body.parts)
```

Functions that describe a single shape, such as [`skin_radius`](@ref), return the value for the root part when
given a composite.
