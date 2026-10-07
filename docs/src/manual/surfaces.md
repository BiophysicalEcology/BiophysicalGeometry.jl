# Surfaces

Parts are joined at their surfaces, so each surface of each shape has a name. `A` and `B` name the two ends of a
shape along its long axis: `A` is the end at the start of the axis for a cylinder or cone, and the end that can
be truncated for an ellipsoid.

```@setup surfaces
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

```@example surfaces
using BiophysicalGeometry, Unitful

mass, density = 2.0u"kg", 1000.0u"kg/m^3"
surface_diagram("Cylinder" => Body(Cylinder(; mass, density, axis_ratio_b = 2.0), Naked()), # hide
    "Cone" => Body(Cone(; mass, density, axis_ratio_b = 2.0, top_ratio = 0.4), Naked()), # hide
    "Ellipsoid" => Body(Ellipsoid(; mass, density, axis_ratio_b = 2.0, axis_ratio_c = 2.0), Naked()), # hide
    "Ellipsoid, truncated" => Body(Ellipsoid(; mass, density, axis_ratio_b = 2.0, axis_ratio_c = 2.0, pole_a_truncation = 0.3), Naked()), # hide
    "HalfCylinder" => Body(HalfCylinder(; mass = mass / 2, density, axis_ratio_b = 2.0), Naked()), # hide
    "HalfEllipsoid" => Body(HalfEllipsoid(; mass = mass / 2, density, axis_ratio_b = 2.0, axis_ratio_c = 4.0), Naked()); ncols = 3) # hide
```

```@example surfaces
surface_diagram("Sphere" => Body(Sphere(; mass, density), Naked()), # hide
    "Plate" => Body(Plate(; mass, density, axis_ratio_b = 2.0, axis_ratio_c = 4.0), Naked()), # hide
    "TriangularPlate" => Body(TriangularPlate(; mass, density, axis_ratio_b = 1.5, axis_ratio_c = 6.0), Naked())) # hide
```

[`attachment_surfaces`](@ref) lists the surfaces of a shape:

```@example surfaces
attachment_surfaces(Cylinder(; mass, density, axis_ratio_b = 2.0))
```

## The surfaces

| Surface | On | Where | Point on it |
|:--|:--|:--|:--|
| [`EndA`](@ref) | cylinder, cone, half-cylinder, half-cone | flat end at ``x = 0``; the base of a cone | `EndA(radius, angle)`: distance from the axis, angle around it |
| [`EndB`](@ref) | cylinder, cone, half-cylinder, half-cone | flat end at ``x`` = length; the top of a cone | `EndB(radius, angle)` |
| [`Lateral`](@ref) | cylinder, cone, half-cylinder, half-cone | the curved side | `Lateral(position, angle)`: distance along the axis, angle around it |
| [`PoleA`](@ref) | ellipsoid | tip of the long axis at ``+x``; a disc if truncated | `PoleA()` |
| [`PoleB`](@ref) | ellipsoid | tip of the long axis at ``-x`` | `PoleB()` |
| [`Equator`](@ref) | ellipsoid | the ring around the middle | `Equator(angle)`: angle around the long axis |
| [`Radial`](@ref) | sphere | anywhere | `Radial(polar, azimuth)`: angle from ``+z``, angle around it from ``+x`` |
| [`Flat`](@ref) | half shapes | the cut face, at ``z = 0`` | `Flat(x, y)`: position on the face |
| [`Dome`](@ref) | half-ellipsoid, half-sphere | the curved side | `Dome(polar, azimuth)`: angle from ``+x``, angle around it |
| [`Top`](@ref), [`Bottom`](@ref) | plate, triangular plate | faces at ``\pm z`` | `Top(x, y)` |
| [`SideA`](@ref), [`SideB`](@ref) | plate; `SideB` on a triangular plate | faces at ``\pm x`` | `SideA(y, z)` |
| [`SideC`](@ref), [`SideD`](@ref) | plate; `SideD` on a triangular plate | faces at ``\pm y`` | `SideC(x, z)` |
| [`Diagonal`](@ref) | triangular plate | the long edge | `Diagonal(position, z)`: distance along it from its ``+x`` end, height |

Distances have units and angles are in radians. Angles around the axis of a cylinder or cone run from ``+y``
towards ``+z``. On a half-cylinder the angle runs from 0 to ``π``, over the curved side only, with ``π / 2`` at
its top.

## A whole surface or a point on it

A surface is written in two ways. With no arguments it means the whole surface, as when asking for its area:

```@example surfaces
body = Body(Cylinder(; mass, density, axis_ratio_b = 2.0), Naked())
surface_area(body.shape, body, EndA()), surface_area(body.shape, body, Lateral())
```

With arguments it means a point on the surface, which is where the centre of a join is placed. The point halfway
along a cylinder, on its upper side:

```@example surfaces
half_length = body.geometry.length.length_skin / 2
BiophysicalGeometry.surface_point(body.shape, body, Lateral(half_length, π / 2))
```

and the direction the surface faces there:

```@example surfaces
BiophysicalGeometry.surface_normal(body.shape, body, Lateral(half_length, π / 2))
```

Moving a point around a surface is how a limb is moved around a body, see [Joins and poses](joins.md).

## Skin and fur

Points are on the skin, so that joined parts meet flesh to flesh however thick their fur. Areas are those of the
outside of the fur.
