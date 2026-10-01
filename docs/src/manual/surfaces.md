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
surface_diagram("Cylinder" => Body(Cylinder(mass, density, 2.0), Naked()), # hide
    "Cone" => Body(Cone(mass, density, 2.0, 0.4), Naked()), # hide
    "Ellipsoid" => Body(Ellipsoid(mass, density, 2.0, 1.0), Naked()), # hide
    "Ellipsoid, truncated" => Body(Ellipsoid(mass, density, 2.0, 1.0, 0.3), Naked()), # hide
    "HalfCylinder" => Body(HalfCylinder(mass / 2, density, 2.0), Naked()), # hide
    "HalfEllipsoid" => Body(HalfEllipsoid(mass / 2, density, 2.0, 1.0), Naked()); ncols = 3) # hide
```

```@example surfaces
surface_diagram("Sphere" => Body(Sphere(mass, density), Naked()), # hide
    "Plate" => Body(Plate(mass, density, 2.0, 4.0), Naked())) # hide
```

[`attachment_surfaces`](@ref) lists the surfaces of a shape:

```@example surfaces
attachment_surfaces(Cylinder(mass, density, 2.0))
```

## The surfaces

| Surface | On | Where | Point on it |
|:--|:--|:--|:--|
| [`EndA`](@ref) | cylinder, cone, half-cylinder | flat end at ``z = 0``; the base of a cone | `EndA(r, φ)`: distance from the axis, angle around it |
| [`EndB`](@ref) | cylinder, cone, half-cylinder | flat end at ``z`` = length; the top of a cone | `EndB(r, φ)` |
| [`Lateral`](@ref) | cylinder, cone, half-cylinder | the curved side | `Lateral(z, φ)`: distance along the axis, angle around it |
| [`PoleA`](@ref) | ellipsoid | tip of the long axis at ``+x``; a disc if truncated | `PoleA()` |
| [`PoleB`](@ref) | ellipsoid | tip of the long axis at ``-x`` | `PoleB()` |
| [`Equator`](@ref) | ellipsoid | the ring around the middle | `Equator(φ)`: angle around the long axis |
| [`Radial`](@ref) | sphere | anywhere | `Radial(θ, φ)`: angle from ``+z``, angle around it |
| [`Flat`](@ref) | half shapes | the cut face | `Flat(u, v)`: position on the face |
| [`Dome`](@ref) | half-ellipsoid, half-sphere | the curved side | `Dome(α, β)`: angle from ``+x``, angle around it |
| [`Top`](@ref), [`Bottom`](@ref) | plate | faces at ``\pm z`` | `Top(x, y)` |
| [`SideA`](@ref), [`SideB`](@ref) | plate | faces at ``\pm x`` | `SideA(y, z)` |
| [`SideC`](@ref), [`SideD`](@ref) | plate | faces at ``\pm y`` | `SideC(x, z)` |

Distances have units and angles are in radians. On a half-cylinder the angle ``φ`` runs from 0 to ``π``, over the
curved side only.

## A whole surface or a point on it

A surface is written in two ways. With no arguments it means the whole surface, as when asking for its area:

```@example surfaces
body = Body(Cylinder(mass, density, 2.0), Naked())
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
