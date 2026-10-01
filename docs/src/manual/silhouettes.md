# Silhouettes

The *silhouette area* is the area of the shadow a body would cast on a surface facing the sun. It is the area that
intercepts direct sunlight, and it changes with the direction of the sun.

```@setup silhouettes
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## A single shape

For one shape the silhouette is known exactly. [`silhouette`](@ref) gives the largest and smallest: `normal`, with
the long axis at right angles to the sun, and `parallel`, with the long axis pointing at it.

```@example silhouettes
using BiophysicalGeometry, Unitful

body = Body(Cylinder(2.0u"kg", 1000.0u"kg/m^3", 3.0), Naked())
silhouette(body)
```

An orientation picks one of them, or their mean:

```@example silhouettes
silhouette(body, NormalToSun()), silhouette(body, ParallelToSun()), silhouette(body, Intermediate())
```

An angle gives the silhouette with the sun at that angle from the long axis of the shape, for cylinders, cones and
ellipsoids:

```@example silhouettes
angles = (0.0:1.0:90.0) .* u"°"
fig, ax = figure_axis("Angle between the long axis and the sun (°)", "Silhouette area (cm²)")
for (label, shape) in ("Cylinder" => Cylinder(2.0u"kg", 1000.0u"kg/m^3", 3.0),
                       "Ellipsoid" => Ellipsoid(2.0u"kg", 1000.0u"kg/m^3", 3.0, 1.0))
    areas = [silhouette(Body(shape, Naked()), angle) for angle in angles]
    lines!(ax, ustrip.(angles), ustrip.(u"cm^2", areas); linewidth = 2, label)
end
axislegend(ax; position = :rb)
fig
```

For an animal standing upright this angle is the zenith angle of the sun, and for one lying along the ground it
depends on the azimuth too.

## Several parts

Parts shade each other, so the silhouette of a body is less than the sum of those of its parts. With a
[`Beam`](@ref), the direction towards the sun as ``(x, y, z)``, the body is projected as a whole, and each part
gets only the area that is lit:

```@example silhouettes
density = 1000.0u"kg/m^3"
torso = Body(Cylinder(20.0u"kg", density, 3.0), Naked())
head = Body(Sphere(2.0u"kg", density), Naked())
patch = Disc(5.0u"cm")
body = CompositeBody(; parts = (; torso, head), joins = (
    Join(torso = Attachment(EndB(0.0u"m", 0.0), patch), head = Attachment(Radial(π, 0.0), patch)),))

overhead = silhouette(body, Beam(0.0, 0.0, 1.0))
```

With the sun overhead the head shades the middle of the upright torso below it, which is lit only around its
rim. From the side both are in full view:

```@example silhouettes
side = silhouette(body, Beam(1.0, 0.0, 0.0))
```

```@example silhouettes
fig = Figure(size = (640, 330)) # hide
draw_parts!(body_axis(fig[1, 1]; decorations = false), body) # hide
for (i, (label, direction)) in enumerate(("Overhead" => (0.0, 0.0, 1.0), "Side" => (1.0, 0.0, 0.0))) # hide
    ax = Axis(fig[1, i + 1]; xlabel = "cm", ylabel = "cm") # hide
    area = silhouette_panel!(ax, body, direction) # hide
    ax.title = "$label: $(round(Int, ustrip(u"cm^2", area))) cm²" # hide
end # hide
fig # hide
```

[`silhouette_rasterized`](@ref) gives the total for a direction, and the sum over the parts ignoring shade is
`silhouette(body)`:

```@example silhouettes
silhouette_rasterized(body, (0.0, 0.0, 1.0)), sum(overhead), silhouette(body).parallel
```

These are computed by drawing the body on a grid, so they are approximate. `resolution`, the number of cells
across the grid, sets the accuracy.

## Sky, ground and neighbours

Diffuse sunlight and thermal radiation come from every direction. [`silhouette_factors`](@ref) gives, for each
part, the fractions of its view taken up by the sky, the ground and each of the other parts. They sum to one.

```@example silhouettes
factors = silhouette_factors(body, Sky(0.5))
factors.head
```

The second argument says which directions are sky:

| Region | Meaning |
|:--|:--|
| [`Sky`](@ref)`(fraction)` | the fraction of all directions that is sky: 0.5 on flat open ground, more on a summit, 0 in a closed burrow |
| [`Ground`](@ref)`(fraction)` | the same split, given as the fraction that is ground |
| [`Horizon`](@ref)`(angles)` | hills: the angle of the horizon at each compass bearing, as in [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) |

```@example silhouettes
valley = Horizon(fill(30.0u"°", 24))
silhouette_factors(body, valley).head.sky, factors.head.sky
```

## Which to use

| Have | Want | Use |
|:--|:--|:--|
| one shape | largest and smallest silhouette | `silhouette(body)` |
| one shape | silhouette at an angle | `silhouette(body, angle)` |
| several parts | lit area of each part | `silhouette(body, Beam(x, y, z))` |
| several parts | lit area of the whole | `silhouette_rasterized(body, (x, y, z))` |
| several parts | views of sky, ground and other parts | `silhouette_factors(body, Sky(0.5))` |

The functions for several parts need a [`CompositeBody`](@ref). A single body can be given as one with one part,
`CompositeBody(; parts = (; body), joins = ())`.
