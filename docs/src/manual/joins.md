# Joins and poses

A [`Join`](@ref) fixes one part to another. There are no separate commands to move or rotate a part: its
position follows from where on each surface the join is made, and its remaining freedom, spinning about the join,
is the `twist`.

```@setup joins
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## A join

A join names two parts, and gives each an [`Attachment`](@ref): a point on a surface, see
[Surfaces](surfaces.md), and the patch that the join covers. The first part named is the parent and the second
the child.

```@example joins
using BiophysicalGeometry, Unitful

density = 1000.0u"kg/m^3"
torso = Body(Cylinder(; mass = 20.0u"kg", density, axis_ratio_b = 3.0), Naked())
leg = Body(Cylinder(; mass = 1.0u"kg", density, axis_ratio_b = 5.0), Naked())
r = skin_radius(leg)
half_length = torso.geometry.length.length_skin / 2

join = Join(torso = Attachment(Lateral(half_length, 0.0), Disc(r)),
            leg = Attachment(EndA(0.0u"m", 0.0), Disc(r)))
body = CompositeBody(; parts = (; torso, leg), joins = (join,))
composite_views(body; views = (:oblique, :top)) # hide
```

The child is placed so that the two points touch and the two surfaces face each other. Here the end of the leg
faces the side of the torso, so the leg stands out from it at a right angle.

## Patches

The patch is the area hidden by the join, and is taken off the area of both parts.

| Patch | Covers | Used for |
|:--|:--|:--|
| [`Disc`](@ref)`(radius)` | a circle around the point | limbs, necks, heads |
| [`FullCover`](@ref)`()` | the whole surface | the flat faces of two halves |

The two patches of a join must have the same area, and each must fit on its surface. With `FullCover` the surface
is written without a point, as `Flat()`.

```@example joins
join_area(join, body), π * r^2
```

## Moving a part

Changing the point on the parent moves the child over its surface. Along the torso:

```@example joins
place(position, angle) = CompositeBody(; parts = (; torso, leg), joins = (
    Join(torso = Attachment(Lateral(position * 2half_length, angle), Disc(r)),
         leg = Attachment(EndA(0.0u"m", 0.0), Disc(r))),))
fig = Figure(size = (780, 250)) # hide
for (i, position) in enumerate((0.1, 0.5, 0.9)) # hide
    ax = body_axis(fig[1, i]; decorations = false, title = "Lateral($position * length, 0)", titlesize = 12) # hide
    draw_parts!(ax, place(position, 0.0)) # hide
end # hide
fig # hide
```

and around it, from ``+y`` towards ``+z``, seen from the end:

```@example joins
fig = Figure(size = (780, 250)) # hide
for (i, (angle, label)) in enumerate(((0.0, "0"), (π / 4, "π / 4"), (π / 2, "π / 2"))) # hide
    ax = body_axis(fig[1, i]; decorations = false, azimuth = 0.0, elevation = 0.0, # hide
        title = "Lateral(length / 2, $label)", titlesize = 12) # hide
    draw_parts!(ax, place(0.5, angle)) # hide
end # hide
fig # hide
```

Changing the point on the child changes which part of it touches: joined by its side, the leg lies along the
torso.

## Twist

Two surfaces facing each other can still spin about the line through the join. `twist` is that angle, in radians.
It makes no difference to a round leg joined by its end, but it does to a part joined by its side, or to a half
shape:

```@example joins
arm_length = leg.geometry.length.length_skin
twisted(twist) = CompositeBody(; parts = (; torso, leg), joins = (
    Join(torso = Attachment(Lateral(half_length, 0.0), Disc(r / 2)),
         leg = Attachment(Lateral(0.1arm_length, 0.0), Disc(r / 2)); twist),))
fig = Figure(size = (780, 250)) # hide
for (i, (twist, label)) in enumerate(((0.0, "0"), (π / 4, "π / 4"), (π / 2, "π / 2"))) # hide
    ax = body_axis(fig[1, i]; decorations = false, azimuth = -π / 2, elevation = π / 2, title = "twist = $label", # hide
        titlesize = 12) # hide
    draw_parts!(ax, twisted(twist)) # hide
end # hide
fig # hide
```

## The pose of the root

The first part listed is the root, and by default it sits as its shape is defined: every shape lies along ``x``
with its height along ``z``. `root_pose` gives it another position and orientation, as a [`Pose`](@ref) of a
translation and a rotation matrix, and every other part follows. The columns of the matrix are where the ``x``,
``y`` and ``z`` axes of the root end up. To stand a cylinder up:

```@example joins
rotation = [0.0 0.0 -1.0;   # columns: x of the cylinder, its axis, points up along z;
            0.0 1.0 0.0;    # y stays along y;
            1.0 0.0 0.0]    # and z points along -x
standing = CompositeBody(; parts = (; torso, leg), joins = (join,),
                         root_pose = Pose((0.0u"m", 0.0u"m", 0.0u"m"), rotation))
fig = Figure(size = (560, 250)) # hide
draw_parts!(body_axis(fig[1, 1]; title = "default", titlesize = 12), body) # hide
draw_parts!(body_axis(fig[1, 2]; title = "root_pose", titlesize = 12), standing) # hide
fig # hide
```

The solved pose of every part is kept in the body:

```@example joins
standing.poses.leg.translation
```

## Joins and heat

For conduction of heat between parts, [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl)
needs the area of each join and the distance from the centre of each part to it:

```@example joins
join_area(join, body), internal_distance(join, body), join_position(join, body)
```

| Function | Returns |
|:--|:--|
| [`join_partners`](@ref) | names of the parent and child |
| [`join_area`](@ref) | area of the patch |
| [`join_position`](@ref) | position of the centre of the join |
| [`internal_distance`](@ref) | distance from the centre of each part to the join |
| [`flesh_centroid`](@ref) | centre of a part, relative to its own axes |
