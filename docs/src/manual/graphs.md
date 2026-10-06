# Bodies as graphs

A body with several parts is described by what is connected to what. This is a *graph*: a set of *nodes*, here
the parts, and *edges* between pairs of nodes, here the joins.

```@setup graphs
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## Two parts: back and belly

The simplest graph has two nodes and one edge. It is also the commonest body: a shape split into a dorsal and a
ventral half, so that the back and the belly can have different coats and face different surroundings.

```@example graphs
using BiophysicalGeometry, Unitful

density = 1000.0u"kg/m^3"
thick = FibrousLayer(20.0u"mm", 30.0u"μm", 3000u"cm^-2")
thin = FibrousLayer(5.0u"mm", 30.0u"μm", 3000u"cm^-2")
dorsal = Body(HalfEllipsoid(; mass = 1.0u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 6.0), thick)
ventral = Body(HalfEllipsoid(; mass = 1.0u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 6.0), thin)

animal = CompositeBody(;
    parts = (; dorsal, ventral),
    joins = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),),
)
body_graph(animal; size = (760, 240)) # hide
```

[`CompositeBody`](@ref) takes the two halves of the graph:

- `parts` are the nodes: named bodies, written `(; dorsal, ventral)`.
- `joins` are the edges. A [`Join`](@ref) names two parts and, for each, an [`Attachment`](@ref): the surface it
  is joined at and the patch that is covered, here the whole of the [`Flat`](@ref) face of each half. See
  [Joins and poses](joins.md).

## More parts: a dog

Limbs and a head are more nodes, each with an edge to the part it grows from. One body can be used under several
names, as `leg` is for the four legs:

```@example graphs
torso = Body(Cylinder(; mass = 20.0u"kg", density, axis_ratio_b = 3.0), Naked())
head = Body(Ellipsoid(; mass = 2.0u"kg", density, axis_ratio_b = 1.5, axis_ratio_c = 1.5), Naked())
leg = Body(Cylinder(; mass = 1.0u"kg", density, axis_ratio_b = 5.0), Naked())

torso_length = torso.geometry.length.length_skin
r = skin_radius(leg)
hip(z, φ) = Attachment(Lateral(z * torso_length, φ), Disc(r))
leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(r))

dog = CompositeBody(;
    parts = (; torso, head, leg_fl = leg, leg_fr = leg, leg_bl = leg, leg_br = leg),
    joins = (
        Join(torso = Attachment(EndB(0.0u"m", 0.0), Disc(4.0u"cm")), head = Attachment(PoleB(), Disc(4.0u"cm"))),
        Join(torso = hip(0.8, -π / 2 + 0.35), leg_fl = leg_top),
        Join(torso = hip(0.8, -π / 2 - 0.35), leg_fr = leg_top),
        Join(torso = hip(0.2, -π / 2 + 0.35), leg_bl = leg_top),
        Join(torso = hip(0.2, -π / 2 - 0.35), leg_br = leg_top),
    ),
)
body_graph(dog) # hide
```

```@example graphs
keys(dog.parts), join_partners.(dog.joins)
```

## Trees and the root

The graph of a body is a *tree*: it has no loops, so there is one path between any two parts. The first part
listed is the *root*, drawn with a black outline above. It is placed first, and every other part is placed from
the part it is joined to, working outwards along the edges. So each join must connect to a part that has already
been reached: list the joins from the root outwards.

Changing the root does not change the body, only which part stays still and which part is reported by the
functions that describe one shape, such as [`skin_radius`](@ref).

## A longer chain

A join can hang from any part, so a neck can go between the torso and the head. The graph gains a node and an
edge:

```@example graphs
neck = Body(Cone(; mass = 1.5u"kg", density, axis_ratio_b = 1.2, top_ratio = 0.7), Naked())
r_top = 0.7 * skin_radius(neck)
long_neck = CompositeBody(;
    parts = (; torso, neck, head),
    joins = (
        Join(torso = Attachment(EndB(0.0u"m", 0.0), Disc(6.0u"cm")), neck = Attachment(EndA(0.0u"m", 0.0), Disc(6.0u"cm"))),
        Join(neck = Attachment(EndB(0.0u"m", 0.0), Disc(r_top)), head = Attachment(PoleB(), Disc(r_top))),
    ),
)
body_graph(long_neck; size = (760, 240)) # hide
```

## What the graph gives

The position of every part is worked out from the graph when the body is made, and kept in `poses`. Areas of the
whole body are sums over the nodes, less the patches covered at each edge, see [Areas and volumes](areas.md):

```@example graphs
total_area(dog), sum(total_area, values(dog.parts))
```

[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) uses the same graph for heat: each part
has its own heat budget, and heat is conducted between parts along the edges.
