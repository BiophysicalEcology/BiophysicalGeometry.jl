# Back and belly: two halves

The back of an animal faces the sun and the sky, and its belly faces the ground, usually with thinner fur. This
tutorial splits a body into a dorsal and a ventral half, the commonest body of more than one part, and then adds a
head.

```@setup two
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## Why split the body

A single shape has one coat and one set of surroundings. The way around this in NicheMapR is to solve the heat
budget twice, for an animal that is all back under the sky and one that is all belly over the ground, and to
average the two answers. With two halves joined into one body, each half has its own coat and its own view, and
heat can pass between them, as solved by [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl).

## Two halves

A [`HalfEllipsoid`](@ref) is an ellipsoid cut lengthwise. The mass is that of the half:

```@example two
using BiophysicalGeometry, Unitful

mass = 500.0u"g"
density = 1000.0u"kg/m^3"
back_fur = FibrousLayer(10.0u"mm", 20.0u"μm", 8000u"cm^-2")
belly_fur = FibrousLayer(3.0u"mm", 20.0u"μm", 8000u"cm^-2")

dorsal = Body(HalfEllipsoid(mass / 2, density, 3.0, 1.0), back_fur)
ventral = Body(HalfEllipsoid(mass / 2, density, 3.0, 1.0), belly_fur)
nothing # hide
```

They are joined over the whole of their flat faces:

```@example two
animal = CompositeBody(;
    parts = (; dorsal, ventral),
    joins = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),),
)
composite_views(animal; views = (:oblique, :side, :front)) # hide
```

`parts` names the two bodies, and the [`Join`](@ref) says where each is attached: at its [`Flat`](@ref) face, with
the whole face covered. The first part listed, `dorsal`, stays where it is and `ventral` is placed against it.

## Areas

The flat faces are hidden, so the area of the body is that of the two domes:

```@example two
map(x -> uconvert(u"cm^2", x), (; whole = total_area(animal), dorsal = total_area(dorsal), ventral = total_area(ventral),
                                 join = join_area(animal.joins[1], animal)))
```

## Sun, sky and ground

With the sun overhead, only the back is lit:

```@example two
map(x -> uconvert(u"cm^2", x), silhouette(animal, Beam(0.0, 0.0, 1.0)))
```

and with the sun low, both halves are:

```@example two
map(x -> uconvert(u"cm^2", x), silhouette(animal, Beam(1.0, 0.0, 0.3)))
```

Over flat ground the back sees mostly sky and the belly mostly ground:

```@example two
views = silhouette_factors(animal, Sky(0.5))
(; dorsal = views.dorsal, ventral = views.ventral)
```

The rest of the view of each half is the other half, across the join. The split between sky and ground is not
complete, as the sides of each half see past the horizontal. An average of an all-back and an all-belly animal
assumes that it is.

## Adding a head

Other parts are joined to either half. Here a sphere is joined to the front of the dorsal half over a small disc.
A point on the dome is given by two angles, from the long axis and around it:

```@example two
head = Body(Sphere(60.0u"g", density), back_fur)
patch = Disc(1.0u"cm")
with_head = CompositeBody(;
    parts = (; dorsal, ventral, head),
    joins = (
        Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),
        Join(dorsal = Attachment(Dome(0.25, π / 2), patch), head = Attachment(Radial(π / 2, 0.0), patch)),
    ),
)
body_graph(with_head; size = (760, 260)) # hide
```

The body is now a graph of three parts and two joins, see [Bodies as graphs](../manual/graphs.md). The head takes
some of the sun from the back:

```@example two
map(x -> uconvert(u"cm^2", x), silhouette(with_head, Beam(1.0, 0.0, 1.0)))
```

## Half cylinders

A [`HalfCylinder`](@ref) is split in the same way. A cylinder stands along ``z`` with its flat face to the side,
so the body is laid down with a `root_pose`, and the ventral half is turned to lie along the dorsal one with a
`twist`, see [Joins and poses](../manual/joins.md):

```@example two
lying = Pose((0.0u"m", 0.0u"m", 0.0u"m"), [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0])
torso = CompositeBody(;
    parts = (; dorsal = Body(HalfCylinder(mass / 2, density, 3.0), back_fur),
               ventral = Body(HalfCylinder(mass / 2, density, 3.0), belly_fur)),
    joins = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover()); twist = -π / 2),),
    root_pose = lying,
)
composite_views(torso; views = (:oblique, :side, :front)) # hide
```

The [next tutorial](dog.md) adds legs and a head to a torso like this one.
