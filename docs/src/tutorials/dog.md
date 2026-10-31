# A quadruped: the dog

This tutorial builds a 23 kg dog from seven parts: a torso of two halves with different fur, a head, and four
tapering legs. It uses everything in the manual: [layers](../manual/layers.md),
[frustums, truncation and halves](../manual/truncation.md), [surfaces](../manual/surfaces.md) and
[joins](../manual/joins.md).

```@setup dog
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## The parts

The torso is two half cylinders, with thick fur on the back and thin fur on the belly:

```@example dog
using BiophysicalGeometry, Unitful

density = 1000.0u"kg/m^3"
back_fur = FibrousLayer(15.0u"mm", 30.0u"μm", 3000u"cm^-2")
belly_fur = FibrousLayer(5.0u"mm", 30.0u"μm", 3000u"cm^-2")
limb_fur = FibrousLayer(8.0u"mm", 30.0u"μm", 3000u"cm^-2")

dorsal = Body(HalfCylinder(; mass = 10.0u"kg", density, axis_ratio_b = 3.0), back_fur)
ventral = Body(HalfCylinder(; mass = 10.0u"kg", density, axis_ratio_b = 3.0), belly_fur)
nothing # hide
```

Each leg is a frustum, narrowing to 0.4 of its width at the foot:

```@example dog
leg = Body(Cone(; mass = 0.5u"kg", density, axis_ratio_b = 5.0, top_ratio = 0.4), limb_fur)
map(x -> uconvert(u"cm", x), leg.geometry.length)
```

The head is an ellipsoid, truncated to leave a flat disc where it meets the torso. The disc is to have 0.7 of the
radius of the head:

```@example dog
head = Body(Ellipsoid(; mass = 2.0u"kg", density, axis_ratio_b = 1.5, axis_ratio_c = 1.5, pole_a_truncation = 1 - sqrt(1 - 0.7^2)), limb_fur)
shape_gallery("dorsal" => dorsal, "leg" => leg, "head" => head; decorations = false) # hide
```

## The joins

The joins need a few lengths from the parts:

```@example dog
torso_length = dorsal.geometry.length.length_skin
leg_radius = skin_radius(leg)
neck_radius = 0.7 * min(skin_radius(dorsal), skin_radius(head))
nothing # hide
```

The two halves are joined over their flat faces. The head goes on the end of the
dorsal half. The legs hang from the curved side of the ventral half, at 0.2 and 0.8 of its length, and 0.35
radians either side of its midline:

```@example dog
hip(along, around) = Attachment(Lateral(along * torso_length, π / 2 + around), Disc(leg_radius))
leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(leg_radius))

joins = (
    Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),
    Join(dorsal = Attachment(EndB(0.0u"m", 0.0), Disc(neck_radius)), head = Attachment(PoleA(), Disc(neck_radius))),
    Join(ventral = hip(0.8, 0.35), leg_fl = leg_top),
    Join(ventral = hip(0.8, -0.35), leg_fr = leg_top),
    Join(ventral = hip(0.2, 0.35), leg_bl = leg_top),
    Join(ventral = hip(0.2, -0.35), leg_br = leg_top),
)
nothing # hide
```

`EndA` of a cone is its wide end, so the legs narrow towards the ground.

## The dog

The torso lies along ``x`` with its back up, as every half does, so the parts are just put together:

```@example dog
dog = CompositeBody(;
    parts = (; dorsal, ventral, head, leg_fl = leg, leg_fr = leg, leg_bl = leg, leg_br = leg),
    joins,
)
composite_views(dog) # hide
```

```@example dog
body_graph(dog) # hide
```

## Its areas

```@example dog
map(x -> uconvert(u"m^2", x), (; total = total_area(dog), skin = skin_area(dog), evaporation = evaporation_area(dog)))
```

By part, as outer area less the area hidden by joins:

```@example dog
hidden = BiophysicalGeometry.covered_areas(dog.parts, dog.joins)
markdown_table(["Part", "Outer area", "Hidden by joins", "Exposed"],
               [(name, uconvert(u"cm^2", total_area(part)), uconvert(u"cm^2", hidden[name]),
                 uconvert(u"cm^2", total_area(part) - hidden[name])) for (name, part) in pairs(dog.parts)])
```

## In the sun

The silhouette from above, from the side and from the front:

```@example dog
fig = Figure(size = (760, 300))
for (i, (label, direction)) in enumerate(("Above" => (0.0, 0.0, 1.0), "Side" => (0.0, 1.0, 0.0), "Front" => (1.0, 0.0, 0.0)))
    ax = Axis(fig[1, i]; xlabel = "cm", ylabel = "cm")
    area = silhouette_panel!(ax, dog, direction)
    ax.title = "$label: $(round(Int, ustrip(u"cm^2", area))) cm²"
end
fig
```

Here `silhouette_panel!` is a helper of these docs that draws the image returned by
`silhouette_rasterized(dog, direction; return_image = true)`.

With the sun overhead the back takes nearly all of it, and the belly and most of the legs are in shade:

```@example dog
lit = silhouette(dog, Beam(0.0, 0.0, 1.0))
map(x -> round(Int, ustrip(u"cm^2", x)), lit)
```

Adding up the silhouettes of the parts, as if none shaded another, would overestimate the total:

```@example dog
uconvert(u"cm^2", sum(lit)), uconvert(u"cm^2", silhouette(dog).normal)
```

## Sky and ground

The fraction of each part's view that is sky, ground, or another part of the dog, over flat ground:

```@example dog
views = silhouette_factors(dog, Sky(0.5))
markdown_table(["Part", "Sky", "Ground", "Rest of the dog"],
               [(name, v.sky, v.ground, sum(v.neighbours)) for (name, v) in pairs(views)])
```

The back sees the sky, the belly sees the ground and the legs, and the legs see each other. These are the inputs to
the heat budget of each part in [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl).
