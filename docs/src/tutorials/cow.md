# Consider a non-spherical cow

The spherical cow is the physicist's joke about simple models. This tutorial takes it seriously, and then improves
on it: a 682 kg Holstein as a sphere, a cylinder, an ellipsoid, and a body of eight parts, with the surface area of
each compared to that measured on cows.

```@setup cow
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## The cow

The cow is the Holstein of Porter et al. (2023), whose model of the heat and water budgets of dairy cattle used a
three-dimensional image of a cow: 682 kg, 147.3 cm tall at the shoulder, with a coat 3.36 mm deep. The diameter and
density of the hairs are not given there and are estimates.

```@example cow
using BiophysicalGeometry, Unitful

mass = 682.0u"kg"
density = 1000.0u"kg/m^3"
coat = FibrousLayer(3.36u"mm", 60.0u"μm", 1000u"cm^-2")
nothing # hide
```

Brody (1945) measured the surface area of Holsteins from 41 to 617 kg, and fitted ``A = 0.14 M^{0.57}`` m²,
with ``M`` in kg:

```@example cow
measured(mass) = 0.14 * ustrip(u"kg", mass)^0.57 * u"m^2"
measured(mass)
```

## A spherical cow

```@example cow
sphere = Body(Sphere(; mass, density), coat)
uconvert(u"cm", 2 * insulation_radius(sphere)), total_area(sphere)
```

A ball just over a metre across, with two thirds of the surface of a cow. A sphere has the least area of any
shape, so it must fall short.

## A cylindrical cow, and an ellipsoidal one

A torso 2.4 times as long as it is wide is nearer the mark:

```@example cow
cylinder = Body(Cylinder(; mass, density, axis_ratio_b = 2.4), coat)
ellipsoid = Body(Ellipsoid(; mass, density, axis_ratio_b = 2.4, axis_ratio_c = 2.4), coat)
shape_gallery("Sphere" => sphere, "Cylinder" => cylinder, "Ellipsoid" => ellipsoid; decorations = false) # hide
```

```@example cow
total_area(cylinder), total_area(ellipsoid)
```

Still short. Much of the surface of a cow is on its legs, neck and head.

## A cow with parts

The body is given a torso of two halves, a neck and four legs that taper, and a head with a flat end to meet the
neck. The fractions of the mass in each part are estimates, as are the proportions of each part. The ratio of
length to width of the legs is set so that the cow stands 147.3 cm at the shoulder.

```@example cow
fractions = (torso = 0.83, neck = 0.05, head = 0.04, leg = 0.02)

function cow(mass; leg_ratio = 4.33)
    dorsal = Body(HalfCylinder(; mass = fractions.torso * mass / 2, density, axis_ratio_b = 2.4), coat)
    ventral = Body(HalfCylinder(; mass = fractions.torso * mass / 2, density, axis_ratio_b = 2.4), coat)
    neck = Body(Cone(; mass = fractions.neck * mass, density, axis_ratio_b = 1.0, top_ratio = 0.7), coat)
    head = Body(Ellipsoid(; mass = fractions.head * mass, density, axis_ratio_b = 1.8, axis_ratio_c = 1.8, pole_a_truncation = 1 - sqrt(1 - 0.7^2)), coat)
    leg = Body(Cone(; mass = fractions.leg * mass, density, axis_ratio_b = leg_ratio, top_ratio = 0.5), coat)

    torso_length = dorsal.geometry.length.length_skin
    torso_radius = skin_radius(dorsal)
    r_leg, r_neck, r_head = skin_radius(leg), 0.45 * torso_radius, 0.7 * skin_radius(head)
    hip(along, around) = Attachment(Lateral(along * torso_length, π / 2 + around), Disc(r_leg))
    leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))

    CompositeBody(;
        parts = (; dorsal, ventral, neck, head, leg_fl = leg, leg_fr = leg, leg_bl = leg, leg_br = leg),
        joins = (
            Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),
            Join(dorsal = Attachment(EndB(0.5 * torso_radius, π / 2), Disc(r_neck)),
                 neck = Attachment(EndA(0.0u"m", 0.0), Disc(r_neck))),
            Join(neck = Attachment(EndB(0.0u"m", 0.0), Disc(r_head)), head = Attachment(PoleA(), Disc(r_head))),
            Join(ventral = hip(0.88, 0.2), leg_fl = leg_top), Join(ventral = hip(0.88, -0.2), leg_fr = leg_top),
            Join(ventral = hip(0.12, 0.2), leg_bl = leg_top), Join(ventral = hip(0.12, -0.2), leg_br = leg_top),
        ),
    )
end

holstein = cow(mass)
composite_views(holstein; views = (:three_quarter, :side, :front, :top), titles = ["", "side", "front", "top"]) # hide
```

```@example cow
body_graph!(Axis(Figure(size = (640, 330))[1, 1]), holstein); current_figure() # hide
```

The neck leaves the end of the dorsal half, above the middle of the torso, and the head follows from the neck: a
chain of three parts. The height at the shoulder is from the top of the coat to the feet:

```@example cow
function shoulder_height(body)
    leg_length = body.parts.leg_fl.geometry.length.length_skin
    foot = BiophysicalGeometry.apply_pose(body.poses.leg_fl, (0.0u"m", 0.0u"m", leg_length))
    insulation_radius(body.parts.dorsal) - foot[3]
end
uconvert(u"cm", shoulder_height(holstein))
```

## How close?

```@example cow
meeh(area) = ustrip(u"m^2", area) / ustrip(u"kg", mass)^(2 / 3)
models = ("Spherical cow" => sphere, "Ellipsoidal cow" => ellipsoid, "Cylindrical cow" => cylinder,
          "Cow with parts" => holstein)
markdown_table(["Model", "Surface area", "Relative to measured", "Meeh coefficient"],
               vcat([(name, total_area(body), total_area(body) / measured(mass), meeh(total_area(body)))
                     for (name, body) in models],
                    [("Measured (Brody 1945)", measured(mass), 1.0, meeh(measured(mass)))]))
```

The cow with parts has more surface than the measured cow, where the single shapes had less. Its legs are round
and smooth and its torso has square ends. The area by part shows where the surface is:

```@example cow
hidden = BiophysicalGeometry.covered_areas(holstein.parts, holstein.joins)
exposed = map((part, h) -> total_area(part) - h, holstein.parts, hidden)
(; torso = exposed.dorsal + exposed.ventral, neck_and_head = exposed.neck + exposed.head,
   legs = 4 * exposed.leg_fl)
```

## Other sizes

Porter et al. (2023) also modelled cows of 364 kg, the size of a Jersey, and of 818 and 1273 kg. The cow with
parts keeps its shape as it grows, so its area rises with mass to the power 2/3, faster than in real cattle:

```@example cow
markdown_table(["Mass", "Cow with parts", "Measured", "Ratio"],
               [(m, total_area(cow(m)), measured(m), total_area(cow(m)) / measured(m))
                for m in [364.0, 682.0, 818.0, 1273.0] .* u"kg"])
```

Larger cattle are stockier. That change of shape could be added by making the proportions depend on mass, as in
[Geometry and allometry](scaling.md).

## What the sphere cannot say

Surface area is only one of the things a heat budget needs, and the others depend on shape still more. The
silhouette of the cow with the sun overhead, to its side and in front:

```@example cow
fig = Figure(size = (760, 300))
for (i, (label, direction)) in enumerate(("Above" => (0.0, 0.0, 1.0), "Side" => (0.0, 1.0, 0.0), "Front" => (1.0, 0.0, 0.0)))
    ax = Axis(fig[1, i]; xlabel = "cm", ylabel = "cm")
    area = silhouette_panel!(ax, holstein, direction)
    ax.title = "$label: $(round(ustrip(u"m^2", area); digits = 2)) m²"
end
fig
```

A sphere has the same silhouette from every side:

```@example cow
silhouette(sphere).normal
```

A cow can more than halve the sunlight it intercepts by turning to face the sun. And with a back and a belly, its coat and
its surroundings can differ above and below:

```@example cow
views = silhouette_factors(holstein, Sky(0.5))
(; dorsal_sky = views.dorsal.sky, ventral_sky = views.ventral.sky, ventral_legs = sum(views.ventral.neighbours))
```

## Two routes to a cow

Porter et al. (2023) took their areas and heat-transfer coefficients from a detailed three-dimensional image of a
cow, analysed with computational fluid dynamics, following the "non-spherical elephant" of Dudley et al. (2013).
That is the most faithful geometry, for one animal in one posture. A body built from simple shapes is cruder, and
can be reshaped, re-coated and re-posed in a line of code, which is what is needed to ask how size, shape, coat
and behaviour matter.
