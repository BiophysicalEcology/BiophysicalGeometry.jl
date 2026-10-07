# Geometry and allometry

How does surface area change with body mass? An allometric equation answers from measurements, and geometry
answers from shape. This tutorial compares the surface areas from this package with the scaling equations of
[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl), for mammals and for birds.

```@setup scaling
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## Two answers

BiologicalScaling.jl gives the surface area of a mammal as ``0.1 M^{0.667}`` m², with ``M`` in kg
(Schmidt-Nielsen 1975), and the area of its skin as ``1110 M^{0.65}`` cm² (Walsberg and King 1978). Both packages
have functions called `surface_area` and `skin_area`, so BiologicalScaling.jl is loaded with `import`:

```@example scaling
using BiophysicalGeometry, Unitful
import BiologicalScaling
const BS = BiologicalScaling

mass = 20.0u"kg"
density = 1000.0u"kg/m^3"
BS.surface_area(BS.EutherianMammal(), mass), total_area(Body(Cylinder(; mass, density, axis_ratio_b = 3.0), Naked()))
```

## Simple shapes

A shape that keeps its proportions as it grows is *isometric*: its area increases with mass to the power 2/3. On
logarithmic axes every such shape is a straight line with the same slope, and shapes differ only in their height:

```@example scaling
masses = exp10.(range(-2, 3; length = 40)) .* u"kg"
area(shape_of, masses) = [ustrip(u"m^2", total_area(Body(shape_of(m), Naked()))) for m in masses]
shapes = ("Sphere" => m -> Sphere(; mass = m, density),
          "Ellipsoid, 3 : 1" => m -> Ellipsoid(; mass = m, density, axis_ratio_b = 3.0, axis_ratio_c = 3.0),
          "Cylinder, 3 : 1" => m -> Cylinder(; mass = m, density, axis_ratio_b = 3.0))

fig, ax = figure_axis("Body mass (kg)", "Surface area (m²)"; xscale = log10, yscale = log10)
for (label, shape_of) in shapes
    lines!(ax, ustrip.(u"kg", masses), area(shape_of, masses); linewidth = 2, label)
end
lines!(ax, ustrip.(u"kg", masses), ustrip.(u"m^2", BS.surface_area.(Ref(BS.EutherianMammal()), masses));
       linewidth = 2, linestyle = :dash, color = :black, label = "Mammal allometry")
axislegend(ax; position = :lt)
fig
```

## The Meeh coefficient

The height of each line is the *Meeh coefficient* ``k`` in ``A = k M^{2/3}``, a measure of shape. With area in
m² and mass in kg it is 0.1 in the mammal equation:

```@example scaling
meeh(body_area, mass) = ustrip(u"m^2", body_area) / ustrip(u"kg", mass)^(2 / 3)
markdown_table(["Shape", "Meeh coefficient"],
               [(label, meeh(total_area(Body(shape_of(mass), Naked())), mass)) for (label, shape_of) in shapes])
```

A sphere has the least area that a body can have, about half that of a real mammal. No single simple shape reaches
0.1 without being stretched a long way: mammals have legs.

## A body with legs

A quadruped of any mass, with fixed proportions: a torso with 80% of the mass, a head with 8% and four legs with
3% each.

```@example scaling
function quadruped(mass; leg_mass = 0.03mass, leg_ratio = 4.0)
    torso = Body(Cylinder(; mass = 0.92mass - 4leg_mass, density, axis_ratio_b = 2.5), Naked())
    head = Body(Ellipsoid(; mass = 0.08mass, density, axis_ratio_b = 1.5, axis_ratio_c = 1.5), Naked())
    leg = Body(Cylinder(; mass = leg_mass, density, axis_ratio_b = leg_ratio), Naked())
    L = torso.geometry.length.length_skin
    r_leg = skin_radius(leg)
    r_neck = 0.5 * skin_radius(head)
    hip(z, φ) = Attachment(Lateral(z * L, -π / 2 + φ), Disc(r_leg))
    leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))
    CompositeBody(;
        parts = (; torso, head, leg_fl = leg, leg_fr = leg, leg_bl = leg, leg_br = leg),
        joins = (Join(torso = Attachment(EndB(0.0u"m", 0.0), Disc(r_neck)), head = Attachment(PoleB(), Disc(r_neck))),
                 Join(torso = hip(0.85, 0.4), leg_fl = leg_top), Join(torso = hip(0.85, -0.4), leg_fr = leg_top),
                 Join(torso = hip(0.15, 0.4), leg_bl = leg_top), Join(torso = hip(0.15, -0.4), leg_br = leg_top)),
    )
end
composite_views(quadruped(20.0u"kg"); views = (:oblique, :side, :front), legend = false) # hide
```

```@example scaling
meeh(total_area(quadruped(20.0u"kg")), 20.0u"kg")
```

Legs bring the coefficient most of the way to the 0.1 of the allometry, and it is the same at every mass, because
the proportions are.

## Proportions that change with size

Real proportions are not fixed. Under *elastic similarity* (McMahon 1973) the legs of larger animals are
relatively thicker, so that they do not buckle. BiologicalScaling.jl gives the length and diameter of a leg from
the mass of the body:

```@example scaling
markdown_table(["Body mass", "Leg length", "Leg diameter", "Length / diameter"],
               [(m, uconvert(u"cm", BS.limb_length(BS.ElasticSimilarity(), m)),
                 uconvert(u"cm", BS.limb_diameter(BS.ElasticSimilarity(), m)), BS.limb_aspect_ratio(BS.ElasticSimilarity(), m))
                for m in [0.02, 2.0, 200.0] .* u"kg"])
```

A cylinder is sized from its mass and its length over its diameter, so the legs are given the mass of a cylinder
of that length and diameter:

```@example scaling
function elastic(mass)
    leg_length = BS.limb_length(BS.ElasticSimilarity(), mass)
    leg_diameter = BS.limb_diameter(BS.ElasticSimilarity(), mass)
    leg_mass = density * π * (leg_diameter / 2)^2 * leg_length
    quadruped(mass; leg_mass, leg_ratio = leg_length / leg_diameter)
end
fig = Figure(size = (700, 260)) # hide
for (i, m) in enumerate([0.02, 2.0, 200.0] .* u"kg") # hide
    draw_parts!(body_axis(fig[1, i]; decorations = false, azimuth = -π / 2, elevation = 0.0, title = "$m", # hide
        titlesize = 12), elastic(m)) # hide
end # hide
fig # hide
```

The area of this animal no longer scales with an exponent of exactly 2/3:

```@example scaling
elastic_areas = [ustrip(u"m^2", total_area(elastic(m))) for m in masses]
fixed_areas = [ustrip(u"m^2", total_area(quadruped(m))) for m in masses]
slope(areas) = (log(areas[end]) - log(areas[1])) / (log(masses[end] / masses[1]))
(; fixed_proportions = slope(fixed_areas), elastic_similarity = slope(elastic_areas))
```

```@example scaling
fig, ax = figure_axis("Body mass (kg)", "Meeh coefficient (m² kg⁻²ᐟ³)"; xscale = log10)
kg = ustrip.(u"kg", masses)
lines!(ax, kg, fixed_areas ./ kg .^ (2 / 3); linewidth = 2, label = "Quadruped, fixed proportions")
lines!(ax, kg, elastic_areas ./ kg .^ (2 / 3); linewidth = 2, label = "Quadruped, elastic similarity")
lines!(ax, kg, ustrip.(u"m^2", BS.surface_area.(Ref(BS.EutherianMammal()), masses)) ./ kg .^ (2 / 3);
       linewidth = 2, linestyle = :dash, color = :black, label = "Mammal allometry")
lines!(ax, kg, ustrip.(u"m^2", BS.skin_area.(Ref(BS.EutherianMammal()), masses)) ./ kg .^ (2 / 3);
       linewidth = 2, linestyle = :dot, color = :black, label = "Mammal skin allometry")
axislegend(ax; position = :rt)
fig
```

## Birds: less plumage than skin

In birds the outer surface of the plumage is *smaller* than the skin under it. Walsberg and King (1978) found a
skin area of ``10.0 M^{0.667}`` cm² and a plumage area of ``8.11 M^{0.667}`` cm², with ``M`` in g:

```@example scaling
BS.skin_area(BS.PasserineBird(), 20.0u"g"), BS.plumage_area(BS.PasserineBird(), 20.0u"g")
```

A layer of fibres on a single shape cannot do this, as the outside of a layer is always larger than the skin it
covers. The reason is in the shape. The skin of a bird follows its neck, its folded wings and its legs, and the
feathers smooth over all of them into one rounded outline. So the skin is a body of parts:

```@example scaling
function bird_skin(mass)
    mass = uconvert(u"kg", mass)
    body = Body(Ellipsoid(; mass = 0.70mass, density, axis_ratio_b = 1.6, axis_ratio_c = 1.6), Naked())
    head = Body(Sphere(; mass = 0.08mass, density), Naked())
    neck = Body(Cylinder(; mass = 0.04mass, density, axis_ratio_b = 2.0), Naked())
    wing = Body(Plate(; mass = 0.06mass, density, axis_ratio_b = 2.5, axis_ratio_c = 12.0), Naked())    # folded against the body
    leg = Body(Cylinder(; mass = 0.03mass, density, axis_ratio_b = 8.0), Naked())
    r_neck, r_leg = skin_radius(neck), skin_radius(leg)
    shoulder = Disc(0.3 * wing.geometry.length.width_skin)
    under_wing = Attachment(Bottom(0.0u"m", 0.0u"m"), shoulder)
    leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))
    CompositeBody(;
        parts = (; body, neck, head, wing_left = wing, wing_right = wing, leg_left = leg, leg_right = leg),
        joins = (
            Join(body = Attachment(PoleA(), Disc(r_neck)), neck = Attachment(EndA(0.0u"m", 0.0), Disc(r_neck))),
            Join(neck = Attachment(EndB(0.0u"m", 0.0), Disc(r_neck)), head = Attachment(Radial(π / 2, π), Disc(r_neck))),
            Join(body = Attachment(Equator(0.0), shoulder), wing_left = under_wing),
            Join(body = Attachment(Equator(π), shoulder), wing_right = under_wing),
            Join(body = Attachment(Equator(-π / 2 + 0.4), Disc(r_leg)), leg_left = leg_top),
            Join(body = Attachment(Equator(-π / 2 - 0.4), Disc(r_leg)), leg_right = leg_top),
        ),
    )
end
nothing # hide
```

and the plumage is a single ellipsoid of the whole mass, with feathers 4 mm deep on a 20 g bird and in proportion
on others:

```@example scaling
function bird_plumage(mass)
    depth = 4.0u"mm" * (mass / 20.0u"g")^(1 / 3)
    Body(Ellipsoid(; mass = uconvert(u"kg", mass), density, axis_ratio_b = 1.6, axis_ratio_c = 1.6), FibrousLayer(depth, 30.0u"μm", 3000u"cm^-2"))
end
fig = Figure(size = (640, 260)) # hide
draw_parts!(body_axis(fig[1, 1]; decorations = false, azimuth = -0.3π, title = "Skin", titlesize = 12), # hide
    bird_skin(20.0u"g")) # hide
draw_layers!(body_axis(fig[1, 2]; decorations = false, azimuth = -0.3π, title = "Plumage", titlesize = 12), # hide
    bird_plumage(20.0u"g")) # hide
fig # hide
```

```@example scaling
markdown_table(["Mass", "Skin, parts", "Skin, allometry", "Plumage, ellipsoid", "Plumage, allometry"],
               [(m, uconvert(u"cm^2", skin_area(bird_skin(m))), BS.skin_area(BS.PasserineBird(), m),
                 uconvert(u"cm^2", total_area(bird_plumage(m))), BS.plumage_area(BS.PasserineBird(), m))
                for m in [20.0, 200.0, 2000.0] .* u"g"])
```

The proportions here are illustrative, and with them geometry gives what was measured: more skin than plumage.
For a heat budget the two areas do different jobs. Heat leaves the bird from the plumage, and water evaporates
from the skin.

## Which to use

An allometric equation gives a surface area with nothing more than a mass, and is the place to start. Geometry is
needed when the question is about shape: the area of a part, the effect of a posture or a coat, what the sun and
the sky see. The two also work together, as here, with allometry supplying the proportions that geometry turns
into a body. The page [Build an animal](../builder.md) offers legs sized by elastic similarity in the same way.
