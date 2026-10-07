# Frustums, truncation and halves

Limbs taper, heads sit flush against necks, and backs differ from bellies. Three variations on the basic shapes
allow for this: a cone cut short, an ellipsoid with one end sliced off, and a shape cut in half lengthwise.

```@setup truncation
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## Cones and frustums

A *frustum* is a cone with its tip cut off parallel to the base. [`Cone`](@ref) takes a `top_ratio`, the radius of
the top over the radius of the base: 0 is a sharp cone, and values up to 1 are frustums that approach a cylinder.

```@example truncation
using BiophysicalGeometry, Unitful

m, density = 1.0u"kg", 1000.0u"kg/m^3"
shape_gallery(("top_ratio = $t" => Body(Cone(; mass = m, density, axis_ratio_b = 2.5, top_ratio = t), Naked()) for t in (0.0, 0.3, 0.6, 1.0))...; # hide
    ncols = 4) # hide
```

The mass is the same in each, so a sharper cone is wider and longer:

```@example truncation
leg = Body(Cone(; mass = m, density, axis_ratio_b = 2.5, top_ratio = 0.4), Naked())
leg.geometry.length, total_area(leg)
```

`radius_skin` is the radius of the base, and the top has radius `top_ratio * radius_skin`.

## Truncated ellipsoids

[`Ellipsoid`](@ref) takes a keyword `pole_a_truncation`: the fraction of the long semi-axis sliced off the ``+x``
end, leaving a flat disc. At 0 the ellipsoid is whole, and at 1 it is cut through the middle.

```@example truncation
fig = Figure(size = (840, 230)) # hide
for (i, t) in enumerate((0.0, 0.15, 0.4, 0.8)) # hide
    ax = body_axis(fig[1, i]; decorations = false, azimuth = -0.25π, title = "pole_a_truncation = $t", titlesize = 12) # hide
    draw_parts!(ax, Body(Ellipsoid(; mass = m, density, axis_ratio_b = 2.0, axis_ratio_c = 2.0, pole_a_truncation = t), Naked())) # hide
end # hide
fig # hide
```

The disc is where another part is joined, such as a head onto a neck. Its radius follows from the truncation:
for a disc of radius ``f`` times the short semi-axis, the truncation is ``1 - \sqrt{1 - f^2}``.

```@example truncation
fraction = 0.7
head = Body(Ellipsoid(; mass = m, density, axis_ratio_b = 2.0, axis_ratio_c = 2.0, pole_a_truncation = 1 - sqrt(1 - fraction^2)), Naked())
surface_area(head.shape, head, PoleA()), π * (fraction * skin_radius(head))^2
```

The cut body keeps its mass: the ellipsoid is made larger until what is left after the cut holds the volume, and
its area is that of the cut surface — the whole ellipsoid less the cap cut away, plus the disc. `length`, `width`
and `height` describe the ellipsoid before the cut. Its silhouette is that of the cut body too.

```@example truncation
cut = [Body(Ellipsoid(; mass = m, density, axis_ratio_b = 2.0, axis_ratio_c = 2.0, pole_a_truncation = t), Naked())
       for t in (0.0, 0.4, 1.0)]
markdown_table(["pole_a_truncation", "volume", "length before the cut", "total area"],
    [(t, b.geometry.volume, b.geometry.length.length_skin, total_area(b)) for (t, b) in zip((0.0, 0.4, 1.0), cut)])
```

## Halves

[`HalfCylinder`](@ref), [`HalfEllipsoid`](@ref) and [`HalfSphere`](@ref) are a shape cut lengthwise through its
centre, with a flat face where it was cut. Two halves can then have different layers, such as thick fur on the back
and thin fur on the belly.

```@example truncation
shape_gallery("HalfCylinder" => Body(HalfCylinder(; mass = m, density, axis_ratio_b = 2.0), Naked()), # hide
    "HalfEllipsoid" => Body(HalfEllipsoid(; mass = m, density, axis_ratio_b = 2.0, axis_ratio_c = 4.0), Naked()), # hide
    "HalfSphere" => Body(HalfSphere(; mass = m, density), Naked())) # hide
```

The keywords describe the half: its mass and volume are its own, and so is the `height` of a half-ellipsoid, the
height of its dome. A half has the dimensions of the whole shape of twice its mass, half its volume, and half its
area plus the flat face:

```@example truncation
whole = Body(Cylinder(; mass = 2m, density, axis_ratio_b = 2.0), Naked())
half = Body(HalfCylinder(; mass = m, density, axis_ratio_b = 2.0), Naked())
flat_face = surface_area(half.shape, half, Flat())
total_area(half), total_area(whole) / 2 + flat_face
```

Each of these is a [`Half`](@ref) wrapped around the whole shape, so use [`mass`](@ref) to get its mass:

```@example truncation
mass(half.shape), mass(half.shape.parent)
```

Two halves joined on their flat faces make the whole shape again, see [Areas and volumes](areas.md).

A half plate is just a thinner plate, so there is no half plate: make the [`Plate`](@ref) with half the height. For
wings and ears, a [`TriangularPlate`](@ref) is a plate cut along its diagonal instead.
