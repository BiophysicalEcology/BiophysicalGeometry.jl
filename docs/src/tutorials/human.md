# A human: comparison with NicheMapR

The human model of NicheMapR, `HomoTherm`, is made of a head, a trunk, two arms and two legs. This tutorial builds
the same 70 kg person with this package, and compares the geometry with that from NicheMapR 3.3.3, part by part
and as a whole.

```@setup human
using Main.FigureHelpers
using CairoMakie
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Top, Bottom
```

## The parts

`HomoTherm` gives each part a fraction of the body mass and a ratio of length to width. The same numbers are in
[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl):

```@example human
using BiophysicalGeometry, Unitful
import BiologicalScaling

proportions = BiologicalScaling.body_part_proportions(BiologicalScaling.Human())
part_names = (:head, :trunk, :arm, :leg)
part_mass = NamedTuple{part_names}(Tuple(proportions.mass_fraction .* 70.0u"kg"))
ratio = NamedTuple{part_names}(Tuple(proportions.aspect_ratio))
```

The head is an ellipsoid and the other parts are cylinders, all with a density of 1050 kg/m³. Each has its own
fraction of fat, and 6 mm of clothing, treated as a layer of fibres. The head has 10 mm of hair above and none
below, and `HomoTherm` uses the mean of 5 mm for its dimensions:

```@example human
density = 1050.0u"kg/m^3"
fat = (head = 0.035, trunk = 0.252, arm = 0.07, leg = 0.161)
layers(part, depth) = CompositeInsulation(FibrousLayer(depth, 1.0u"μm", 3.0e8u"m^-2"), FatLayer(fat[part], density))

head = Body(Ellipsoid(part_mass.head, density, ratio.head, 1.0), layers(:head, 5.0u"mm"))
trunk = Body(Cylinder(part_mass.trunk, density, ratio.trunk), layers(:trunk, 6.0u"mm"))
arm = Body(Cylinder(part_mass.arm, density, ratio.arm), layers(:arm, 6.0u"mm"))
leg = Body(Cylinder(part_mass.leg, density, ratio.leg), layers(:leg, 6.0u"mm"))
parts = (; head, trunk, arm, leg)
nothing # hide
```

## Part by part

The reference values were written by the script `docs/src/data/homotherm_reference.R`, which calls `GEOM_ENDO`,
`HomoTherm` and `human_silhouette_area` with their defaults:

```@example human
using DelimitedFiles

data(file) = readdlm(joinpath(pkgdir(BiophysicalGeometry), "docs", "src", "data", file), ','; header = true)
table, header = data("homotherm_parts.csv")
column(name) = Float64.(table[:, findfirst(==(name), vec(header))])
nothing # hide
```

Every dimension and area of every part agrees:

```@example human
compare(label, name, f, unit) = [(label, part, uconvert(unit, f(parts[part])), uconvert(unit, column(name)[i] * u"m" ^ (unit == u"cm^2" ? 2 : 1)))
                                 for (i, part) in enumerate(part_names)]
markdown_table(["Quantity", "Part", "BiophysicalGeometry.jl", "NicheMapR"],
               vcat(compare("Skin radius", "skin_radius_m", skin_radius, u"mm"),
                    compare("Fat thickness", "fat_thickness_m", b -> b.geometry.length.fat, u"mm"),
                    compare("Skin area", "skin_area_m2", skin_area, u"cm^2"),
                    compare("Outer area", "total_area_m2", total_area, u"cm^2")))
```

This is expected, as the calculations for a single shape came from NicheMapR, see
[Origins](../manual/introduction.md#Origins).

## Putting them together

NicheMapR joins the parts too: each gives up a fixed fraction of its area to its joins, the `PJOINs` argument,
heat is conducted across them, and the area of the body is the sum of what is left. What its parts do not have is
a place in space. Here they are given one: the head on top of the trunk, the legs under it, and the arms hung from
its sides by a join between two curved surfaces, turned with a `twist` to point down.

```@example human
trunk_length = trunk.geometry.length.length_skin
r_arm, r_leg = skin_radius(arm), skin_radius(leg)
shoulder(side) = Attachment(Lateral(trunk_length - r_arm, side), Disc(r_arm))
hip(side) = Attachment(EndA(1.1r_leg, side), Disc(r_leg))
leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))

human = CompositeBody(;
    parts = (; trunk, head, arm_left = arm, arm_right = arm, leg_left = leg, leg_right = leg),
    joins = (
        Join(trunk = Attachment(EndB(0.0u"m", 0.0), Disc(r_arm)), head = Attachment(PoleB(), Disc(r_arm))),
        Join(trunk = shoulder(0.0), arm_left = Attachment(Lateral(r_arm, π), Disc(r_arm)); twist = π),
        Join(trunk = shoulder(π), arm_right = Attachment(Lateral(r_arm, 0.0), Disc(r_arm)); twist = π),
        Join(trunk = hip(0.0), leg_left = leg_top),
        Join(trunk = hip(π), leg_right = leg_top),
    ),
)
composite_views(human; views = (:oblique, :side, :front), titles = ["", "front", "side"], size = (700, 340)) # hide
```

```@example human
body_graph!(Axis(Figure(size = (520, 300))[1, 1]), human); current_figure() # hide
```

A trunk is upright by default, so no `root_pose` is needed. The patch under the head has the radius of an arm, a
neck of sorts.

## The whole body

```@example human
totals, _ = data("homotherm_totals.csv")
reference = Dict(totals[i, 1] => Float64(totals[i, 2]) for i in axes(totals, 1))

height = leg.geometry.length.length_skin + trunk_length + 2 * head.geometry.length.a_semi_major_skin
joined = sum(2 * join_area(join, human) for join in human.joins)
markdown_table(["Quantity", "BiophysicalGeometry.jl", "NicheMapR"], [
    ("Height", height, reference["height_m"] * u"m"),
    ("Outer area of all parts", sum(total_area, values(human.parts)), reference["area_gross_m2"] * u"m^2"),
    ("Area hidden by joins", joined, reference["join_area_m2"] * u"m^2"),
    ("Outer area of the body", total_area(human), reference["area_net_m2"] * u"m^2"),
])
```

The parts add up to the same area. The difference is in the joins. The fractions that NicheMapR takes off each
part were chosen so that the shares of area in the head, trunk, arms and legs match those measured on people. Here the hidden area follows from the patches of the joins that were made, and is about
a fifth less. NicheMapR's height also includes 6 mm of shoe.

For comparison, the equation of DuBois and DuBois (1916) for the skin area of a person of this mass and height
gives

```@example human
dubois = 0.00718 * 70.0^0.425 * ustrip(u"cm", height)^0.725 * u"m^2"
dubois, skin_area(human)
```

## In the sun

NicheMapR finds the silhouette of each part with the sun at a given zenith angle, and adds them up. The same sum
here is `silhouette(human, angle)`, and it agrees except for the hair, which NicheMapR takes at its full 10 mm
depth for this purpose:

```@example human
sun, _ = data("homotherm_silhouette.csv")
zenith = Float64.(sun[:, 1])
parts_sum = [ustrip(u"m^2", silhouette(human, z * u"°")) for z in zenith]
maximum(abs.(parts_sum .- sun[:, 2]))
```

Adding the parts ignores the shade of one part on another: an arm on the trunk with the sun to the side, and the
head and shoulders on everything else with the sun overhead. NicheMapR corrects for this with a fit to the
silhouettes of people measured by Underwood and Ward (1966). Here the shade is computed, for a person facing the
sun or side-on to it:

```@example human
direction(z, azimuth) = (sind(z) * cosd(azimuth), sind(z) * sind(azimuth), cosd(z))
facing = [ustrip(u"m^2", silhouette_rasterized(human, direction(z, 90.0))) for z in zenith]
side_on = [ustrip(u"m^2", silhouette_rasterized(human, direction(z, 0.0))) for z in zenith]

fig, ax = figure_axis("Zenith angle of the sun (°)", "Silhouette area (m²)")
lines!(ax, zenith, Float64.(sun[:, 2]); linewidth = 2, color = :grey60, label = "NicheMapR, sum of parts")
scatter!(ax, zenith, parts_sum; color = :grey30, markersize = 7, label = "Sum of parts")
lines!(ax, zenith, facing; linewidth = 2, label = "With shading, facing the sun")
lines!(ax, zenith, side_on; linewidth = 2, label = "With shading, side-on")
lines!(ax, zenith, Float64.(sun[:, 3]); linewidth = 2, linestyle = :dash, color = :black,
       label = "Underwood and Ward (1966), facing the sun")
axislegend(ax; position = :lt, labelsize = 11)
fig
```

```@example human
fig = Figure(size = (620, 330)) # hide
for (i, (label, d)) in enumerate(("Facing, sun at 60°" => direction(60.0, 90.0), "Side-on, sun at 60°" => direction(60.0, 0.0), # hide
                                  "Sun overhead" => direction(0.0, 0.0))) # hide
    ax = Axis(fig[1, i]; title = label, titlesize = 12) # hide
    silhouette_panel!(ax, human, d) # hide
    hidedecorations!(ax); hidespines!(ax) # hide
end # hide
fig # hide
```

With the shading computed, a person side-on to the sun is distinguished from one facing it, which a sum of parts
cannot do. The measured silhouettes are smaller still than the computed ones facing the sun, which says that six
round parts are wider than a person: a real trunk is flatter than a cylinder.

## Sky and ground

`HomoTherm` is given the fraction of each part's view that is sky and ground as arguments. Here they are computed:

```@example human
views = silhouette_factors(human, Sky(0.5))
markdown_table(["Part", "Sky", "Ground", "Rest of the body"],
               [(name, v.sky, v.ground, sum(v.neighbours)) for (name, v) in pairs(views)])
```

The defaults of `HomoTherm` for the view of the sky are 0.50 for the head, 0.42 for the trunk and 0.35 for the arms
and legs, close to those computed.
