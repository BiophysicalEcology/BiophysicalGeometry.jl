# Plotting

The plotting functions are available once a [Makie](https://docs.makie.org) backend is loaded: CairoMakie for
figures in files and documents, GLMakie for windows that can be rotated.

```@setup plotting
using Main.FigureHelpers
```

```@example plotting
using BiophysicalGeometry, Unitful
using CairoMakie
import BiophysicalGeometry: Sphere, Top, Bottom
```

The `import` is needed because Makie exports `Sphere`, `Top` and `Bottom` too.

## A body

[`plot_body`](@ref) draws a body with a quarter of its fat and fibres cut away:

```@example plotting
fur = FibrousLayer(20.0u"mm", 30.0u"μm", 3000u"cm^-2")
fat = FatLayer(0.2, 901.0u"kg/m^3")
body = Body(Ellipsoid(; mass = 2.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 2.0, axis_ratio_c = 2.0), CompositeInsulation(fur, fat))
plot_body(body)
```

[`plot_cross_sections`](@ref) cuts it along and across, for every shape:

```@example plotting
plot_cross_sections(body)
```

## Several parts

A [`CompositeBody`](@ref) is drawn with each part in its place. Parts are not cut away.

```@example plotting
dorsal = Body(HalfEllipsoid(; mass = 1.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 3.0, axis_ratio_c = 6.0), fur)
ventral = Body(HalfEllipsoid(; mass = 1.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 3.0, axis_ratio_c = 6.0), Naked())
animal = CompositeBody(; parts = (; dorsal, ventral), joins = (
    Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),))
plot_body(animal)
```

[`plot_body_silhouette`](@ref) shows the body beside its silhouette, with sliders for the zenith and azimuth
of the sun. The sliders work in GLMakie and WGLMakie.

## Your own figures

[`draw_cutaway!`](@ref) and [`draw_cross_sections!`](@ref) draw into axes of your own, to arrange several bodies
in one figure. `sc` is the scale from metres to the units of the axes, 100 for cm.

```@example plotting
fig = Figure(size = (700, 320))
for (i, shape) in enumerate((Cylinder(; mass = 2.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 2.0), Sphere(; mass = 2.0u"kg", density = 1000.0u"kg/m^3")))
    ax = Axis3(fig[1, i]; aspect = :data, title = string(nameof(typeof(shape))))
    draw_cutaway!(ax, Body(shape, CompositeInsulation(fur, fat)))
end
fig
```

## Fibres

[`plot_insulation_properties`](@ref) draws a [`FibrousLayer`](@ref) to scale beside the fraction of the skin that
its fibres cover; [`draw_insulation_schematic!`](@ref) and [`draw_insulation_coverage!`](@ref) draw each panel
alone:

```@example plotting
plot_insulation_properties(fur)
```

## Summary

| Function | Draws | For |
|:--|:--|:--|
| [`plot_body`](@ref), [`draw_cutaway!`](@ref) | body in three dimensions | all shapes and composites |
| [`plot_cross_sections`](@ref), [`draw_cross_sections!`](@ref) | sections along and across | cylinder, sphere, ellipsoid, plate |
| [`plot_body_silhouette`](@ref) | body and its silhouette, with sliders | composites |
| [`plot_insulation_properties`](@ref) | fibres and their cover of the skin | a [`FibrousLayer`](@ref) |
