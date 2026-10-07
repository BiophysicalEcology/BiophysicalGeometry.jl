# Build an animal

Build an animal from simple shapes with sliders and see its areas change. The torso is a dorsal and a ventral half
with their own coats, see [Back and belly](tutorials/two_parts.md). To it can be joined a head, a neck, a nose or
beak, ears, legs, wings and a tail. The body can lie level, be tipped up as in a kangaroo, or stand upright with
arms, as in the [human tutorial](tutorials/human.md).

```@example builder
using Markdown # hide
url = get(ENV, "BUILDER_URL", "") # hide
if isempty(url) # hide
    Markdown.parse("The app runs on your own computer, see below.") # hide
else # hide
    HTML("""<iframe src="$url" title="Build an animal" style="width: 100%; height: 1100px; border: 1px solid var(--vp-c-divider); border-radius: 8px;"></iframe>""") # hide
end # hide
```

## Run it

The app is a small web page served by Julia:

```julia
using BiophysicalGeometry, Bonito, WGLMakie
server = BiophysicalGeometry.app()   # opens http://127.0.0.1:8080 in your browser
# ...
close(server)
```

Legs sized by similarity also need [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl)
installed. The `app` folder of the repository has the same as a script, and a container to serve it.

## The code is the model

Every change in the app writes the Julia code for the animal, runs it with BiophysicalGeometry.jl, and shows what the
package makes of it: the body, its areas part by part, and its silhouette to the sun. The code is shown under the
animal, to copy into Julia and go further: change a join, add a part, or ask for each part's view of the sky, see the
[tutorials](tutorials/dog.md).

The same code can be had without the app, from the settings of the page:

```@example builder
using BiophysicalGeometry
const AB = BiophysicalGeometry.AnimalBuilder
settings = (; mass = 20.0, neck = true, nose = "Nose", ears = "Cone", tail = true)
code = AB.animal_code(settings)
print(code)
```

```@example builder
animal = AB.build_animal(code)
total_area(animal), skin_area(animal)
```

## The animals it starts from

```@example builder
using CairoMakie, Main.FigureHelpers
fig = Figure(size = (1100, 480))
for (i, (name, preset)) in enumerate(AB.PRESETS)
    ax = body_axis(fig[fld1(i, 5), mod1(i, 5)]; decorations = false, azimuth = 1.25π, elevation = π / 7,
        title = name, titlesize = 13)
    draw_parts!(ax, AB.build_animal(AB.animal_code(preset)))
end
fig
```

## Controls

| Control | What it sets |
|:--|:--|
| mass, posture, pitch | the whole animal; the torso lies level, tips head-up, or stands upright |
| torso shape, length / width, fat | the two halves of the torso, a `HalfCylinder` or `HalfEllipsoid`, with a `FatLayer` |
| coat depths | `FibrousLayer`s on the back, the belly, and the head and limbs |
| head, neck, nose or beak, ears | each a fraction of the mass, with its shape and posture |
| legs | none, two or four; set by hand, or sized from body mass by elastic or geometric similarity, see [Geometry and allometry](tutorials/scaling.md); hind legs can differ |
| arms | for an upright animal |
| wings, tail | for a level animal |
| sun | the direction of the silhouette |
