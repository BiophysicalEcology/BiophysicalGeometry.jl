# Build an animal

Move the sliders to build an animal from simple shapes and see its areas change. The torso is a dorsal and a ventral
half with their own coats, see [Back and belly](tutorials/two_parts.md). Each animal has its own parts, among a
head, a neck, a nose or beak, ears, legs, arms, wings and a tail, and each can be resized and angled.

```@raw html
<script setup>
import AnimalBuilder from './components/AnimalBuilder.vue'
</script>

<ClientOnly>
  <AnimalBuilder />
</ClientOnly>
```

The numbers are BiophysicalGeometry.jl's own. The page runs the package in your browser, compiled to WebAssembly
with [Whisk.jl](https://github.com/SimonDanisch/Whisk.jl): `AnimalBuilder.build` makes the animal from the sliders
with the package's `Unchecked()` constructors, and the areas, the triangles drawn and the silhouette to the sun come
from the same functions you would call in Julia. No Julia server is involved.

The same animals can be built in Julia:

```julia
using BiophysicalGeometry
using BiophysicalGeometry.AnimalBuilder: build, animal, settings

dog = build(animal("Dog"), merge(settings("Dog"), (; mass = 30.0, neckAngle = 20.0)))
total_area(dog), silhouette(dog, Beam(0.0, 0.0, 1.0))
```
