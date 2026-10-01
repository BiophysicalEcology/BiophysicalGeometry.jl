# Build an animal

Move the sliders to build an animal from simple shapes and see its areas change. The torso is a dorsal and a ventral
half with their own coats, see [Back and belly](tutorials/two_parts.md), with a head and legs joined to it, and a neck, a nose or beak, ears, wings and a tail if wanted.

```@raw html
<script setup>
import AnimalBuilder from './components/AnimalBuilder.vue'
</script>

<ClientOnly>
  <AnimalBuilder />
</ClientOnly>
```

Under the animal is the Julia code that builds the same body with BiophysicalGeometry.jl, updated as the sliders
move. Copy it into Julia to go further:
change a join, add a neck or a tail, or ask for each part's view of the sky, see the [tutorials](tutorials/dog.md).

This page computes its numbers in the browser, with a copy of the package's calculations for these shapes. The copy
is checked against the package each time the documentation is built.
