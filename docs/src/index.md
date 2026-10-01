```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "BiophysicalGeometry.jl"
  text: "Consider a non-spherical cow"
  tagline: "areas, lengths and volumes of organisms from simple shapes, with layers of fat and fur, joined into bodies of many parts, with units."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/BiophysicalGeometry.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: 🥚 Shapes
    details: <a class="highlight-link">Spheres, cylinders, ellipsoids, plates and cones</a> sized from mass and density, with the same functions for every shape.
    link: /manual/shapes
  - title: 🧥 Layers
    details: <a class="highlight-link">Fat inside the skin and fur or feathers outside it</a>, giving the radii and areas that heat passes through.
    link: /manual/layers
  - title: 🐕 Bodies of many parts
    details: A back and a belly, or a torso, head and legs, <a class="highlight-link">joined at named surfaces</a> as a graph of parts.
    link: /manual/graphs
  - title: ☀️ Silhouettes
    details: The area that intercepts sunlight from any direction, with <a class="highlight-link">parts shading each other</a>, and each part's view of sky and ground.
    link: /manual/silhouettes
  - title: 🛠️ Build an animal
    details: An interactive page to <a class="highlight-link">assemble an animal with sliders</a>, watch its areas and silhouette change, and copy the code that builds it.
    link: /builder
  - title: 🐄 Tutorials
    details: From a first sphere to <a class="highlight-link">a cow and a human</a>, compared with allometry and with NicheMapR.
    link: /tutorials/first_body
  - title: 📏 Units
    details: Every mass, length, area and volume is a <a class="highlight-link">Unitful.jl</a> quantity.
    link: /manual/introduction
---
```

## How to install BiophysicalGeometry.jl?

BiophysicalGeometry.jl can be installed from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/BiophysicalGeometry.jl")
```

## Manual

BiophysicalGeometry.jl computes areas, lengths and volumes of organisms for biophysical calculations, from
geometry: simple shapes, alone or joined into bodies of many parts. The other way to get them, from empirical scaling
equations, is the domain of [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl).

It is part of the [BiophysicalEcology](https://github.com/BiophysicalEcology) ecosystem for mechanistic niche
modelling, where it provides the bodies whose heat budgets are solved by
[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) and whose postures are changed by
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl). See the
[Introduction](manual/introduction.md) for the design of the package and its origin in
[NicheMapR](https://github.com/mrke/NicheMapR).
