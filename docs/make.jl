using Documenter
using DocumenterVitepress
using BiophysicalGeometry
using CairoMakie
using Unitful
import BiophysicalGeometry: Sphere, Top, Bottom  # also exported by Makie

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# Helpers for the figures, loaded in the examples with `using Main.FigureHelpers`
include("figure_helpers.jl")

# The "Build an animal" page runs the package in the browser: its animals, compiled to wasm.
using Whisk
const AB = BiophysicalGeometry.AnimalBuilder
compile_wasm(joinpath(@__DIR__, "src", "components", "builder"),
             NamedTuple{Symbol.(AB.ANIMAL_NAMES)}(map(name -> (AB.animal(name), AB.settings(name)), AB.ANIMAL_NAMES));
             name = "animals")

makedocs(
    modules = [BiophysicalGeometry, Base.get_extension(BiophysicalGeometry, :BiophysicalGeometryMakieExt)],
    sitename = "BiophysicalGeometry.jl",
    authors = "Michael Kearney, Rafael Schouten et al.",
    clean = true,
    doctest = false,
    checkdocs = :exports,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/BiophysicalGeometry.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/BiophysicalEcology/BiophysicalGeometry.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
