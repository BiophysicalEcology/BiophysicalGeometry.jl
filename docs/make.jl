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
