module BiophysicalGeometry

using LinearAlgebra: det, lu
using StaticArrays: SVector, SMatrix, @SMatrix
using Unitful

export AbstractGeometryPars, AbstractBody, Body
export AbstractShape, Cylinder, Sphere, Ellipsoid, Plate, TriangularPlate, Cone, Unchecked
export Half, HalfCylinder, HalfCone, HalfEllipsoid, HalfSphere
export AbstractCylindrical, AbstractSpherical, AbstractEllipsoidal, AbstractSlab
export AbstractInsulationLayer, AbstractSolidLayer, AbstractPorousLayer
export CompositeInsulation, Naked, FibrousLayer, FatLayer
export SolarOrientation, Intermediate, ParallelToSun, NormalToSun
export SurfaceAreas
export CompositeBody, Join, Attachment, Disc, FullCover, AbstractAttachmentShape, Pose
export AbstractSurface, EndA, EndB, Lateral, Flat, Dome, PoleA, PoleB, Equator, Radial
export Top, Bottom, SideA, SideB, SideC, SideD, Diagonal
export SolarOrientation, Intermediate, ParallelToSun, NormalToSun, ZenithAngleVarying
export Beam, Sky, Ground, Horizon
export SilhouetteResult

export attachment_surfaces
export join_area, join_position, join_partners, internal_distance, flesh_centroid
export geometry, shape, mass, insulation, outer_dims
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius, flesh_volume
export surface_area, silhouette, silhouette!, silhouette_factors, silhouette_factors!
export silhouette_rasterized, silhouette_rasterized!
export outer_insulation
export plot_body, draw_cutaway!, plot_cross_sections, draw_cross_sections!
export plot_body_silhouette
export draw_insulation_schematic!, draw_insulation_coverage!, plot_insulation_properties
export compile_wasm

# Stubs — implemented in BiophysicalGeometryMakieExt when Makie is loaded.
# Fallbacks give a clear message if no Makie backend has been loaded.
const REQUIRES_MAKIE = " requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first."

draw_cutaway!(args...; kwargs...) = error("draw_cutaway! $REQUIRES_MAKIE")
plot_body(args...; kwargs...) = error("plot_body $REQUIRES_MAKIE")
draw_cross_sections!(args...; kwargs...) = error("draw_cross_sections! $REQUIRES_MAKIE")
plot_cross_sections(args...; kwargs...) = error("plot_cross_sections $REQUIRES_MAKIE")
draw_insulation_schematic!(args...; kwargs...) = error("draw_insulation_schematic! $REQUIRES_MAKIE")
draw_insulation_coverage!(args...; kwargs...) = error("draw_insulation_coverage! $REQUIRES_MAKIE")
plot_insulation_properties(args...; kwargs...) = error("plot_insulation_properties $REQUIRES_MAKIE")
plot_body_silhouette(args...; kwargs...) = error("plot_body_silhouette $REQUIRES_MAKIE")

"""
    compile_wasm(model, defaults) -> (; wasm, spec)
    compile_wasm(models) -> (; wasm, spec)
    compile_wasm(path, model, defaults; name = "model")
    compile_wasm(path, models; name = "model")

Compile a model to WebAssembly, to run in a web page. A model is a function from its settings to a
`CompositeBody`, built with `Unchecked` constructors so that nothing can throw: a top-level function, or any other
callable without fields. `defaults` is a `NamedTuple` of its settings, all numbers; their names and order are the
settings' in the page. Several models compile together as a `NamedTuple` of `(model, defaults)` pairs, one module
running any of them.

The first two return the module's bytes, `wasm`, and `spec`, which describes it: each model's settings, defaults
and parts, and the layout of what it writes. With a `path`, they write `name.wasm`, `name.json` and the JavaScript
that loads them, `biophysical.mjs`, into the folder `path`.

In the page, `load` from `biophysical.mjs` takes the module (a URL, an `ArrayBuffer` or a `Response`) and its spec,
and gives `run(model, settings, sun)`: the body posed in triangles to draw, each part's exposed area and mass, the
totals, and the shadow toward the sun. It needs no framework, so it works in Vue, Bonito or plain HTML.

Needs [Whisk.jl](https://github.com/SimonDanisch/Whisk.jl): `using Whisk`.
"""
function compile_wasm end

include("geometry.jl")
include("construction.jl")
include("composition.jl")
include("shapes/plate.jl")
include("shapes/triangular_plate.jl")
include("shapes/cylinder.jl")
include("shapes/sphere.jl")
include("shapes/ellipsoid.jl")
include("shapes/cone.jl")
include("shapes/half.jl")
include("meshes.jl")
include("silhouette.jl")
include("joins.jl")
include("display.jl")

end
