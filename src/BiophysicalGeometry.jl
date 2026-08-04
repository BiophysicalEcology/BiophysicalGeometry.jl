module BiophysicalGeometry

using Unitful

export AbstractGeometryModel, AbstractGeometryPars, AbstractBody, Body
export AbstractShape, Cylinder, Sphere, Ellipsoid, Plate, Cone, LeopardFrog, DesertIguana
export HalfCylinder, HalfEllipsoid, TriMesh
export AbstractCylindrical, AbstractSpherical, AbstractEllipsoidal, AbstractSlab
export AbstractInsulationLayer, CompositeInsulation, Naked, FibrousLayer, FatLayer
export AbstractPorousLayer, AbstractSolidLayer
export SolarOrientation, Intermediate, ParallelToSun, NormalToSun
export SurfaceAreas
export CompositeBody, Join, Attachment, Disc, FullCover, AbstractAttachmentShape, Pose
export AbstractSurface, EndA, EndB, Lateral, Flat, Dome, PoleA, PoleB, Equator, Radial
export Top, Bottom, SideA, SideB, SideC, SideD
export attachment_surfaces
export join_area, join_position, join_partners, internal_distance, flesh_centroid

export geometry, shape, insulation, outer_dims
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius, flesh_volume
export surface_area, silhouette_area, silhouette_rasterized, SilhouetteResult
export SolarOrientation, Intermediate, ParallelToSun, NormalToSun, ZenithAngleVarying
export plot_body, draw_cutaway!, plot_cross_sections, draw_cross_sections!
export plot_body_silhouette
export draw_insulation_schematic!, draw_insulation_coverage!, plot_insulation_properties

# Stubs — implemented in BiophysicalGeometryMakieExt when Makie is loaded.
# Fallbacks give a clear message if no Makie backend has been loaded.
function draw_cutaway!(args...; kwargs...)
    error("draw_cutaway! requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end
function plot_body(args...; kwargs...)
    error("plot_body requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end
function draw_cross_sections!(args...; kwargs...)
    error("draw_cross_sections! requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end
function plot_cross_sections(args...; kwargs...)
    error("plot_cross_sections requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end
function draw_insulation_schematic!(args...; kwargs...)
    error("draw_insulation_schematic! requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end
function draw_insulation_coverage!(args...; kwargs...)
    error("draw_insulation_coverage! requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end
function plot_insulation_properties(args...; kwargs...)
    error("plot_insulation_properties requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end
function plot_body_silhouette(args...; kwargs...)
    error("plot_body_silhouette requires a Makie backend — add `using GLMakie` (or CairoMakie / WGLMakie) first.")
end

include("geometry.jl")
include("composition.jl")
include("shapes/plate.jl")
include("shapes/cylinder.jl")
include("shapes/sphere.jl")
include("shapes/ellipsoid.jl")
include("shapes/cone.jl")
include("shapes/half_cylinder.jl")
include("shapes/half_ellipsoid.jl")
include("shapes/trimesh.jl")
include("shapes/desert_iguana.jl")
include("shapes/leopard_frog.jl")
include("meshes.jl")
include("silhouette.jl")
include("joins.jl")

end
