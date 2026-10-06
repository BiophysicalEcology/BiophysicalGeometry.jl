module BiophysicalGeometry

using LinearAlgebra: svd
using Unitful

export AbstractGeometryModel, AbstractGeometryPars, AbstractBody, Body
export AbstractShape, Cylinder, Sphere, Ellipsoid, Plate, Cone
export Half, HalfCylinder, HalfCone, HalfEllipsoid, HalfSphere
export AbstractCylindrical, AbstractSpherical, AbstractEllipsoidal, AbstractSlab
export AbstractInsulationLayer, AbstractSolidLayer, AbstractPorousLayer
export CompositeInsulation, Naked, FibrousLayer, FatLayer
export SolarOrientation, Intermediate, ParallelToSun, NormalToSun
export SurfaceAreas
export CompositeBody, Join, Attachment, Disc, FullCover, AbstractAttachmentShape, Pose
export AbstractSurface, EndA, EndB, Lateral, Flat, Dome, PoleA, PoleB, Equator, Radial
export Top, Bottom, SideA, SideB, SideC, SideD
export SolarOrientation, Intermediate, ParallelToSun, NormalToSun, ZenithAngleVarying
export Beam, Sky, Ground, Horizon
export SilhouetteResult

export attachment_surfaces
export join_area, join_position, join_partners, internal_distance, flesh_centroid
export geometry, shape, mass, insulation, outer_dims
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius, flesh_volume
export surface_area, silhouette, silhouette_factors
export silhouette_rasterized
export outer_insulation
export plot_body, draw_cutaway!, plot_cross_sections, draw_cross_sections!
export plot_body_silhouette
export draw_insulation_schematic!, draw_insulation_coverage!, plot_insulation_properties

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

include("geometry.jl")
include("construction.jl")
include("composition.jl")
include("shapes/plate.jl")
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
