module BiophysicalGeometry

using Unitful

export AbstractGeometryModel, AbstractGeometryPars, AbstractShape, AbstractInsulation, AbstractBody, Body
export Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana
export CompositeInsulation, Naked, Fur, Fat
export SurfaceAreas
export geometry, shape, insulation
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius, flesh_volume, total_volume
export surface_area, silhouette_area, SolarOrientation, Intermediate, ParallelToSun, NormalToSun, ZenithAngleVarying
export plot_body, draw_cutaway!, plot_cross_sections, draw_cross_sections!, draw_organism!, plot_organism
export draw_organism_prone!, draw_organism_upright!
# Multi-part organisms
export Organism, BodyPart, Join, OrganismGeometry
export OrganismPosture, Upright, Prone
export junction_area, junction_radius, part_positions, PartPosition
export butterfly_config_factors

# Fallback stubs — real implementations live in BiophysicalGeometryMakieExt.
# These throw a helpful error when Makie has not been loaded.
const _MAKIE_MSG = "load GLMakie or CairoMakie first, e.g. `using CairoMakie`"
draw_cutaway!(args...; kwargs...)      = error("draw_cutaway!: $_MAKIE_MSG")
plot_body(args...; kwargs...)          = error("plot_body: $_MAKIE_MSG")
draw_cross_sections!(args...; kwargs...) = error("draw_cross_sections!: $_MAKIE_MSG")
plot_cross_sections(args...; kwargs...)  = error("plot_cross_sections: $_MAKIE_MSG")
draw_organism!(args...; kwargs...)        = error("draw_organism!: $_MAKIE_MSG")
plot_organism(args...; kwargs...)         = error("plot_organism: $_MAKIE_MSG")
draw_organism_prone!(args...; kwargs...)  = error("draw_organism_prone!: $_MAKIE_MSG")
draw_organism_upright!(args...; kwargs...) = error("draw_organism_upright!: $_MAKIE_MSG")

include("geometry.jl")
include("shapes/plate.jl")
include("shapes/cylinder.jl")
include("shapes/sphere.jl")
include("shapes/ellipsoid.jl")
include("shapes/desert_iguana.jl")
include("shapes/leopard_frog.jl")
include("organism.jl")

end
