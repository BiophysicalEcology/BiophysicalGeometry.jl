module BiophysicalGeometry

using Unitful

export AbstractGeometryModel, AbstractGeometryPars, AbstractShape, AbstractInsulation, AbstractBody, Body
export Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana
export CompositeInsulation, Naked, Fur, Fat
export SurfaceAreas
export geometry, shape, insulation
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius, flesh_volume
export surface_area, silhouette_area, SolarOrientation, Intermediate, ParallelToSun, NormalToSun, ZenithAngleVarying
export plot_body, draw_cutaway!, plot_cross_sections, draw_cross_sections!

# Stubs — implemented in BiophysicalGeometryMakieExt when Makie is loaded
function draw_cutaway! end
function plot_body end
function draw_cross_sections! end
function plot_cross_sections end

include("geometry.jl")
include("shapes/plate.jl")
include("shapes/cylinder.jl")
include("shapes/sphere.jl")
include("shapes/ellipsoid.jl")
include("shapes/desert_iguana.jl")
include("shapes/leopard_frog.jl")

end
