module BiophysicalGeometry

using Unitful

export AbstractGeometryModel, AbstractGeometryPars, AbstractShape, AbstractInsulation, AbstractBody, Body
export Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana
export CompositeInsulation, Naked, Fur, Fat
export SurfaceAreas
export geometry, shape, insulation
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius, flesh_volume
export surface_area, silhouette_area, SolarOrientation, Intermediate, ParallelToSun, NormalToSun, ZenithAngleVarying
export outer_insulation

include("geometry.jl")
include("shapes/plate.jl")
include("shapes/cylinder.jl")
include("shapes/sphere.jl")
include("shapes/ellipsoid.jl")
include("shapes/desert_iguana.jl")
include("shapes/leopard_frog.jl")

end
