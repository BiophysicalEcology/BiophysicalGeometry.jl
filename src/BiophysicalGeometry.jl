module BiophysicalGeometry

using Unitful

export AbstractGeometryModel, AbstractGeometryPars, AbstractBody, Body
export AbstractShape, Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana
export AbstractInsulation, CompositeInsulation, Naked, Fur, Fat
export SolarOrientation, Intermediate, ParallelToSun, NormalToSun
export SurfaceAreas

export geometry, shape, insulation
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius, flesh_volume
export surface_area, silhouette_area, SolarOrientation, Intermediate, ParallelToSun, NormalToSun, ZenithAngleVarying

include("geometry.jl")
include("shapes/plate.jl")
include("shapes/cylinder.jl")
include("shapes/sphere.jl")
include("shapes/ellipsoid.jl")
include("shapes/desert_iguana.jl")
include("shapes/leopard_frog.jl")

end
