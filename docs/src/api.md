# API

## Bodies

```@docs
Body
CompositeBody
AbstractBody
geometry
BiophysicalGeometry.Geometry
shape
insulation
mass
SurfaceAreas
```

## Shapes

```@docs
Sphere
Cylinder
Cone
Ellipsoid
Plate
TriangularPlate
Half
HalfCylinder
HalfCone
HalfEllipsoid
HalfSphere
Unchecked
AbstractShape
AbstractCylindrical
AbstractSpherical
AbstractEllipsoidal
AbstractSlab
```

## Layers

```@docs
Naked
FibrousLayer
FatLayer
CompositeInsulation
outer_insulation
AbstractInsulationLayer
AbstractPorousLayer
AbstractSolidLayer
```

## Areas, radii and volumes

```@docs
total_area
skin_area
evaporation_area
surface_area
flesh_radius
skin_radius
insulation_radius
flesh_volume
outer_dims
```

## Joins

```@docs
Join
Attachment
Disc
FullCover
AbstractAttachmentShape
Pose
join_partners
join_area
join_position
internal_distance
flesh_centroid
```

## Surfaces

```@docs
AbstractSurface
attachment_surfaces
EndA
EndB
Lateral
PoleA
PoleB
Equator
Radial
Flat
Dome
Top
Bottom
SideA
SideB
SideC
SideD
Diagonal
```

## Silhouettes

```@docs
silhouette
silhouette!
silhouette_rasterized
silhouette_rasterized!
SilhouetteResult
silhouette_factors
silhouette_factors!
Beam
Sky
Ground
Horizon
SolarOrientation
NormalToSun
ParallelToSun
Intermediate
ZenithAngleVarying
```

## Plotting

```@docs
plot_body
draw_cutaway!
plot_cross_sections
draw_cross_sections!
plot_body_silhouette
plot_insulation_properties
draw_insulation_schematic!
draw_insulation_coverage!
```

## WebAssembly

```@docs
compile_wasm
```
