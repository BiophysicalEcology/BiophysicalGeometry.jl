"""
    HalfCylinder <: AbstractShape

A half-cylinder organism shape — a cylinder cut along a plane through its
axis. Dimensions are sized so that two half-cylinders of equal mass joined
at their flat faces reconstruct the full `Cylinder(2*mass, density, b)`.

In local coordinates the axis lies along `+z`, the curved surface occupies
`y ≥ 0` half-space, and the flat face is at `y = 0`.

The flat-face area is reported at *skin* level even when fur is present,
so that two halves (mixed insulation) join cleanly with `FullCover`.
"""
mutable struct HalfCylinder{M,D,B} <: AbstractCylindrical
    mass::M
    density::D
    axis_ratio_b::B
end

# A half of mass m is the full `Cylinder` of mass 2m split along its axis, so
# all radius/length/fat/fur math is inherited from `Cylinder` (single source of
# truth) via `full`; the half only differs in area (half the lateral + the two
# semicircular ends, plus a flat rectangular face at skin level) and volume.
_full_cylinder(sh::HalfCylinder) = Cylinder(2 * sh.mass, sh.density, sh.axis_ratio_b)

geometry(sh::HalfCylinder, ins::Naked)        = _half_geometry(sh, geometry(_full_cylinder(sh), ins))
geometry(sh::HalfCylinder, fat::FatLayer)     = _half_geometry(sh, geometry(_full_cylinder(sh), fat))
geometry(sh::HalfCylinder, fur::FibrousLayer) = _half_geometry(sh, geometry(_full_cylinder(sh), fur), fur)
geometry(sh::HalfCylinder, fur::FibrousLayer, fat::FatLayer) =
    _half_geometry(sh, geometry(_full_cylinder(sh), fur, fat), fur)

_half_cyl_dome(r, L) = π * r * L + π * r^2  # half the full lateral + two half-discs

function _half_geometry(sh::HalfCylinder, full)
    len = full.length
    flat = 2 * len.radius_skin * len.length_skin
    total = _half_cyl_dome(len.radius_skin, len.length_skin) + flat
    vol = sh.mass / sh.density
    Geometry(vol, vol^(1 / 3), len, SurfaceAreas(; total))
end
function _half_geometry(sh::HalfCylinder, full, fur::FibrousLayer)
    len = full.length
    flat = 2 * len.radius_skin * len.length_skin
    total = _half_cyl_dome(len.radius_fibrous, len.length_fibrous) + flat
    skin = _half_cyl_dome(len.radius_skin, len.length_skin) + flat
    convection = skin - insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    vol = sh.mass / sh.density
    Geometry(vol, vol^(1 / 3) + fur.thickness, len, SurfaceAreas(; total, skin, convection))
end

# Silhouette area: half of the equivalent full cylinder's, so two halves
# joined back together reproduce the full cylinder's silhouettes. `outer_dims`
# picks skin- vs fibrous-level (r, L) so one wrapper per arity suffices.
function silhouette_area(::HalfCylinder, r, L, θ)
    r * L * sin(θ) + (π * r^2 / 2) * cos(θ)
end
function silhouette_area(sh::HalfCylinder, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette_area(sh, d.r, d.L, θ)
end
function silhouette_area(sh::HalfCylinder, ::AbstractInsulationLayer, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = d.r * d.L, parallel = π * d.r^2 / 2)
end

# Radii come from the shared `AbstractCylindrical` dispatch in cylinder.jl.

# Composition

attachment_surfaces(::HalfCylinder) = (EndA, EndB, Lateral, Flat)

# Outer (insulation-aware) dimensions. `r` matches insulation_radius(body).
outer_dims(sh::HalfCylinder, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::HalfCylinder, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (r = body.geometry.length.radius_skin, L = body.geometry.length.length_skin)
outer_dims(::HalfCylinder, ::FibrousLayer, body::AbstractBody) =
    (r = body.geometry.length.radius_fibrous, L = body.geometry.length.length_fibrous)

# FibrousLayer axial half-thickness (0 for naked); fur extends past the flesh by this on each end.
_halfcyl_fur_pad(body) =
    (outer_dims(shape(body), body).L - body.geometry.length.length_skin) / 2

# Per-surface areas. Outer (insulation-aware) for the curved/end faces;
# skin level for the flat face so two halves join cleanly under FullCover.
function surface_area(::HalfCylinder, body::AbstractBody, ::EndA)
    R = insulation_radius(body)
    π * R^2 / 2
end
surface_area(sh::HalfCylinder, body::AbstractBody, ::EndB) =
    surface_area(sh, body, EndA())
function surface_area(sh::HalfCylinder, body::AbstractBody, ::Lateral)
    d = outer_dims(sh, body)
    π * d.r * d.L
end
function surface_area(::HalfCylinder, body::AbstractBody, ::Flat)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    2 * r * L
end

# Attachment positions are anchored at flesh (skin) level so joins are
# flesh-to-flesh: local z=0 is the flesh end at -z, z=length_skin is the
# flesh end at +z. Lateral attachment radius = radius_skin. The fur axial
# overhang past z=0 / z=length_skin is part of the outer mesh and counted
# in surface_area, but bears no attachments.

function validate_range(::HalfCylinder, body::AbstractBody, loc::EndA)
    R = body.geometry.length.radius_skin
    loc.r ≥ zero(loc.r) && loc.r ≤ R || error("EndA r out of range [0, $R]: $(loc.r)")
    0 ≤ loc.φ ≤ π || error("EndA φ out of range [0, π]: $(loc.φ)")
end
validate_range(sh::HalfCylinder, body::AbstractBody, loc::EndB) =
    validate_range(sh, body, EndA(loc.r, loc.φ))
function validate_range(::HalfCylinder, body::AbstractBody, loc::Lateral)
    L = body.geometry.length.length_skin
    loc.z ≥ zero(loc.z) && loc.z ≤ L || error("Lateral z out of range [0, $L]: $(loc.z)")
    0 ≤ loc.φ ≤ π || error("Lateral φ out of range [0, π]: $(loc.φ)")
end
# For HalfCylinder, Flat uses (u=z, v=x).
function validate_range(::HalfCylinder, body::AbstractBody, loc::Flat)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    loc.u ≥ zero(loc.u) && loc.u ≤ L || error("Flat z out of range [0, $L]: $(loc.u)")
    abs(loc.v) ≤ r || error("Flat x out of range ±$r: $(loc.v)")
end

surface_point(::HalfCylinder, body::AbstractBody, loc::EndA) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), zero(loc.r))
surface_point(::HalfCylinder, body::AbstractBody, loc::EndB) =
    (loc.r * cos(loc.φ), loc.r * sin(loc.φ), body.geometry.length.length_skin)
function surface_point(::HalfCylinder, body::AbstractBody, loc::Lateral)
    R = body.geometry.length.radius_skin
    (R * cos(loc.φ), R * sin(loc.φ), loc.z)
end
function surface_point(::HalfCylinder, body::AbstractBody, loc::Flat)
    (loc.v, zero(loc.v), loc.u)
end

surface_normal(::HalfCylinder, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_normal(::HalfCylinder, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_normal(::HalfCylinder, ::AbstractBody, loc::Lateral) =
    (cos(loc.φ), sin(loc.φ), 0.0)
surface_normal(::HalfCylinder, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)

# Centroids — also at flesh level.
function surface_centroid(::HalfCylinder, body::AbstractBody, ::EndA)
    R = body.geometry.length.radius_skin; (zero(R), R/2, zero(R))
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::EndB)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), R/2, L)
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), R, L/2)
end
function surface_centroid(::HalfCylinder, body::AbstractBody, ::Flat)
    L = body.geometry.length.length_skin
    (zero(L), zero(L), L/2)
end
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Lateral) = (0.0, 1.0, 0.0)
surface_centroid_normal(::HalfCylinder, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)
