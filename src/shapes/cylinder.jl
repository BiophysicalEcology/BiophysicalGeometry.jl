"""
    Cylinder(; mass, density, volume, length, radius, axis_ratio_b) <: AbstractShape

A cylindrical organism shape, lying along `+x` from `x = 0` to `x = length`.
`axis_ratio_b` is length / diameter. Give any sufficient set of keywords — e.g.
`mass` and `density` with `axis_ratio_b`, or `length` and `radius` with one of
`mass` / `density` — and the rest is solved for. Dimensions are at skin level.
"""
struct Cylinder{M,D,B} <: AbstractCylindrical
    mass::M
    density::D
    axis_ratio_b::B
    Cylinder(::_Resolved, mass::M, density::D, axis_ratio_b::B) where {M,D,B} =
        new{M,D,B}(mass, density, axis_ratio_b)
end

# volume = π·radius²·length; axis_ratio_b = length / (2·radius)
const _CYLINDER_SPEC = _ShapeSpec((:length, :radius), (1, 2), log(π), (:axis_ratio_b => (1, 2, 2.0),))

function Cylinder(; kw...)
    s = _resolve_shape("Cylinder", _CYLINDER_SPEC, NamedTuple(kw))
    Cylinder(_RESOLVED, s.mass, s.density, s.ratios...)
end

# Radial dimension from an enclosed volume; used for both skin and flesh radii.
_cylinder_radius(shape::Cylinder, volume) = cbrt(volume / (shape.axis_ratio_b * π * 2))

function _skin_level(shape::Cylinder, volume)
    radius_skin = _cylinder_radius(shape, volume)
    length_skin = shape.axis_ratio_b * radius_skin * 2
    (; dims = (; radius_skin, length_skin), area = surface_area(shape, radius_skin, length_skin))
end
function _fibrous_level(shape::Cylinder, skin, thickness)
    radius_fibrous = skin.radius_skin + thickness
    length_fibrous = skin.length_skin + thickness * 2
    (; dims = (; radius_fibrous, length_fibrous), area = surface_area(shape, radius_fibrous, length_fibrous))
end
_fat_thickness(shape::Cylinder, skin, flesh_volume, fat_volume) =
    skin.radius_skin - _cylinder_radius(shape, flesh_volume)

# Surface area

surface_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2

# Silhouette area. `outer_dims` selects skin- vs fibrous-level (r, L) by
# insulation, so a single insulation-dispatched wrapper per arity covers
# all four insulation kinds.
silhouette(shape::Cylinder, r, l, θ) = 2 * r * l * abs(sin(θ)) + π * r^2 * abs(cos(θ))
function silhouette(sh::Cylinder, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette(sh, d.radius, d.length, θ)
end
function silhouette(sh::Cylinder, ::AbstractInsulationLayer, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = 2 * d.radius * d.length, parallel = π * d.radius^2)
end

# Radius accessors — shared by every cylindrical shape (`Cylinder`, `HalfCylinder`,
# `Cone`); all store the same `radius_skin` / `radius_fibrous` fields, so the
# dispatch lives once on the family type.

_skin_radius(::AbstractCylindrical, length) = length.radius_skin
_fibrous_radius(::AbstractCylindrical, length) = length.radius_fibrous

# Ratio of the end-B (top) radius to the end-A (base) radius: 1 for a
# cylinder, `top_ratio` for a cone. Lets cylindrical halves, meshes and plots
# treat every cylindrical shape as a frustum.
_top_ratio(::Cylinder) = 1.0

# Composition

attachment_surfaces(::Cylinder) = (EndA, EndB, Lateral)

# Outer (insulation-aware) dimensions. Insulation-dispatched; no runtime
# field-lookup. `radius` matches insulation_radius(body); `length` is the axial extent
# of the outer surface (skin + fur padding at each end).
outer_dims(sh::Cylinder, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Cylinder, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (radius = body.geometry.length.radius_skin, length = body.geometry.length.length_skin)
outer_dims(::Cylinder, ::FibrousLayer, body::AbstractBody) =
    (radius = body.geometry.length.radius_fibrous, length = body.geometry.length.length_fibrous)

# Surface AREAS report the actual outer (insulation-aware) area. This drives
# composition's patch-fits-surface validation and FullCover lookup.
surface_area(::Cylinder, body::AbstractBody, ::EndA) =
    π * insulation_radius(body)^2
surface_area(::Cylinder, body::AbstractBody, ::EndB) =
    π * insulation_radius(body)^2
function surface_area(sh::Cylinder, body::AbstractBody, ::Lateral)
    d = outer_dims(sh, body)
    2 * π * d.radius * d.length
end

# Attachment POSITIONS are anchored at flesh (skin) level so joins meet
# flesh-to-flesh. Local frame: the axis is +x, with the flesh ends at x = 0
# (EndA) and x = length_skin (EndB); the lateral attachment surface has radius
# radius_skin. Angles run around the axis from +y towards +z, so a point at
# `angle` sits at (y, z) = radius·(cos angle, sin angle). The fur overhang past
# either end is part of the outer mesh (and counted in surface_area) but never
# bears attachments.

function validate_range(::Cylinder, body::AbstractBody, loc::EndA)
    R = skin_radius(body)
    loc.radius ≥ zero(loc.radius) && loc.radius ≤ R ||
        error("EndA radius out of range [0, $R]: got $(loc.radius)")
end
validate_range(sh::Cylinder, body::AbstractBody, loc::EndB) =
    validate_range(sh, body, EndA(loc.radius, loc.angle))

function validate_range(::Cylinder, body::AbstractBody, loc::Lateral)
    L = body.geometry.length.length_skin
    loc.position ≥ zero(loc.position) && loc.position ≤ L ||
        error("Lateral position out of range [0, $L]: got $(loc.position)")
end

surface_point(::Cylinder, body::AbstractBody, loc::EndA) =
    (zero(loc.radius), loc.radius * cos(loc.angle), loc.radius * sin(loc.angle))
surface_point(::Cylinder, body::AbstractBody, loc::EndB) =
    (body.geometry.length.length_skin, loc.radius * cos(loc.angle), loc.radius * sin(loc.angle))
function surface_point(::Cylinder, body::AbstractBody, loc::Lateral)
    R = skin_radius(body)
    (loc.position, R * cos(loc.angle), R * sin(loc.angle))
end

surface_normal(::Cylinder, ::AbstractBody, ::EndA) = (-1.0, 0.0, 0.0)
surface_normal(::Cylinder, ::AbstractBody, ::EndB) = (1.0, 0.0, 0.0)
surface_normal(::Cylinder, ::AbstractBody, loc::Lateral) =
    (0.0, cos(loc.angle), sin(loc.angle))

# Centroids (used by FullCover) — also at flesh level.
function surface_centroid(::Cylinder, body::AbstractBody, ::EndA)
    R = skin_radius(body); (zero(R), zero(R), zero(R))
end
function surface_centroid(::Cylinder, body::AbstractBody, ::EndB)
    L = body.geometry.length.length_skin; (L, zero(L), zero(L))
end
function surface_centroid(::Cylinder, body::AbstractBody, ::Lateral)
    R = skin_radius(body); L = body.geometry.length.length_skin
    (L/2, R, zero(R))
end
surface_centroid_normal(::Cylinder, ::AbstractBody, ::EndA) = (-1.0, 0.0, 0.0)
surface_centroid_normal(::Cylinder, ::AbstractBody, ::EndB) = (1.0, 0.0, 0.0)
surface_centroid_normal(::Cylinder, ::AbstractBody, ::Lateral) = (0.0, 1.0, 0.0)
