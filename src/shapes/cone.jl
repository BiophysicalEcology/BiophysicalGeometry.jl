"""
    Cone(; mass, density, volume, length, radius, axis_ratio_b, top_ratio=0.0) <: AbstractShape

A right circular cone or truncated cone (frustum), lying along `+x` with its
base disc (radius `radius`) at `x = 0` and its top disc at `x = length`.
`axis_ratio_b` is length / base diameter; `top_ratio` is top / base radius —
`0` makes a sharp cone, values in `(0, 1)` a frustum. Give any sufficient set
of the other keywords, as for [`Cylinder`](@ref). Insulation expands radii and
length as for `Cylinder`; attachment positions stay at flesh level.
"""
struct Cone{M,D,B,T} <: AbstractCylindrical
    mass::M
    density::D
    axis_ratio_b::B
    top_ratio::T
    Cone(::Resolved, mass::M, density::D, axis_ratio_b::B, top_ratio::T) where {M,D,B,T} =
        new{M,D,B,T}(mass, density, axis_ratio_b, top_ratio)
end

# volume = (π/3)(1 + t + t²)·radius²·length; axis_ratio_b = length / (2·radius)
_cone_spec(t) = ShapeSpec((:length, :radius), (1, 2), log(π / 3 * _cone_volume_factor(t)),
                           (:axis_ratio_b => (1, 2, 2.0),))

function Cone(; top_ratio = 0.0, kw...)
    0 <= top_ratio <= 1 || throw(ArgumentError("Cone `top_ratio` must be in [0, 1], got $top_ratio"))
    s = _resolve_shape("Cone", _cone_spec(top_ratio), NamedTuple(kw))
    Cone(RESOLVED, s.mass, s.density, s.axis_ratio_b, top_ratio)
end

# Volume of a frustum = (π/3) · L · (R² + R·r + r²) with r = top_ratio·R.
# With L = 2·b·R: V = (2π/3) · b · R³ · (1 + t + t²).
_cone_volume_factor(t) = 1 + t + t^2

function _cone_radius(volume, b, t)
    cbrt(3 * volume / (2π * b * _cone_volume_factor(t)))
end

function _skin_level(shape::Cone, volume)
    radius_skin = _cone_radius(volume, shape.axis_ratio_b, shape.top_ratio)
    length_skin = 2 * shape.axis_ratio_b * radius_skin
    (; dims = (; radius_skin, length_skin), area = surface_area(shape, radius_skin, length_skin))
end
function _fibrous_level(shape::Cone, skin, thickness)
    radius_fibrous = skin.radius_skin + thickness
    length_fibrous = skin.length_skin + 2 * thickness
    (; dims = (; radius_fibrous, length_fibrous), area = surface_area(shape, radius_fibrous, length_fibrous))
end
_fat_thickness(shape::Cone, skin, flesh_volume, fat_volume) =
    skin.radius_skin - _cone_radius(flesh_volume, shape.axis_ratio_b, shape.top_ratio)

# Aggregate surface area for a frustum: base disc + top disc + slant surface.
# Slant length s = sqrt((R - r)² + L²) where r = t·R.
function surface_area(shape::Cone, R, L)
    t = shape.top_ratio
    r = t * R
    s = sqrt((R - r)^2 + L^2)
    π * R^2 + π * r^2 + π * (R + r) * s
end

# Silhouette. The shadow of a frustum is the convex hull of the shadows of its
# two end discs: homothetic ellipses (minor/major = |cos θ|) whose centres sit
# L·|sin θ| apart along the minor axis. Stretching the minor axis by 1/|cos θ|
# turns them into circles of radii R ≥ r a distance D = L·|tan θ| apart, whose
# hull is two arcs joined by external tangents at angle α, sin α = (R - r)/D:
#     (π/2 + α)·R² + (π/2 - α)·r² + (R + r)·D·cos α
# (just π·R² once one disc's shadow lies inside the other's). Shrinking back
# by |cos θ| gives the area below. A cylinder (r = R) reduces to
# 2·R·L·sin θ + π·R²·cos θ. `outer_dims` picks skin- vs fibrous-level (r, L).
function _frustum_silhouette(R1, R2, L, θ)
    R, r = max(R1, R2), min(R1, R2)
    c, s = abs(cos(θ)), abs(sin(θ))
    Ls = L * s # projected axis length
    Ls <= (R - r) * c && return π * R^2 * c
    sinα = (R - r) * c / Ls
    α = asin(sinα)
    c * ((π / 2 + α) * R^2 + (π / 2 - α) * r^2) + (R + r) * sqrt(1 - sinα^2) * Ls
end

silhouette(sh::Cone, r, L, θ) = _frustum_silhouette(r, sh.top_ratio * r, L, θ)
function silhouette(sh::Cone, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette(sh, d.radius, d.length, θ)
end
function silhouette(sh::Cone, ::AbstractInsulationLayer, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = (1 + sh.top_ratio) * d.radius * d.length, parallel = π * max(1, sh.top_ratio)^2 * d.radius^2)
end

# Radii come from the shared `AbstractCylindrical` dispatch in cylinder.jl.

top_ratio(sh::Cone) = sh.top_ratio

# Composition
#
# `EndA` is the base disc (z=0), `EndB` is the top disc (z=length_skin,
# radius = top_ratio * radius_skin), `Lateral` is the slant surface.

attachment_surfaces(::Cone) = (EndA, EndB, Lateral)

# Outer (insulation-aware) dimensions.
outer_dims(sh::Cone, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Cone, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (radius = body.geometry.length.radius_skin, length = body.geometry.length.length_skin)
outer_dims(::Cone, ::FibrousLayer, body::AbstractBody) =
    (radius = body.geometry.length.radius_fibrous, length = body.geometry.length.length_fibrous)

function surface_area(sh::Cone, body::AbstractBody, ::EndA)
    d = outer_dims(sh, body); π * d.radius^2
end
function surface_area(sh::Cone, body::AbstractBody, ::EndB)
    d = outer_dims(sh, body); π * (sh.top_ratio * d.radius)^2
end
function surface_area(sh::Cone, body::AbstractBody, ::Lateral)
    d = outer_dims(sh, body)
    r = sh.top_ratio * d.radius
    s = sqrt((d.radius - r)^2 + d.length^2)
    π * (d.radius + r) * s
end

# Local frame as for `Cylinder`: axis +x, base (EndA) at x = 0, top (EndB) at
# x = length_skin, angles around the axis from +y towards +z.

function validate_range(::Cone, body::AbstractBody, loc::EndA)
    R = body.geometry.length.radius_skin
    loc.radius ≥ zero(loc.radius) && loc.radius ≤ R ||
        error("EndA radius out of range [0, $R]: $(loc.radius)")
end
function validate_range(shape::Cone, body::AbstractBody, loc::EndB)
    Rt = shape.top_ratio * body.geometry.length.radius_skin
    Rt > zero(Rt) || error("EndB has zero radius (top_ratio=0); use a Disc(0) only")
    loc.radius ≥ zero(loc.radius) && loc.radius ≤ Rt ||
        error("EndB radius out of range [0, $Rt]: $(loc.radius)")
end
function validate_range(::Cone, body::AbstractBody, loc::Lateral)
    L = body.geometry.length.length_skin
    loc.position ≥ zero(loc.position) && loc.position ≤ L ||
        error("Lateral position out of range [0, $L]: $(loc.position)")
end

# Attachment positions at flesh (skin) level.
surface_point(::Cone, body::AbstractBody, loc::EndA) =
    (zero(loc.radius), loc.radius * cos(loc.angle), loc.radius * sin(loc.angle))
surface_point(::Cone, body::AbstractBody, loc::EndB) =
    (body.geometry.length.length_skin, loc.radius * cos(loc.angle), loc.radius * sin(loc.angle))
function surface_point(shape::Cone, body::AbstractBody, loc::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    rx = R * (1 - (1 - shape.top_ratio) * loc.position / L)
    (loc.position, rx * cos(loc.angle), rx * sin(loc.angle))
end

surface_normal(::Cone, ::AbstractBody, ::EndA) = (-1.0, 0.0, 0.0)
surface_normal(::Cone, ::AbstractBody, ::EndB) = (1.0, 0.0, 0.0)
# The slant narrows by Δr = R - r over the length L, so the outward normal is
# (radial)·(L/s) + (axial)·(Δr/s) with s the slant length. Each ratio is a
# ratio of lengths, so unitless.
function _cone_slant_normal(shape::Cone, body, angle)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    Δr = R * (1 - shape.top_ratio)
    s = sqrt(Δr^2 + L^2)
    (Δr/s, L/s * cos(angle), L/s * sin(angle))
end
surface_normal(shape::Cone, body::AbstractBody, loc::Lateral) =
    _cone_slant_normal(shape, body, loc.angle)

# Centroids (FullCover).
function surface_centroid(::Cone, body::AbstractBody, ::EndA)
    R = body.geometry.length.radius_skin; (zero(R), zero(R), zero(R))
end
function surface_centroid(::Cone, body::AbstractBody, ::EndB)
    L = body.geometry.length.length_skin; (L, zero(L), zero(L))
end
function surface_centroid(shape::Cone, body::AbstractBody, ::Lateral)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (L/2, R * (1 + shape.top_ratio) / 2, zero(R)) # midpoint radius
end
surface_centroid_normal(::Cone, ::AbstractBody, ::EndA) = (-1.0, 0.0, 0.0)
surface_centroid_normal(::Cone, ::AbstractBody, ::EndB) = (1.0, 0.0, 0.0)
surface_centroid_normal(shape::Cone, body::AbstractBody, ::Lateral) =
    _cone_slant_normal(shape, body, 0.0)
