"""
    Cone(mass, density, b, top_ratio=0.0) <: AbstractShape

A right circular cone or truncated cone (frustum). `b` is the length to
base-diameter ratio. `top_ratio` is the ratio of apex (top) radius to base
radius — `0` makes a sharp cone, values in `(0, 1)` make a frustum.

Local frame: axis along `+z`, base disc at `z = 0` with radius `radius_skin`,
top disc at `z = length_skin` with radius `top_ratio * radius_skin`.
Insulation expands radii and length as for `Cylinder`; attachment positions
stay at flesh level.
"""
mutable struct Cone{M,D,B,T} <: AbstractShape
    mass::M
    density::D
    b::B
    top_ratio::T
end
Cone(mass, density, b) = Cone(mass, density, b, 0.0)

# Volume of a frustum = (π/3) · L · (R² + R·r + r²) with r = top_ratio·R.
# With L = 2·b·R: V = (2π/3) · b · R³ · (1 + t + t²).
_cone_top_ratio(s::Cone) = s.top_ratio
_cone_volume_factor(t) = 1 + t + t^2

function _cone_radius(volume, b, t)
    (3 * volume / (2π * b * _cone_volume_factor(t)))^(1/3)
end

function geometry(shape::Cone, ::Naked)
    volume = shape.mass / shape.density
    t = _cone_top_ratio(shape)
    radius_skin = _cone_radius(volume, shape.b, t)
    length_skin = 2 * shape.b * radius_skin
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1/3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin),
                    SurfaceAreas(; total))
end
function geometry(shape::Cone, fur::Fur)
    volume = shape.mass / shape.density
    t = _cone_top_ratio(shape)
    radius_skin = _cone_radius(volume, shape.b, t)
    radius_fur = radius_skin + fur.thickness
    length_skin = 2 * shape.b * radius_skin
    length_fur = length_skin + 2 * fur.thickness
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1/3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fur, length_skin, length_fur),
                    SurfaceAreas(; total, skin, convection))
end
function geometry(shape::Cone, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    t = _cone_top_ratio(shape)
    radius_skin = _cone_radius(volume, shape.b, t)
    length_skin = 2 * shape.b * radius_skin
    radius_flesh = _cone_radius(flesh_volume, shape.b, t)
    fat_thickness = radius_skin - radius_flesh
    total = surface_area(shape, radius_skin, length_skin)
    characteristic_dimension = volume^(1/3)
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, length_skin, fat=fat_thickness),
                    SurfaceAreas(; total))
end
function geometry(shape::Cone, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    t = _cone_top_ratio(shape)
    radius_skin = _cone_radius(volume, shape.b, t)
    radius_fur = radius_skin + fur.thickness
    length_skin = 2 * shape.b * radius_skin
    length_fur = length_skin + 2 * fur.thickness
    radius_flesh = _cone_radius(flesh_volume, shape.b, t)
    fat_thickness = radius_skin - radius_flesh
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = volume^(1/3) + fur.thickness
    return Geometry(volume, characteristic_dimension,
                    (; radius_skin, radius_fur, length_skin, length_fur, fat=fat_thickness),
                    SurfaceAreas(; total, skin, convection))
end

# Aggregate surface area for a frustum: base disc + top disc + slant surface.
# Slant length s = sqrt((R - r)² + L²) where r = t·R.
function surface_area(shape::Cone, R, L)
    t = _cone_top_ratio(shape)
    r = t * R
    s = sqrt((R - r)^2 + L^2)
    π * R^2 + π * r^2 + π * (R + r) * s
end

# Silhouette: same shape as cylinder pattern, but the lateral profile is a
# triangle (r·L) rather than a rectangle (2·r·L).
function silhouette_area(::Cone, r, L, θ)
    r * L * sin(θ) + π * r^2 * cos(θ)
end
function silhouette_area(::Cone, ::Naked, body::AbstractBody, θ)
    silhouette_area(shape(body), body.geometry.length.radius_skin,
                    body.geometry.length.length_skin, θ)
end
function silhouette_area(::Cone, ::Fat, body::AbstractBody, θ)
    silhouette_area(shape(body), body.geometry.length.radius_skin,
                    body.geometry.length.length_skin, θ)
end
function silhouette_area(::Cone, ::Fur, body::AbstractBody, θ)
    silhouette_area(shape(body), body.geometry.length.radius_fur,
                    body.geometry.length.length_fur, θ)
end
function silhouette_area(::Cone, ::CompositeInsulation, body::AbstractBody, θ)
    silhouette_area(shape(body), body.geometry.length.radius_fur,
                    body.geometry.length.length_fur, θ)
end
function silhouette_area(::Cone, ::Union{Naked,Fat}, body::AbstractBody)
    r = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (; normal = r * L, parallel = π * r^2)
end
function silhouette_area(::Cone, ::Union{Fur,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fur
    L = body.geometry.length.length_fur
    (; normal = r * L, parallel = π * r^2)
end

# Radii — mirror Cylinder.
skin_radius(::Cone, ::AbstractInsulation, body) = body.geometry.length.radius_skin

insulation_radius(::Cone, ::Naked, body) = body.geometry.length.radius_skin
flesh_radius(::Cone, ::Naked, body) = body.geometry.length.radius_skin

insulation_radius(::Cone, ::Fur, body) = body.geometry.length.radius_fur
flesh_radius(::Cone, ::Fur, body) = body.geometry.length.radius_skin

insulation_radius(::Cone, ::Fat, body) = body.geometry.length.radius_skin
flesh_radius(::Cone, ::Fat, body) = body.geometry.length.radius_skin - body.geometry.length.fat

insulation_radius(::Cone, ::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(::Cone, ::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# Composition

attachment_surfaces(::Cone) = (:base, :top, :lateral)

# Per-surface area (outer / fur-aware).
_cone_outer_radius(body) = haskey(body.geometry.length, :radius_fur) ?
    body.geometry.length.radius_fur : body.geometry.length.radius_skin
_cone_outer_length(body) = haskey(body.geometry.length, :length_fur) ?
    body.geometry.length.length_fur : body.geometry.length.length_skin

function surface_area(::Cone, body::AbstractBody, ::Val{:base})
    R = _cone_outer_radius(body); π * R^2
end
function surface_area(shape::Cone, body::AbstractBody, ::Val{:top})
    R = _cone_outer_radius(body); π * (shape.top_ratio * R)^2
end
function surface_area(shape::Cone, body::AbstractBody, ::Val{:lateral})
    R = _cone_outer_radius(body); L = _cone_outer_length(body)
    r = shape.top_ratio * R
    s = sqrt((R - r)^2 + L^2)
    π * (R + r) * s
end

function validate_position(::Cone, body::AbstractBody, ::Val{:base}, pos)
    issetequal(keys(pos), (:r, :φ)) || error(":base needs (r, φ); got $(keys(pos))")
    R = body.geometry.length.radius_skin
    pos.r ≥ zero(pos.r) && pos.r ≤ R || error(":base r out of range [0, $R]: $(pos.r)")
end
function validate_position(shape::Cone, body::AbstractBody, ::Val{:top}, pos)
    issetequal(keys(pos), (:r, :φ)) || error(":top needs (r, φ); got $(keys(pos))")
    Rt = shape.top_ratio * body.geometry.length.radius_skin
    Rt > zero(Rt) || error(":top has zero radius (top_ratio=0); use a Disc(0) only")
    pos.r ≥ zero(pos.r) && pos.r ≤ Rt || error(":top r out of range [0, $Rt]: $(pos.r)")
end
function validate_position(::Cone, body::AbstractBody, ::Val{:lateral}, pos)
    issetequal(keys(pos), (:z, :φ)) || error(":lateral needs (z, φ); got $(keys(pos))")
    L = body.geometry.length.length_skin
    pos.z ≥ zero(pos.z) && pos.z ≤ L || error(":lateral z out of range [0, $L]: $(pos.z)")
end

# Attachment positions at flesh (skin) level.
surface_point(::Cone, body::AbstractBody, ::Val{:base}, pos) =
    (pos.r * cos(pos.φ), pos.r * sin(pos.φ), zero(pos.r))
surface_point(::Cone, body::AbstractBody, ::Val{:top}, pos) =
    (pos.r * cos(pos.φ), pos.r * sin(pos.φ), body.geometry.length.length_skin)
function surface_point(shape::Cone, body::AbstractBody, ::Val{:lateral}, pos)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    t = shape.top_ratio
    rz = R * (1 - (1 - t) * pos.z / L)
    (rz * cos(pos.φ), rz * sin(pos.φ), pos.z)
end

surface_normal(::Cone, ::AbstractBody, ::Val{:base}, _) = (0.0, 0.0, -1.0)
surface_normal(::Cone, ::AbstractBody, ::Val{:top}, _)  = (0.0, 0.0,  1.0)
function surface_normal(shape::Cone, body::AbstractBody, ::Val{:lateral}, pos)
    R = ustrip(u"m", body.geometry.length.radius_skin)
    L = ustrip(u"m", body.geometry.length.length_skin)
    t = shape.top_ratio
    # Slant surface tangent direction has Δr = -(R - r) over Δz = L → outward
    # normal is (radial)·(L/s) + (axial)·((R-r)/s).
    Δr = R - t * R
    s = sqrt(Δr^2 + L^2)
    ((L/s) * cos(pos.φ), (L/s) * sin(pos.φ), Δr/s)
end

# Centroids (FullCover).
function surface_centroid(::Cone, body::AbstractBody, ::Val{:base})
    R = body.geometry.length.radius_skin; (zero(R), zero(R), zero(R))
end
function surface_centroid(::Cone, body::AbstractBody, ::Val{:top})
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), zero(R), L)
end
function surface_centroid(shape::Cone, body::AbstractBody, ::Val{:lateral})
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    rz = R * (1 + shape.top_ratio) / 2  # midpoint radius
    (rz, zero(R), L/2)
end
surface_centroid_normal(::Cone, ::AbstractBody, ::Val{:base}) = (0.0, 0.0, -1.0)
surface_centroid_normal(::Cone, ::AbstractBody, ::Val{:top})  = (0.0, 0.0,  1.0)
function surface_centroid_normal(shape::Cone, body::AbstractBody, ::Val{:lateral})
    R = ustrip(u"m", body.geometry.length.radius_skin)
    L = ustrip(u"m", body.geometry.length.length_skin)
    Δr = R - shape.top_ratio * R
    s = sqrt(Δr^2 + L^2)
    (L/s, 0.0, Δr/s)
end
