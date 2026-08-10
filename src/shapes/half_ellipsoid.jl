"""
    HalfEllipsoid <: AbstractShape

A half-ellipsoid organism shape — a prolate ellipsoid (b = c) cut by the
plane perpendicular to the minor axis through the centre. Dimensions are
sized so that two half-ellipsoids of equal mass joined at their flat faces
reconstruct the full `Ellipsoid(2*mass, density, b, c)`.

In local coordinates the long axis lies along `+x`, and the flat face is
the elliptical disc at `z = 0` with the dome in `z ≥ 0`.

The flat-face area is reported at *skin* level even when fur is present,
so two halves with mixed insulation can join cleanly with `FullCover`.
"""
mutable struct HalfEllipsoid{M,D,B,C} <: AbstractEllipsoidal
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
end

"""
    HalfSphere(mass, density) -> HalfEllipsoid

A hemisphere — the dome half of a sphere cut through its centre. A half-sphere is
exactly the equal-axis half-ellipsoid, so this is a convenience constructor for
`HalfEllipsoid(mass, density, 1.0, 1.0)` (the dome-area formula already carries an
exact half-sphere branch), reusing all of `HalfEllipsoid`'s geometry, silhouette,
meshing and composition machinery rather than duplicating it.

Sized so that two `HalfSphere(mass, density)` joined at their flat faces reconstruct
`Sphere(2*mass, density)`: a `HalfSphere(m/2, ρ)` has the same skin radius as
`Sphere(m, ρ)`, and its two domes sum to the full sphere's surface area. This is the
sphere counterpart of `HalfCylinder`, letting a spherical body take the genuine
two-part (dorsal/ventral) decomposition.
"""
HalfSphere(mass, density) = HalfEllipsoid(mass, density, 1.0, 1.0)

# Dome area for one half of a prolate spheroid (b = c assumed). A dome is
# exactly half the full spheroid's surface, so this reuses `_prolate_area`
# (which already carries the ustrip / asin dance and the sphere-limit branch)
# rather than repeating the eccentricity math.
_half_dome_area(a, b, c) = _prolate_area(a, b, c) / 2

# A half of mass m is the full `Ellipsoid` of mass 2m cut through the centre,
# so all radius/flesh/fat/fur math is inherited from `Ellipsoid` (single source
# of truth) via `full`; the half only differs in area (half dome + flat face,
# at skin level so mixed-insulation halves join under `FullCover`) and volume.
_full_ellipsoid(sh::HalfEllipsoid) = Ellipsoid(2 * sh.mass, sh.density, sh.axis_ratio_b, sh.axis_ratio_c)

geometry(sh::HalfEllipsoid, ins::Naked)        = _half_geometry(sh, geometry(_full_ellipsoid(sh), ins))
geometry(sh::HalfEllipsoid, fat::FatLayer)     = _half_geometry(sh, geometry(_full_ellipsoid(sh), fat))
geometry(sh::HalfEllipsoid, fur::FibrousLayer) = _half_geometry(sh, geometry(_full_ellipsoid(sh), fur), fur)
geometry(sh::HalfEllipsoid, fur::FibrousLayer, fat::FatLayer) =
    _half_geometry(sh, geometry(_full_ellipsoid(sh), fur, fat), fur)

function _half_geometry(sh::HalfEllipsoid, full)
    len = full.length
    flat = π * len.b_semi_minor_skin * len.c_semi_minor_skin
    total = _half_dome_area(len.a_semi_major_skin, len.b_semi_minor_skin, len.c_semi_minor_skin) + flat
    vol = sh.mass / sh.density
    Geometry(vol, vol^(1/3), len, SurfaceAreas(; total))
end
function _half_geometry(sh::HalfEllipsoid, full, fur::FibrousLayer)
    len = full.length
    flat = π * len.b_semi_minor_skin * len.c_semi_minor_skin
    total = _half_dome_area(len.a_semi_major_fibrous, len.b_semi_minor_fibrous, len.c_semi_minor_fibrous) + flat
    skin = _half_dome_area(len.a_semi_major_skin, len.b_semi_minor_skin, len.c_semi_minor_skin) + flat
    convection = skin - insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    vol = sh.mass / sh.density
    Geometry(vol, vol^(1/3) + fur.thickness, len, SurfaceAreas(; total, skin, convection))
end

# Silhouette area: half of full ellipsoid's, so two halves reproduce the
# full. For an ellipsoid (a, b, c) the silhouette projected along direction
# d is an ellipse of area π·sqrt(b²c²·d_x² + a²c²·d_y² + a²b²·d_z²). With
# the sun in the equatorial plane at angle θ from the long (a) axis,
# d = (cos θ, sin θ, 0) gives π·c·sqrt(b²·cos²θ + a²·sin²θ).
function silhouette_area(::HalfEllipsoid, a, b, c, θ)
    π * c * sqrt(b^2 * cos(θ)^2 + a^2 * sin(θ)^2) / 2
end

function silhouette_area(sh::HalfEllipsoid, ::AbstractInsulationLayer, body::AbstractBody)
    d = outer_dims(sh, body)
    (; normal = π * d.a * d.b / 2, parallel = π * d.b * d.c / 2)
end
function silhouette_area(sh::HalfEllipsoid, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(sh, body)
    silhouette_area(sh, d.a, d.b, d.c, θ)
end

# Radii come from the shared `AbstractEllipsoidal` dispatch in ellipsoid.jl.

# Composition

attachment_surfaces(::HalfEllipsoid) = (Dome, Flat)

# Outer (insulation-aware) semi-axes of the dome. Insulation-dispatched.
outer_dims(sh::HalfEllipsoid, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::HalfEllipsoid, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_skin,
     b = body.geometry.length.b_semi_minor_skin,
     c = body.geometry.length.c_semi_minor_skin)
outer_dims(::HalfEllipsoid, ::FibrousLayer, body::AbstractBody) =
    (a = body.geometry.length.a_semi_major_fibrous,
     b = body.geometry.length.b_semi_minor_fibrous,
     c = body.geometry.length.c_semi_minor_fibrous)

# Skin-level semi-axes — drives flesh-anchored attachment positions.
function _halfellip_skin(body::AbstractBody)
    gl = body.geometry.length
    (gl.a_semi_major_skin, gl.b_semi_minor_skin, gl.c_semi_minor_skin)
end

# Loose bound for dome attachments.
surface_area(::HalfEllipsoid, body::AbstractBody, ::Dome) =
    body.geometry.area.total

# Flat at skin level so two halves with mixed insulation join cleanly.
function surface_area(::HalfEllipsoid, body::AbstractBody, ::Flat)
    π * body.geometry.length.b_semi_minor_skin * body.geometry.length.c_semi_minor_skin
end

function validate_range(::HalfEllipsoid, ::AbstractBody, loc::Dome)
    0 ≤ loc.α ≤ π || error("Dome α out of range [0, π]: $(loc.α)")
    0 ≤ loc.β ≤ π || error("Dome β out of range [0, π]: $(loc.β)")
end
# For HalfEllipsoid, Flat uses (u=x, v=y).
function validate_range(::HalfEllipsoid, body::AbstractBody, loc::Flat)
    a = body.geometry.length.a_semi_major_skin
    b = body.geometry.length.b_semi_minor_skin
    abs(loc.u) ≤ a || error("Flat x out of range ±$a: $(loc.u)")
    abs(loc.v) ≤ b || error("Flat y out of range ±$b: $(loc.v)")
    rel = (loc.u / a)^2 + (loc.v / b)^2
    rel ≤ 1 + 1e-9 || error("Flat (x,y) outside ellipse boundary")
end

function surface_point(::HalfEllipsoid, body::AbstractBody, loc::Dome)
    a, b, c = _halfellip_skin(body)
    (a * cos(loc.α), b * sin(loc.α) * cos(loc.β), c * sin(loc.α) * sin(loc.β))
end
function surface_point(::HalfEllipsoid, body::AbstractBody, loc::Flat)
    (loc.u, loc.v, zero(loc.u))
end

function surface_normal(::HalfEllipsoid, body::AbstractBody, loc::Dome)
    a, b, c = _halfellip_skin(body)
    # Gradient of (x/a)² + (y/b)² + (z/c)² = 1 is (2x/a², 2y/b², 2z/c²).
    # Scale each component by a common length (b·c, a·c, a·b) so units
    # cancel: multiply top and bottom by a·b·c to get a unitless vector.
    nx = cos(loc.α) * (b * c)
    ny = sin(loc.α) * cos(loc.β) * (a * c)
    nz = sin(loc.α) * sin(loc.β) * (a * b)
    n = sqrt(nx^2 + ny^2 + nz^2)
    (nx / n, ny / n, nz / n)
end
surface_normal(::HalfEllipsoid, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)

# Centroids — flesh level.
function surface_centroid(::HalfEllipsoid, body::AbstractBody, ::Dome)
    a, _, c = _halfellip_skin(body)
    (zero(a), zero(a), c)  # top of the dome
end
function surface_centroid(::HalfEllipsoid, body::AbstractBody, ::Flat)
    a = body.geometry.length.a_semi_major_skin
    (zero(a), zero(a), zero(a))
end
surface_centroid_normal(::HalfEllipsoid, ::AbstractBody, ::Dome) = (0.0, 0.0, 1.0)
surface_centroid_normal(::HalfEllipsoid, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)
