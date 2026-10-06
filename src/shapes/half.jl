# ── Half shapes ──────────────────────────────────────────────────────────
#
# A `Half{S}` (see the type in geometry.jl) is the full shape `S` of double
# mass, cut through the centre and mirror-joined at a flat face. It inherits
# the parent's dimension math and thermal family; only surface area, the flat
# join face, and the mesh are half-specific. Two `Half`s of equal mass joined
# at their flat faces reconstruct the full shape of double mass.

"""
    HalfCylinder(mass, density, axis_ratio_b) -> Half{<:Cylinder}

Dorsal/ventral half-cylinder: `Half(Cylinder(2mass, density, axis_ratio_b))`.
Axis along `+z`, curved surface in `y ≥ 0`, flat face at `y = 0`.
"""
HalfCylinder(mass, density, axis_ratio_b) = Half(Cylinder(2mass, density, axis_ratio_b))

"""
    HalfCone(mass, density, axis_ratio_b, top_ratio=0.0) -> Half{<:Cone}

Dorsal/ventral half-cone (or half-frustum): `Half(Cone(2mass, density,
axis_ratio_b, top_ratio))`. Axis along `+z` with the base at `z = 0`, curved
surface in `y ≥ 0`, flat trapezoidal face at `y = 0`.
"""
HalfCone(mass, density, axis_ratio_b, top_ratio=0.0) =
    Half(Cone(2mass, density, axis_ratio_b, top_ratio))

"""
    HalfEllipsoid(mass, density, axis_ratio_b, axis_ratio_c) -> Half{<:Ellipsoid}

Dorsal/ventral half-ellipsoid: `Half(Ellipsoid(2mass, density, b, c))`. Long
axis along `+x`, dome in `z ≥ 0`, flat elliptical face at `z = 0`.
"""
HalfEllipsoid(mass, density, axis_ratio_b, axis_ratio_c) =
    Half(Ellipsoid(2mass, density, axis_ratio_b, axis_ratio_c))

"""
    HalfSphere(mass, density) -> Half{<:Sphere}

Hemisphere: `Half(Sphere(2mass, density))`. Wrapping `Sphere` (not an equal-axis
ellipsoid) means the geometry uses the fast spherical closed forms rather than
the eccentricity/asin and fat-cubic of the ellipsoidal path.
"""
HalfSphere(mass, density) = Half(Sphere(2mass, density))

# Halving a box through its centre just gives a thinner box, so a half plate is
# a `Plate` — there's no `Half{<:Plate}`.
Half(p::Plate) = error("a half plate is a plate: use Plate(mass / 2, density, axis_ratio_b, " *
                       "2 * axis_ratio_c) (halved along its height) instead of Half(Plate(...))")

# ── Geometry: delegate dimensions to the parent, override area + volume ───

geometry(h::Half, ins::Naked)        = _halve(h, geometry(h.parent, ins))
geometry(h::Half, fat::FatLayer)     = _halve(h, geometry(h.parent, fat))
geometry(h::Half, fur::FibrousLayer) = _halve(h, geometry(h.parent, fur), fur)
geometry(h::Half, fur::FibrousLayer, fat::FatLayer) =
    _halve(h, geometry(h.parent, fur, fat), fur)

# Assemble a half Geometry from the parent's: same `length`, half the volume.
function _halfgeom(full, total)
    vol = full.volume / 2
    Geometry(vol, full.length, SurfaceAreas(; total))
end
function _halfgeom(full, total, skin, fur::FibrousLayer)
    vol = full.volume / 2
    convection = skin - insulation_area(fur.fibre_diameter, fur.fibre_density, skin)
    Geometry(vol, full.length, SurfaceAreas(; total, skin, convection))
end

# A half's surface is the parent's, halved, plus the cut face — the mirror plane
# through the centre where the two halves join. The parent area is already
# computed (`full.area.*`); only this cut face is family-specific.
# The cylindrical cut is the axial section: a trapezoid with parallel sides
# 2R and 2tR (a rectangle for a cylinder). The domed cut contains the long
# axis: an ellipse with semi-axes a and b.
_cut_face_area(h::Half{<:AbstractCylindrical}, l) = (1 + _top_ratio(h)) * l.radius_skin * l.length_skin
_cut_face_area(::Half{<:AbstractEllipsoidal}, l) = π * l.a_semi_major_skin * l.b_semi_minor_skin
_cut_face_area(::Half{<:AbstractSpherical},   l) = π * l.radius_skin^2

function _halve(h::Half, full)
    cut_face = _cut_face_area(h, full.length)
    _halfgeom(full, full.area.total / 2 + cut_face)
end
function _halve(h::Half, full, fur::FibrousLayer)
    cut_face = _cut_face_area(h, full.length)
    _halfgeom(full, full.area.total / 2 + cut_face, full.area.skin / 2 + cut_face, fur)
end

# ── Route-around-the-wrapper forwards ────────────────────────────────────
# The half's `length` NamedTuple is the parent's, so these read identically.

_skin_radius(h::Half, length)    = _skin_radius(h.parent, length)
_fibrous_radius(h::Half, length) = _fibrous_radius(h.parent, length)
outer_dims(h::Half, body::AbstractBody) = outer_dims(h.parent, body)
_top_ratio(h::Half) = _top_ratio(h.parent)

# ── Silhouette ────────────────────────────────────────────────────────────
#
# A half's shadow depends on where the sun sits relative to the flat face,
# not just on θ. The convention here is the dorsal/ventral one: the sun lies
# in the plane holding the long axis and the flat-face normal, on the dome
# side, θ from the long axis as for the parent. θ = π/2 looks straight down
# on the dome (`normal`, the same shadow as the full shape); θ = 0 looks
# along the axis (`parallel`). In the local frame the sun direction is
# (0, sin θ, cos θ) for a cylindrical half and (cos θ, 0, sin θ) for a domed
# half.
#
# Cauchy's projection formula, A = ½∮|n·d| dA, splits the half's surface into
# its dome and its flat face. For a centrally symmetric parent (cylinder,
# ellipsoid, sphere) the dome carries exactly the parent's share, so
#     A_half = A_parent / 2 + A_flat · |sin θ| / 2.
# A frustum isn't centrally symmetric; its half is handled below.

_outer_cut_face_area(h::Half{<:AbstractCylindrical}, body) =
    (d = outer_dims(h, body); (1 + _top_ratio(h)) * d.r * d.L)
_outer_cut_face_area(h::Half{<:AbstractEllipsoidal}, body) =
    (d = outer_dims(h, body); π * d.a * d.b)
_outer_cut_face_area(h::Half{<:AbstractSpherical}, body) =
    π * insulation_radius(body)^2

silhouette(h::Half, ins::AbstractInsulationLayer, body::AbstractBody, θ) =
    silhouette(h.parent, ins, body, θ) / 2 + _outer_cut_face_area(h, body) * abs(sin(θ)) / 2

# Half-frustum (base radius R at z = 0, top r = tR at z = L). Project onto the
# plane ⊥ the sun and stretch by 1/|cos θ| as for the full frustum: the end
# half-discs become half-circles a distance D = L·|tan θ| apart, bulging the
# same way. Seen from the top end (cos θ ≥ 0) the base half-disc caps a
# trapezoid of the two diameters and the top half-disc sits inside it:
#     π·R²/2 + (R + r)·D.
# Seen from the base end the shadow is the full frustum's hull less the
# base's far half-disc, π·R²/2. Shrinking back by |cos θ| gives the two
# branches below; a cylinder (r = R) gives the same value either way.
function silhouette(h::Half{<:AbstractCylindrical}, ::AbstractInsulationLayer, body::AbstractBody, θ)
    d = outer_dims(h, body)
    s, c = sin(θ), cos(θ)
    if s < 0                               # sun on the flat side: same shadow as -d
        s, c = -s, -c
    end
    R, r = d.r, _top_ratio(h) * d.r
    half_base = π * R^2 * abs(c) / 2
    c >= 0 ? half_base + (R + r) * d.L * s :
             _frustum_silhouette(R, r, d.L, θ) - half_base
end

function silhouette(h::Half, ins::AbstractInsulationLayer, body::AbstractBody)
    (; normal = silhouette(h, ins, body, π / 2), parallel = silhouette(h, ins, body, 0.0))
end

# ── Composition: surfaces are half-specific (extra flat join face) ────────
#
# Attachment positions are anchored at flesh (skin) level so joins meet
# flesh-to-flesh; the flat face is reported at skin level so mixed-insulation
# halves join cleanly under FullCover.

# Cylindrical half (cylinder or cone): EndA, EndB (half discs), Lateral (half
# tube), Flat (axial section). The curved surfaces are the parent's restricted
# to y ≥ 0 (φ ∈ [0, π]), so points, normals and ranges forward to the parent
# and areas are half the parent's; only Flat is half-specific.
attachment_surfaces(::Half{<:AbstractCylindrical}) = (EndA, EndB, Lateral, Flat)

const _CurvedEnd = Union{EndA,EndB,Lateral}

surface_area(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd) =
    surface_area(h.parent, body, loc) / 2
surface_area(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::Flat) =
    _cut_face_area(h, body.geometry.length)

function validate_range(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd)
    validate_range(h.parent, body, loc)
    0 ≤ loc.φ ≤ π || error("$(nameof(typeof(loc))) φ out of range [0, π]: $(loc.φ)")
end
function validate_range(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Flat)  # (u=z, v=x)
    L = body.geometry.length.length_skin
    loc.u ≥ zero(loc.u) && loc.u ≤ L || error("Flat z out of range [0, $L]: $(loc.u)")
    r = _half_radius_at(h, body, loc.u)
    abs(loc.v) ≤ r || error("Flat x out of range ±$r at z = $(loc.u): $(loc.v)")
end

# Skin radius of the axial section at height z (linear from R to tR).
function _half_radius_at(h::Half{<:AbstractCylindrical}, body, z)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    R * (1 - (1 - _top_ratio(h)) * z / L)
end

surface_point(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd) =
    surface_point(h.parent, body, loc)
surface_point(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Flat) =
    (loc.v, zero(loc.v), loc.u)

surface_normal(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd) =
    surface_normal(h.parent, body, loc)
surface_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)

# Centroids sit on the half's symmetry plane x = 0: halfway out the end half
# discs, on the crest of the lateral surface, and mid-length on the flat face.
function surface_centroid(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndA)
    R = body.geometry.length.radius_skin; (zero(R), R / 2, zero(R))
end
function surface_centroid(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndB)
    L = body.geometry.length.length_skin
    (zero(L), _half_radius_at(h, body, L) / 2, L)
end
surface_centroid(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::Lateral) =
    surface_point(h.parent, body, Lateral(body.geometry.length.length_skin / 2, π / 2))
function surface_centroid(::Half{<:AbstractCylindrical}, body::AbstractBody, ::Flat)
    L = body.geometry.length.length_skin; (zero(L), zero(L), L / 2)
end
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndA) = (0.0, 0.0, -1.0)
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndB) = (0.0, 0.0, 1.0)
surface_centroid_normal(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::Lateral) =
    surface_normal(h.parent, body, Lateral(body.geometry.length.length_skin / 2, π / 2))
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::Flat) = (0.0, -1.0, 0.0)

# Domed half (ellipsoidal or spherical): Dome + Flat (elliptical disc at z = 0).
# A sphere is the a=b=c case, so both families share one parametrization; only
# the skin semi-axes are read differently.
const HalfDomed = Half{<:Union{AbstractEllipsoidal,AbstractSpherical}}

_domed_semiaxes(::Half{<:AbstractEllipsoidal}, body) =
    (body.geometry.length.a_semi_major_skin, body.geometry.length.b_semi_minor_skin, body.geometry.length.c_semi_minor_skin)
function _domed_semiaxes(::Half{<:AbstractSpherical}, body)
    r = body.geometry.length.radius_skin; (r, r, r)
end

attachment_surfaces(::HalfDomed) = (Dome, Flat)

# The dome is everything but the (skin-level) flat face.
surface_area(sh::HalfDomed, body::AbstractBody, ::Dome) =
    body.geometry.area.total - surface_area(sh, body, Flat())
surface_area(sh::HalfDomed, body::AbstractBody, ::Flat) =
    _cut_face_area(sh, body.geometry.length)

function validate_range(::HalfDomed, ::AbstractBody, loc::Dome)
    0 ≤ loc.α ≤ π || error("Dome α out of range [0, π]: $(loc.α)")
    0 ≤ loc.β ≤ π || error("Dome β out of range [0, π]: $(loc.β)")
end
function validate_range(sh::HalfDomed, body::AbstractBody, loc::Flat)  # (u=x, v=y)
    a, b, _ = _domed_semiaxes(sh, body)
    (loc.u / a)^2 + (loc.v / b)^2 ≤ 1 + 1e-9 || error("Flat (x,y) outside boundary ellipse")
end

function surface_point(sh::HalfDomed, body::AbstractBody, loc::Dome)
    a, b, c = _domed_semiaxes(sh, body)
    (a * cos(loc.α), b * sin(loc.α) * cos(loc.β), c * sin(loc.α) * sin(loc.β))
end
surface_point(::HalfDomed, body::AbstractBody, loc::Flat) =
    (loc.u, loc.v, zero(loc.u))

function surface_normal(sh::HalfDomed, body::AbstractBody, loc::Dome)
    a, b, c = _domed_semiaxes(sh, body)
    nx = cos(loc.α) * (b * c)
    ny = sin(loc.α) * cos(loc.β) * (a * c)
    nz = sin(loc.α) * sin(loc.β) * (a * b)
    n = sqrt(nx^2 + ny^2 + nz^2)
    (nx / n, ny / n, nz / n)
end
surface_normal(::HalfDomed, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)

function surface_centroid(sh::HalfDomed, body::AbstractBody, ::Dome)
    _, _, c = _domed_semiaxes(sh, body); (zero(c), zero(c), c)
end
function surface_centroid(sh::HalfDomed, body::AbstractBody, ::Flat)
    a, = _domed_semiaxes(sh, body); (zero(a), zero(a), zero(a))
end
surface_centroid_normal(::HalfDomed, ::AbstractBody, ::Dome) = (0.0, 0.0, 1.0)
surface_centroid_normal(::HalfDomed, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)
