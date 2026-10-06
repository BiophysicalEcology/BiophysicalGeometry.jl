# ── Half shapes ──────────────────────────────────────────────────────────
#
# A `Half{S}` (see the type in geometry.jl) is the full shape `S` of double
# mass, cut through the centre and mirror-joined at a flat face. It inherits
# the parent's dimension math and thermal family; only surface area, the flat
# join face, and the mesh are half-specific. Two `Half`s of equal mass joined
# at their flat faces reconstruct the full shape of double mass.

# The constructors take the same keywords as the full shape, but every value
# describes the half itself: its own mass and volume, and for a HalfEllipsoid
# its own height (the dome, half the full ellipsoid's). They are translated to
# the full parent shape here. The stored geometry of a `Half` is the parent's.

_double(kw, key) = haskey(kw, key) ? merge(kw, NamedTuple{(key,)}((2 * kw[key],))) : kw
_full_shape_keywords(kw) = _double(_double(NamedTuple(kw), :mass), :volume)

"""
    HalfCylinder(; mass, density, volume, length, radius, axis_ratio_b) -> Half{<:Cylinder}

Dorsal/ventral half-cylinder, lying along `+x` with its curved surface in
`z ≥ 0` and its flat face at `z = 0`. Keywords are those of [`Cylinder`](@ref),
describing the half: `mass` and `volume` are the half's own.
"""
HalfCylinder(; kw...) = Half(Cylinder(; _full_shape_keywords(kw)...))

"""
    HalfCone(; mass, density, volume, length, radius, axis_ratio_b, top_ratio=0.0) -> Half{<:Cone}

Dorsal/ventral half-cone (or half-frustum), lying along `+x` with its base at
`x = 0`, its curved surface in `z ≥ 0` and its flat trapezoidal face at `z = 0`.
Keywords are those of [`Cone`](@ref), describing the half.
"""
HalfCone(; kw...) = Half(Cone(; _full_shape_keywords(kw)...))

"""
    HalfEllipsoid(; mass, density, volume, length, width, height, axis_ratio_b, axis_ratio_c) -> Half{<:Ellipsoid}

Dorsal/ventral half-ellipsoid: `length` along `x`, `width` along `y`, a dome of
`height` in `z ≥ 0` over a flat elliptical face at `z = 0`. Keywords describe the
half: `height` is the dome's and `axis_ratio_c` is the half's length / height
(each half the full ellipsoid's).
"""
function HalfEllipsoid(; kw...)
    full = _double(_full_shape_keywords(kw), :height)
    if haskey(full, :axis_ratio_c)
        full = merge(full, (; axis_ratio_c = full.axis_ratio_c / 2))
    end
    Half(Ellipsoid(; full...))
end

"""
    HalfSphere(; mass, density, volume, radius) -> Half{<:Sphere}

Hemisphere: a dome of `radius` in `z ≥ 0` over a flat disc at `z = 0`. Wrapping
`Sphere` (not an equal-axis ellipsoid) means the geometry uses the fast
spherical closed forms. Keywords describe the half.
"""
HalfSphere(; kw...) = Half(Sphere(; _full_shape_keywords(kw)...))

# Halving a box through its centre just gives a thinner box, so a half plate is
# a `Plate` — there's no `Half{<:Plate}`.
Half(p::Plate) = error("a half plate is a plate: build the half directly with Plate(; ...) " *
                       "(e.g. half the height) instead of Half(Plate(...))")
Half(p::TriangularPlate) = error("a half triangular plate is a triangular plate: build it " *
                                 "directly with TriangularPlate(; ...) instead of Half(...)")

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
_cut_face_area(::Half{<:AbstractEllipsoidal}, l) = π * (l.length_skin / 2) * (l.width_skin / 2)
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
# along the axis (`parallel`). Every half lies along x with its dome up, so in
# the local frame the sun direction is (cos θ, 0, sin θ).
#
# Cauchy's projection formula, A = ½∮|n·d| dA, splits the half's surface into
# its dome and its flat face. For a centrally symmetric parent (cylinder,
# ellipsoid, sphere) the dome carries exactly the parent's share, so
#     A_half = A_parent / 2 + A_flat · |sin θ| / 2.
# A frustum isn't centrally symmetric; its half is handled below.

_outer_cut_face_area(h::Half{<:AbstractCylindrical}, body) =
    (d = outer_dims(h, body); (1 + _top_ratio(h)) * d.radius * d.length)
_outer_cut_face_area(h::Half{<:AbstractEllipsoidal}, body) =
    (d = outer_dims(h, body); π * (d.length / 2) * (d.width / 2))
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
    if s < 0 # sun on the flat side: same shadow as -d
        s, c = -s, -c
    end
    R, r = d.radius, _top_ratio(h) * d.radius
    half_base = π * R^2 * abs(c) / 2
    c >= 0 ? half_base + (R + r) * d.length * s :
             _frustum_silhouette(R, r, d.length, θ) - half_base
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
# to z ≥ 0 (angle ∈ [0, π]), so points, normals and ranges forward to the parent
# and areas are half the parent's; only Flat is half-specific.
attachment_surfaces(::Half{<:AbstractCylindrical}) = (EndA, EndB, Lateral, Flat)

const _CurvedEnd = Union{EndA,EndB,Lateral}

surface_area(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd) =
    surface_area(h.parent, body, loc) / 2
surface_area(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::Flat) =
    _cut_face_area(h, body.geometry.length)

function validate_range(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd)
    validate_range(h.parent, body, loc)
    0 ≤ loc.angle ≤ π || error("$(nameof(typeof(loc))) angle out of range [0, π]: $(loc.angle)")
end
function validate_range(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Flat)
    L = body.geometry.length.length_skin
    loc.x ≥ zero(loc.x) && loc.x ≤ L || error("Flat x out of range [0, $L]: $(loc.x)")
    r = _half_radius_at(h, body, loc.x)
    abs(loc.y) ≤ r || error("Flat y out of range ±$r at x = $(loc.x): $(loc.y)")
end

# Skin radius of the axial section at position x (linear from R to tR).
function _half_radius_at(h::Half{<:AbstractCylindrical}, body, x)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    R * (1 - (1 - _top_ratio(h)) * x / L)
end

surface_point(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd) =
    surface_point(h.parent, body, loc)
surface_point(::Half{<:AbstractCylindrical}, body::AbstractBody, loc::Flat) =
    (loc.x, loc.y, zero(loc.x))

surface_normal(h::Half{<:AbstractCylindrical}, body::AbstractBody, loc::_CurvedEnd) =
    surface_normal(h.parent, body, loc)
surface_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)

# Centroids sit on the half's symmetry plane y = 0: halfway up the end half
# discs, on the crest of the lateral surface, and mid-length on the flat face.
function surface_centroid(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndA)
    R = body.geometry.length.radius_skin; (zero(R), zero(R), R / 2)
end
function surface_centroid(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::EndB)
    L = body.geometry.length.length_skin
    (L, zero(L), _half_radius_at(h, body, L) / 2)
end
surface_centroid(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::Lateral) =
    surface_point(h.parent, body, Lateral(body.geometry.length.length_skin / 2, π / 2))
function surface_centroid(::Half{<:AbstractCylindrical}, body::AbstractBody, ::Flat)
    L = body.geometry.length.length_skin; (L / 2, zero(L), zero(L))
end
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndA) = (-1.0, 0.0, 0.0)
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::EndB) = (1.0, 0.0, 0.0)
surface_centroid_normal(h::Half{<:AbstractCylindrical}, body::AbstractBody, ::Lateral) =
    surface_normal(h.parent, body, Lateral(body.geometry.length.length_skin / 2, π / 2))
surface_centroid_normal(::Half{<:AbstractCylindrical}, ::AbstractBody, ::Flat) = (0.0, 0.0, -1.0)

# Domed half (ellipsoidal or spherical): Dome + Flat (elliptical disc at z = 0).
# A sphere is the a=b=c case, so both families share one parametrization; only
# the skin semi-axes are read differently.
const HalfDomed = Half{<:Union{AbstractEllipsoidal,AbstractSpherical}}

_domed_semiaxes(::Half{<:AbstractEllipsoidal}, body) =
    _skin_semiaxes(body.geometry.length)
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
    0 ≤ loc.polar ≤ π || error("Dome polar angle out of range [0, π]: $(loc.polar)")
    0 ≤ loc.azimuth ≤ π || error("Dome azimuth out of range [0, π]: $(loc.azimuth)")
end
function validate_range(sh::HalfDomed, body::AbstractBody, loc::Flat)
    a, b, _ = _domed_semiaxes(sh, body)
    (loc.x / a)^2 + (loc.y / b)^2 ≤ 1 + 1e-9 || error("Flat (x, y) outside boundary ellipse")
end

function surface_point(sh::HalfDomed, body::AbstractBody, loc::Dome)
    a, b, c = _domed_semiaxes(sh, body)
    (a * cos(loc.polar), b * sin(loc.polar) * cos(loc.azimuth), c * sin(loc.polar) * sin(loc.azimuth))
end
surface_point(::HalfDomed, body::AbstractBody, loc::Flat) =
    (loc.x, loc.y, zero(loc.x))

function surface_normal(sh::HalfDomed, body::AbstractBody, loc::Dome)
    a, b, c = _domed_semiaxes(sh, body)
    nx = cos(loc.polar) * (b * c)
    ny = sin(loc.polar) * cos(loc.azimuth) * (a * c)
    nz = sin(loc.polar) * sin(loc.azimuth) * (a * b)
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
