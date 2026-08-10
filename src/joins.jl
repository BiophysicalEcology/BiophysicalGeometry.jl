# Join geometry accessors and internal conduction-path lengths.
#
# These expose, for thermal consumers (HeatExchange's per-join couplings),
# the three facts a coupling needs about a `Join` — its interface area, its
# world-frame location, and the two part names it connects — plus the
# centroid-to-interface conduction path length within each part
# (`internal_distance`), used to derive series conductive resistance across
# a join from each part's own flesh conductivity.
#
# All are pure geometry: they read the poses and dimensions already computed
# at `CompositeBody` construction and never touch physics.

# ── Join accessors ────────────────────────────────────────────────────────

"""
    join_partners(join) -> (Symbol, Symbol)

The `(parent, child)` part names a `Join` connects. Lifted out of the
`Join`'s type parameters, so constant-folded.
"""
join_partners(j::Join) = (_parent(j), _child(j))

"""
    join_area(join, body) -> Area

Interface area of the join — the contact-patch area shared by its two
parts. `π·radius²` for a `Disc` attachment, the full surface area for a
`FullCover`. Parent and child patch areas are equal by construction
(validated when the `CompositeBody` is built), so the parent side is
returned.
"""
join_area(j::Join, body::CompositeBody) =
    patch_area(getfield(body.parts, _parent(j)), j.parent_attachment)

"""
    join_position(join, body) -> NTuple{3,Length}

World-frame position of the join interface centre, obtained by applying the
parent part's pose to the parent attachment point.
"""
function join_position(j::Join, body::CompositeBody)
    part = getfield(body.parts, _parent(j))
    pose = getfield(body.poses, _parent(j))
    p_local = _attach_point(shape(part), part, j.parent_attachment)
    return apply_pose(pose, p_local)
end

# ── Internal conduction-path length ───────────────────────────────────────

"""
    flesh_centroid(shape, body) -> NTuple{3,Length}

Volumetric (flesh) centroid of a part in its local frame. Closed-form per
shape family. Shapes with no method (the animal shapes) error here rather
than silently using the wrong path length for a lumped-resistance solve.
"""
function flesh_centroid(sh::AbstractShape, body::AbstractBody)
    error("flesh_centroid not defined for $(typeof(sh)) — no lumped conduction path")
end

# Axial shapes: centroid on the axis. Cylinder is centred at mid-length.
function flesh_centroid(::Cylinder, body::AbstractBody)
    z = body.geometry.length.length_skin
    (zero(z), zero(z), z / 2)
end
# Cone (frustum): axial centroid z̄ = (L/4)(1 + 2t + 3t²)/(1 + t + t²),
# with t the top/base radius ratio (base at z=0).
function flesh_centroid(sh::Cone, body::AbstractBody)
    L = body.geometry.length.length_skin
    t = sh.top_ratio
    (zero(L), zero(L), (L / 4) * (1 + 2t + 3t^2) / (1 + t + t^2))
end
# Sphere / Ellipsoid: centred at the origin.
function flesh_centroid(::Sphere, body::AbstractBody)
    z = body.geometry.length.radius_skin
    (zero(z), zero(z), zero(z))
end
function flesh_centroid(::Ellipsoid, body::AbstractBody)
    z = body.geometry.length.a_semi_major_skin
    (zero(z), zero(z), zero(z))
end
# Cylindrical half: half-disc cross-section, bulge along +y; centroid of a
# half-disc sits at 4R/(3π) off the flat plane, at mid-length.
function flesh_centroid(::Half{<:AbstractCylindrical}, body::AbstractBody)
    R = body.geometry.length.radius_skin
    L = body.geometry.length.length_skin
    (zero(R), 4R / (3π), L / 2)
end
# Domed half (ellipsoidal/spherical): dome along +z from the flat face; centroid
# sits at 3c/8 off the flat plane (hemisphere centroid generalised).
function flesh_centroid(sh::HalfDomed, body::AbstractBody)
    _, _, c = _domed_semiaxes(sh, body)
    (zero(c), zero(c), 3c / 8)
end
# Plate: centred at the origin.
function flesh_centroid(::Plate, body::AbstractBody)
    z = body.geometry.length.height_skin
    (zero(z), zero(z), zero(z))
end

"""
    internal_distance(body, location) -> Length
    internal_distance(shape, body, location) -> Length

Straight-line path length from a part's flesh centroid to the centroid of
the named surface `location`. Reduces to the expected closed forms —
`radius` for a lateral surface, half-length for an end cap, a semi-axis for
a pole — and is the conduction length used when deriving a
`ConductiveCoupling`'s series resistance.
"""
internal_distance(body::AbstractBody, loc::AbstractSurface) =
    internal_distance(shape(body), body, loc)

function internal_distance(sh::AbstractShape, body::AbstractBody, loc::AbstractSurface)
    c = flesh_centroid(sh, body)
    s = surface_centroid(sh, body, loc)
    return sqrt((c[1] - s[1])^2 + (c[2] - s[2])^2 + (c[3] - s[3])^2)
end

"""
    internal_distance(join, body) -> (; parent, child)

The two centroid-to-interface conduction path lengths for a join: the
parent part's centroid to the parent attachment surface, and likewise for
the child. These feed the series-resistance derivation of a
`ConductiveCoupling`.
"""
function internal_distance(j::Join, body::CompositeBody)
    p = getfield(body.parts, _parent(j))
    c = getfield(body.parts, _child(j))
    (; parent = internal_distance(p, j.parent_attachment.location),
       child  = internal_distance(c, j.child_attachment.location))
end
