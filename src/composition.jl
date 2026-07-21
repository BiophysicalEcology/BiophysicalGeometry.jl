# Geometry composition for multi-part organisms.
#
# A `CompositeBody` is a tree of `Body` parts joined at named surfaces.
# Each `Join` carries an `Attachment` on each side: a surface location
# (a singleton subtype of `AbstractSurface`, optionally carrying a
# parametric point on that surface) and the patch shape hidden by the join.
# Surface-area accounting subtracts each patch from the part it's on.
# World-frame poses for every part are derived once at construction by
# walking the join tree from `root`, so plotting falls out automatically.
#
# No Symbols and no Dicts appear in the composition machinery: part
# identifiers are user-defined singleton types, surface identifiers are
# `AbstractSurface` singletons, and the pose/coverage accumulators are
# Tuples of Pairs — all fully type-stable and autodiff-friendly.

# ── Surface singletons ────────────────────────────────────────────────────
#
# Each surface identifier is a singleton subtype of `AbstractSurface`.
# The "point on surface" flavour parametrises the same type with concrete
# coordinate types; the "whole surface" flavour uses `Nothing` fields
# (used by `FullCover` attachments, which don't need a point).

"""
    AbstractSurface

Supertype for surface identifiers used in `Attachment`. Every concrete
subtype has zero-arg constructors that produce the *bare* (whole-surface)
form used with `FullCover`, and positional constructors that produce the
*located* form used with `Disc`.
"""
abstract type AbstractSurface end

"""
    EndA(r=nothing, φ=nothing) <: AbstractSurface

Disc-shaped end cap on the negative-axis side of an axial shape
(Cylinder, HalfCylinder, Cone). Located form uses polar coordinates
`(r, φ)` in the disc's local frame.
"""
struct EndA{R,P} <: AbstractSurface
    r::R
    φ::P
end
EndA() = EndA(nothing, nothing)

"""
    EndB(r=nothing, φ=nothing) <: AbstractSurface

Disc-shaped end cap on the positive-axis side of an axial shape.
"""
struct EndB{R,P} <: AbstractSurface
    r::R
    φ::P
end
EndB() = EndB(nothing, nothing)

"""
    Lateral(z=nothing, φ=nothing) <: AbstractSurface

Curved side surface of an axial shape (Cylinder, HalfCylinder, Cone).
`(z, φ)` are cylindrical coordinates along and around the axis.
"""
struct Lateral{Z,P} <: AbstractSurface
    z::Z
    φ::P
end
Lateral() = Lateral(nothing, nothing)

"""
    Flat(u=nothing, v=nothing) <: AbstractSurface

Flat (cut-plane) face of a half-shape. Coordinate meaning is shape-dependent:
for `HalfCylinder`, `(u, v) = (z, x)`; for `HalfEllipsoid`, `(u, v) = (x, y)`.
"""
struct Flat{U,V} <: AbstractSurface
    u::U
    v::V
end
Flat() = Flat(nothing, nothing)

"""
    Dome(α=nothing, β=nothing) <: AbstractSurface

Curved (dome) surface of a `HalfEllipsoid`. `(α, β)` are spheroidal
angles: `α` measures from the long-axis pole, `β` around the axis.
"""
struct Dome{A,B} <: AbstractSurface
    α::A
    β::B
end
Dome() = Dome(nothing, nothing)

"""
    PoleA(r=nothing, φ=nothing) <: AbstractSurface

Pole on the positive-x end of an `Ellipsoid`. For an untruncated
ellipsoid, the located form is not used (the pole is a point). For a
truncated pole (`pole_a_truncation > 0`), the pole becomes a disc and
`(r, φ)` selects a point on it.
"""
struct PoleA{R,P} <: AbstractSurface
    r::R
    φ::P
end
PoleA() = PoleA(nothing, nothing)

"""
    PoleB() <: AbstractSurface

Pole on the negative-x end of an `Ellipsoid`. Untruncated; the located
form is not meaningful.
"""
struct PoleB <: AbstractSurface end

"""
    Equator(φ=nothing) <: AbstractSurface

Ring at the equator of an `Ellipsoid`. Located form gives the angular
coordinate `φ` around the long axis.
"""
struct Equator{P} <: AbstractSurface
    φ::P
end
Equator() = Equator(nothing)

"""
    Radial(θ=nothing, φ=nothing) <: AbstractSurface

Any point on the surface of a `Sphere`, given as spherical angles
`(θ, φ)`.
"""
struct Radial{T,P} <: AbstractSurface
    θ::T
    φ::P
end
Radial() = Radial(nothing, nothing)

"""
    Top(x=nothing, y=nothing) <: AbstractSurface

Top face of a `Plate` (z = +H/2), located by `(x, y)`.
"""
struct Top{X,Y} <: AbstractSurface
    x::X
    y::Y
end
Top() = Top(nothing, nothing)

"""
    Bottom(x=nothing, y=nothing) <: AbstractSurface

Bottom face of a `Plate` (z = -H/2), located by `(x, y)`.
"""
struct Bottom{X,Y} <: AbstractSurface
    x::X
    y::Y
end
Bottom() = Bottom(nothing, nothing)

"""
    SideA(y=nothing, z=nothing) <: AbstractSurface

+x face of a `Plate`, located by `(y, z)`.
"""
struct SideA{Y,Z} <: AbstractSurface
    y::Y
    z::Z
end
SideA() = SideA(nothing, nothing)

"""
    SideB(y=nothing, z=nothing) <: AbstractSurface

-x face of a `Plate`, located by `(y, z)`.
"""
struct SideB{Y,Z} <: AbstractSurface
    y::Y
    z::Z
end
SideB() = SideB(nothing, nothing)

"""
    SideC(x=nothing, z=nothing) <: AbstractSurface

+y face of a `Plate`, located by `(x, z)`.
"""
struct SideC{X,Z} <: AbstractSurface
    x::X
    z::Z
end
SideC() = SideC(nothing, nothing)

"""
    SideD(x=nothing, z=nothing) <: AbstractSurface

-y face of a `Plate`, located by `(x, z)`.
"""
struct SideD{X,Z} <: AbstractSurface
    x::X
    z::Z
end
SideD() = SideD(nothing, nothing)

# ── Attachment shapes ─────────────────────────────────────────────────────

abstract type AbstractAttachmentShape end

"""
    Disc{R} <: AbstractAttachmentShape

A flat circular contact patch of radius `R`. Used for legs into bodies,
head poles into cylinder ends, etc.
"""
struct Disc{R} <: AbstractAttachmentShape
    radius::R
end

"""
    FullCover <: AbstractAttachmentShape

Covers the whole named surface. Used for half-shape dorsal/ventral joins
where two flat faces match exactly.
"""
struct FullCover <: AbstractAttachmentShape end

"""
    Attachment(location, shape)

One side of a `Join`. `location` is an `AbstractSurface` instance (either
the bare form like `Lateral()` used with `FullCover`, or a located form
like `Lateral(z, φ)` used with `Disc`). `shape` is the patch shape
(`Disc` or `FullCover`).
"""
struct Attachment{L<:AbstractSurface, S<:AbstractAttachmentShape}
    location::L
    shape::S
end

# ── Join ──────────────────────────────────────────────────────────────────

"""
    Join(parent_tag, parent_attachment, child_tag, child_attachment; twist=0.0)

A connection between two parts of a `CompositeBody`, identified by
singleton tag types provided by the user (e.g. user defines
`struct Torso end; struct Head end` and passes `Torso()`, `Head()` here).
`twist` (radians) sets the rotation about the joint axis — the 6th DOF
that two anti-aligned surface normals don't fix.
"""
struct Join{PT, CT, A1<:Attachment, A2<:Attachment, T<:Real}
    parent::PT
    parent_attachment::A1
    child::CT
    child_attachment::A2
    twist::T
end
Join(parent, pa::Attachment, child, ca::Attachment; twist=0.0) =
    Join(parent, pa, child, ca, twist)

# Reverse a join (parent/child swapped, twist negated).
_reverse_join(j::Join) =
    Join(j.child, j.child_attachment, j.parent, j.parent_attachment, -j.twist)

# ── Pose ──────────────────────────────────────────────────────────────────

"""
    Pose{T,M}

World-frame pose of a part: a translation (3-tuple of length quantities)
and a 3×3 rotation matrix (dimensionless). Rotation matrices are never
mutated after construction, so the identity matrix is shared.
"""
struct Pose{T,M<:AbstractMatrix}
    translation::NTuple{3,T}
    rotation::M
end

const _IDENTITY_ROTATION = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

identity_pose(::Type{T}) where {T} =
    Pose((zero(T), zero(T), zero(T)), _IDENTITY_ROTATION)

"""
    apply_pose(pose, point) -> NTuple{3,Length}

Transform a local point `(x,y,z)` to world coordinates.
"""
function apply_pose(p::Pose, v::NTuple{3})
    R = p.rotation
    x, y, z = v
    (R[1,1]*x + R[1,2]*y + R[1,3]*z + p.translation[1],
     R[2,1]*x + R[2,2]*y + R[2,3]*z + p.translation[2],
     R[3,1]*x + R[3,2]*y + R[3,3]*z + p.translation[3])
end

"""
    apply_rotation(R, v) -> NTuple{3}

Rotate a (dimensional or dimensionless) 3-vector by matrix `R`.
"""
function apply_rotation(R::AbstractMatrix, v::NTuple{3})
    (R[1,1]*v[1] + R[1,2]*v[2] + R[1,3]*v[3],
     R[2,1]*v[1] + R[2,2]*v[2] + R[2,3]*v[3],
     R[3,1]*v[1] + R[3,2]*v[2] + R[3,3]*v[3])
end

# Rodrigues' formula for an axis-angle rotation matrix.
function rotation_axis_angle(axis::NTuple{3,<:Real}, θ::Real)
    c = cos(θ); s = sin(θ); t = 1 - c
    x, y, z = axis
    [t*x*x + c     t*x*y - s*z   t*x*z + s*y;
     t*x*y + s*z   t*y*y + c     t*y*z - s*x;
     t*x*z - s*y   t*y*z + s*x   t*z*z + c]
end

# Rotation matrix that takes unit vector `a` to unit vector `b`.
function rotation_align(a::NTuple{3,<:Real}, b::NTuple{3,<:Real})
    d = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
    if d > 1.0 - 1e-12
        return _IDENTITY_ROTATION
    elseif d < -1.0 + 1e-12
        # 180° rotation; pick any axis ⟂ a
        ax = abs(a[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
        proj = a[1]*ax[1] + a[2]*ax[2] + a[3]*ax[3]
        u = (ax[1] - proj*a[1], ax[2] - proj*a[2], ax[3] - proj*a[3])
        n = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        return rotation_axis_angle((u[1]/n, u[2]/n, u[3]/n), π)
    else
        c = (a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
        cn = sqrt(c[1]^2 + c[2]^2 + c[3]^2)
        return rotation_axis_angle((c[1]/cn, c[2]/cn, c[3]/cn), acos(d))
    end
end

# ── Per-shape interface (defaults) ────────────────────────────────────────
#
# Shapes opt into composition by overriding these for each surface type
# they support. The default `attachment_surfaces` is empty, meaning a
# shape cannot be joined into a `CompositeBody`. The other methods throw a
# clear error if composition code reaches them on an unsupported
# (shape, surface) pair.

"""
    attachment_surfaces(shape) -> Tuple of AbstractSurface subtypes

Types (not instances) of surfaces on `shape` that can be used in a
`Join`. Default `()` means the shape cannot be joined.
"""
attachment_surfaces(::AbstractShape) = ()

"""
    surface_area(shape, body, location::AbstractSurface) -> Area

Area of the named surface alone (one face / one side of `shape`).
The `location`'s parametric fields (if any) are ignored — only its type
matters.
"""
function surface_area(shape::AbstractShape, body::AbstractBody, loc::AbstractSurface)
    error("surface_area not defined for $(typeof(shape)) surface $(typeof(loc))")
end

"""
    surface_point(shape, body, location) -> NTuple{3,Length}

Local 3D point on the named surface at the given located `location`.
"""
function surface_point(shape::AbstractShape, body::AbstractBody, loc::AbstractSurface)
    error("surface_point not defined for $(typeof(shape)) surface $(typeof(loc))")
end

"""
    surface_normal(shape, body, location) -> NTuple{3,Float64}

Local outward unit normal at the located `location`.
"""
function surface_normal(shape::AbstractShape, body::AbstractBody, loc::AbstractSurface)
    error("surface_normal not defined for $(typeof(shape)) surface $(typeof(loc))")
end

"""
    validate_range(shape, body, location)

Throw if the location's parametric coordinates are out of range for the
(shape, surface) pair. Default: accept anything (override per shape).
Field-set validation is unnecessary — the surface type fixes the fields.
"""
validate_range(::AbstractShape, ::AbstractBody, ::AbstractSurface) = nothing

"""
    surface_centroid(shape, body, ::S) where {S<:AbstractSurface} -> NTuple{3,Length}

Local 3D centroid of the surface `S`. Used by `FullCover` attachments,
which have no parametric point.
"""
function surface_centroid(shape::AbstractShape, body::AbstractBody, loc::AbstractSurface)
    error("surface_centroid not defined for $(typeof(shape)) surface $(typeof(loc))")
end

"""
    surface_centroid_normal(shape, body, ::S) where {S<:AbstractSurface} -> NTuple{3,Float64}

Local outward unit normal at the surface centroid. Used by `FullCover`
attachments.
"""
function surface_centroid_normal(shape::AbstractShape, body::AbstractBody, loc::AbstractSurface)
    error("surface_centroid_normal not defined for $(typeof(shape)) surface $(typeof(loc))")
end

# ── Patch area dispatch ───────────────────────────────────────────────────

patch_area(body::AbstractBody, att::Attachment{<:AbstractSurface, <:Disc}) =
    π * att.shape.radius^2

patch_area(body::AbstractBody, att::Attachment{<:AbstractSurface, FullCover}) =
    surface_area(shape(body), body, att.location)

# ── Validation ────────────────────────────────────────────────────────────

# A shape supports a surface type if any element of `attachment_surfaces`
# is that surface type (compared by `isa`, so `Lateral(z, φ) isa Lateral`
# works). Tuple iteration is unrolled by the compiler for the small tuples
# used here, so this reduces to a compile-time boolean.
_supports_surface(::Tuple{}, ::AbstractSurface) = false
_supports_surface(surfaces::Tuple, loc::AbstractSurface) =
    (loc isa surfaces[1]) || _supports_surface(Base.tail(surfaces), loc)

function validate_attachment(body::AbstractBody, att::Attachment)
    sh = shape(body)
    surfaces = attachment_surfaces(sh)
    isempty(surfaces) &&
        error("$(typeof(sh)) does not support being joined (attachment_surfaces is empty)")
    _supports_surface(surfaces, att.location) ||
        error("$(typeof(sh)) has no surface $(typeof(att.location)); valid: $surfaces")
    # FullCover has no parametric point; skip range validation.
    if !(att.shape isa FullCover)
        validate_range(sh, body, att.location)
    end
    Asurface = surface_area(sh, body, att.location)
    Apatch = patch_area(body, att)
    if Apatch > Asurface * (1 + 1e-9)
        error("attachment patch area ($Apatch) exceeds surface area of " *
              "$(typeof(att.location)) on $(typeof(sh)) ($Asurface)")
    end
    return nothing
end

# ── Typed parts container helpers ─────────────────────────────────────────
#
# `parts` is a Tuple of `Pair{Tag, Body}`; `poses` is a Tuple of
# `Pair{Tag, Pose}`. Tags are user-defined singleton types. Lookup by tag
# instance uses a generated function so the position is resolved at compile
# time — no Symbol, no Dict, no runtime search.

"""
    _lookup(entries::Tuple, tag) -> value

Return `entries[i].second` where `entries[i].first isa typeof(tag)`.
Type-stable via `@generated`.
"""
@generated function _lookup(entries::Tuple, ::T) where {T}
    for (i, ptype) in enumerate(entries.parameters)
        if ptype <: Pair && T <: ptype.parameters[1]
            return :(entries[$i].second)
        end
    end
    :(error("no entry for tag ", $T))
end

@generated function _has_tag(entries::Tuple, ::Type{T}) where {T}
    for ptype in entries.parameters
        if ptype <: Pair && T === ptype.parameters[1]
            return :(true)
        end
    end
    :(false)
end

# Extract all first-elements of a Tuple of Pairs.
_tags(::Tuple{}) = ()
_tags(entries::Tuple) = (entries[1].first, _tags(Base.tail(entries))...)

function validate_parts(parts::Tuple)
    isempty(parts) && error("CompositeBody must have at least one part")
    for pair in parts
        pair isa Pair || error("parts must be Pairs (tag => Body); got $(typeof(pair))")
        pair.second isa AbstractBody ||
            error("part value must be AbstractBody; got $(typeof(pair.second))")
        pair.second isa CompositeBody &&
            error("nested CompositeBody is not supported")
    end
    # No duplicate tag types.
    _check_unique_tags(map(p -> typeof(p.first), parts))
    return nothing
end

function _check_unique_tags(tag_types::Tuple)
    seen = Set{DataType}()
    for T in tag_types
        T in seen && error("duplicate part tag type $T in parts")
        push!(seen, T)
    end
end

function validate_join(parts::Tuple, j::Join)
    _has_tag(parts, typeof(j.parent)) ||
        error("CompositeBody has no part with tag $(typeof(j.parent))")
    _has_tag(parts, typeof(j.child)) ||
        error("CompositeBody has no part with tag $(typeof(j.child))")
    typeof(j.parent) === typeof(j.child) &&
        error("Join cannot connect a part to itself ($(typeof(j.parent)))")
    p = _lookup(parts, j.parent); c = _lookup(parts, j.child)
    validate_attachment(p, j.parent_attachment)
    validate_attachment(c, j.child_attachment)
    Ap = patch_area(p, j.parent_attachment)
    Ac = patch_area(c, j.child_attachment)
    rel = abs(Ap - Ac) / max(Ap, Ac)
    if rel > 1e-6
        ps, cs = j.parent_attachment.shape, j.child_attachment.shape
        if ps isa Disc && cs isa Disc
            error("Join Disc radii must match: parent r=$(ps.radius), child r=$(cs.radius)")
        elseif ps isa FullCover && cs isa FullCover
            error("Join FullCover surfaces have unequal area: " *
                  "$(typeof(j.parent_attachment.location))=$Ap vs " *
                  "$(typeof(j.child_attachment.location))=$Ac")
        else
            error("Join patch areas differ: parent $(typeof(ps))=$Ap vs child $(typeof(cs))=$Ac")
        end
    end
    return nothing
end

# ── Covered area accumulator (Dict-free) ──────────────────────────────────

covered_areas(parts::Tuple, joins::Tuple) =
    map(parts) do pair
        tag = pair.first
        body = pair.second
        A0 = zero(body.geometry.area.total)
        pair.first => _sum_patches_for_tag(A0, body, tag, joins)
    end

_sum_patches_for_tag(acc, body, tag, ::Tuple{}) = acc
function _sum_patches_for_tag(acc, body, tag, joins::Tuple)
    j = joins[1]
    acc2 = acc
    if typeof(j.parent) === typeof(tag)
        acc2 += patch_area(body, j.parent_attachment)
    end
    if typeof(j.child) === typeof(tag)
        acc2 += patch_area(body, j.child_attachment)
    end
    _sum_patches_for_tag(acc2, body, tag, Base.tail(joins))
end

# ── Pose tree solver (Dict-free, tuple-recursive) ─────────────────────────
#
# Each child's world pose is determined by its parent's world pose and the
# join: position the child so the two attachment points coincide, orient
# it so the two outward surface normals are anti-aligned, then apply
# `twist` about the joint axis.
#
# Joins must be given in an order such that, when processed left-to-right,
# at least one endpoint of each join has already been visited (starting
# from `root`). Any depth-first or breadth-first ordering from `root`
# satisfies this.

_attach_point(sh, body, att::Attachment) =
    att.shape isa FullCover ?
        surface_centroid(sh, body, att.location) :
        surface_point(sh, body, att.location)

_attach_normal(sh, body, att::Attachment) =
    att.shape isa FullCover ?
        surface_centroid_normal(sh, body, att.location) :
        surface_normal(sh, body, att.location)

function _child_pose(parent_body, parent_pose::Pose, child_body, j::Join)
    sh_p = shape(parent_body); sh_c = shape(child_body)
    pa = j.parent_attachment;  ca = j.child_attachment

    p_local = _attach_point(sh_p, parent_body, pa)
    n_local = _attach_normal(sh_p, parent_body, pa)
    p_world = apply_pose(parent_pose, p_local)
    n_world = apply_rotation(parent_pose.rotation, n_local)

    c_point = _attach_point(sh_c, child_body, ca)
    c_normal = _attach_normal(sh_c, child_body, ca)
    target_normal = (-n_world[1], -n_world[2], -n_world[3])

    R0 = rotation_align(c_normal, target_normal)
    Rtwist = rotation_axis_angle(target_normal, j.twist)
    R = Rtwist * R0

    Rc = apply_rotation(R, c_point)
    t = (p_world[1] - Rc[1], p_world[2] - Rc[2], p_world[3] - Rc[3])
    return Pose(t, R)
end

# Apply one join, given the current `poses` accumulator (Tuple of Pairs).
# Requires exactly one endpoint of `j` to already be in `poses`.
function _apply_join(parts::Tuple, j::Join{PT,CT}, poses::Tuple) where {PT,CT}
    parent_known = _has_tag(poses, PT)
    child_known  = _has_tag(poses, CT)
    if parent_known && !child_known
        pb = _lookup(parts, j.parent)
        cb = _lookup(parts, j.child)
        parent_pose = _lookup(poses, j.parent)
        cp = _child_pose(pb, parent_pose, cb, j)
        return (poses..., j.child => cp)
    elseif child_known && !parent_known
        pb = _lookup(parts, j.parent)
        cb = _lookup(parts, j.child)
        child_pose = _lookup(poses, j.child)
        rj = _reverse_join(j)
        pp = _child_pose(cb, child_pose, pb, rj)
        return (poses..., j.parent => pp)
    elseif parent_known && child_known
        # Cycle — the extra join is a constraint we don't validate here;
        # skip to keep tree structure.
        return poses
    else
        error("Join $(typeof(j.parent)) ↔ $(typeof(j.child)) has neither " *
              "endpoint reachable from root yet — reorder joins so a " *
              "connected endpoint comes first.")
    end
end

# Recursively fold joins into the poses tuple. Type-stable via tuple recursion.
_fold_joins(parts::Tuple, ::Tuple{}, poses::Tuple) = poses
_fold_joins(parts::Tuple, joins::Tuple, poses::Tuple) =
    _fold_joins(parts, Base.tail(joins), _apply_join(parts, joins[1], poses))

function solve_poses(parts::Tuple, joins::Tuple, root_tag, root_pose::Pose)
    _has_tag(parts, typeof(root_tag)) ||
        error("root $(typeof(root_tag)) is not among the parts")
    poses = (root_tag => root_pose,)
    poses = _fold_joins(parts, joins, poses)
    if length(poses) != length(parts)
        visited_types = Set{DataType}(map(p -> typeof(p.first), poses))
        missing_types = filter(t -> !(t in visited_types),
                               collect(map(p -> typeof(p.first), parts)))
        error("parts not reachable from root $(typeof(root_tag)) via joins: $missing_types")
    end
    return poses
end

# Pull a length zero out of a part for the pose translation type.
_length_unit(b::AbstractBody) = zero(b.geometry.characteristic_dimension)

# ── CompositeBody ─────────────────────────────────────────────────────────

"""
    CompositeBody(; parts, joins, root)

A multi-part organism: a Tuple of `tag => Body` pairs joined by `Join`s.
Tags are user-defined singleton types (e.g. `struct Torso end`, then
pass `Torso()`). `root` is the kinematic root tag and serves as the
"primary" part for scalar accessors (`skin_radius`,
`characteristic_dimension`, …) that aren't defined for a composite as a
whole.

The constructor validates each `Join` (surface types, coordinate ranges,
patch sizes) and derives world-frame `poses` for every part by walking
the join tree from `root`.
"""
struct CompositeBody{P<:Tuple, J<:Tuple, R, RP<:Tuple} <: AbstractBody
    parts::P
    joins::J
    root::R
    poses::RP
end

function CompositeBody(; parts, joins, root,
                         root_pose::Union{Pose,Nothing} = nothing)
    parts_t = parts isa Tuple ? parts : Tuple(parts)
    joins_t = joins isa Tuple ? joins : Tuple(joins)
    validate_parts(parts_t)
    _has_tag(parts_t, typeof(root)) ||
        error("root $(typeof(root)) is not among parts $(map(p -> typeof(p.first), parts_t))")
    for j in joins_t
        validate_join(parts_t, j)
    end
    rp = root_pose === nothing ?
         identity_pose(typeof(_length_unit(_lookup(parts_t, root)))) : root_pose
    poses = solve_poses(parts_t, joins_t, root, rp)
    CompositeBody(parts_t, joins_t, root, poses)
end

# ── Accessors that delegate to root ───────────────────────────────────────

shape(b::CompositeBody)      = shape(_lookup(b.parts, b.root))
insulation(b::CompositeBody) = insulation(_lookup(b.parts, b.root))
geometry(b::CompositeBody)   = geometry(_lookup(b.parts, b.root))

skin_radius(b::CompositeBody)       = skin_radius(_lookup(b.parts, b.root))
insulation_radius(b::CompositeBody) = insulation_radius(_lookup(b.parts, b.root))
flesh_radius(b::CompositeBody)      = flesh_radius(_lookup(b.parts, b.root))

# ── Aggregate accessors over parts ────────────────────────────────────────

# For each part in `parts`, subtract the covered area for that tag, then sum.
function _sum_over_parts(f, b::CompositeBody)
    cov = covered_areas(b.parts, b.joins)
    _sum_parts_impl(f, b.parts, cov)
end

_sum_parts_impl(f, ::Tuple{}, ::Tuple{}) =
    error("cannot sum over zero parts")
_sum_parts_impl(f, parts::Tuple{<:Pair}, cov::Tuple{<:Pair}) =
    f(parts[1].second) - cov[1].second
_sum_parts_impl(f, parts::Tuple, cov::Tuple) =
    (f(parts[1].second) - cov[1].second) +
    _sum_parts_impl(f, Base.tail(parts), Base.tail(cov))

total_area(b::CompositeBody)       = _sum_over_parts(total_area, b)
skin_area(b::CompositeBody)        = _sum_over_parts(skin_area, b)
evaporation_area(b::CompositeBody) = _sum_over_parts(evaporation_area, b)

# Simple sum (no covered-area accounting): flesh_volume.
_sum_map(f, ::Tuple{}) = error("cannot sum over zero parts")
_sum_map(f, parts::Tuple{<:Pair}) = f(parts[1].second)
_sum_map(f, parts::Tuple) = f(parts[1].second) + _sum_map(f, Base.tail(parts))

flesh_volume(b::CompositeBody) = _sum_map(flesh_volume, b.parts)

# Silhouette for a composite is the per-part sum — an *upper bound* that
# ignores part-on-part shadowing. For accurate values, project all parts
# together and rasterise: see `silhouette_rasterized`.
_sil_normal(parts::Tuple{<:Pair})   = silhouette_area(parts[1].second).normal
_sil_normal(parts::Tuple)           = silhouette_area(parts[1].second).normal + _sil_normal(Base.tail(parts))
_sil_parallel(parts::Tuple{<:Pair}) = silhouette_area(parts[1].second).parallel
_sil_parallel(parts::Tuple)         = silhouette_area(parts[1].second).parallel + _sil_parallel(Base.tail(parts))

silhouette_area(b::CompositeBody) =
    (; normal = _sil_normal(b.parts), parallel = _sil_parallel(b.parts))

silhouette_area(b::CompositeBody, θ) = _sum_map(p -> silhouette_area(p, θ), b.parts)
silhouette_area(b::CompositeBody, ::NormalToSun)   = silhouette_area(b).normal
silhouette_area(b::CompositeBody, ::ParallelToSun) = silhouette_area(b).parallel
silhouette_area(b::CompositeBody, ::Intermediate)  =
    (silhouette_area(b).normal + silhouette_area(b).parallel) * 0.5
