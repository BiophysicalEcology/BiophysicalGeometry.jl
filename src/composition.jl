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
    Join(; twist=0.0, <parent_name>=parent_attachment, <child_name>=child_attachment)

A connection between two parts of a `CompositeBody`. The two keyword
argument names are the names of the parent and child parts (matching keys
of the `parts` NamedTuple passed to `CompositeBody`); their values are
`Attachment`s. Order matters — the first-listed part is the parent, the
second is the child.

`twist` (radians) sets the rotation about the joint axis — the 6th DOF
that two anti-aligned surface normals don't fix. `twist` is a reserved
kwarg name; a part cannot be named `twist`.

The names are lifted into `Join`'s type parameters (`Parent`, `Child`),
so no `Symbol` value ever appears in the composition machinery at runtime.

# Example

    Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r)),
         head = Attachment(PoleA(), Disc(r)))
"""
struct Join{Parent, Child, A1<:Attachment, A2<:Attachment, T<:Real}
    parent_attachment::A1
    child_attachment::A2
    twist::T
end

# Kwarg constructor. `twist` is reserved; the remaining two kwargs are the
# named parent/child attachments. Their names are lifted into the type
# parameters `Parent` and `Child` — Symbols only ever exist in types.
function Join(; twist::Real=0.0, kwargs...)
    nt = NamedTuple(kwargs)
    _make_join(nt, twist)
end

function _make_join(nt::NamedTuple{Names, Vals}, twist) where {Names, Vals}
    length(Names) == 2 ||
        error("Join needs exactly two named attachments (parent then child); got $(length(Names)): $Names")
    Vals <: NTuple{2, Attachment} ||
        error("Join arguments must be `Attachment`s (got $(Vals.parameters[1]) and $(Vals.parameters[2]))")
    P = Names[1]; C = Names[2]
    P === C && error("Join cannot connect a part to itself (`$P`)")
    Join{P, C, Vals.parameters[1], Vals.parameters[2], typeof(twist)}(nt[1], nt[2], twist)
end

# Reverse a join (parent/child swapped, twist negated). Names swap in the
# type parameters; attachments and twist swap in the fields.
_reverse_join(j::Join{P, C, A1, A2}) where {P, C, A1, A2} =
    Join{C, P, A2, A1, typeof(-j.twist)}(j.child_attachment, j.parent_attachment, -j.twist)

# Part name accessors — pull from the type parameters, so constant-folded.
_parent(::Join{P}) where {P} = P
_child(::Join{P, C}) where {P, C} = C

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

# ── Parts / poses machinery ───────────────────────────────────────────────
#
# `parts` is a `NamedTuple` — user writes `(; torso=body, head=head_body,
# ...)`. `poses` mirrors it. Part names live only as `Symbol` type
# parameters of the NamedTuple and of `Join{Parent, Child}`; no `Symbol`
# ever appears as a runtime value in the composition machinery.
#
# Lookups are plain `getfield(nt, P)` where `P` comes from a
# `where`-bound type parameter — Julia constant-folds those.

function validate_parts(parts::NamedTuple)
    isempty(parts) && error("CompositeBody must have at least one part")
    for (name, body) in pairs(parts)
        body isa AbstractBody ||
            error("part `$name` must be an AbstractBody; got $(typeof(body))")
        body isa CompositeBody &&
            error("part `$name` is a nested CompositeBody, which is not supported")
    end
    return nothing
end

function validate_join(parts::NamedTuple, j::Join{P, C}) where {P, C}
    haskey(parts, P) || error("CompositeBody has no part named `$P`")
    haskey(parts, C) || error("CompositeBody has no part named `$C`")
    pb = getfield(parts, P); cb = getfield(parts, C)
    validate_attachment(pb, j.parent_attachment)
    validate_attachment(cb, j.child_attachment)
    Ap = patch_area(pb, j.parent_attachment)
    Ac = patch_area(cb, j.child_attachment)
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

# ── Covered area accumulator ──────────────────────────────────────────────
#
# Folds over `joins` — each join already knows its two endpoint names in
# its type parameters. Start with zero areas per part, and for each join
# add the two patch areas into the corresponding entries.

covered_areas(parts::NamedTuple, joins::Tuple) =
    _fold_cov(parts, _zero_cov(parts), joins)

_zero_cov(parts::NamedTuple) = map(b -> zero(b.geometry.area.total), parts)

_fold_cov(parts, cov, ::Tuple{}) = cov
_fold_cov(parts, cov, joins::Tuple) =
    _fold_cov(parts, _add_join_cov(parts, cov, joins[1]), Base.tail(joins))

function _add_join_cov(parts::NamedTuple, cov::NamedTuple, j::Join{P, C}) where {P, C}
    p_area = patch_area(getfield(parts, P), j.parent_attachment)
    c_area = patch_area(getfield(parts, C), j.child_attachment)
    cov = merge(cov, NamedTuple{(P,)}((getfield(cov, P) + p_area,)))
    cov = merge(cov, NamedTuple{(C,)}((getfield(cov, C) + c_area,)))
    cov
end

# ── Pose tree solver ──────────────────────────────────────────────────────
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

# Apply one join, given the current `poses` NamedTuple accumulator.
# Requires exactly one endpoint of `j` to already be in `poses`.
function _apply_join(parts::NamedTuple, j::Join{P, C}, poses::NamedTuple) where {P, C}
    parent_known = haskey(poses, P)
    child_known = haskey(poses, C)
    if parent_known && !child_known
        cp = _child_pose(getfield(parts, P), getfield(poses, P),
                         getfield(parts, C), j)
        return merge(poses, NamedTuple{(C,)}((cp,)))
    elseif child_known && !parent_known
        rj = _reverse_join(j)
        pp = _child_pose(getfield(parts, C), getfield(poses, C),
                         getfield(parts, P), rj)
        return merge(poses, NamedTuple{(P,)}((pp,)))
    elseif parent_known && child_known
        # Cycle — extra join carries a constraint we don't validate here.
        return poses
    else
        error("Join `$P` ↔ `$C` has neither endpoint reachable from root " *
              "yet — reorder joins so a connected endpoint comes first.")
    end
end

# Tuple-recursive fold of joins into a growing poses NamedTuple.
_fold_joins(parts, ::Tuple{}, poses) = poses
_fold_joins(parts, joins::Tuple, poses) =
    _fold_joins(parts, Base.tail(joins), _apply_join(parts, joins[1], poses))

function solve_poses(parts::NamedTuple, joins::Tuple, root::Symbol, root_pose::Pose)
    poses = NamedTuple{(root,)}((root_pose,))
    poses = _fold_joins(parts, joins, poses)
    if length(poses) != length(parts)
        missing_names = filter(n -> !(n in propertynames(poses)), propertynames(parts))
        error("parts not reachable from root `$root` via joins: $missing_names")
    end
    return poses
end

# Pull a length zero out of a part for the pose translation type.
_length_unit(b::AbstractBody) = zero(b.geometry.characteristic_dimension)

# ── CompositeBody ─────────────────────────────────────────────────────────

"""
    CompositeBody(; parts, joins)

A multi-part organism: a `NamedTuple` of `Body` parts connected by
`Join`s.

`parts` is a NamedTuple like `(; torso=body, head=head_body, leg_fl=leg,
leg_fr=leg, ...)` — the keys are ordinary Julia identifiers, never
`Symbol` literals with a colon. The first-listed part is the kinematic
`root` and serves as the "primary" part for scalar accessors
(`skin_radius`, `characteristic_dimension`, …) that aren't defined for a
composite as a whole. Reorder `parts` to change the root.

`joins` is a Tuple of `Join`s. Each `Join(<parent_name>=..., <child_name>=...)`
takes exactly two attachment kwargs whose names match keys of `parts`.
Joins must be given in an order such that, when processed left-to-right,
at least one endpoint of each join has already been reached from `root`
(any DFS/BFS ordering works).

The constructor validates each `Join` (surface types, coordinate ranges,
patch sizes) and derives world-frame `poses` for every part.
"""
struct CompositeBody{Root, P<:NamedTuple, J<:Tuple, RP<:NamedTuple} <: AbstractBody
    parts::P
    joins::J
    poses::RP
end

function CompositeBody(; parts::NamedTuple, joins,
                         root_pose::Union{Pose,Nothing} = nothing)
    validate_parts(parts)
    joins_t = joins isa Tuple ? joins : Tuple(joins)
    for j in joins_t
        validate_join(parts, j)
    end
    root = first(propertynames(parts))
    root_body = getfield(parts, root)
    rp = root_pose === nothing ? identity_pose(typeof(_length_unit(root_body))) : root_pose
    poses = solve_poses(parts, joins_t, root, rp)
    CompositeBody{root, typeof(parts), typeof(joins_t), typeof(poses)}(parts, joins_t, poses)
end

# ── Accessors that delegate to root ───────────────────────────────────────

_root_part(b::CompositeBody{Root}) where {Root} = getfield(b.parts, Root)

shape(b::CompositeBody) = shape(_root_part(b))
insulation(b::CompositeBody) = insulation(_root_part(b))
geometry(b::CompositeBody) = geometry(_root_part(b))

skin_radius(b::CompositeBody) = skin_radius(_root_part(b))
insulation_radius(b::CompositeBody) = insulation_radius(_root_part(b))
flesh_radius(b::CompositeBody) = flesh_radius(_root_part(b))

# ── Aggregate accessors over parts ────────────────────────────────────────

# Sum an area-like accessor `f(body)` over all parts, subtracting the
# per-part covered patch area (accounts for attached joints).
function _sum_over_parts(f, b::CompositeBody)
    cov = covered_areas(b.parts, b.joins)
    sum(map((body, c) -> f(body) - c, values(b.parts), values(cov)))
end

total_area(b::CompositeBody) = _sum_over_parts(total_area, b)
skin_area(b::CompositeBody) = _sum_over_parts(skin_area, b)
evaporation_area(b::CompositeBody) = _sum_over_parts(evaporation_area, b)

flesh_volume(b::CompositeBody) = sum(map(flesh_volume, values(b.parts)))

# Silhouette for a composite is the per-part sum — an *upper bound* that
# ignores part-on-part shadowing. For accurate values, project all parts
# together and rasterise: see `silhouette_rasterized`.
function silhouette_area(b::CompositeBody)
    sils = map(silhouette_area, values(b.parts))
    (; normal = sum(map(s -> s.normal, sils)),
       parallel = sum(map(s -> s.parallel, sils)))
end

silhouette_area(b::CompositeBody, θ) = sum(map(p -> silhouette_area(p, θ), values(b.parts)))
silhouette_area(b::CompositeBody, ::NormalToSun) = silhouette_area(b).normal
silhouette_area(b::CompositeBody, ::ParallelToSun) = silhouette_area(b).parallel
silhouette_area(b::CompositeBody, ::Intermediate) =
    (silhouette_area(b).normal + silhouette_area(b).parallel) * 0.5
