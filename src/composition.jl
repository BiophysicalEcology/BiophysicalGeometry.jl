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
    EndA(radius=nothing, angle=nothing) <: AbstractSurface

Disc-shaped end cap at the start (`x = 0`) of an axial shape (Cylinder, Cone
and their halves). The located form gives polar coordinates on the disc:
`radius` from the axis and `angle` around it, from `+y` towards `+z`.
"""
struct EndA{R,A} <: AbstractSurface
    radius::R
    angle::A
end
EndA() = EndA(nothing, nothing)

"""
    EndB(radius=nothing, angle=nothing) <: AbstractSurface

Disc-shaped end cap at the far end (`x = length`) of an axial shape; located
as for [`EndA`](@ref).
"""
struct EndB{R,A} <: AbstractSurface
    radius::R
    angle::A
end
EndB() = EndB(nothing, nothing)

"""
    Lateral(position=nothing, angle=nothing) <: AbstractSurface

Curved side surface of an axial shape (Cylinder, Cone and their halves). The
located form gives `position` along the axis from `EndA` and `angle` around it,
from `+y` towards `+z`.
"""
struct Lateral{P,A} <: AbstractSurface
    position::P
    angle::A
end
Lateral() = Lateral(nothing, nothing)

"""
    Flat(x=nothing, y=nothing) <: AbstractSurface

Flat cut face of a half shape, lying in its `z = 0` plane; located by the
`(x, y)` coordinates of a point on it.
"""
struct Flat{X,Y} <: AbstractSurface
    x::X
    y::Y
end
Flat() = Flat(nothing, nothing)

"""
    Dome(polar=nothing, azimuth=nothing) <: AbstractSurface

Curved (dome) surface of a `HalfEllipsoid` or `HalfSphere`. The located form
gives ellipsoidal angles: `polar` from the `+x` pole, `azimuth` around the
long axis from `+y` towards `+z`.
"""
struct Dome{P,A} <: AbstractSurface
    polar::P
    azimuth::A
end
Dome() = Dome(nothing, nothing)

"""
    PoleA(radius=nothing, angle=nothing) <: AbstractSurface

Pole on the positive-x end of an `Ellipsoid`. For an untruncated ellipsoid the
located form is not used (the pole is a point). For a truncated pole
(`pole_a_truncation > 0`) the pole is a disc, located as for [`EndA`](@ref).
"""
struct PoleA{R,A} <: AbstractSurface
    radius::R
    angle::A
end
PoleA() = PoleA(nothing, nothing)

"""
    PoleB() <: AbstractSurface

Pole on the negative-x end of an `Ellipsoid`. Untruncated; the located
form is not meaningful.
"""
struct PoleB <: AbstractSurface end

"""
    Equator(angle=nothing) <: AbstractSurface

Ring at the equator of an `Ellipsoid`; the located form gives the `angle`
around the long axis, from `+y` towards `+z`.
"""
struct Equator{A} <: AbstractSurface
    angle::A
end
Equator() = Equator(nothing)

"""
    Radial(polar=nothing, azimuth=nothing) <: AbstractSurface

Any point on the surface of a `Sphere`, given as spherical angles: `polar`
from `+z`, `azimuth` around it from `+x`.
"""
struct Radial{P,A} <: AbstractSurface
    polar::P
    azimuth::A
end
Radial() = Radial(nothing, nothing)

"""
    Diagonal(position=nothing, z=nothing) <: AbstractSurface

The long diagonal face of a `TriangularPlate`, located by `position` along the
diagonal from its `+x` end and the height `z`.
"""
struct Diagonal{P,Z} <: AbstractSurface
    position::P
    z::Z
end
Diagonal() = Diagonal(nothing, nothing)

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
like `Lateral(position, angle)` used with `Disc`). `shape` is the patch shape
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
        error("Join needs exactly two named attachments, parent then child")
    Vals <: NTuple{2, Attachment} ||
        error("Join arguments must be `Attachment`s")
    P = Names[1]; C = Names[2]
    P === C && error("Join cannot connect a part to itself")
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
    Pose(translation, rotation)

World-frame pose of a part: a translation (3-tuple of length quantities)
and a 3×3 rotation matrix (dimensionless), whose columns are where the
part's x, y and z axes point.
"""
struct Pose{T}
    translation::NTuple{3,T}
    rotation::SMatrix{3,3,Float64,9}
end
Pose(translation, rotation::AbstractMatrix) = Pose(translation, SMatrix{3,3,Float64}(rotation))

const IDENTITY_ROTATION = SMatrix{3,3,Float64}(1, 0, 0, 0, 1, 0, 0, 0, 1)

identity_pose(::Type{T}) where {T} =
    Pose((zero(T), zero(T), zero(T)), IDENTITY_ROTATION)

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
    @SMatrix [t*x*x + c     t*x*y - s*z   t*x*z + s*y;
              t*x*y + s*z   t*y*y + c     t*y*z - s*x;
              t*x*z - s*y   t*y*z + s*x   t*z*z + c]
end

# Rotation matrix that takes unit vector `a` to unit vector `b`.
function rotation_align(a::NTuple{3,<:Real}, b::NTuple{3,<:Real})
    d = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
    if d > 1.0 - 1e-12
        return IDENTITY_ROTATION
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
function surface_area end

"""
    surface_point(shape, body, location) -> NTuple{3,Length}

Local 3D point on the named surface at the given located `location`.
"""
function surface_point end

"""
    surface_normal(shape, body, location) -> NTuple{3,Float64}

Local outward unit normal at the located `location`.
"""
function surface_normal end

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
function surface_centroid end

"""
    surface_centroid_normal(shape, body, ::S) where {S<:AbstractSurface} -> NTuple{3,Float64}

Local outward unit normal at the surface centroid. Used by `FullCover`
attachments.
"""
function surface_centroid_normal end

# ── Patch area dispatch ───────────────────────────────────────────────────

patch_area(body::AbstractBody, att::Attachment{<:AbstractSurface, <:Disc}) =
    π * att.shape.radius^2

patch_area(body::AbstractBody, att::Attachment{<:AbstractSurface, FullCover}) =
    surface_area(shape(body), body, att.location)

# ── Validation ────────────────────────────────────────────────────────────

# A shape supports a surface type if any element of `attachment_surfaces`
# is that surface type (compared by `isa`, so `Lateral(position, angle) isa Lateral`
# works). Tuple iteration is unrolled by the compiler for the small tuples
# used here, so this reduces to a compile-time boolean.
_supports_surface(::Tuple{}, ::AbstractSurface) = false
_supports_surface(surfaces::Tuple, loc::AbstractSurface) =
    (loc isa surfaces[1]) || _supports_surface(Base.tail(surfaces), loc)

function validate_attachment(body::AbstractBody, att::Attachment)
    sh = shape(body)
    surfaces = attachment_surfaces(sh)
    isempty(surfaces) &&
        error("this shape does not support being joined")
    _supports_surface(surfaces, att.location) ||
        error("the shape has no such surface; see `attachment_surfaces`")
    # FullCover has no parametric point; skip range validation.
    if !(att.shape isa FullCover)
        validate_range(sh, body, att.location)
    end
    Asurface = surface_area(sh, body, att.location)
    Apatch = patch_area(body, att)
    if Apatch > Asurface * (1 + 1e-9)
        error("attachment patch area exceeds the area of its surface")
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

validate_parts(::NamedTuple{<:Any,<:Tuple{Body,Vararg{Body}}}) = nothing
validate_parts(parts) = error("CompositeBody parts must be a non-empty NamedTuple of `Body`s")

function validate_join(parts::NamedTuple, j::Join{P, C}) where {P, C}
    pb = getfield(parts, P); cb = getfield(parts, C)
    validate_attachment(pb, j.parent_attachment)
    validate_attachment(cb, j.child_attachment)
    Ap = patch_area(pb, j.parent_attachment)
    Ac = patch_area(cb, j.child_attachment)
    rel = abs(Ap - Ac) / max(Ap, Ac)
    if rel > 1e-6
        ps, cs = j.parent_attachment.shape, j.child_attachment.shape
        if ps isa Disc && cs isa Disc
            error("Join Disc radii must match")
        elseif ps isa FullCover && cs isa FullCover
            error("Join FullCover surfaces must have equal areas")
        else
            error("Join patch areas must match")
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

attach_normal(sh, body, att::Attachment) =
    att.shape isa FullCover ?
        surface_centroid_normal(sh, body, att.location) :
        surface_normal(sh, body, att.location)

function child_pose(parent_body, parent_pose::Pose, child_body, j::Join)
    sh_p = shape(parent_body); sh_c = shape(child_body)
    pa = j.parent_attachment;  ca = j.child_attachment

    p_local = _attach_point(sh_p, parent_body, pa)
    n_local = attach_normal(sh_p, parent_body, pa)
    p_world = apply_pose(parent_pose, p_local)
    n_world = apply_rotation(parent_pose.rotation, n_local)

    c_point = _attach_point(sh_c, child_body, ca)
    c_normal = attach_normal(sh_c, child_body, ca)
    target_normal = (-n_world[1], -n_world[2], -n_world[3])

    R0 = rotation_align(c_normal, target_normal)
    Rtwist = rotation_axis_angle(target_normal, j.twist)
    R = Rtwist * R0

    Rc = apply_rotation(R, c_point)
    t = (p_world[1] - Rc[1], p_world[2] - Rc[2], p_world[3] - Rc[3])
    return Pose(t, R)
end

# Apply one join to the poses found so far, where a part not yet placed has pose
# `nothing`. Which ends are placed is known from the types.
_apply_join(parts, j::Join{P,C}, poses) where {P,C} =
    _apply_join(parts, j, poses, getfield(poses, P), getfield(poses, C))
function _apply_join(parts, j::Join{P,C}, poses, parent::Pose, ::Nothing) where {P,C}
    pose = child_pose(getfield(parts, P), parent, getfield(parts, C), j)
    merge(poses, NamedTuple{(C,)}((pose,)))
end
function _apply_join(parts, j::Join{P,C}, poses, ::Nothing, child::Pose) where {P,C}
    pose = child_pose(getfield(parts, C), child, getfield(parts, P), _reverse_join(j))
    merge(poses, NamedTuple{(P,)}((pose,)))
end
# Both placed: a cycle, whose extra join carries a constraint not checked here.
_apply_join(parts, j, poses, ::Pose, ::Pose) = poses
_apply_join(parts, j, poses, ::Nothing, ::Nothing) =
    error("a Join has neither end joined to the root yet; put a join to one of them first")

# Tuple-recursive fold of joins into the poses.
_fold_joins(parts, ::Tuple{}, poses) = poses
_fold_joins(parts, joins::Tuple, poses) =
    _fold_joins(parts, Base.tail(joins), _apply_join(parts, joins[1], poses))

# The root is the first part.
function solve_poses(parts::NamedTuple{K}, joins::Tuple, root_pose::Pose) where {K}
    unplaced = map(_ -> nothing, Base.tail(values(parts)))
    _placed(_fold_joins(parts, joins, NamedTuple{K}((root_pose, unplaced...))))
end
_placed(poses::NamedTuple{<:Any,<:Tuple{Vararg{Pose}}}) = poses
_placed(poses) = error("some parts are not joined to the root")

# Pull a length zero out of a part for the pose translation type.
_length_unit(b::AbstractBody) = zero(cbrt(b.geometry.volume))

# ── CompositeBody ─────────────────────────────────────────────────────────

"""
    CompositeBody(; parts, joins)

A multi-part organism: a `NamedTuple` of `Body` parts connected by
`Join`s.

`parts` is a NamedTuple like `(; torso=body, head=head_body, leg_fl=leg,
leg_fr=leg, ...)` — the keys are ordinary Julia identifiers, never
`Symbol` literals with a colon. The first-listed part is the kinematic
`root` and serves as the "primary" part for scalar accessors
(`skin_radius`, `insulation_radius`, …) that aren't defined for a
composite as a whole. Reorder `parts` to change the root.

`joins` is a Tuple of `Join`s. Each `Join(<parent_name>=..., <child_name>=...)`
takes exactly two attachment kwargs whose names match keys of `parts`.
Joins must be given in an order such that, when processed left-to-right,
at least one endpoint of each join has already been reached from `root`
(any DFS/BFS ordering works).

The constructor validates each `Join` (surface types, coordinate ranges,
patch sizes) and derives world-frame `poses` for every part.
`CompositeBody(Unchecked(); parts, joins)` only derives the poses: `joins` must
then be a `Tuple`, and nothing is validated.
"""
struct CompositeBody{P<:NamedTuple, J<:Tuple, RP<:NamedTuple} <: AbstractBody
    parts::P
    joins::J
    poses::RP
end

function CompositeBody(; parts::NamedTuple, joins, root_pose::Union{Pose,Nothing} = nothing)
    validate_parts(parts)
    joins_t = joins isa Tuple ? joins : Tuple(joins)
    map(j -> validate_join(parts, j), joins_t)
    CompositeBody(Unchecked(); parts, joins = joins_t, root_pose)
end
function CompositeBody(::Unchecked; parts::NamedTuple, joins::Tuple, root_pose = nothing)
    rp = root_pose === nothing ? identity_pose(typeof(_length_unit(first(parts)))) : root_pose
    poses = solve_poses(parts, joins, rp)
    CompositeBody{typeof(parts), typeof(joins), typeof(poses)}(parts, joins, poses)
end

# ── Accessors that delegate to root ───────────────────────────────────────

_root_part(b::CompositeBody) = first(b.parts)

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
function silhouette(b::CompositeBody)
    sils = map(silhouette, values(b.parts))
    (; normal = sum(map(s -> s.normal, sils)),
       parallel = sum(map(s -> s.parallel, sils)))
end

silhouette(b::CompositeBody, θ) = sum(map(p -> silhouette(p, θ), values(b.parts)))
silhouette(b::CompositeBody, ::NormalToSun) = silhouette(b).normal
silhouette(b::CompositeBody, ::ParallelToSun) = silhouette(b).parallel
silhouette(b::CompositeBody, ::Intermediate) =
    (silhouette(b).normal + silhouette(b).parallel) * 0.5
