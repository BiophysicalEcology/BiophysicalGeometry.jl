# Geometry composition for multi-part organisms.
#
# A `CompositeBody` is a tree of `Body` parts joined at named surfaces.
# Each `Join` carries an `Attachment` on each side: the surface name, a
# parametric position on it, and the patch shape that is hidden by the join.
# Surface area accounting subtracts each patch from the part it's on.
# World-frame poses for every part are derived once at construction by
# walking the join tree from `root`, so plotting falls out automatically.

# ── Attachment shapes ─────────────────────────────────────────────────────

abstract type AbstractAttachmentShape end

"""
    Disc{R} <: AbstractAttachmentShape

A flat circular contact patch of radius `R`. Used for legs into bodies, head
poles into cylinder ends, etc.
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
    Attachment{P,S}

One side of a `Join`. Names the attachment surface (e.g. `:end_a`, `:lateral`,
`:flat`), a parametric position on it (`NamedTuple`; `(;)` when there's no
choice), and the patch shape (`Disc` or `FullCover`).
"""
struct Attachment{P<:NamedTuple, S<:AbstractAttachmentShape}
    surface::Symbol
    position::P
    shape::S
end

"""
    Join(parent, parent_attachment, child, child_attachment; twist=0.0)

A connection between two parts of a `CompositeBody`, identified by the keys
`parent` and `child` in `CompositeBody.parts`. `twist` (radians) sets the
rotation about the joint axis — the 6th DOF that two anti-aligned surface
normals don't fix.
"""
struct Join{A1<:Attachment, A2<:Attachment, T<:Real}
    parent::Symbol
    parent_attachment::A1
    child::Symbol
    child_attachment::A2
    twist::T
end
Join(parent::Symbol, pa::Attachment, child::Symbol, ca::Attachment; twist=0.0) =
    Join(parent, pa, child, ca, twist)

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
# Shapes opt into composition by overriding these for each named surface.
# The default `attachment_surfaces` is empty, meaning a shape cannot be
# joined into a `CompositeBody`. The other methods throw a clear error if
# composition code reaches them on an unsupported (shape, surface) pair.

"""
    attachment_surfaces(shape) -> Tuple{Symbol,...}

Names of surfaces on `shape` that can be used in a `Join`. Default `()`
means the shape cannot be joined.
"""
attachment_surfaces(::AbstractShape) = ()

"""
    surface_area(shape, body, ::Val{name}) -> Area

Area of the named surface alone (one face / one side of `shape`).
"""
function surface_area(shape::AbstractShape, body::AbstractBody, ::Val{name}) where {name}
    error("surface_area not defined for $(typeof(shape)) surface :$name")
end

"""
    surface_point(shape, body, ::Val{name}, position::NamedTuple) -> NTuple{3,Length}

Local 3D point on the named surface at the given parametric `position`.
"""
function surface_point(shape::AbstractShape, body::AbstractBody, ::Val{name}, position) where {name}
    error("surface_point not defined for $(typeof(shape)) surface :$name")
end

"""
    surface_normal(shape, body, ::Val{name}, position::NamedTuple) -> NTuple{3,Float64}

Local outward unit normal at `position` on the named surface.
"""
function surface_normal(shape::AbstractShape, body::AbstractBody, ::Val{name}, position) where {name}
    error("surface_normal not defined for $(typeof(shape)) surface :$name")
end

"""
    validate_position(shape, body, ::Val{name}, position::NamedTuple)

Throw if `position` has the wrong field set or is out of range for the
(shape, surface) pair. Default: accept anything (override per shape).
"""
validate_position(::AbstractShape, ::AbstractBody, ::Val, _) = nothing

"""
    surface_centroid(shape, body, ::Val{name}) -> NTuple{3,Length}

Local 3D centroid of the named surface. Used by `FullCover` attachments,
which have no parametric `position`. Falls back to `surface_point` with
an empty NamedTuple if not overridden.
"""
function surface_centroid(shape::AbstractShape, body::AbstractBody, ::Val{name}) where {name}
    error("surface_centroid not defined for $(typeof(shape)) surface :$name")
end

"""
    surface_centroid_normal(shape, body, ::Val{name}) -> NTuple{3,Float64}

Local outward unit normal at the surface centroid. Used by `FullCover`
attachments.
"""
function surface_centroid_normal(shape::AbstractShape, body::AbstractBody, ::Val{name}) where {name}
    error("surface_centroid_normal not defined for $(typeof(shape)) surface :$name")
end

# ── Patch area dispatch ───────────────────────────────────────────────────

patch_area(body::AbstractBody, att::Attachment{<:NamedTuple, <:Disc}) =
    π * att.shape.radius^2

patch_area(body::AbstractBody, att::Attachment{<:NamedTuple, FullCover}) =
    surface_area(shape(body), body, Val(att.surface))

# ── Validation ────────────────────────────────────────────────────────────

function validate_attachment(body::AbstractBody, att::Attachment)
    sh = shape(body)
    surfaces = attachment_surfaces(sh)
    isempty(surfaces) &&
        error("$(typeof(sh)) does not support being joined (attachment_surfaces is empty)")
    att.surface in surfaces ||
        error("$(typeof(sh)) has no surface :$(att.surface); valid: $surfaces")
    # FullCover has no parametric position; skip position validation.
    if !(att.shape isa FullCover)
        validate_position(sh, body, Val(att.surface), att.position)
    end
    Asurface = surface_area(sh, body, Val(att.surface))
    Apatch = patch_area(body, att)
    if Apatch > Asurface * (1 + 1e-9)
        error("attachment patch area ($Apatch) exceeds surface area of " *
              ":$(att.surface) on $(typeof(sh)) ($Asurface)")
    end
    return nothing
end

function validate_join(parts::NamedTuple, j::Join)
    haskey(parts, j.parent) || error("CompositeBody has no part :$(j.parent)")
    haskey(parts, j.child)  || error("CompositeBody has no part :$(j.child)")
    j.parent === j.child &&
        error("Join cannot connect part :$(j.parent) to itself")
    p = parts[j.parent]; c = parts[j.child]
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
                  ":$(j.parent_attachment.surface)=$Ap vs :$(j.child_attachment.surface)=$Ac")
        else
            error("Join patch areas differ: parent $(typeof(ps))=$Ap vs child $(typeof(cs))=$Ac")
        end
    end
    return nothing
end

# ── Covered area accumulator ──────────────────────────────────────────────

function covered_areas(parts::NamedTuple, joins)
    ks = keys(parts)
    A = typeof(parts[first(ks)].geometry.area.total)
    cov = Dict{Symbol, A}(k => zero(parts[k].geometry.area.total) for k in ks)
    for j in joins
        cov[j.parent] += patch_area(parts[j.parent], j.parent_attachment)
        cov[j.child]  += patch_area(parts[j.child],  j.child_attachment)
    end
    NamedTuple{ks}(Tuple(cov[k] for k in ks))
end

# ── Pose tree solver ──────────────────────────────────────────────────────
#
# Each child's world pose is determined by its parent's world pose and the
# join: position the child so the two attachment points coincide, orient it
# so the two outward surface normals are anti-aligned, then apply `twist`
# about the joint axis.

_attach_point(sh, body, att::Attachment) =
    att.shape isa FullCover ?
        surface_centroid(sh, body, Val(att.surface)) :
        surface_point(sh, body, Val(att.surface), att.position)

_attach_normal(sh, body, att::Attachment) =
    att.shape isa FullCover ?
        surface_centroid_normal(sh, body, Val(att.surface)) :
        surface_normal(sh, body, Val(att.surface), att.position)

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

function solve_poses(parts::NamedTuple, joins, root::Symbol, root_pose::Pose)
    haskey(parts, root) || error("root :$root is not a part")
    # Adjacency: which joins involve which part
    neighbours = Dict{Symbol, Vector{Tuple{Symbol, Join, Bool}}}()
    for j in joins
        push!(get!(() -> [], neighbours, j.parent), (j.child, j, true))
        push!(get!(() -> [], neighbours, j.child),  (j.parent, j, false))
    end

    poses = Dict{Symbol, Pose}()
    poses[root] = root_pose
    visited = Set{Symbol}([root])
    stack = [root]
    while !isempty(stack)
        cur = pop!(stack)
        for (other, j, parent_is_cur) in get(neighbours, cur, ())
            other in visited && continue
            push!(visited, other)
            poses[other] = if parent_is_cur
                _child_pose(parts[cur], poses[cur], parts[other], j)
            else
                # Reverse the join: child becomes parent. Swap and recompute.
                rj = Join(j.child, j.child_attachment, j.parent, j.parent_attachment, -j.twist)
                _child_pose(parts[cur], poses[cur], parts[other], rj)
            end
            push!(stack, other)
        end
    end
    if length(visited) != length(parts)
        missing_parts = setdiff(collect(keys(parts)), collect(visited))
        error("parts not reachable from root :$root via joins: $missing_parts")
    end
    ks = keys(parts)
    NamedTuple{ks}(Tuple(poses[k] for k in ks))
end

# Pull a length zero out of a part for the pose translation type.
_length_unit(b::AbstractBody) = zero(b.geometry.characteristic_dimension)

# ── CompositeBody ─────────────────────────────────────────────────────────

"""
    CompositeBody(; parts, joins, root)

A multi-part organism: a `NamedTuple` of `Body` parts joined by `Join`s.
`root` is the kinematic root and serves as the "primary" part for scalar
accessors (`skin_radius`, `characteristic_dimension`, …) that aren't
defined for a composite as a whole.

The constructor validates each `Join` (surface names, positions, patch
sizes) and derives world-frame `poses` for every part by walking the join
tree from `root`.
"""
struct CompositeBody{P<:NamedTuple, J<:Tuple, R<:NamedTuple} <: AbstractBody
    parts::P
    joins::J
    root::Symbol
    poses::R
end

function CompositeBody(; parts::NamedTuple, joins, root::Symbol,
                         root_pose::Union{Pose,Nothing} = nothing)
    haskey(parts, root) || error("root :$root is not a part of $(keys(parts))")
    js = joins isa Tuple ? joins : Tuple(joins)
    for p in values(parts)
        p isa AbstractBody || error("parts must be AbstractBody, got $(typeof(p))")
        p isa CompositeBody && error("nested CompositeBody is not supported")
    end
    for j in js
        validate_join(parts, j)
    end
    rp = root_pose === nothing ? identity_pose(typeof(_length_unit(parts[root]))) : root_pose
    poses = solve_poses(parts, js, root, rp)
    CompositeBody(parts, js, root, poses)
end

# ── Accessors that delegate to root ───────────────────────────────────────

shape(b::CompositeBody)      = shape(b.parts[b.root])
insulation(b::CompositeBody) = insulation(b.parts[b.root])
geometry(b::CompositeBody)   = geometry(b.parts[b.root])

skin_radius(b::CompositeBody)       = skin_radius(b.parts[b.root])
insulation_radius(b::CompositeBody) = insulation_radius(b.parts[b.root])
flesh_radius(b::CompositeBody)      = flesh_radius(b.parts[b.root])

# ── Aggregate accessors over parts ────────────────────────────────────────

function total_area(b::CompositeBody)
    cov = covered_areas(b.parts, b.joins)
    sum(total_area(p) - cov[k] for (k, p) in pairs(b.parts))
end

function skin_area(b::CompositeBody)
    cov = covered_areas(b.parts, b.joins)
    sum(skin_area(p) - cov[k] for (k, p) in pairs(b.parts))
end

function evaporation_area(b::CompositeBody)
    cov = covered_areas(b.parts, b.joins)
    sum(evaporation_area(p) - cov[k] for (k, p) in pairs(b.parts))
end

flesh_volume(b::CompositeBody) = sum(flesh_volume, values(b.parts))

# Silhouette for a composite is the per-part sum — an *upper bound* that
# ignores part-on-part shadowing. For accurate values, project all parts
# together and rasterise: see `silhouette_rasterized`.
function silhouette_area(b::CompositeBody)
    s1 = silhouette_area(first(values(b.parts)))
    total_normal   = zero(s1.normal)
    total_parallel = zero(s1.parallel)
    for p in values(b.parts)
        s = silhouette_area(p)
        total_normal   += s.normal
        total_parallel += s.parallel
    end
    (; normal = total_normal, parallel = total_parallel)
end

silhouette_area(b::CompositeBody, θ) =
    sum(silhouette_area(p, θ) for p in values(b.parts))
silhouette_area(b::CompositeBody, ::NormalToSun)   = silhouette_area(b).normal
silhouette_area(b::CompositeBody, ::ParallelToSun) = silhouette_area(b).parallel
silhouette_area(b::CompositeBody, ::Intermediate)  =
    (silhouette_area(b).normal + silhouette_area(b).parallel) * 0.5
