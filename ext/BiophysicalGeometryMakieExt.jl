module BiophysicalGeometryMakieExt

using Makie
using Unitful
using BiophysicalGeometry
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Plate, TriangularPlate, Cone, Half
import BiophysicalGeometry: Naked
import BiophysicalGeometry: CompositeBody, Pose, apply_pose, apply_rotation, silhouette_rasterized
import BiophysicalGeometry: AbstractCylindrical, AbstractEllipsoidal, AbstractSpherical
# Mesh helpers now live in core (src/meshes.jl); reuse them here.
# TODO: these should not have leading underscores and be imported in an extension
import BiophysicalGeometry: _cylinder_tube, _cylinder_cap, _ellipsoid_mesh, _cone_tube,
    _half_cylinder_flat, _half_ellipsoid_flat_mesh,
    _box_face_x, _box_face_y, _box_face_z,
    _part_outer_meshes, _transform_mesh, outer_dims, _top_ratio, HalfDomed,
    _triangle_face, _prism_side, _inradius, _ellipsoid_mesh_truncated, _ellipsoid_pole_a_cap

# ══════════════════════════════════════════════════════════════════════════════
# GENERIC HELPERS
# ══════════════════════════════════════════════════════════════════════════════

# Two ways units leave this module:
#   `_pu(x)`      → Float32 cm, for Point2f coordinates in 2-D recipes.
#   `_m(x, sc)`   → Float64 m·sc, for 3-D drawing that uses the caller's
#                   `sc` scale factor.
_pu(x) = Float32(ustrip(u"cm", x))
_m(x, sc) = ustrip(u"m", x) * sc

_radii(body) = (flesh=flesh_radius(body), skin=skin_radius(body), ins=insulation_radius(body))

_scaled_radii(body, sc) = map(x -> _m(x, sc), _radii(body))

_layer_flags(r) = (fat_layer=r.skin > r.flesh + 1e-9, fibrous_layer=r.ins > r.skin + 1e-9)

# Semi-axes (a, b, c) along x, y, z of each layer of an ellipsoidal or spherical
# body: the skin, the flesh inside it (one fat thickness in on every axis) and
# the outer insulation (the fibrous shell, or the skin when there is none).
_domed_layers(::Union{Sphere,Half{<:AbstractSpherical}}, body) =
    map(r -> (r, r, r), _radii(body))
function _domed_layers(::Union{Ellipsoid,Half{<:AbstractEllipsoidal}}, body)
    gl = body.geometry.length
    skin = (gl.length_skin, gl.width_skin, gl.height_skin) ./ 2
    fat = haskey(gl, :fat) ? gl.fat : zero(skin[1])
    ins = haskey(gl, :length_fibrous) ?
        (gl.length_fibrous, gl.width_fibrous, gl.height_fibrous) ./ 2 : skin
    (flesh = skin .- fat, skin = skin, ins = ins)
end

_colors(p) = (flesh=p[:flesh_col][], fat_layer=p[:fat_layer_col][], fibrous_layer=p[:fibrous_layer_col][])
_root_name(::CompositeBody{Root}) where {Root} = Root
# Title name; a Half reads as its constructor (`HalfCone`, `HalfCylinder`, …).
_shape_name(s) = string(nameof(typeof(s)))
_shape_name(h::Half) = "Half" * _shape_name(h.parent)

_limits(lx, ly, tx, ty) = (
    long_x=(-lx, lx), long_y=(-ly, ly),
    tran_x=(-tx, tx), tran_y=(-ty, ty)
)

function _draw_layers!(target, layers)
    for (cond, geom, col) in layers
        cond && _layer!(target, geom(), col)
    end
end

_x_ratio(s::Ellipsoid) = 1 - s.pole_a_truncation
_x_ratio(::Any) = 1.0

# Each layer of a triangular plate is the skin triangle scaled about its
# incentre to that layer's inradius: the fur grows every face by its thickness,
# the flesh is the skin prism scaled down. Returned in the units of `f`.
function _triangle_layers(body, f)
    gl = body.geometry.length
    L, W, H = gl.length_skin, gl.width_skin, gl.height_skin
    r = _inradius(L, W)
    map(_radii(body)) do ri
        t = ri - r
        height = t >= zero(t) ? H + 2t : H * ri / r
        (corner = f(-t), length = f(L * ri / r), width = f(W * ri / r), height = f(height))
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 3-D DRAWING
# ══════════════════════════════════════════════════════════════════════════════
#
# Everything drawn for one body goes into a single mesh. A backend without a depth
# buffer (CairoMakie) sorts the faces of one mesh by depth but not separate plots,
# so separate surfaces for flesh, fat and fibres would hide one another. Colours
# are shaded here from the surface normals, for a light beside the viewer.
#
# A tile is `(X, Y, Z, colour)`: a parametric grid of vertices, as in `src/meshes.jl`.

_view_direction(azimuth, elevation) =
    (cos(elevation) * cos(azimuth), cos(elevation) * sin(azimuth), sin(elevation))

_opaque(c) = (c = RGBAf(c); RGBf(c.r, c.g, c.b))

# Subdivide a coarse grid by linear interpolation, so that no face is long enough
# to be sorted wrongly.
_refine(A::AbstractMatrix, n=16) = _refine_along(_refine_along(A, 1, n), 2, n)

function _refine_along(A::AbstractMatrix{T}, dim, n) where {T}
    m = size(A, dim)
    m >= n && return A
    B = Matrix{T}(undef, dim == 1 ? (n, size(A, 2)) : (size(A, 1), n))
    for k in 1:n
        t = 1 + (k - 1) * (m - 1) / (n - 1)
        lo = clamp(floor(Int, t), 1, m - 1); w = t - lo
        if dim == 1
            for j in axes(A, 2)
                B[k, j] = (1 - w) * A[lo, j] + w * A[lo + 1, j]
            end
        else
            for i in axes(A, 1)
                B[i, k] = (1 - w) * A[i, lo] + w * A[i, lo + 1]
            end
        end
    end
    return B
end

# Shading at grid point (i, j): the cosine between the surface normal and the light, NaN where the grid is
# degenerate (a cone's tip, the centre of a cap).
function _shade_at(X, Y, Z, i, j, light)
    n1, n2 = size(X)
    ip, im, jp, jm = min(i + 1, n1), max(i - 1, 1), min(j + 1, n2), max(j - 1, 1)
    du = (X[ip, j] - X[im, j], Y[ip, j] - Y[im, j], Z[ip, j] - Z[im, j])
    dv = (X[i, jp] - X[i, jm], Y[i, jp] - Y[i, jm], Z[i, jp] - Z[i, jm])
    n = (du[2] * dv[3] - du[3] * dv[2], du[3] * dv[1] - du[1] * dv[3], du[1] * dv[2] - du[2] * dv[1])
    len = sqrt(n[1]^2 + n[2]^2 + n[3]^2)
    return len > 1e-10 ? abs(n[1] * light[1] + n[2] * light[2] + n[3] * light[3]) / len : NaN
end

function _mesh_tiles!(target, tiles, azimuth, elevation)
    d = _view_direction(azimuth, elevation)
    light = (d[1] - 0.35 * d[2], d[2] + 0.35 * d[1], d[3] + 0.6)
    light = light ./ sqrt(sum(abs2, light))
    points = Point3f[]; colours = RGBAf[]; faces = Int[]
    npoints = sum(((X, _, _, _),) -> max(size(X, 1), 16) * max(size(X, 2), 16), tiles; init=0)
    sizehint!(points, npoints); sizehint!(colours, npoints); sizehint!(faces, 6 * npoints)
    for (X, Y, Z, col) in tiles
        _append_tile!(points, colours, faces, _refine(X), _refine(Y), _refine(Z), col, light)
    end
    mesh!(target, points, permutedims(reshape(faces, 3, :)); color=colours, shading=NoShading)
end

function _append_tile!(points, colours, faces, X::Matrix{Float64}, Y::Matrix{Float64}, Z::Matrix{Float64}, col, light)
    n1, n2 = size(X)
    offset = length(points)
    # Degenerate points take the tile's mean shade.
    total, valid = 0.0, 0
    for j in 1:n2, i in 1:n1
        x = _shade_at(X, Y, Z, i, j, light)
        isfinite(x) && (total += x; valid += 1)
    end
    fallback = valid == 0 ? 1.0 : total / valid
    c = RGBAf(col)
    for j in 1:n2, i in 1:n1
        x = _shade_at(X, Y, Z, i, j, light)
        f = Float32(0.55 + 0.45 * (isfinite(x) ? x : fallback))
        push!(points, Point3f(X[i, j], Y[i, j], Z[i, j]))
        push!(colours, RGBAf(c.r * f, c.g * f, c.b * f, c.alpha))
    end
    for j in 1:n2-1, i in 1:n1-1
        k = offset + (j - 1) * n1 + i
        append!(faces, (k, k + 1, k + n1, k + 1, k + n1 + 1, k + n1))
    end
    return nothing
end

# A tile of a drawing: a parametric grid of points and its colour.
const _Tile = Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64},RGBf}

_grid(f, us, vs) = ([f(u, v)[1] for u in us, v in vs], [f(u, v)[2] for u in us, v in vs],
                    [f(u, v)[3] for u in us, v in vs])

# ══════════════════════════════════════════════════════════════════════════════
# CUTAWAYS
# ══════════════════════════════════════════════════════════════════════════════
#
# Flesh is drawn whole. Fat and fibres are drawn over the angles `θs` around the
# long axis, which leave out the part facing the viewer, and the faces of the cut
# are filled in. Every shape lies along x: an axial shape's angles run around x
# from +y towards +z, the others' around z, in the x–y plane.

# The angle around a shape's axis that a view direction `d` (in its frame) comes from.
_view_angle(::Union{AbstractCylindrical,Half{<:AbstractCylindrical}}, d) = atan(d[3], d[2])
_view_angle(::Any, d) = atan(d[2], d[1])
_cut_angles(shape, d, cut) =
    (a = _view_angle(shape, d); range(a + cut/2, a + 2π - cut/2; length=73))

# Radii and lengths of an axial shape; radius `r * taper(x)` at position x along it.
function _axial_layers(body, sc)
    t = Float64(_top_ratio(body.shape))
    L = _m(body.geometry.length.length_skin, sc)
    pad = _m(outer_dims(body.shape, body).length, sc) / 2 - L / 2
    taper(x) = 1 - (1 - t) * clamp(x, 0, L) / L
    r = _scaled_radii(body, sc)
    return (; L, pad, taper, rf=r.flesh, rs=r.skin, ri=r.ins)
end

# A point at radius `r` from the x axis, at position `x` and angle `a`.
_ring(r, x, a) = (x, r * cos(a), r * sin(a))

# Flesh whole over the angles `full`; fat and fibres over each run of angles in `runs`, with the faces of the cut
# filled in, except where a run ends at one of the `edges` (the flat face of a half, which shows the layers itself).
function _axial_shells!(tiles, (; L, pad, taper, rf, rs, ri), cols, full, runs, edges)
    fl = _layer_flags((flesh=rf, skin=rs, ins=ri))
    function shell(r, x0, x1, θ, col)
        push!(tiles, (_grid((a, x) -> _ring(r * taper(x), x, a), θ, range(x0, x1; length=2))..., col))
        for x in (x0, x1)
            taper(x) > 0 && push!(tiles, (_grid((a, ρ) -> _ring(ρ * r * taper(x), x, a),
                θ, range(0, 1; length=2))..., col))
        end
    end
    shell(rf, 0.0, L, full, cols.flesh)
    for θs in runs
        fl.fat_layer && shell(rs, 0.0, L, θs, cols.fat_layer)
        fl.fibrous_layer && shell(ri, -pad, L + pad, θs, cols.fibrous_layer)
        for a in (first(θs), last(θs))
            any(e -> abs(rem2pi(a - e, RoundNearest)) < 1e-9, edges) && continue
            face(r0, r1, x0, x1, col) = push!(tiles, (_grid((ρ, x) -> _ring((r0 + ρ * (r1 - r0)) * taper(x), x, a),
                range(0, 1; length=2), range(x0, x1; length=2))..., col))
            fl.fat_layer && face(rf, rs, 0.0, L, cols.fat_layer)
            if fl.fibrous_layer
                face(rs, ri, 0.0, L, cols.fibrous_layer)
                face(0.0, ri, -pad, 0.0, cols.fibrous_layer); face(0.0, ri, L, L + pad, cols.fibrous_layer)
            end
        end
    end
    return tiles
end

# The parts of the range of angles `θs` (less than a turn) that lie within [lo, hi].
function _runs_within(θs, lo, hi)
    s, e = first(θs), last(θs)
    runs = StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int}[]
    for k in -2:2
        a, b = max(s + 2π * k, lo), min(e + 2π * k, hi)
        b - a > 1e-9 && push!(runs, range(a, b; length=max(2, ceil(Int, 36 * (b - a) / π) + 1)))
    end
    return runs
end

_cutaway_tiles(::Union{Cylinder,Cone}, body, sc, cols, θs) =
    _axial_shells!(_Tile[], _axial_layers(body, sc), cols, range(0, 2π; length=73), [θs], ())

_cutaway_tiles(sh::Union{Sphere,Ellipsoid}, body, sc, cols, θs) = _ellipsoidal_tiles(sh, body, sc, cols, θs, π)
_cutaway_tiles(sh::HalfDomed, body, sc, cols, θs) = _ellipsoidal_tiles(sh.parent, body, sc, cols, θs, π / 2)

# `φ_max = π` is the whole shape; `π / 2` the upper half, closed by its flat face.
# Each layer has its own semi-axes (a, b, c). On a cut ellipsoid, points beyond the
# cut plane are pressed onto it, which draws the flat face where the cap was.
function _ellipsoidal_tiles(sh, body, sc, cols, θs, φ_max)
    fl = _layer_flags(_scaled_radii(body, sc))
    l = map(axes -> map(x -> _m(x, sc), axes), _domed_layers(body.shape, body))
    x_ratio = _x_ratio(sh)
    φs = range(0, φ_max; length=49)
    full = range(0, 2π; length=73)
    point((a, b, c), θ, φ) = (min(a * sin(φ) * cos(θ), x_ratio * a), b * sin(φ) * sin(θ), c * cos(φ))
    between(l0, l1, ρ) = l0 .+ ρ .* (l1 .- l0)
    tiles = _Tile[(_grid((θ, φ) -> point(l.flesh, θ, φ), full, φs)..., cols.flesh)]
    fl.fat_layer && push!(tiles, (_grid((θ, φ) -> point(l.skin, θ, φ), θs, φs)..., cols.fat_layer))
    fl.fibrous_layer && push!(tiles, (_grid((θ, φ) -> point(l.ins, θ, φ), θs, φs)..., cols.fibrous_layer))
    for θ in (first(θs), last(θs))
        face(l0, l1, col) = push!(tiles, (_grid((ρ, φ) -> point(between(l0, l1, ρ), θ, φ),
            range(0, 1; length=2), φs)..., col))
        fl.fat_layer && face(l.flesh, l.skin, cols.fat_layer)
        fl.fibrous_layer && face(l.skin, l.ins, cols.fibrous_layer)
    end
    if φ_max < π
        ring(l0, l1, θ, col) = push!(tiles, (_grid((ρ, t) -> point(between(l0, l1, ρ), t, π / 2),
            range(0, 1; length=2), θ)..., col))
        ring(map(zero, l.flesh), l.flesh, full, cols.flesh)
        fl.fat_layer && ring(l.flesh, l.skin, θs, cols.fat_layer)
        fl.fibrous_layer && ring(l.skin, l.ins, θs, cols.fibrous_layer)
    end
    return tiles
end

# A half cylinder or half cone shows its layers on its flat face (z = 0), which is a
# section along it, and its dome is cut away like a cylinder's.
function _cutaway_tiles(::Half{<:AbstractCylindrical}, body, sc, cols, θs)
    layers = _axial_layers(body, sc)
    (; L, pad, taper, rf, rs, ri) = layers
    fl = _layer_flags((flesh=rf, skin=rs, ins=ri))
    tiles = _axial_shells!(_Tile[], layers, cols, range(0, π; length=37), _runs_within(θs, 0.0, π), (0.0, π))
    flat(y0, y1, xa, xb, col) = push!(tiles, (_grid((y, x) -> (x, y * taper(x), 0.0), range(y0, y1; length=2),
        range(xa, xb; length=2))..., col))
    flat(-rf, rf, 0.0, L, cols.flesh)
    if fl.fat_layer
        flat(rf, rs, 0.0, L, cols.fat_layer); flat(-rs, -rf, 0.0, L, cols.fat_layer)
    end
    if fl.fibrous_layer
        flat(rs, ri, -pad, L + pad, cols.fibrous_layer); flat(-ri, -rs, -pad, L + pad, cols.fibrous_layer)
        flat(-rs, rs, -pad, 0.0, cols.fibrous_layer); flat(-rs, rs, L, L + pad, cols.fibrous_layer)
    end
    return tiles
end

# A plate: nested boxes, the outer two without their top and the two sides nearest the viewer.
function _cutaway_tiles(sh::Plate, body, sc, cols, θs)
    r = _scaled_radii(body, sc)
    fl = _layer_flags(r)
    gl = body.geometry.length
    di = outer_dims(sh, body)
    hw_f = r.flesh
    hl_f = hw_f * Float64(sh.axis_ratio_b)
    hh_f = hl_f / Float64(sh.axis_ratio_c)
    az = (first(θs) + last(θs)) / 2 - π          # towards the viewer
    sx, sy = cos(az) >= 0 ? 1 : -1, sin(az) >= 0 ? 1 : -1
    tiles = _Tile[]
    function box(hl, hw, hh, col; open=false)
        push!(tiles, (_box_face_z(-hl, hl, -hw, hw, -hh)..., col))
        open || push!(tiles, (_box_face_z(-hl, hl, -hw, hw, hh)..., col))
        for s in (-1, 1)
            (open && s == sx) || push!(tiles, (_box_face_x(s * hl, -hw, hw, -hh, hh)..., col))
            (open && s == sy) || push!(tiles, (_box_face_y(-hl, hl, s * hw, -hh, hh)..., col))
        end
    end
    box(hl_f, hw_f, hh_f, cols.flesh)
    fl.fat_layer && box(_m(gl.length_skin, sc) / 2, _m(gl.width_skin, sc) / 2, _m(gl.height_skin, sc) / 2,
                        cols.fat_layer; open=true)
    fl.fibrous_layer && box(_m(di.length, sc) / 2, _m(di.width, sc) / 2, _m(di.height, sc) / 2,
                            cols.fibrous_layer; open=true)
    return tiles
end

# A triangular plate: nested prisms, the outer two without their top.
function _cutaway_tiles(::TriangularPlate, body, sc, cols, θs)
    fl = _layer_flags(_scaled_radii(body, sc))
    l = _triangle_layers(body, x -> _m(x, sc))
    tiles = _Tile[]
    function prism(l, col; open=false)
        p1 = (l.corner, l.corner); p2 = (l.corner + l.length, l.corner); p3 = (l.corner, l.corner + l.width)
        push!(tiles, (_triangle_face(p1, p2, p3, -l.height / 2)..., col))
        open || push!(tiles, (_triangle_face(p1, p2, p3, l.height / 2)..., col))
        for (p, q) in ((p1, p2), (p2, p3), (p3, p1))
            push!(tiles, (_prism_side(p, q, l.height)..., col))
        end
    end
    prism(l.flesh, cols.flesh)
    fl.fat_layer && prism(l.skin, cols.fat_layer; open=true)
    fl.fibrous_layer && prism(l.ins, cols.fibrous_layer; open=true)
    return tiles
end

# ══════════════════════════════════════════════════════════════════════════════
# COMPOSITE BODY DRAWING
# ══════════════════════════════════════════════════════════════════════════════

# Pick a fill colour per part — outer insulation if any, else flesh.
_part_color(body, cols) = body.insulation isa Naked ? cols.flesh : cols.fibrous_layer

function _composite_tiles(b::CompositeBody, sc, cols)
    tiles = _Tile[]
    for name in propertynames(b.parts)
        part = getfield(b.parts, name)
        pose = getfield(b.poses, name)
        col = _part_color(part, cols)
        for mesh in _part_outer_meshes(part.shape, part, sc)
            push!(tiles, (_transform_mesh(mesh..., pose, sc)..., col))
        end
    end
    return tiles
end

# Each part of a composite (those in `parts`, or all) cut open, in place: the cut faces the viewer, whose direction
# is turned into the part's own frame.
function _composite_cutaway_tiles(b::CompositeBody, sc, cols, azimuth, elevation, cut, parts)
    d = _view_direction(azimuth, elevation)
    tiles = _Tile[]
    for name in propertynames(b.parts)
        parts === nothing || name in parts || continue
        _append_cutaway!(tiles, getfield(b.parts, name), getfield(b.poses, name), sc, cols, d, cut)
    end
    return tiles
end

function _append_cutaway!(tiles, part, pose, sc, cols, d, cut)
    l = apply_rotation(transpose(pose.rotation), d)
    for (X, Y, Z, col) in _cutaway_tiles(part.shape, part, sc, cols, _cut_angles(part.shape, l, cut))
        push!(tiles, (_transform_mesh(X, Y, Z, pose, sc)..., col))
    end
    return tiles
end

# Posed world-frame bounding box of a composite (in `sc` units): (mins, maxs).
function _composite_bbox(b::CompositeBody, sc)
    xmin, xmax = Inf, -Inf
    ymin, ymax = Inf, -Inf
    zmin, zmax = Inf, -Inf
    for name in propertynames(b.parts)
        part = getfield(b.parts, name)
        pose = getfield(b.poses, name)
        for grid in _part_outer_meshes(part.shape, part, sc)
            X, Y, Z = _transform_mesh(grid..., pose, sc)
            xmin = min(xmin, minimum(X)); xmax = max(xmax, maximum(X))
            ymin = min(ymin, minimum(Y)); ymax = max(ymax, maximum(Y))
            zmin = min(zmin, minimum(Z)); zmax = max(zmax, maximum(Z))
        end
    end
    ((xmin, ymin, zmin), (xmax, ymax, zmax))
end

_composite_extent(b::CompositeBody, sc) =
    let (mins, maxs) = _composite_bbox(b, sc)
        max(maxs[1]-mins[1], maxs[2]-mins[2], maxs[3]-mins[3])
    end

# Posed bounding box for a single part (in `sc` units).
function _part_bbox(shape, body, pose::Pose, sc)
    xmin, xmax = Inf, -Inf
    ymin, ymax = Inf, -Inf
    zmin, zmax = Inf, -Inf
    for grid in _part_outer_meshes(shape, body, sc)
        X, Y, Z = _transform_mesh(grid..., pose, sc)
        xmin = min(xmin, minimum(X)); xmax = max(xmax, maximum(X))
        ymin = min(ymin, minimum(Y)); ymax = max(ymax, maximum(Y))
        zmin = min(zmin, minimum(Z)); zmax = max(zmax, maximum(Z))
    end
    ((xmin, ymin, zmin), (xmax, ymax, zmax))
end

# ══════════════════════════════════════════════════════════════════════════════
# 2-D PRIMITIVES
# ══════════════════════════════════════════════════════════════════════════════

_circle_pts(r; n=300) =
    [Point2f(_pu(r)*cos(t), _pu(r)*sin(t)) for t in LinRange(0, 2π, n+1)]

_ellipse_pts(xr, yr; n=300) =
    [Point2f(_pu(xr)*cos(t), _pu(yr)*sin(t)) for t in LinRange(0, 2π, n+1)]

_rect_pts(hw, hh) =
    Point2f[(-_pu(hw), -_pu(hh)), (_pu(hw), -_pu(hh)),
            ( _pu(hw), _pu(hh)), (-_pu(hw), _pu(hh))]

# Plan view of a triangular plate layer: the right angle at (corner, corner).
_tri_pts(l) = Point2f[(_pu(l.corner), _pu(l.corner)), (_pu(l.corner + l.length), _pu(l.corner)),
                      (_pu(l.corner), _pu(l.corner + l.width))]

# Frustum long section: base half-width `rb` at y = -hl, top `rt` at y = +hl.
_trap_pts(rb, rt, hl) =
    Point2f[(-_pu(rb), -_pu(hl)), (_pu(rb), -_pu(hl)),
            ( _pu(rt), _pu(hl)), (-_pu(rt), _pu(hl))]

# Half shapes are drawn bulging toward +x from the cut plane x = 0.
_half_trap_pts(rb, rt, hl) =
    Point2f[(0, -_pu(hl)), (_pu(rb), -_pu(hl)), (_pu(rt), _pu(hl)), (0, _pu(hl))]
_half_ellipse_pts(xr, yr; n=150) =
    [Point2f(_pu(xr)*cos(t), _pu(yr)*sin(t)) for t in LinRange(-π/2, π/2, n+1)]

_layer!(target, pts, col) = poly!(target, pts; color=col, strokecolor=col, strokewidth=1)

# ══════════════════════════════════════════════════════════════════════════════
# 2-D SECTION LAYER DISPATCH
# ══════════════════════════════════════════════════════════════════════════════

function _section_layers(sh::Cylinder, body, mode, r, cols)
    r_f, r_s, r_i = r.flesh, r.skin, r.ins
    hl_s = body.geometry.length.length_skin / 2
    hl_i = outer_dims(sh, body).length / 2
    if mode === :long
        [
            (r_i > r_s, () -> _rect_pts(r_i, hl_i), cols.fibrous_layer),
            (r_s > r_f, () -> _rect_pts(r_s, hl_s), cols.fat_layer),
            (true, () -> _rect_pts(r_f, hl_s), cols.flesh),
        ]
    else
        [
            (r_i > r_s, () -> _circle_pts(r_i), cols.fibrous_layer),
            (r_s > r_f, () -> _circle_pts(r_s), cols.fat_layer),
            (true, () -> _circle_pts(r_f), cols.flesh),
        ]
    end
end

# Cone and the cylindrical halves are frustums with top/base ratio t (t = 1
# for a half cylinder); the long section runs base (bottom) to top.
function _section_layers(sh::Union{Cone,Half{<:AbstractCylindrical}}, body, mode, r, cols)
    t = _top_ratio(sh)
    half = sh isa Half
    hl_s = body.geometry.length.length_skin / 2
    hl_i = outer_dims(sh, body).length / 2
    long(x, hl) = half ? _half_trap_pts(x, t * x, hl) : _trap_pts(x, t * x, hl)
    tran(x) = half ? _half_ellipse_pts(x, x) : _circle_pts(x)
    if mode === :long
        [
            (r.ins > r.skin, () -> long(r.ins, hl_i), cols.fibrous_layer),
            (r.skin > r.flesh, () -> long(r.skin, hl_s), cols.fat_layer),
            (true, () -> long(r.flesh, hl_s), cols.flesh),
        ]
    else
        [
            (r.ins > r.skin, () -> tran(r.ins), cols.fibrous_layer),
            (r.skin > r.flesh, () -> tran(r.skin), cols.fat_layer),
            (true, () -> tran(r.flesh), cols.flesh),
        ]
    end
end

# Ellipsoidal sections: the long section shows height (across) × length (up),
# the transverse section height × width; halves bulge toward +x.
function _section_layers(sh::Union{Sphere,Ellipsoid,HalfDomed}, body, mode, r, cols)
    pts = sh isa Half ? _half_ellipse_pts : _ellipse_pts
    # A cut ellipsoid's long section stops at the cut plane, up the plot.
    cut(points, a) = filter(q -> q[2] <= _pu(a) * _x_ratio(sh), points)
    geom((a, b, c)) = mode === :long ? cut(pts(c, a), a) : pts(c, b)
    l = _domed_layers(sh, body)
    [
        (r.ins > r.skin, () -> geom(l.ins), cols.fibrous_layer),
        (r.skin > r.flesh, () -> geom(l.skin), cols.fat_layer),
        (true, () -> geom(l.flesh), cols.flesh),
    ]
end

# Triangular plate: the long section is the side view (length × height), the
# transverse section the plan view of the triangle.
function _section_layers(sh::TriangularPlate, body, mode, r, cols)
    l = _triangle_layers(body, identity)
    geom(x) = mode === :long ? _rect_pts(x.length / 2, x.height / 2) : _tri_pts(x)
    [
        (r.ins > r.skin, () -> geom(l.ins), cols.fibrous_layer),
        (r.skin > r.flesh, () -> geom(l.skin), cols.fat_layer),
        (true, () -> geom(l.flesh), cols.flesh),
    ]
end

function _section_layers(sh::Plate, body, mode, r, cols)
    gl = body.geometry.length
    di = outer_dims(sh, body)
    r_f, r_s, r_i = r.flesh, r.skin, r.ins
    if mode === :long
        d_s = (gl.length_skin / 2, gl.height_skin / 2)
        d_i = (di.length / 2, di.height / 2)
        d_f = (r_f * sh.axis_ratio_b, r_f * sh.axis_ratio_b / sh.axis_ratio_c)
    else
        d_s = (gl.width_skin / 2, gl.height_skin / 2)
        d_i = (di.width / 2, di.height / 2)
        d_f = (r_f, (r_f * sh.axis_ratio_b) / sh.axis_ratio_c)
    end
    [
        (r_i > r_s, () -> _rect_pts(d_i...), cols.fibrous_layer),
        (r_s > r_f, () -> _rect_pts(d_s...), cols.fat_layer),
        (true, () -> _rect_pts(d_f...), cols.flesh),
    ]
end


# ══════════════════════════════════════════════════════════════════════════════
# 2-D AXIS LIMIT DISPATCH
# ══════════════════════════════════════════════════════════════════════════════

function _section_limits(sh::Union{AbstractCylindrical,Half{<:AbstractCylindrical}}, body, r, pad)
    hl_i = outer_dims(sh, body).length / 2
    ri = _pu(r.ins) * (1 + pad)
    li = _pu(hl_i) * (1 + pad)
    _limits(ri, li, ri, li)
end

function _section_limits(sh::TriangularPlate, body, r, pad)
    o = _triangle_layers(body, _pu).ins
    span = max(o.length, o.width) * (1 + pad)
    margin = span * pad
    (long_x = (-o.length / 2, o.length / 2) .* (1 + pad), long_y = (-span, span) ./ 2,
     tran_x = (o.corner - margin, o.corner + span), tran_y = (o.corner - margin, o.corner + span))
end

function _section_limits(sh::Plate, body, r, pad)
    di = outer_dims(sh, body)
    hl_i = di.length / 2
    hh_i = di.height / 2
    hw_i = di.width / 2
    _limits(_pu(hl_i) * (1 + pad), _pu(hh_i) * (1 + pad),
            _pu(hw_i) * (1 + pad), _pu(hh_i) * (1 + pad))
end

function _section_limits(shape::Union{Sphere,Ellipsoid,HalfDomed}, body, r, pad)
    (a, b, c) = map(x -> _pu(x) * (1 + pad), _domed_layers(shape, body).ins)
    _limits(c, a, c, b)
end

# ══════════════════════════════════════════════════════════════════════════════
# RECIPES
# ══════════════════════════════════════════════════════════════════════════════

@recipe(BodyCutaway, body) do scene
    Theme(
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_layer_col = RGBf(1.00, 0.93, 0.55),
        fibrous_layer_col = RGBf(0.76, 0.62, 0.42),
        sc = 100.0,
        azimuth = 5π/4,      # direction of the viewer: the cut faces it, and so does the light
        elevation = π/7,
        cut = π/2,           # angle of fat and fibres cut away
        parts = nothing,     # of a composite, the names of the parts to draw; nothing for all
    )
end

function Makie.plot!(p::BodyCutaway)
    body = p[:body][];  sc = p[:sc][]
    az = p[:azimuth][];  cut = p[:cut][]
    cols = map(_opaque, _colors(p))
    tiles = if body isa CompositeBody
        _composite_cutaway_tiles(body, sc, cols, az, p[:elevation][], cut, p[:parts][])
    else
        _cutaway_tiles(body.shape, body, sc, cols, _cut_angles(body.shape, _view_direction(az, p[:elevation][]), cut))
    end
    _mesh_tiles!(p, tiles, az, p[:elevation][])
    p
end

@recipe(BodyLongSection, body) do scene
    Theme(
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_layer_col = RGBf(1.00, 0.97, 0.60),
        fibrous_layer_col = RGBf(0.76, 0.62, 0.42),
    )
end

function Makie.plot!(p::BodyLongSection)
    body = p[:body][];  r = _radii(body)
    _draw_layers!(p, _section_layers(body.shape, body, :long, r, _colors(p)))
    p
end

@recipe(BodyTransSection, body) do scene
    Theme(
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_layer_col = RGBf(1.00, 0.97, 0.60),
        fibrous_layer_col = RGBf(0.76, 0.62, 0.42),
    )
end

function Makie.plot!(p::BodyTransSection)
    body = p[:body][];  r = _radii(body)
    _draw_layers!(p, _section_layers(body.shape, body, :trans, r, _colors(p)))
    p
end

# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC API — extends the stubs declared in BiophysicalGeometry
# ══════════════════════════════════════════════════════════════════════════════

# The direction an `Axis3` is viewed from, so that the cut faces the viewer.
_view_angles(ax::Axis3) = (ax.azimuth[], ax.elevation[])
_view_angles(ax) = (5π/4, π/7)

"""
    draw_cutaway!(ax::Axis3, body; sc=100.0, cut=π/2, parts=nothing, flesh_col=…, fat_layer_col=…, fibrous_layer_col=…)

Draw `body` into an existing `Axis3`, with the fat and fibres facing the viewer cut away
over the angle `cut`. `sc` converts metres to axis units (default 100 → cm labels).
Each part of a `CompositeBody` is cut open where it faces the viewer; `parts`, a
collection of part names, draws only those.
"""
function BiophysicalGeometry.draw_cutaway!(ax, body;
        sc = 100.0, cut = π/2, parts = nothing,
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_layer_col = RGBf(1.00, 0.93, 0.55),
        fibrous_layer_col = RGBf(0.76, 0.62, 0.42))
    azimuth, elevation = _view_angles(ax)
    bodycutaway!(ax, body; sc, cut, azimuth, elevation, parts, flesh_col, fat_layer_col, fibrous_layer_col)
    ax.xlabel = "x (cm)"; ax.ylabel = "y (cm)"; ax.zlabel = "z (cm)"
end

"""
    plot_body(body; sc=100.0, flesh_col=…, fat_layer_col=…, fibrous_layer_col=…) -> Figure

Create a labelled `Figure` with `Axis3`, draw a quarter-cutaway of `body`,
and attach a legend. Returns the `Figure`.
"""
function BiophysicalGeometry.plot_body(body;
        sc = 100.0,
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_layer_col = RGBf(1.00, 0.93, 0.55),
        fibrous_layer_col = RGBf(0.76, 0.62, 0.42))

    if body isa CompositeBody
        shape_name = "Composite($(length(body.parts)) parts)"
        ins_name = ""
    else
        shape_name = _shape_name(body.shape)
        ins_name = string(nameof(typeof(body.insulation)))
    end

    fig = Figure(size=(600, 480), backgroundcolor=:white)
    Label(fig[0, 1], isempty(ins_name) ? shape_name : "$(shape_name) · $(ins_name)";
          fontsize=13, font=:bold, padding=(0, 0, 10, 0))
    ax = Axis3(fig[1, 1];
               perspectiveness=0.3, viewmode=:fitzoom, aspect=:data,
               elevation=π/7, azimuth=5π/4)
    draw_cutaway!(ax, body; sc, flesh_col, fat_layer_col, fibrous_layer_col)
    Legend(fig[2, 1],
        [PolyElement(polycolor=flesh_col, strokecolor=:saddlebrown, strokewidth=1),
         PolyElement(polycolor=fat_layer_col, strokecolor=:saddlebrown, strokewidth=1),
         PolyElement(polycolor=fibrous_layer_col, strokecolor=:black, strokewidth=1)],
        ["Flesh", "Fat", "Fibres"];
        orientation=:horizontal, framevisible=false)
    return fig
end

"""
    plot_body_silhouette(body::CompositeBody; resolution=128, sc=100.0, kwargs...) -> Figure

Interactive Makie figure: 3-D cutaway of `body` on the left, rasterised
silhouette projection on the right, with `zenith θ` and `azimuth φ`
sliders below. The silhouette area updates live as you drag the sliders.

The `resolution` kwarg controls the silhouette bitmap (≈ accuracy /
performance trade-off). `sc` scales metres to plot units (default 100 → cm).
"""
function BiophysicalGeometry.plot_body_silhouette(body::CompositeBody;
        resolution::Integer = 128,
        sc = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_layer_col = RGBAf(1.00, 0.97, 0.60, 0.75),
        fibrous_layer_col = RGBAf(0.76, 0.62, 0.42, 0.45))

    fig = Figure(size=(1100, 640), backgroundcolor=:white)
    Label(fig[0, 1:2],
          "BiophysicalGeometry.jl — silhouette projection (interactive)";
          fontsize=13, font=:bold, padding=(0, 0, 8, 0))

    ax3 = Axis3(fig[1, 1];
                perspectiveness=0.3, viewmode=:fitzoom, aspect=:data,
                elevation=π/7, azimuth=5π/4,
                xlabel="x (cm)", ylabel="y (cm)", zlabel="z (cm)",
                title="Body + sun direction")
    ax2 = Axis(fig[1, 2]; aspect=DataAspect(), title="Silhouette projection",
               xlabel="u (cm)", ylabel="v (cm)")

    cols = map(_opaque, (flesh=flesh_col, fat_layer=fat_layer_col, fibrous_layer=fibrous_layer_col))
    _mesh_tiles!(ax3, _composite_tiles(body, sc, cols), _view_angles(ax3)...)

    sg = SliderGrid(fig[2, 1:2],
        (label="zenith θ (0=overhead, π/2=horizon)", range=range(0.0, π/2, length=91),
         startvalue=π/4, format="{:.2f} rad"),
        (label="azimuth φ", range=range(0.0, 2π, length=361),
         startvalue=0.0, format="{:.2f} rad"))
    θobs = sg.sliders[1].value
    φobs = sg.sliders[2].value

    sun_dir = lift(θobs, φobs) do θ, φ
        (sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ))
    end

    # Silhouette image + area (single rasteriser call per slider event).
    result = lift(sun_dir) do d
        silhouette_rasterized(body, d; resolution=resolution, return_image=true)
    end

    # Display bitmap as a heatmap in cm-coords on ax2.
    img_x = lift(result) do r
        LinRange(r.x_range[1] * sc, r.x_range[2] * sc, size(r.bitmap, 1))
    end
    img_y = lift(result) do r
        LinRange(r.y_range[1] * sc, r.y_range[2] * sc, size(r.bitmap, 2))
    end
    img = lift(r -> Float32.(r.bitmap), result)
    heatmap!(ax2, img_x, img_y, img; colormap=[:white, RGBAf(0.2, 0.2, 0.25, 1.0)],
             colorrange=(0.0f0, 1.0f0))
    # Refit ax2 limits to the projected bbox each time the sun direction
    # changes (Makie won't auto-shrink the limits otherwise).
    on(result) do r
        xlims!(ax2, r.x_range[1] * sc, r.x_range[2] * sc)
        ylims!(ax2, r.y_range[1] * sc, r.y_range[2] * sc)
    end
    notify(result)

    # Sun-direction indicator: line from a sun marker to the root body's
    # centre (so the target stays at the main body, not pulled around by
    # legs/head when the bbox centre moves).
    root = first(propertynames(body.parts))
    root_part = body.parts[root]
    root_pose = body.poses[root]
    (rb_min, rb_max) = _part_bbox(root_part.shape, root_part, root_pose, sc)
    body_centre = Point3f((rb_min[1] + rb_max[1]) / 2,
                          (rb_min[2] + rb_max[2]) / 2,
                          (rb_min[3] + rb_max[3]) / 2)
    # Whole-composite bbox just used to pick a sensible arrow length and
    # axis limits.
    (bb_min, bb_max) = _composite_bbox(body, sc)
    extent = max(bb_max[1]-bb_min[1], bb_max[2]-bb_min[2], bb_max[3]-bb_min[3])
    arrow_len = 0.7 * extent

    sun_marker_pos = lift(d -> Point3f(body_centre[1] + arrow_len * d[1],
                                        body_centre[2] + arrow_len * d[2],
                                        body_centre[3] + arrow_len * d[3]), sun_dir)
    sun_line = lift(p -> [p, body_centre], sun_marker_pos)
    lines!(ax3, sun_line; color=:goldenrod, linewidth=3)
    scatter!(ax3, lift(p -> [p], sun_marker_pos);
             color=:goldenrod, markersize=20, marker=:circle,
             strokecolor=:black, strokewidth=1)

    # Expand axis limits so the sun marker is always visible.
    pad = arrow_len * 1.05
    xlims!(ax3, body_centre[1] - pad, body_centre[1] + pad)
    ylims!(ax3, body_centre[2] - pad, body_centre[2] + pad)
    zlims!(ax3, body_centre[3] - pad, body_centre[3] + pad)

    Label(fig[3, 1:2],
          lift(r -> "silhouette area = " * string(round(ustrip(u"cm^2", r.area), digits=1)) * " cm²", result);
          fontsize=14, font=:bold)

    return fig
end

"""
    draw_cross_sections!(ax_long, ax_tran, body; flesh_col=…, fat_layer_col=…, fibrous_layer_col=…)

Draw longitudinal and transverse cross-section polygons into a pair of `Axis`
objects.  Layers are drawn outer-to-inner so each inner layer paints over the
outer one.
"""
function BiophysicalGeometry.draw_cross_sections!(ax_long, ax_tran, body;
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_layer_col = RGBf(1.00, 0.97, 0.60),
        fibrous_layer_col = RGBf(0.76, 0.62, 0.42))

    bodylongsection!(ax_long, body; flesh_col, fat_layer_col, fibrous_layer_col)
    bodytranssection!(ax_tran, body; flesh_col, fat_layer_col, fibrous_layer_col)

    r = _radii(body)
    lims = _section_limits(body.shape, body, r, 0.12)
    xlims!(ax_long, lims.long_x...); ylims!(ax_long, lims.long_y...)
    xlims!(ax_tran, lims.tran_x...); ylims!(ax_tran, lims.tran_y...)
end

"""
    plot_cross_sections(body; flesh_col=…, fat_layer_col=…, fibrous_layer_col=…) -> Figure

Create a two-panel `Figure` (longitudinal + transverse cross-sections) with
legend. Returns the `Figure`.
"""
function BiophysicalGeometry.plot_cross_sections(body;
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_layer_col = RGBf(1.00, 0.97, 0.60),
        fibrous_layer_col = RGBf(0.76, 0.62, 0.42))

    shape_name = _shape_name(body.shape)
    ins_name = string(nameof(typeof(body.insulation)))

    fig = Figure(backgroundcolor=:white)
    Label(fig[0, 1:2],
          "$(shape_name) · $(ins_name)  (mass = $(BiophysicalGeometry.mass(body.shape)))";
          fontsize=13, font=:bold, padding=(0, 0, 8, 0))
    ax_long = Axis(fig[1, 1];
                   title="Longitudinal section",
                   xlabel="x (cm)", ylabel="y (cm)", aspect=DataAspect())
    ax_tran = Axis(fig[1, 2];
                   title="Transverse section",
                   xlabel="x (cm)", ylabel="y (cm)", aspect=DataAspect())
    draw_cross_sections!(ax_long, ax_tran, body; flesh_col, fat_layer_col, fibrous_layer_col)
    Legend(fig[2, 1:2],
        [PolyElement(polycolor=flesh_col, strokecolor=flesh_col, strokewidth=1),
         PolyElement(polycolor=fat_layer_col, strokecolor=fat_layer_col, strokewidth=1),
         PolyElement(polycolor=fibrous_layer_col, strokecolor=fibrous_layer_col, strokewidth=1)],
        ["Flesh", "Fat", "Fibres"];
        orientation=:horizontal, framevisible=false)
    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 30)
    return fig
end

"""
    draw_insulation_schematic!(ax, fibrous_layer::FibrousLayer; fibre_length=fibrous_layer.thickness)

Draw a side-view schematic of a `FibrousLayer` insulation layer into `ax`.  Fibre width
is exaggerated for clarity.  When `fibre_length > fibrous_layer.thickness` the fibres are
drawn as tilted parallelograms.
"""
function BiophysicalGeometry.draw_insulation_schematic!(ax, fibrous_layer::FibrousLayer;
        fibre_length = fibrous_layer.thickness)

    thick_mm = ustrip(u"mm", fibrous_layer.thickness)
    fibre_len_mm = ustrip(u"mm", fibre_length)
    d_μm = ustrip(u"μm", fibrous_layer.fibre_diameter)
    n_cm2 = ustrip(u"cm^-2", fibrous_layer.fibre_density)

    spacing_mm = 1.0 / sqrt(n_cm2 / 100.0)
    d_display = spacing_mm * 0.40
    dx = sqrt(max(0.0, fibre_len_mm^2 - thick_mm^2))
    n_show = 8
    skin_h = thick_mm * 0.12
    W = n_show * spacing_mm

    poly!(ax,
          [Point2f(0, 0), Point2f(W, 0),
           Point2f(W, -skin_h), Point2f(0, -skin_h)],
          color=RGBf(0.88, 0.48, 0.42),
          strokecolor=:black, strokewidth=0.5)
    text!(ax, W / 2, -skin_h / 2;
          text="Skin", fontsize=12, align=(:center, :center))

    for i in 0:(n_show - 1)
        xc = (i + 0.5) * spacing_mm
        poly!(ax,
              [Point2f(xc - d_display/2, 0.0),
               Point2f(xc + d_display/2, 0.0),
               Point2f(xc + d_display/2 + dx, thick_mm),
               Point2f(xc - d_display/2 + dx, thick_mm)],
              color=RGBf(0.76, 0.62, 0.42),
              strokecolor=(:saddlebrown, 0.6), strokewidth=0.4)
    end

    lines!(ax, [0.0, W, W + dx, dx, 0.0],
               [0.0, 0.0, thick_mm, thick_mm, 0.0];
           color=(:grey40, 0.5), linewidth=0.8, linestyle=:dash)

    x_ann = W + dx + 0.15 * W
    lines!(ax, [x_ann, x_ann], [0.0, thick_mm]; color=:black, linewidth=1)
    scatter!(ax, [x_ann, x_ann], [0.0, thick_mm];
             color=:black, markersize=6, marker=:rect)
    text!(ax, x_ann + 0.04*W, thick_mm / 2;
          text="thickness\n$(round(Int, thick_mm)) mm",
          fontsize=12, align=(:left, :center))

    if dx > 1e-6
        i_ann = n_show ÷ 2
        xc_ann = (i_ann + 0.5) * spacing_mm
        xfl0 = xc_ann + d_display / 2
        xfl1 = xfl0 + dx
        tilt_deg = round(Int, atand(dx, thick_mm))
        lines!(ax, [xfl0, xfl1], [0.0, thick_mm]; color=:purple, linewidth=1.4)
        scatter!(ax, [xfl0, xfl1], [0.0, thick_mm];
                 color=:purple, markersize=5, marker=:vline)
        text!(ax, (xfl0 + xfl1)/2, thick_mm * 1.02;
              text="L = $(round(fibre_len_mm, digits=1)) mm  ($(tilt_deg)° from vertical)",
              fontsize=11, align=(:center, :bottom), color=:purple)
    end

    y_d = -skin_h * 2.4
    x0 = 0.5 * spacing_mm - d_display / 2
    x1 = 0.5 * spacing_mm + d_display / 2
    lines!(ax, [x0, x1], [y_d, y_d]; color=:darkblue, linewidth=1.2)
    scatter!(ax, [x0, x1], [y_d, y_d]; color=:darkblue, markersize=5, marker=:vline)
    text!(ax, x0, y_d - 0.3;
          text="diameter $(round(Int, d_μm)) μm (drawn wider)",
          fontsize=11, align=(:left, :top), color=:darkblue)

    x_sp0 = dx + 0.5 * spacing_mm
    x_sp1 = dx + 1.5 * spacing_mm
    y_sp = thick_mm * 1.22
    lines!(ax, [x_sp0, x_sp1], [y_sp, y_sp]; color=:darkgreen, linewidth=1.2)
    scatter!(ax, [x_sp0, x_sp1], [y_sp, y_sp]; color=:darkgreen, markersize=5, marker=:vline)
    text!(ax, x_sp0, y_sp + 0.2;
          text="spacing $(round(spacing_mm, digits=2)) mm ($(round(Int, n_cm2)) cm⁻²)",
          fontsize=11, align=(:left, :bottom), color=:darkgreen)

    ax.title = "Fibres"
    ax.xlabel = "position (mm)"
    ax.ylabel = "height above skin (mm)"
    ylims!(ax, -skin_h * 6, thick_mm * 1.55)
    xlims!(ax, -0.1 * W, x_ann + 0.75 * W)
    hidespines!(ax, :t, :r)
end

"""
    draw_insulation_coverage!(ax, fibrous_layer::FibrousLayer; d_range=LinRange(10,120,200), N_range=LinRange(200,9000,200))

Draw a coverage-fraction heatmap (plasma colormap) with contour lines at f = 0.25,
0.50, 0.75, 1.0, and mark the reference point for `fibrous_layer`.  Returns the `Heatmap`
object so the caller can attach a `Colorbar`.
"""
function BiophysicalGeometry.draw_insulation_coverage!(ax, fibrous_layer::FibrousLayer;
        d_range = LinRange(10.0, 120.0, 200),
        N_range = LinRange(200.0, 9000.0, 200))

    cov = [π * ((d_μm * 1e-6) / 2)^2 * (N_cm2 * 1e4)
           for d_μm in d_range, N_cm2 in N_range]

    hm = heatmap!(ax, d_range, N_range, cov;
                  colormap=:plasma, colorrange=(0.0, 1.0), highclip=:white)

    for (level, clr) in [(0.25, :white), (0.50, :white),
                          (0.75, :white), (1.00, :cyan)]
        contour!(ax, d_range, N_range, cov;
                 levels=[level], color=clr, linewidth=1.2, linestyle=:dash)
        N_label_cm2 = level / (π * ((d_range[end] * 1e-6) / 2)^2 * 1e4) * 1e-4
        if 200 < N_label_cm2 < 9000
            text!(ax, d_range[end] - 3, N_label_cm2;
                  text="f=$(level)", fontsize=11, align=(:right, :bottom), color=clr)
        end
    end

    d_ref = ustrip(u"μm", fibrous_layer.fibre_diameter)
    N_ref = ustrip(u"cm^-2", fibrous_layer.fibre_density)
    cov_ref = π * ((d_ref * 1e-6) / 2)^2 * (N_ref * 1e4)
    scatter!(ax, [d_ref], [N_ref]; color=:lime, markersize=11,
             strokecolor=:black, strokewidth=1)
    text!(ax, d_ref + 2, N_ref + 150;
          text="$(d_ref) μm, $(N_ref) cm⁻²\nf ≈ $(round(cov_ref, digits=3))",
          fontsize=11, color=:lime)

    ax.title = "Fraction of skin covered, f = π(d/2)² N"
    ax.xlabel = "Fibre diameter d (μm)"
    ax.ylabel = "Fibre density N (cm⁻²)"
    return hm
end

"""
    plot_insulation_properties(fibrous_layer::FibrousLayer; fibre_length=fibrous_layer.thickness, kwargs...) → Figure

Create a two-panel `Figure` showing a fibrous layer schematic and coverage heatmap for
the supplied `FibrousLayer` object.  `fibre_length` may exceed `fibrous_layer.thickness` to show
tilted fibres.  `d_range` and `N_range` control the heatmap axes (μm / cm⁻²).
"""
function BiophysicalGeometry.plot_insulation_properties(fibrous_layer::FibrousLayer;
        fibre_length = fibrous_layer.thickness,
        d_range = LinRange(10.0, 120.0, 200),
        N_range = LinRange(200.0, 9000.0, 200))

    fig = Figure(size=(1000, 460), backgroundcolor=:white)
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    draw_insulation_schematic!(ax1, fibrous_layer; fibre_length)
    hm = draw_insulation_coverage!(ax2, fibrous_layer; d_range, N_range)
    Colorbar(fig[1, 3], hm; label="Fraction covered, f", width=14, labelsize=13)
    colgap!(fig.layout, 12)
    colsize!(fig.layout, 3, Auto(0.05))
    return fig
end

end  # module BiophysicalGeometryMakieExt
