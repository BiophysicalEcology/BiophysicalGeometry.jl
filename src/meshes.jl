# Parametric surface mesh helpers.
#
# Each `_*_mesh` / `_*_tube` / `_*_cap` returns three Float64 matrices
# `(X, Y, Z)` representing a parametric grid of vertices. Consumers can
# either feed the triplet directly to Makie's `surface!` or iterate
# triangles via `_each_triangle` for silhouette rasterisation, etc.
#
# All helpers are unit-agnostic: caller passes raw numbers in whatever
# units they want (typically metres or cm).
#
# Internal API: these underscore-prefixed names are consumed by
# `silhouette.jl` and the Makie extension (`BiophysicalGeometryMakieExt`).
# Not part of the public package API.

# ── Cylinder (axis +x; angle θ around it from +y towards +z) ─────────────
#
# The axial meshes take `x0` / `L` along the axis and sweep `θ` from 0 to
# `θ_end`, so `θ_end = π` is the upper (z ≥ 0) half.

function cylinder_tube(r, L; nθ=72, nx=2, θ_end=2π, x0=0.0)
    θ = LinRange(0.0, θ_end, nθ)
    x = LinRange(x0, x0 + Float64(L), nx)

    [xi        for _  in θ, xi in x],
    [r*cos(θi) for θi in θ, _ in x],
    [r*sin(θi) for θi in θ, _ in x]
end

function cylinder_cap(r, x0; nθ=72, nr=12, θ_end=2π)
    θ = LinRange(0.0, θ_end, nθ)
    rv = LinRange(0.0, r, nr)

    fill(Float64(x0), nθ, nr),
    [rr*cos(θi) for θi in θ, rr in rv],
    [rr*sin(θi) for θi in θ, rr in rv]
end

# ── Ellipsoid (sphere is the special case ellipsoid_mesh(r, r)) ─────────

# `φ_end=π/2` gives the upper (z ≥ 0) dome — the half-ellipsoid / hemisphere.
# Semi-axes a, b, c along x, y, z; `c` defaults to `b` (a spheroid).
function ellipsoid_mesh(a, b, c=b; n=60, θ_end=2π, φ_end=π)
    θ = LinRange(0.0, θ_end, n)
    φ = LinRange(0.0, φ_end, n)

    [a*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [c*cos(φj)         for _  in θ, φj in φ]
end

# Prolate-ellipsoid mesh truncated at the +x pole, parameterised in (α, β):
# x = a·cos(α), y = b·sin(α)·cos(β), z = c·sin(α)·sin(β). α ∈ [α_min, π],
# β ∈ [0, 2π]. `α_min = acos(x_ratio)` with `x_ratio` ∈ (-1, 1] is the
# fraction-of-a where the cut sits (1 → no cut).
function ellipsoid_mesh_truncated(a, b, c, x_ratio; n=60)
    α_min = acos(clamp(x_ratio, -1.0, 1.0))
    α = LinRange(α_min, π, n)
    β = LinRange(0.0, 2π, n)

    [a*cos(αj)            for βi in β, αj in α],
    [b*sin(αj)*cos(βi)    for βi in β, αj in α],
    [c*sin(αj)*sin(βi)    for βi in β, αj in α]
end

# Flat elliptical disc at x = x_cut perpendicular to the long axis, with
# y/z extents = scale·b and scale·c where scale = sqrt(1 - (x_cut/a)²).
function ellipsoid_pole_a_cap(a, b, c, x_ratio; n=40)
    scale = sqrt(max(0.0, 1 - x_ratio^2))
    rs = LinRange(0.0, 1.0, n)
    βs = LinRange(0.0, 2π, n)

    [a * x_ratio          for βi in βs, ri in rs],
    [scale * b * ri * cos(βi) for βi in βs, ri in rs],
    [scale * c * ri * sin(βi) for βi in βs, ri in rs]
end

# ── Half-cylinder (axis +x, dome z ≥ 0, flat at z = 0) ───────────────────
# The dome and end caps are `cylinder_tube` / `cylinder_cap` at `θ_end=π`;
# only the flat cut face at z = 0 is unique to the half.

# Flat axial face of a half cylinder in z = 0, |y| ≤ radius. With `r_top` it is
# the trapezoidal face of a half cone, radius linear from `r` to `r_top`.
function half_cylinder_flat(r, L; r_top=r, nx=2, ny=12, x0=0.0)
    xs = LinRange(x0, x0 + Float64(L), nx)
    us = LinRange(-1.0, 1.0, ny)
    [xi for _ in us, xi in xs],
    [u * (r + (r_top - r) * (xi - x0) / L) for u in us, xi in xs],
    fill(0.0, ny, nx)
end

# ── Half-ellipsoid (long axis +x, dome z ≥ 0, flat at z = 0) ─────────────
# The dome is `ellipsoid_mesh(a, b; φ_end=π/2)`; only the flat cut face at
# z = 0 is unique to the half.

function half_ellipsoid_flat_mesh(a, b; n=60)
    rs = LinRange(0.0, 1.0, n)
    ts = LinRange(0.0, 2π, n)
    [a*r*cos(t) for t in ts, r in rs],
    [b*r*sin(t) for t in ts, r in rs],
    fill(0.0, n, n)
end

# ── Cone / frustum (axis +x, base at x0 of radius r_base, top at x0 + L of
#                   radius r_top = top_ratio*r_base; r_top = 0 for sharp cone)

function cone_tube(r_base, r_top, L; nθ=72, nx=2, θ_end=2π, x0=0.0)
    θ = LinRange(0.0, θ_end, nθ)
    x = LinRange(x0, x0 + Float64(L), nx)
    [xi                                                  for _  in θ, xi in x],
    [(r_base + (r_top - r_base) * (xi - x0)/L) * cos(θi) for θi in θ, xi in x],
    [(r_base + (r_top - r_base) * (xi - x0)/L) * sin(θi) for θi in θ, xi in x]
end

# ── Box faces (Plate) ────────────────────────────────────────────────────

box_face_z(x1, x2, y1, y2, z) =
    ([xi for xi in (x1, x2), _ in (y1, y2)],
     [yi for _ in (x1, x2), yi in (y1, y2)],
     fill(Float64(z), 2, 2))

box_face_y(x1, x2, y, z1, z2) =
    ([xi for xi in (x1, x2), _ in (z1, z2)],
     fill(Float64(y), 2, 2),
     [zi for _ in (x1, x2), zi in (z1, z2)])

box_face_x(x, y1, y2, z1, z2) =
    (fill(Float64(x), 2, 2),
     [yi for yi in (y1, y2), _ in (z1, z2)],
     [zi for _ in (y1, y2), zi in (z1, z2)])

# ── Per-shape outer mesh tiles ────────────────────────────────────────────
#
# Each `part_outer_meshes` returns a Vector of (X, Y, Z) matrix triplets
# in part-local metres scaled by `sc`. Used by the Makie ext for plotting
# and by silhouette.jl for rasterisation (with sc=1 to stay in metres).
#
# Unit handling is confined to one boundary per shape: `_mesh_dims(shape,
# body, sc)` reads dimensional fields off `body.geometry`, ustrips them
# once, and returns a NamedTuple of plain `Float64` scaled by `sc`. The
# mesh helpers (`cylinder_tube`, `ellipsoid_mesh`, etc.) know nothing
# about units.

_ustrip_m(x, sc) = ustrip(u"m", x) * sc

function _mesh_dims(sh::Cylinder, body, sc)
    d = outer_dims(sh, body)
    Lo = _ustrip_m(d.length, sc)
    Ls = _ustrip_m(body.geometry.length.length_skin, sc)
    (r = _ustrip_m(d.radius, sc), Lo = Lo, Ls = Ls, pad = (Lo - Ls) / 2)
end

function _mesh_dims(sh::Cone, body, sc)
    d = outer_dims(sh, body)
    Lo = _ustrip_m(d.length, sc)
    Ls = _ustrip_m(body.geometry.length.length_skin, sc)
    (r = _ustrip_m(d.radius, sc), Lo = Lo, Ls = Ls, pad = (Lo - Ls) / 2)
end

_mesh_dims(::Sphere, body, sc) = (r = _ustrip_m(insulation_radius(body), sc),)

function _mesh_dims(sh::Ellipsoid, body, sc)
    d = outer_dims(sh, body)
    (a = _ustrip_m(d.length, sc) / 2, b = _ustrip_m(d.width, sc) / 2, c = _ustrip_m(d.height, sc) / 2)
end

function _mesh_dims(sh::Plate, body, sc)
    d = outer_dims(sh, body)
    (hL = _ustrip_m(d.length, sc) / 2,
     hW = _ustrip_m(d.width, sc) / 2,
     hH = _ustrip_m(d.height, sc) / 2)
end

function _mesh_dims(sh::Half{<:AbstractCylindrical}, body, sc)
    d = outer_dims(sh, body)
    Lo = _ustrip_m(d.length, sc)
    Ls = _ustrip_m(body.geometry.length.length_skin, sc)
    (r = _ustrip_m(d.radius, sc),
     Lo = Lo, Ls = Ls, pad = (Lo - Ls) / 2,
     r_skin = _ustrip_m(body.geometry.length.radius_skin, sc),
     t = Float64(top_ratio(sh)))
end

function _mesh_dims(sh::Half{<:AbstractEllipsoidal}, body, sc)
    d = outer_dims(sh, body)
    (a = _ustrip_m(d.length, sc) / 2,
     b = _ustrip_m(d.width, sc) / 2,
     c = _ustrip_m(d.height, sc) / 2,
     a_skin = _ustrip_m(body.geometry.length.length_skin, sc) / 2,
     b_skin = _ustrip_m(body.geometry.length.width_skin, sc) / 2)
end

function _mesh_dims(sh::Half{<:AbstractSpherical}, body, sc)
    r = _ustrip_m(insulation_radius(body), sc)
    rs = _ustrip_m(body.geometry.length.radius_skin, sc)
    (a = r, b = r, c = r, a_skin = rs, b_skin = rs)
end

function part_outer_meshes(sh::Cylinder, body, sc)
    d = _mesh_dims(sh, body, sc)
    x0 = -d.pad
    [cylinder_tube(d.r, d.Lo; x0=x0), cylinder_cap(d.r, x0), cylinder_cap(d.r, x0 + d.Lo)]
end

function part_outer_meshes(sh::Cone, body, sc)
    d = _mesh_dims(sh, body, sc)
    x0 = -d.pad
    Rtop = sh.top_ratio * d.r
    meshes = Any[cone_tube(d.r, Rtop, d.Lo; x0=x0), cylinder_cap(d.r, x0)]
    if Rtop > 0
        push!(meshes, cylinder_cap(Rtop, x0 + d.Lo))
    end
    meshes
end

function part_outer_meshes(sh::Sphere, body, sc)
    d = _mesh_dims(sh, body, sc)
    [ellipsoid_mesh(d.r, d.r)]
end

function part_outer_meshes(sh::Ellipsoid, body, sc)
    d = _mesh_dims(sh, body, sc)
    if sh.pole_a_truncation == 0
        [ellipsoid_mesh(d.a, d.b, d.c)]
    else
        x_ratio = 1 - sh.pole_a_truncation
        [ellipsoid_mesh_truncated(d.a, d.b, d.c, x_ratio),
         ellipsoid_pole_a_cap(d.a, d.b, d.c, x_ratio)]
    end
end

function part_outer_meshes(sh::Plate, body, sc)
    d = _mesh_dims(sh, body, sc)
    [box_face_z(-d.hL, d.hL, -d.hW, d.hW, -d.hH),
     box_face_z(-d.hL, d.hL, -d.hW, d.hW, d.hH),
     box_face_y(-d.hL, d.hL, -d.hW, -d.hH, d.hH),
     box_face_y(-d.hL, d.hL, d.hW, -d.hH, d.hH),
     box_face_x(-d.hL, -d.hW, d.hW, -d.hH, d.hH),
     box_face_x( d.hL, -d.hW, d.hW, -d.hH, d.hH)]
end

function part_outer_meshes(sh::Half{<:AbstractCylindrical}, body, sc)
    d = _mesh_dims(sh, body, sc)
    x0 = -d.pad
    meshes = Any[cone_tube(d.r, d.t * d.r, d.Lo; θ_end=π, x0=x0),
                 cylinder_cap(d.r, x0; θ_end=π),
                 half_cylinder_flat(d.r_skin, d.Ls; r_top=d.t * d.r_skin, x0=0.0)]
    d.t > 0 && push!(meshes, cylinder_cap(d.t * d.r, x0 + d.Lo; θ_end=π))
    meshes
end

function part_outer_meshes(sh::Union{Half{<:AbstractEllipsoidal},Half{<:AbstractSpherical}}, body, sc)
    d = _mesh_dims(sh, body, sc)
    [ellipsoid_mesh(d.a, d.b, d.c; φ_end=π/2), half_ellipsoid_flat_mesh(d.a_skin, d.b_skin)]
end

# ── Triangular prism (TriangularPlate) ───────────────────────────────────
#
# A triangle is a 2×2 grid with one corner repeated (the second triangle of the
# cell is degenerate); a side is a 2×2 quad from edge p→q over z ∈ [-h/2, h/2].

triangle_face(p1, p2, p3, z) =
    ([p1[1] p3[1]; p2[1] p3[1]], [p1[2] p3[2]; p2[2] p3[2]], fill(Float64(z), 2, 2))
prism_side(p, q, h) =
    ([p[1] p[1]; q[1] q[1]], [p[2] p[2]; q[2] q[2]], [-h/2 h/2; -h/2 h/2])

function _mesh_dims(sh::TriangularPlate, body, sc)
    o = outer_dims(sh, body)
    gl = body.geometry.length
    # The fibrous triangle's right angle sits at (-t, -t), t the inradius growth.
    t = inradius(o.length, o.width) - inradius(gl.length_skin, gl.width_skin)
    (corner = -_ustrip_m(t, sc), L = _ustrip_m(o.length, sc),
     W = _ustrip_m(o.width, sc), H = _ustrip_m(o.height, sc))
end

function part_outer_meshes(sh::TriangularPlate, body, sc)
    d = _mesh_dims(sh, body, sc)
    p1 = (d.corner, d.corner); p2 = (d.corner + d.L, d.corner); p3 = (d.corner, d.corner + d.W)
    [triangle_face(p1, p2, p3, -d.H / 2), triangle_face(p1, p2, p3, d.H / 2),
     prism_side(p1, p2, d.H), prism_side(p2, p3, d.H), prism_side(p3, p1, d.H)]
end

# ── Pose application ─────────────────────────────────────────────────────

# Apply a Pose (translation in m, dimensionless rotation) to a triple of
# mesh coordinate arrays expressed in (m * sc) units; output in same units.
function transform_mesh(X, Y, Z, pose::Pose, sc)
    R = pose.rotation
    tx = _ustrip_m(pose.translation[1], sc)
    ty = _ustrip_m(pose.translation[2], sc)
    tz = _ustrip_m(pose.translation[3], sc)
    Xp = similar(X, Float64); Yp = similar(Y, Float64); Zp = similar(Z, Float64)
    @inbounds for i in eachindex(X)
        x, y, z = X[i], Y[i], Z[i]
        Xp[i] = R[1,1]*x + R[1,2]*y + R[1,3]*z + tx
        Yp[i] = R[2,1]*x + R[2,2]*y + R[2,3]*z + ty
        Zp[i] = R[3,1]*x + R[3,2]*y + R[3,3]*z + tz
    end
    return (Xp, Yp, Zp)
end

# ── Triangle iterator ─────────────────────────────────────────────────────

# Each (X, Y, Z) is an `nrow × ncol` parametric grid. Each interior cell
# becomes two triangles. Returns a `Vector{NTuple{3, NTuple{3,Float64}}}`.
function _each_triangle(X, Y, Z)
    nrow, ncol = size(X)
    out = Vector{NTuple{3, NTuple{3, Float64}}}(undef, 2 * (nrow - 1) * (ncol - 1))
    k = 0
    @inbounds for j in 1:ncol-1, i in 1:nrow-1
        v00 = (X[i,j], Y[i,j], Z[i,j])
        v10 = (X[i+1,j], Y[i+1,j], Z[i+1,j])
        v01 = (X[i,j+1], Y[i,j+1], Z[i,j+1])
        v11 = (X[i+1,j+1], Y[i+1,j+1], Z[i+1,j+1])
        k += 1; out[k] = (v00, v10, v01)
        k += 1; out[k] = (v10, v11, v01)
    end
    out
end
