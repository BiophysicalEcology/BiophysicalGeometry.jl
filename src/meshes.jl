# Parametric surface mesh helpers.
#
# Each `_*_mesh` / `_*_tube` / `_*_cap` returns three Float64 matrices
# `(X, Y, Z)` representing a parametric grid of vertices. Consumers can
# either feed the triplet directly to Makie's `surface!` or iterate
# triangles via `_each_triangle` for silhouette rasterisation, etc.
#
# All helpers are unit-agnostic: caller passes raw numbers in whatever
# units they want (typically metres or cm).

# ── Cylinder ──────────────────────────────────────────────────────────────

function _cylinder_tube(r, L; nθ=72, nz=2, θ_end=2π, z0=0.0)
    θ = LinRange(0.0, θ_end, nθ);  z = LinRange(z0, z0 + Float64(L), nz)
    [r*cos(θi) for θi in θ, _ in z],
    [r*sin(θi) for θi in θ, _ in z],
    [zi        for _  in θ, zi in z]
end

function _cylinder_cap(r, z0; nθ=72, nr=12, θ_end=2π)
    θ = LinRange(0.0, θ_end, nθ);  rv = LinRange(0.0, r, nr)
    [rr*cos(θi) for θi in θ, rr in rv],
    [rr*sin(θi) for θi in θ, rr in rv],
    fill(Float64(z0), nθ, nr)
end

# ── Ellipsoid (sphere is the special case _ellipsoid_mesh(r, r)) ─────────

function _ellipsoid_mesh(a, b; n=60, θ_end=2π)
    θ = LinRange(0.0, θ_end, n);  φ = LinRange(0.0, π, n)
    [a*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [b*cos(φj)         for _  in θ, φj in φ]
end

# Prolate-ellipsoid mesh truncated at the +x pole, parameterised in (α, β):
# x = a·cos(α), y = b·sin(α)·cos(β), z = c·sin(α)·sin(β). α ∈ [α_min, π],
# β ∈ [0, 2π]. `α_min = acos(x_ratio)` with `x_ratio` ∈ (-1, 1] is the
# fraction-of-a where the cut sits (1 → no cut).
function _ellipsoid_mesh_truncated(a, b, c, x_ratio; n=60)
    α_min = acos(clamp(x_ratio, -1.0, 1.0))
    α = LinRange(α_min, π, n);  β = LinRange(0.0, 2π, n)
    [a*cos(αj)            for βi in β, αj in α],
    [b*sin(αj)*cos(βi)    for βi in β, αj in α],
    [c*sin(αj)*sin(βi)    for βi in β, αj in α]
end

# Flat elliptical disc at x = x_cut perpendicular to the long axis, with
# y/z extents = scale·b and scale·c where scale = sqrt(1 - (x_cut/a)²).
function _ellipsoid_pole_a_cap(a, b, c, x_ratio; n=40)
    scale = sqrt(max(0.0, 1 - x_ratio^2))
    rs = LinRange(0.0, 1.0, n);  βs = LinRange(0.0, 2π, n)
    [a * x_ratio          for βi in βs, ri in rs],
    [scale * b * ri * cos(βi) for βi in βs, ri in rs],
    [scale * c * ri * sin(βi) for βi in βs, ri in rs]
end

# ── Half-cylinder (axis +z, dome y ≥ 0, flat at y = 0) ───────────────────

function _half_cylinder_tube(r, L; nθ=72, nz=2, z0=0.0)
    θ = LinRange(0.0, π, nθ);  z = LinRange(z0, z0 + Float64(L), nz)
    [r*cos(θi) for θi in θ, _ in z],
    [r*sin(θi) for θi in θ, _ in z],
    [zi        for _  in θ, zi in z]
end
function _half_cylinder_cap(r, z0; nθ=72, nr=12)
    θ = LinRange(0.0, π, nθ);  rv = LinRange(0.0, r, nr)
    [rr*cos(θi) for θi in θ, rr in rv],
    [rr*sin(θi) for θi in θ, rr in rv],
    fill(Float64(z0), nθ, nr)
end
function _half_cylinder_flat(r, L; nx=12, nz=2, z0=0.0)
    xs = LinRange(-r, r, nx);  zs = LinRange(z0, z0 + Float64(L), nz)
    [xi for xi in xs, _ in zs],
    fill(0.0, nx, nz),
    [zi for _  in xs, zi in zs]
end

# ── Half-ellipsoid (long axis +x, dome z ≥ 0, flat at z = 0) ─────────────

function _half_ellipsoid_dome_mesh(a, b; n=60)
    θ = LinRange(0.0, 2π, n);  φ = LinRange(0.0, π/2, n)
    [a*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [b*cos(φj)         for _  in θ, φj in φ]
end
function _half_ellipsoid_flat_mesh(a, b; n=60)
    rs = LinRange(0.0, 1.0, n);  ts = LinRange(0.0, 2π, n)
    [a*r*cos(t) for t in ts, r in rs],
    [b*r*sin(t) for t in ts, r in rs],
    fill(0.0, n, n)
end

# ── Cone / frustum (axis +z, base at z=0 of radius r_base, top at z=L of
#                   radius r_top = top_ratio*r_base; r_top = 0 for sharp cone)

function _cone_tube(r_base, r_top, L; nθ=72, nz=2, z0=0.0)
    θ = LinRange(0.0, 2π, nθ);  z = LinRange(z0, z0 + Float64(L), nz)
    [(r_base + (r_top - r_base) * (zi - z0)/L) * cos(θi) for θi in θ, zi in z],
    [(r_base + (r_top - r_base) * (zi - z0)/L) * sin(θi) for θi in θ, zi in z],
    [zi                                                  for _  in θ, zi in z]
end

# ── Box faces (Plate) ────────────────────────────────────────────────────

_box_face_z(x1, x2, y1, y2, z) =
    ([xi for xi in (x1, x2), _ in (y1, y2)],
     [yi for _ in (x1, x2), yi in (y1, y2)],
     fill(Float64(z), 2, 2))

_box_face_y(x1, x2, y, z1, z2) =
    ([xi for xi in (x1, x2), _ in (z1, z2)],
     fill(Float64(y), 2, 2),
     [zi for _ in (x1, x2), zi in (z1, z2)])

_box_face_x(x, y1, y2, z1, z2) =
    (fill(Float64(x), 2, 2),
     [yi for yi in (y1, y2), _ in (z1, z2)],
     [zi for _ in (y1, y2), zi in (z1, z2)])

# ── Per-shape outer mesh tiles ────────────────────────────────────────────
#
# Each returns a Vector of (X, Y, Z) matrix triplets in part-local cm
# (multiplied by `sc`). Used by the Makie ext for plotting and by
# silhouette.jl for rasterisation (with sc=1 to stay in metres).

_outer_length_cyl(body) =
    haskey(body.geometry.length, :length_fur) ? body.geometry.length.length_fur : body.geometry.length.length_skin

function _part_outer_meshes(::Cylinder, body, sc)
    R = ustrip(u"m", insulation_radius(body)) * sc
    Lo = ustrip(u"m", _outer_length_cyl(body)) * sc
    Ls = ustrip(u"m", body.geometry.length.length_skin) * sc
    pad = (Lo - Ls) / 2  # fur axial overhang per end
    z0 = -pad
    [_cylinder_tube(R, Lo; z0=z0), _cylinder_cap(R, z0), _cylinder_cap(R, z0 + Lo)]
end

function _part_outer_meshes(shape::Cone, body, sc)
    R = ustrip(u"m", insulation_radius(body)) * sc
    Lo = ustrip(u"m", _outer_length_cyl(body)) * sc
    Ls = ustrip(u"m", body.geometry.length.length_skin) * sc
    pad = (Lo - Ls) / 2
    z0 = -pad
    Rtop = shape.top_ratio * R
    meshes = Any[_cone_tube(R, Rtop, Lo; z0=z0), _cylinder_cap(R, z0)]
    if Rtop > 0
        push!(meshes, _cylinder_cap(Rtop, z0 + Lo))
    end
    meshes
end

function _part_outer_meshes(::Sphere, body, sc)
    R = ustrip(u"m", insulation_radius(body)) * sc
    [_ellipsoid_mesh(R, R)]
end

function _part_outer_meshes(shape::Ellipsoid, body, sc)
    gl = body.geometry.length
    a = ustrip(u"m", haskey(gl, :a_semi_major_fur) ? gl.a_semi_major_fur : gl.a_semi_major_skin) * sc
    b = ustrip(u"m", haskey(gl, :b_semi_minor_fur) ? gl.b_semi_minor_fur : gl.b_semi_minor_skin) * sc
    c = b  # prolate
    if shape.pole_a_truncation == 0
        [_ellipsoid_mesh(a, b)]
    else
        x_ratio = 1 - shape.pole_a_truncation
        [_ellipsoid_mesh_truncated(a, b, c, x_ratio),
         _ellipsoid_pole_a_cap(a, b, c, x_ratio)]
    end
end

function _part_outer_meshes(::Plate, body, sc)
    gl = body.geometry.length
    L = ustrip(u"m", haskey(gl, :length_fur) ? gl.length_fur : gl.length_skin) * sc
    W = ustrip(u"m", haskey(gl, :width_fur)  ? gl.width_fur  : gl.width_skin)  * sc
    H = ustrip(u"m", haskey(gl, :height_fur) ? gl.height_fur : gl.height_skin) * sc
    hL, hW, hH = L/2, W/2, H/2
    [_box_face_z(-hL, hL, -hW, hW, -hH),
     _box_face_z(-hL, hL, -hW, hW,  hH),
     _box_face_y(-hL, hL, -hW, -hH, hH),
     _box_face_y(-hL, hL,  hW, -hH, hH),
     _box_face_x(-hL, -hW, hW, -hH, hH),
     _box_face_x( hL, -hW, hW, -hH, hH)]
end

function _part_outer_meshes(::HalfCylinder, body, sc)
    R = ustrip(u"m", insulation_radius(body)) * sc
    Lo = ustrip(u"m", _outer_length_cyl(body)) * sc
    r_skin = ustrip(u"m", body.geometry.length.radius_skin) * sc
    L_skin = ustrip(u"m", body.geometry.length.length_skin) * sc
    pad = (Lo - L_skin) / 2
    z0 = -pad
    [_half_cylinder_tube(R, Lo; z0=z0), _half_cylinder_cap(R, z0), _half_cylinder_cap(R, z0 + Lo),
     _half_cylinder_flat(r_skin, L_skin; z0=0.0)]
end

function _part_outer_meshes(::HalfEllipsoid, body, sc)
    gl = body.geometry.length
    a = ustrip(u"m", haskey(gl, :a_semi_major_fur) ? gl.a_semi_major_fur : gl.a_semi_major_skin) * sc
    b = ustrip(u"m", haskey(gl, :b_semi_minor_fur) ? gl.b_semi_minor_fur : gl.b_semi_minor_skin) * sc
    a_skin = ustrip(u"m", gl.a_semi_major_skin) * sc
    b_skin = ustrip(u"m", gl.b_semi_minor_skin) * sc
    [_half_ellipsoid_dome_mesh(a, b), _half_ellipsoid_flat_mesh(a_skin, b_skin)]
end

# Each TriMesh face becomes one degenerate 2×2 grid (v1, v2, v3, v3) so it
# slots into the existing parametric-grid pipeline (plotting + silhouette).
function _part_outer_meshes(shape::TriMesh, body, sc)
    out = Vector{Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}()
    for (i, j, k) in shape.faces
        v1 = shape.vertices[i]
        v2 = shape.vertices[j]
        v3 = shape.vertices[k]
        x1 = ustrip(u"m", v1[1]) * sc; y1 = ustrip(u"m", v1[2]) * sc; z1 = ustrip(u"m", v1[3]) * sc
        x2 = ustrip(u"m", v2[1]) * sc; y2 = ustrip(u"m", v2[2]) * sc; z2 = ustrip(u"m", v2[3]) * sc
        x3 = ustrip(u"m", v3[1]) * sc; y3 = ustrip(u"m", v3[2]) * sc; z3 = ustrip(u"m", v3[3]) * sc
        push!(out, ([x1 x3; x2 x3],
                    [y1 y3; y2 y3],
                    [z1 z3; z2 z3]))
    end
    out
end

# ── Pose application ─────────────────────────────────────────────────────

# Apply a Pose (translation in m, dimensionless rotation) to a triple of
# mesh coordinate arrays expressed in (m * sc) units; output in same units.
function _transform_mesh(X, Y, Z, pose::Pose, sc)
    R = pose.rotation
    tx = ustrip(u"m", pose.translation[1]) * sc
    ty = ustrip(u"m", pose.translation[2]) * sc
    tz = ustrip(u"m", pose.translation[3]) * sc
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
        v00 = (X[i,j],     Y[i,j],     Z[i,j])
        v10 = (X[i+1,j],   Y[i+1,j],   Z[i+1,j])
        v01 = (X[i,j+1],   Y[i,j+1],   Z[i,j+1])
        v11 = (X[i+1,j+1], Y[i+1,j+1], Z[i+1,j+1])
        k += 1; out[k] = (v00, v10, v01)
        k += 1; out[k] = (v10, v11, v01)
    end
    out
end
