# Parametric surface tiles.
#
# A `Tile` is a surface given by a point function over (s, t) ∈ [0, 1]², sampled
# on an `nu × nv` grid. Its triangles are generated from the grid on the fly, so
# meshing and rasterising need no arrays. `tile_matrices` gives the `(X, Y, Z)`
# grids for plotting.
#
# All helpers are unit-agnostic: caller passes raw numbers in whatever
# units they want (typically metres or cm).

struct Tile{F}
    point::F
    nu::Int
    nv::Int
end

vertex(tile::Tile, i, j) = tile.point((i - 1) / (tile.nu - 1), (j - 1) / (tile.nv - 1))

function tile_matrices(tile::Tile)
    points = [vertex(tile, i, j) for i in 1:tile.nu, j in 1:tile.nv]
    (map(p -> p[1], points), map(p -> p[2], points), map(p -> p[3], points))
end

# Fold `f(acc, p1, p2, p3)` over the triangles of a tile. Each grid cell is two
# triangles. Tuples of tiles are combined with `map`, not recursion, which
# inference would give up on for a body of many parts.
function fold_triangles(f, acc, tile::Tile)
    for j in 1:tile.nv - 1, i in 1:tile.nu - 1
        v00, v10 = vertex(tile, i, j), vertex(tile, i + 1, j)
        v01, v11 = vertex(tile, i, j + 1), vertex(tile, i + 1, j + 1)
        acc = f(acc, v00, v10, v01)
        acc = f(acc, v10, v11, v01)
    end
    acc
end

# Call `f(p1, p2, p3)` on every triangle of a tile or of a tuple of tiles.
foreach_triangle(f::F, tile::Tile) where {F} = fold_triangles((_, p1, p2, p3) -> (f(p1, p2, p3); nothing), nothing, tile)
foreach_triangle(f::F, tiles::Tuple) where {F} = (map(tile -> foreach_triangle(f, tile), tiles); nothing)

# ── Cylinder (axis +x; angle θ around it from +y towards +z) ─────────────
#
# The axial meshes take `x0` / `L` along the axis and sweep `θ` from 0 to
# `θ_end`, so `θ_end = π` is the upper (z ≥ 0) half.

cylinder_tube(r, L; nθ=72, nx=2, θ_end=2π, x0=0.0) = cone_tube(r, r, L; nθ, nx, θ_end, x0)

function cylinder_cap(r, x0; nθ=72, nr=12, θ_end=2π)
    Tile(nθ, nr) do s, t
        sθ, cθ = sincos(s * θ_end)
        (Float64(x0), t * r * cθ, t * r * sθ)
    end
end

# ── Ellipsoid (sphere is the special case ellipsoid_mesh(r, r)) ─────────

# `φ_end=π/2` gives the upper (z ≥ 0) dome — the half-ellipsoid / hemisphere.
# Semi-axes a, b, c along x, y, z; `c` defaults to `b` (a spheroid).
function ellipsoid_mesh(a, b, c=b; n=60, θ_end=2π, φ_end=π)
    Tile(n, n) do s, t
        sθ, cθ = sincos(s * θ_end)
        sφ, cφ = sincos(t * φ_end)
        (a * sφ * cθ, b * sφ * sθ, c * cφ)
    end
end

# Ellipsoid truncated at the +x pole, parameterised in (α, β):
# x = a·cos(α), y = b·sin(α)·cos(β), z = c·sin(α)·sin(β). α ∈ [α_min, π],
# β ∈ [0, 2π]. `α_min = acos(x_ratio)` with `x_ratio` ∈ (-1, 1] is the
# fraction-of-a where the cut sits (1 → no cut).
function ellipsoid_mesh_truncated(a, b, c, x_ratio; n=60)
    α_min = acos(clamp(x_ratio, -1.0, 1.0))
    Tile(n, n) do s, t
        sβ, cβ = sincos(2π * s)
        sα, cα = sincos(α_min + (π - α_min) * t)
        (a * cα, b * sα * cβ, c * sα * sβ)
    end
end

# Flat elliptical disc at x = x_cut perpendicular to the long axis, with
# y/z extents = scale·b and scale·c where scale = sqrt(1 - (x_cut/a)²).
function ellipsoid_pole_a_cap(a, b, c, x_ratio; n=40)
    scale = sqrt(max(0.0, 1 - x_ratio^2))
    Tile(n, n) do s, t
        sβ, cβ = sincos(2π * s)
        (a * x_ratio, scale * b * t * cβ, scale * c * t * sβ)
    end
end

# ── Half-cylinder (axis +x, dome z ≥ 0, flat at z = 0) ───────────────────
# The dome and end caps are `cylinder_tube` / `cylinder_cap` at `θ_end=π`;
# only the flat cut face at z = 0 is unique to the half.

# Flat axial face of a half cylinder in z = 0, |y| ≤ radius. With `r_top` it is
# the trapezoidal face of a half cone, radius linear from `r` to `r_top`.
half_cylinder_flat(r, L; r_top=r, nx=2, ny=12, x0=0.0) =
    Tile((s, t) -> (x0 + t * L, (2s - 1) * (r + (r_top - r) * t), 0.0), ny, nx)

# ── Half-ellipsoid (long axis +x, dome z ≥ 0, flat at z = 0) ─────────────
# The dome is `ellipsoid_mesh(a, b; φ_end=π/2)`; only the flat cut face at
# z = 0 is unique to the half.

function half_ellipsoid_flat_mesh(a, b; n=60)
    Tile(n, n) do s, t
        sθ, cθ = sincos(2π * s)
        (a * t * cθ, b * t * sθ, 0.0)
    end
end

# ── Cone / frustum (axis +x, base at x0 of radius r_base, top at x0 + L of
#                   radius r_top = top_ratio*r_base; r_top = 0 for sharp cone)

function cone_tube(r_base, r_top, L; nθ=72, nx=2, θ_end=2π, x0=0.0)
    Tile(nθ, nx) do s, t
        sθ, cθ = sincos(s * θ_end)
        r = r_base + (r_top - r_base) * t
        (x0 + t * L, r * cθ, r * sθ)
    end
end

# ── Box faces (Plate) ────────────────────────────────────────────────────

box_face_z(x1, x2, y1, y2, z) = Tile((s, t) -> (x1 + s * (x2 - x1), y1 + t * (y2 - y1), Float64(z)), 2, 2)
box_face_y(x1, x2, y, z1, z2) = Tile((s, t) -> (x1 + s * (x2 - x1), Float64(y), z1 + t * (z2 - z1)), 2, 2)
box_face_x(x, y1, y2, z1, z2) = Tile((s, t) -> (Float64(x), y1 + s * (y2 - y1), z1 + t * (z2 - z1)), 2, 2)

# ── Per-shape outer mesh tiles ────────────────────────────────────────────
#
# Each `part_outer_meshes` returns a tuple of tiles in part-local metres scaled
# by `sc`. Used by the Makie ext for plotting and by silhouette.jl for
# rasterisation (with sc=1 to stay in metres). The tuple is the same for every
# shape of a type: a cap that a shape lacks is there with zero size.
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
    (cylinder_tube(d.r, d.Lo; x0), cylinder_cap(d.r, x0), cylinder_cap(d.r, x0 + d.Lo))
end

function part_outer_meshes(sh::Cone, body, sc)
    d = _mesh_dims(sh, body, sc)
    x0 = -d.pad
    Rtop = sh.top_ratio * d.r
    (cone_tube(d.r, Rtop, d.Lo; x0), cylinder_cap(d.r, x0), cylinder_cap(Rtop, x0 + d.Lo))
end

function part_outer_meshes(sh::Sphere, body, sc)
    r = _mesh_dims(sh, body, sc).r
    (ellipsoid_mesh(r, r),)
end

function part_outer_meshes(sh::Ellipsoid, body, sc)
    d = _mesh_dims(sh, body, sc)
    x_ratio = 1 - sh.pole_a_truncation
    (ellipsoid_mesh_truncated(d.a, d.b, d.c, x_ratio), ellipsoid_pole_a_cap(d.a, d.b, d.c, x_ratio))
end

function part_outer_meshes(sh::Plate, body, sc)
    d = _mesh_dims(sh, body, sc)
    (box_face_z(-d.hL, d.hL, -d.hW, d.hW, -d.hH),
     box_face_z(-d.hL, d.hL, -d.hW, d.hW, d.hH),
     box_face_y(-d.hL, d.hL, -d.hW, -d.hH, d.hH),
     box_face_y(-d.hL, d.hL, d.hW, -d.hH, d.hH),
     box_face_x(-d.hL, -d.hW, d.hW, -d.hH, d.hH),
     box_face_x( d.hL, -d.hW, d.hW, -d.hH, d.hH))
end

function part_outer_meshes(sh::Half{<:AbstractCylindrical}, body, sc)
    d = _mesh_dims(sh, body, sc)
    x0 = -d.pad
    (cone_tube(d.r, d.t * d.r, d.Lo; θ_end=π, x0),
     cylinder_cap(d.r, x0; θ_end=π),
     half_cylinder_flat(d.r_skin, d.Ls; r_top=d.t * d.r_skin, x0=0.0),
     cylinder_cap(d.t * d.r, x0 + d.Lo; θ_end=π))
end

function part_outer_meshes(sh::Union{Half{<:AbstractEllipsoidal},Half{<:AbstractSpherical}}, body, sc)
    d = _mesh_dims(sh, body, sc)
    (ellipsoid_mesh(d.a, d.b, d.c; φ_end=π/2), half_ellipsoid_flat_mesh(d.a_skin, d.b_skin))
end

# ── Triangular prism (TriangularPlate) ───────────────────────────────────
#
# A triangle is a 2×2 grid with one corner repeated (the second triangle of the
# cell is degenerate); a side is a 2×2 quad from edge p→q over z ∈ [-h/2, h/2].

function triangle_face(p1, p2, p3, z)
    corners = ((p1, p3), (p2, p3))
    Tile(2, 2) do s, t
        p = corners[s < 0.5 ? 1 : 2][t < 0.5 ? 1 : 2]
        (Float64(p[1]), Float64(p[2]), Float64(z))
    end
end
function prism_side(p, q, h)
    Tile(2, 2) do s, t
        e = s < 0.5 ? p : q
        (Float64(e[1]), Float64(e[2]), (t - 0.5) * h)
    end
end

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
    (triangle_face(p1, p2, p3, -d.H / 2), triangle_face(p1, p2, p3, d.H / 2),
     prism_side(p1, p2, d.H), prism_side(p2, p3, d.H), prism_side(p3, p1, d.H))
end

# ── Pose application ─────────────────────────────────────────────────────

# A tile moved by a Pose (translation in m, dimensionless rotation), in (m * sc) units.
function transform_mesh(tile::Tile, pose::Pose, sc)
    R = pose.rotation
    t = map(x -> _ustrip_m(x, sc), pose.translation)
    Tile(tile.nu, tile.nv) do u, v
        x, y, z = tile.point(u, v)
        (R[1,1]*x + R[1,2]*y + R[1,3]*z + t[1],
         R[2,1]*x + R[2,2]*y + R[2,3]*z + t[2],
         R[3,1]*x + R[3,2]*y + R[3,3]*z + t[3])
    end
end
