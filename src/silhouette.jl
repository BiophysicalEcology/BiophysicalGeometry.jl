# Rasterised silhouette area for a CompositeBody.
#
# Project every part's outer mesh triangles onto a plane perpendicular to
# the sun direction, rasterise them into a pixel grid, and count the pixels in shadow.
# This handles part-on-part overlap correctly (unlike the per-part summed
# fallback in composition.jl, which double-counts shadowed regions).

# Pick two orthonormal basis vectors perpendicular to a unit vector `d`, with the
# screen-`v` axis aligned to world `+z` where possible so image axes stay
# meaningful. A `+x` reference is used instead when `d` is nearly vertical, to
# avoid the singularity at the zenith.
function _ortho_basis(d::NTuple{3,<:Real})
    up = abs(d[3]) > 0.999 ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)
    proj = up[1]*d[1] + up[2]*d[2] + up[3]*d[3]
    v = (up[1] - proj*d[1], up[2] - proj*d[2], up[3] - proj*d[3])
    nv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    v = (v[1]/nv, v[2]/nv, v[3]/nv)
    # u = v × d, so (u, v, d) is right-handed.
    u = (v[2]*d[3] - v[3]*d[2], v[3]*d[1] - v[1]*d[3], v[1]*d[2] - v[2]*d[1])
    (u, v)
end

function _normalize3(d::NTuple{3,<:Real})
    n = sqrt(d[1]^2 + d[2]^2 + d[3]^2)
    (d[1]/n, d[2]/n, d[3]/n)
end

function _check_direction(d)
    all(isfinite, d) || throw(ArgumentError("a direction must be finite, got $d"))
    d[1]^2 + d[2]^2 + d[3]^2 > 0 || throw(ArgumentError("a direction must be non-zero, got $d"))
    nothing
end

"""
    Beam(direction)

A collimated beam from one `direction` — parallel rays from a source at infinity,
i.e. the direct (beam) component of radiation, as opposed to the diffuse `Sky` /
`Ground` regions. `silhouette(body, Beam(dir))` gives each part's projected area
facing the beam. `Beam(Unchecked(), direction)` doesn't check the direction.
"""
struct Beam
    direction::NTuple{3,Float64}
    Beam(::Unchecked, d::NTuple{3,<:Real}) = new(_normalize3(d))
end
Beam(d::NTuple{3,<:Real}) = (_check_direction(d); Beam(Unchecked(), d))
Beam(x::Real, y::Real, z::Real) = Beam((x, y, z))

# Edge-function rasteriser shared by the coverage and depth passes. Calls
# `f(i, j, w1, w2, w3)` for every pixel whose centre lies inside the 2D triangle
# (p1, p2, p3), with (w1, w2, w3) the pixel centre's barycentric weights. (x0, y0)
# is the grid origin, (dx, dy) the cell size and `n` the side length.
# A whole-numbered pixel coordinate as an index clamped to 1:n. Clamped first,
# `unsafe_trunc` is exact and, unlike `floor(Int, x)`, has no error to throw.
_pixel(x, n) = unsafe_trunc(Int, clamp(x, 1, n))

@inline function _each_pixel(f, p1, p2, p3, x0, y0, dx, dy, n)
    s = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    abs(s) < 1e-18 && return  # degenerate
    sgn = sign(s)
    inv_area = 1 / abs(s)

    imin = _pixel(floor((min(p1[1], p2[1], p3[1]) - x0) / dx) + 1, n)
    imax = _pixel(ceil((max(p1[1], p2[1], p3[1]) - x0) / dx), n)
    jmin = _pixel(floor((min(p1[2], p2[2], p3[2]) - y0) / dy) + 1, n)
    jmax = _pixel(ceil((max(p1[2], p2[2], p3[2]) - y0) / dy), n)

    @inbounds for j in jmin:jmax
        y = y0 + (j - 0.5) * dy
        for i in imin:imax
            x = x0 + (i - 0.5) * dx
            e1 = sgn * ((p2[1] - p1[1]) * (y - p1[2]) - (p2[2] - p1[2]) * (x - p1[1]))
            e1 < 0 && continue
            e2 = sgn * ((p3[1] - p2[1]) * (y - p2[2]) - (p3[2] - p2[2]) * (x - p2[1]))
            e2 < 0 && continue
            e3 = sgn * ((p1[1] - p3[1]) * (y - p3[2]) - (p1[2] - p3[2]) * (x - p3[1]))
            e3 < 0 && continue
            # e2, e3, e1 are twice the sub-triangle areas opposite p1, p2, p3.
            f(i, j, e2 * inv_area, e3 * inv_area, e1 * inv_area)
        end
    end
end

# Flag the pixels a triangle covers.
_rasterize_triangle!(shadow, p1, p2, p3, x0, y0, dx, dy, n) =
    _each_pixel((i, j, _, _, _) -> (@inbounds shadow[i, j] = true), p1, p2, p3, x0, y0, dx, dy, n)

# Rasterise a triangle into a depth buffer, keeping the nearest depth per pixel.
# Depth is interpolated across the triangle from its vertex depths `z`
# (distance along the view direction; larger = nearer the source), so long
# triangles — a tube's lateral strips span the whole part — order correctly
# against other parts at every pixel, not just at their centroids.
function _rasterize_depth!(depth_buf, k, p1, p2, p3, z, x0, y0, dx, dy, n)
    _each_pixel(p1, p2, p3, x0, y0, dx, dy, n) do i, j, w1, w2, w3
        depth = w1 * z[1] + w2 * z[2] + w3 * z[3]
        @inbounds depth > depth_buf[i, j, k] && (depth_buf[i, j, k] = depth)
    end
end

# ── Shared projection ─────────────────────────────────────────────────────

# Every part's outer mesh tiles, posed in world metres.
_part_tiles(body::CompositeBody) =
    map(values(body.parts), values(body.poses)) do part, pose
        map(tile -> transform_mesh(tile, pose, 1.0), part_outer_meshes(part.shape, part, 1.0))
    end

# The view along `d`: a basis (u, v) of the plane ⟂ `d`, and a `resolution ×
# resolution` pixel grid on it fitted around every triangle, with a 2% margin so
# triangles don't clip the grid edge.
function _view(part_tiles, d, resolution)
    u, v = _ortho_basis(d)
    xmin, xmax, ymin, ymax = _bounds(part_tiles, u, v)
    pad = 0.02 * max(xmax - xmin, ymax - ymin)
    x0, x1 = xmin - pad, xmax + pad
    y0, y1 = ymin - pad, ymax + pad
    (; u, v, d, x0, x1, y0, y1, dx = (x1 - x0) / resolution, dy = (y1 - y0) / resolution, resolution)
end

# (xmin, xmax, ymin, ymax) of a tile, or of tuples of tiles, projected onto (u, v).
_onto(p, a) = p[1]*a[1] + p[2]*a[2] + p[3]*a[3]
function _bounds(tile::Tile, u, v)
    fold_triangles((Inf, -Inf, Inf, -Inf), tile) do b, p1, p2, p3
        xs = (_onto(p1, u), _onto(p2, u), _onto(p3, u))
        ys = (_onto(p1, v), _onto(p2, v), _onto(p3, v))
        (min(b[1], xs...), max(b[2], xs...), min(b[3], ys...), max(b[4], ys...))
    end
end
_bounds(tiles::Tuple, u, v) =
    reduce(map(t -> _bounds(t, u, v), tiles); init = (Inf, -Inf, Inf, -Inf)) do a, b
        (min(a[1], b[1]), max(a[2], b[2]), min(a[3], b[3]), max(a[4], b[4]))
    end

_project(g, p) = (_onto(p, g.u), _onto(p, g.v))
_depth(g, p) = _onto(p, g.d)

# Rasterise each part `k` into `depth[:, :, k]` (-Inf where the part doesn't
# cover a pixel). Keeping parts separate lets each pixel be attributed to the
# part directly in front of another, not just to the frontmost overall.
function _part_depths!(depth, g, part_tiles)
    fill!(depth, -Inf)
    map(part_tiles, ntuple(identity, Val(length(part_tiles)))) do tiles, k
        foreach_triangle(tiles) do p1, p2, p3
            z = (_depth(g, p1), _depth(g, p2), _depth(g, p3))
            _rasterize_depth!(depth, k, _project(g, p1), _project(g, p2), _project(g, p3), z,
                              g.x0, g.y0, g.dx, g.dy, g.resolution)
        end
    end
    depth
end

# A NamedTuple over the parts of `body` of `f(k)` for each part index `k`. `F` makes
# Julia specialise on `f`, which is only passed on.
_per_part(f::F, body::CompositeBody) where {F} = _per_part(f, body.parts)
_per_part(f::F, ::NamedTuple{K,<:NTuple{N,Any}}) where {F,K,N} = NamedTuple{K}(ntuple(f, Val(N)))

_depth_buffer(body::CompositeBody, resolution) = Array{Float64}(undef, resolution, resolution, length(body.parts))

# The part in front at pixel (i, j), or 0 where no part covers it.
function _front(depth, i, j)
    front, best = 0, -Inf
    @inbounds for k in axes(depth, 3)
        depth[i, j, k] > best && ((front, best) = (k, depth[i, j, k]))
    end
    front
end

"""
    SilhouetteResult

Result of a rasterised silhouette projection.

- `shadow` : `Matrix{Bool}` of size `resolution × resolution`; `true` at each pixel the body's shadow falls on.
- `x_range`, `y_range` : `(min, max)` axis extents (m) of the projection plane.
- `area` : silhouette area as a `Quantity{Float64, m²}`.
"""
struct SilhouetteResult
    shadow::Matrix{Bool}
    x_range::NTuple{2,Float64}
    y_range::NTuple{2,Float64}
    area::typeof(1.0u"m^2")
end

"""
    silhouette_rasterized(body::CompositeBody, sun_direction; resolution=256, return_image=false)

Compute the silhouette area of `body` projected onto a plane perpendicular
to `sun_direction` (a 3-tuple of any non-zero numbers — automatically
normalised). Unlike the per-part summed default, this correctly accounts
for parts shadowing each other: each part's posed mesh is triangulated,
projected, and rasterised into a `resolution × resolution` image of its
shadow whose pixels are counted.

Returns a length² `Quantity` (m²) by default, or a `SilhouetteResult`
including the shadow image and projection ranges if `return_image=true`.
Increase `resolution` for more accuracy. [`silhouette_rasterized!`](@ref)
draws into an image you pass in instead.
"""
function silhouette_rasterized(body::CompositeBody, sun_direction::NTuple{3,<:Real};
                                resolution::Integer = 256,
                                return_image::Bool = false)
    shadow = Matrix{Bool}(undef, resolution, resolution)
    area, g = _silhouette_rasterized!(shadow, body, sun_direction)
    return return_image ? SilhouetteResult(shadow, (g.x0, g.x1), (g.y0, g.y1), area) : area
end

"""
    silhouette_rasterized!(shadow, body::CompositeBody, sun_direction)

[`silhouette_rasterized`](@ref) drawn into `shadow`, a square `Matrix{Bool}`
whose size is the resolution. Returns the area, and leaves `shadow` as the image
of the body's shadow: `true` at each pixel the shadow falls on.
"""
silhouette_rasterized!(shadow::AbstractMatrix{Bool}, body::CompositeBody, sun_direction::NTuple{3,<:Real}) =
    first(_silhouette_rasterized!(shadow, body, sun_direction))

function _silhouette_rasterized!(shadow, body, sun_direction)
    resolution = size(shadow, 1)
    part_tiles = _part_tiles(body)
    g = _view(part_tiles, _normalize3(sun_direction), resolution)
    fill!(shadow, false)
    foreach_triangle(part_tiles) do p1, p2, p3
        _rasterize_triangle!(shadow, _project(g, p1), _project(g, p2), _project(g, p3),
                             g.x0, g.y0, g.dx, g.dy, resolution)
    end
    (count(shadow) * g.dx * g.dy * u"m^2", g)
end

"""
    silhouette(body::CompositeBody, beam::Beam; resolution=256)

Per-part lit (unshadowed) silhouette area facing `beam` — a `NamedTuple` keyed
like `body.parts`, each a `Quantity` (m²). Every part's posed mesh is projected
along the beam direction and depth-buffered, so at each pixel only the frontmost
part (nearest the source) is counted. A part occluded by another — a
ground-facing half under a sky-facing half toward an overhead sun — therefore
reports (near) zero, and the parts' areas sum to the composite silhouette (no
double counting), unlike the per-part analytic `silhouette`. [`silhouette!`](@ref)
uses a depth buffer you pass in instead.
"""
silhouette(body::CompositeBody, src::Beam; resolution::Integer = 256) =
    silhouette!(_depth_buffer(body, resolution), body, src)

"""
    silhouette!(depth, body::CompositeBody, beam::Beam)

[`silhouette`](@ref) toward `beam` using `depth`, a `resolution × resolution ×
nparts` array of `Float64`, as the depth buffer. It is left holding each part's
depth toward the beam, `-Inf` where the part doesn't cover a pixel.
"""
function silhouette!(depth::AbstractArray{Float64,3}, body::CompositeBody, src::Beam)
    resolution = size(depth, 1)
    part_tiles = _part_tiles(body)
    g = _view(part_tiles, src.direction, resolution)
    _part_depths!(depth, g, part_tiles)
    cell = g.dx * g.dy * u"m^2"
    _per_part(body) do k
        count(_front(depth, i, j) == k for i in 1:resolution, j in 1:resolution) * cell
    end
end


# Direction `i` of `n` near-uniform directions on the unit sphere (Fibonacci
# spiral). Deterministic (no RNG) so view partitions are reproducible.
function _fibonacci_direction(i, n)
    z = 1 - 2 * (i - 0.5) / n
    r = sqrt(max(0.0, 1 - z * z))
    sθ, cθ = sincos(π * (3 - sqrt(5.0)) * (i - 1))
    (r * cθ, r * sθ, z)
end

_check_size(size) = 0 <= size <= 1 ||
    throw(ArgumentError("a sky or ground size must be a fraction of the sphere in [0, 1], got $size"))

"""
    Sky(size)
    Sky(size, tilt)

The sky as a region of the direction sphere: a spherical cap covering `size` (a
fraction `0`–`1` of the whole sphere), centred straight up. `Sky(0.5)` is a
flat-ground hemisphere, `Sky(0.7)` a mountaintop (sky bulging past horizontal),
`Sky(0.0)` a sealed burrow. Pass a `tilt` direction (a 3-tuple) to point the cap off
vertical — a slope whose sky leans downhill. `Ground` is its complement; hand either
to [`silhouette_factors`](@ref) to split every direction into sky and ground.
`Sky(Unchecked(), size, tilt)` doesn't check its arguments.
"""
struct Sky{S,A}
    size::S
    axis::A
    function Sky(::Unchecked, size::S, tilt::NTuple{3,<:Real}) where {S}
        axis = _normalize3(tilt)
        new{S,typeof(axis)}(size, axis)
    end
end
Sky(size, tilt::NTuple{3,<:Real}) = (_check_size(size); _check_direction(tilt); Sky(Unchecked(), size, tilt))
Sky(size) = Sky(size, (0.0, 0.0, 1.0))

"""
    Ground(size)
    Ground(size, tilt)

The ground as a region of the direction sphere: a spherical cap covering `size` (a
fraction `0`–`1` of the whole sphere), centred straight down. `Ground(g)` is the
complement of `Sky(1 - g)`, so `Ground(0.5)` is flat ground and `Ground(1.0)` a
sealed burrow. `tilt` points the cap off vertical for a slope.
`Ground(Unchecked(), size, tilt)` doesn't check its arguments.
"""
struct Ground{S,A}
    size::S
    axis::A
    function Ground(::Unchecked, size::S, tilt::NTuple{3,<:Real}) where {S}
        axis = _normalize3(tilt)
        new{S,typeof(axis)}(size, axis)
    end
end
Ground(size, tilt::NTuple{3,<:Real}) = (_check_size(size); _check_direction(tilt); Ground(Unchecked(), size, tilt))
Ground(size) = Ground(size, (0.0, 0.0, -1.0))

"""
    Horizon(angles)

The measured horizon as a per-azimuth elevation profile — the real, possibly
asymmetric sky/ground boundary that a symmetric `Sky`/`Ground` cap cannot express.
`angles` is a vector of `n` `Unitful` elevations above horizontal, one per compass
bearing: bin `i` is azimuth `(i-1)·360°/n`, from north (`+y`) clockwise, matched to a
direction by nearest bin — exactly the `horizon_angles` a Microclimate.jl `Site`
carries (`fill(0.0u"°", 24)` for flat ground). Sky is every direction above the
profile, ground every direction below.
"""
struct Horizon{A}
    angles::A
    function Horizon(angles::A) where {A}
        isempty(angles) && throw(ArgumentError("Horizon needs at least one elevation angle"))
        new{A}(angles)
    end
end

# Sky share of a unit direction `d` (0 or 1) under a region; the ground share is its
# complement `1 - sky`. A spherical cap of sphere-fraction `f` has half-angle `θ` with
# `cos θ = 1 - 2f`, so cap membership is the plain dot-product test `d · axis ≥ 1 - 2·size`.
_sky_share(r::Sky, d)    = (d[1]*r.axis[1] + d[2]*r.axis[2] + d[3]*r.axis[3]) >= 1 - 2*r.size ? 1.0 : 0.0
_sky_share(r::Ground, d) = (d[1]*r.axis[1] + d[2]*r.axis[2] + d[3]*r.axis[3]) >= 1 - 2*r.size ? 0.0 : 1.0
function _sky_share(r::Horizon, d)
    n = length(r.angles)
    elevation = asin(clamp(d[3], -1.0, 1.0))                # radians above horizontal
    azimuth = atan(d[1], d[2])                              # from +y (north), clockwise
    azimuth < 0 && (azimuth += 2π)
    idx = mod(round(Int, azimuth / (2π / n)), n) + 1        # nearest bin (Microclimate rule)
    elevation >= ustrip(u"rad", @inbounds r.angles[idx]) ? 1.0 : 0.0
end

"""
    silhouette_factors(body::CompositeBody, region; ndirections=256, resolution=96)

Occlusion-aware partition of every part's radiative view into `sky`, `ground`, and
per-`neighbour` fractions that **sum to 1**. A `NamedTuple` keyed like `body.parts`;
each entry is `(; sky, ground, neighbours)` where `neighbours` is a `NamedTuple` keyed
like `body.parts`, with a part's share of itself 0.

Integrates the depth-buffered per-part silhouette over a Fibonacci-sphere set of
directions: toward each direction a part's *unoccluded* projected area counts as sky
or ground per `region`, while the projected area it loses accrues to the neighbour
directly in front of it there. So the blocked solid angle is never lost — it becomes the
part-to-part term — and the shares exhaust the sphere. General for any parts at any
pose (no shape/orientation assumption); internal mated faces score zero because the
neighbour buries them in the depth buffer.

`region` fixes the sky/ground split of each direction and is one of [`Sky`](@ref),
[`Ground`](@ref), or [`Horizon`](@ref) — an idealised cap from either side, or the
measured horizon profile. The two are complementary, so a single object suffices.

The fractions are true view factors only for convex parts: a part is treated as one
depth layer per direction, so a concave part's view of itself is not resolved.

`ndirections` sets the angular quadrature, `resolution` the raster grid; both trade
accuracy for cost. Runs once per pose/solar configuration (an `init!`-time quantity).
[`silhouette_factors!`](@ref) uses buffers you pass in instead.
"""
function silhouette_factors(body::CompositeBody, region;
                            ndirections::Integer = 256, resolution::Integer = 96)
    npart = length(body.parts)
    silhouette_factors!(_depth_buffer(body, resolution), Matrix{Float64}(undef, npart, npart + 2),
                        body, region; ndirections)
end

"""
    silhouette_factors!(depth, shares, body::CompositeBody, region; ndirections=256)

[`silhouette_factors`](@ref) using `depth`, a `resolution × resolution × nparts`
array of `Float64`, as the depth buffer, and `shares`, an `nparts × (nparts + 2)`
matrix, to sum each part's unnormalised sky, ground and neighbour shares in its
columns: sky, ground, then one per part.
"""
function silhouette_factors!(depth::AbstractArray{Float64,3}, shares::AbstractMatrix{Float64},
                             body::CompositeBody, region; ndirections::Integer = 256)
    resolution = size(depth, 1)
    part_tiles = _part_tiles(body) # only the projection changes per direction
    npart = length(body.parts)
    fill!(shares, 0.0)
    sky = view(shares, :, 1)
    ground = view(shares, :, 2)
    neighbour = view(shares, :, 3:npart + 2)

    for n in 1:ndirections
        d = _fibonacci_direction(n, ndirections)
        g = _view(part_tiles, d, resolution)
        _part_depths!(depth, g, part_tiles)
        # Pixel area; the shared dω cancels in the per-part normalisation, so it is omitted.
        cell = g.dx * g.dy
        sky_fraction = _sky_share(region, d) # 0–1; split each exposed pixel sky/ground

        # At each pixel every covering part looks toward the source: the frontmost
        # sees sky/ground, and every other part sees the part *directly* in front of
        # it — the covering part with the smallest depth beyond its own. (For
        # source → A → B → C, C sees B, not A.)
        @inbounds for j in 1:resolution, i in 1:resolution, p in 1:npart
            dp = depth[i, j, p]
            dp == -Inf && continue
            blocker, nearest = 0, Inf
            for k in 1:npart
                dk = depth[i, j, k]
                k != p && dp < dk < nearest && ((blocker, nearest) = (k, dk))
            end
            if blocker == 0
                sky[p] += cell * sky_fraction
                ground[p] += cell * (1 - sky_fraction)
            else
                neighbour[p, blocker] += cell
            end
        end
    end

    # Normalise each part's shares to sum to 1. A part is its own neighbour with share 0.
    _per_part(body) do p
        total = sky[p] + ground[p]
        for k in 1:npart
            total += neighbour[p, k]
        end
        scale = total > 0 ? 1 / total : 0.0
        (; sky = sky[p] * scale, ground = ground[p] * scale,
           neighbours = _per_part(k -> neighbour[p, k] * scale, body))
    end
end
