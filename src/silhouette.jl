# Rasterised silhouette area for a CompositeBody.
#
# Project every part's outer mesh triangles onto a plane perpendicular to
# the sun direction, rasterise into a 2D bitmap, and count covered pixels.
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
    all(isfinite, d) || throw(ArgumentError("direction must be finite, got $d"))
    n = sqrt(d[1]^2 + d[2]^2 + d[3]^2)
    n > 0 || throw(ArgumentError("direction must be non-zero, got $d"))
    (d[1]/n, d[2]/n, d[3]/n)
end

_check_count(name, n) =
    n > 0 || throw(ArgumentError("`$name` must be positive, got $n"))

"""
    Beam(direction)

A collimated beam from one `direction` — parallel rays from a source at infinity,
i.e. the direct (beam) component of radiation, as opposed to the diffuse `Sky` /
`Ground` regions. `silhouette(body, Beam(dir))` gives each part's projected area
facing the beam.
"""
struct Beam
    direction::NTuple{3,Float64}
    Beam(d::NTuple{3,<:Real}) = new(_normalize3(d))
end
Beam(x::Real, y::Real, z::Real) = Beam((x, y, z))

# Edge-function rasteriser shared by the coverage and depth passes. Calls
# `f(i, j, w1, w2, w3)` for every pixel whose centre lies inside the 2D triangle
# (p1, p2, p3), with (w1, w2, w3) the pixel centre's barycentric weights. (x0, y0)
# is the grid origin, (dx, dy) the cell size and `n` the side length.
@inline function _each_pixel(f, p1, p2, p3, x0, y0, dx, dy, n)
    s = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    abs(s) < 1e-18 && return  # degenerate
    sgn = sign(s)
    inv_area = 1 / abs(s)

    imin = max(1, floor(Int, (min(p1[1], p2[1], p3[1]) - x0) / dx) + 1)
    imax = min(n, ceil(Int, (max(p1[1], p2[1], p3[1]) - x0) / dx))
    jmin = max(1, floor(Int, (min(p1[2], p2[2], p3[2]) - y0) / dy) + 1)
    jmax = min(n, ceil(Int, (max(p1[2], p2[2], p3[2]) - y0) / dy))

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
_rasterize_triangle!(covered, p1, p2, p3, x0, y0, dx, dy, n) =
    _each_pixel((i, j, _, _, _) -> (@inbounds covered[i, j] = true), p1, p2, p3, x0, y0, dx, dy, n)

# Rasterise a triangle into a depth buffer, keeping the nearest depth per pixel.
# Depth is interpolated across the triangle from its vertex depths `z`
# (distance along the view direction; larger = nearer the source), so long
# triangles — a tube's lateral strips span the whole part — order correctly
# against other parts at every pixel, not just at their centroids.
function _rasterize_depth!(depth_buf, p1, p2, p3, z, x0, y0, dx, dy, n)
    _each_pixel(p1, p2, p3, x0, y0, dx, dy, n) do i, j, w1, w2, w3
        depth = w1 * z[1] + w2 * z[2] + w3 * z[3]
        @inbounds depth > depth_buf[i, j] && (depth_buf[i, j] = depth)
    end
end

# ── Shared projection ─────────────────────────────────────────────────────

# Every part's posed outer-mesh triangles in world metres, collected once.
function _part_triangles(body::CompositeBody)
    map(propertynames(body.parts)) do name
        part = getfield(body.parts, name)
        pose = getfield(body.poses, name)
        tris = NTuple{3,NTuple{3,Float64}}[]
        for grid in _part_outer_meshes(part.shape, part, 1.0)  # sc=1 → metres
            X, Y, Z = _transform_mesh(grid..., pose, 1.0)
            append!(tris, _each_triangle(X, Y, Z))
        end
        tris
    end
end

# Project each part's triangles onto the plane ⟂ `d` and fit a
# `resolution × resolution` pixel grid around them all (2% margin so triangles
# don't clip the grid edge). Each projected triangle keeps its vertex depths
# along `d`. Returns `nothing` if there are no triangles.
function _project(part_tris, d, resolution)
    u, v = _ortho_basis(d)
    onto(p, a) = p[1]*a[1] + p[2]*a[2] + p[3]*a[3]
    proj = map(part_tris) do tris
        map(tris) do (p1, p2, p3)
            ((onto(p1, u), onto(p1, v)), (onto(p2, u), onto(p2, v)), (onto(p3, u), onto(p3, v)),
             (onto(p1, d), onto(p2, d), onto(p3, d)))
        end
    end
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    for tris in proj, (q1, q2, q3, _) in tris
        xmin = min(xmin, q1[1], q2[1], q3[1]); xmax = max(xmax, q1[1], q2[1], q3[1])
        ymin = min(ymin, q1[2], q2[2], q3[2]); ymax = max(ymax, q1[2], q2[2], q3[2])
    end
    isfinite(xmin) || return nothing
    pad = 0.02 * max(xmax - xmin, ymax - ymin)
    x0, x1 = xmin - pad, xmax + pad
    y0, y1 = ymin - pad, ymax + pad
    (; proj, x0, x1, y0, y1, dx = (x1 - x0) / resolution, dy = (y1 - y0) / resolution)
end

# One depth buffer per part (-Inf where the part doesn't cover a pixel). Keeping
# them separate lets each pixel be attributed to the part directly in front of
# another, not just to the frontmost overall.
function _part_depths(g, resolution)
    map(g.proj) do tris
        buf = fill(-Inf, resolution, resolution)
        for (q1, q2, q3, z) in tris
            _rasterize_depth!(buf, q1, q2, q3, z, g.x0, g.y0, g.dx, g.dy, resolution)
        end
        buf
    end
end

"""
    SilhouetteResult

Result of a rasterised silhouette projection.

- `bitmap` : `BitMatrix` of size `resolution × resolution`; `true` where the body covers that pixel.
- `x_range`, `y_range` : `(min, max)` axis extents (m) of the projection plane.
- `area` : silhouette area as a `Quantity{Float64, m²}`.
"""
struct SilhouetteResult
    bitmap::BitMatrix
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
projected, and rasterised into a `resolution × resolution` bitmap whose
covered pixels are summed.

Returns a length² `Quantity` (m²) by default, or a `SilhouetteResult`
including the bitmap and projection ranges if `return_image=true`.
Increase `resolution` for more accuracy.
"""
function silhouette_rasterized(body::CompositeBody, sun_direction::NTuple{3,<:Real};
                                resolution::Integer = 256,
                                return_image::Bool = false)
    _check_count("resolution", resolution)
    d = _normalize3(sun_direction)
    g = _project(_part_triangles(body), d, resolution)
    g === nothing && return return_image ?
        SilhouetteResult(falses(resolution, resolution), (0.0, 0.0), (0.0, 0.0), 0.0u"m^2") :
        0.0u"m^2"

    covered = falses(resolution, resolution)
    for tris in g.proj, (q1, q2, q3, _) in tris
        _rasterize_triangle!(covered, q1, q2, q3, g.x0, g.y0, g.dx, g.dy, resolution)
    end
    area = count(covered) * g.dx * g.dy * u"m^2"
    (x0, x1, y0, y1) = (g.x0, g.x1, g.y0, g.y1)
    return return_image ? SilhouetteResult(covered, (x0, x1), (y0, y1), area) : area
end

"""
    silhouette(body::CompositeBody, beam::Beam; resolution=256)

Per-part lit (unshadowed) silhouette area facing `beam` — a `NamedTuple` keyed
like `body.parts`, each a `Quantity` (m²). Every part's posed mesh is projected
along the beam direction and depth-buffered, so at each pixel only the frontmost
part (nearest the source) is counted. A part occluded by another — a
ground-facing half under a sky-facing half toward an overhead sun — therefore
reports (near) zero, and the parts' areas sum to the composite silhouette (no
double counting), unlike the per-part analytic `silhouette`.
"""
function silhouette(body::CompositeBody, src::Beam; resolution::Integer = 256)
    _check_count("resolution", resolution)
    names = propertynames(body.parts)
    g = _project(_part_triangles(body), src.direction, resolution)
    g === nothing && return NamedTuple{names}(ntuple(_ -> 0.0u"m^2", length(names)))

    depths = _part_depths(g, resolution)
    counts = zeros(Int, length(names))
    @inbounds for j in 1:resolution, i in 1:resolution
        front, best = 0, -Inf
        for p in eachindex(depths)
            depths[p][i, j] > best && ((front, best) = (p, depths[p][i, j]))
        end
        front > 0 && (counts[front] += 1)
    end
    cell = g.dx * g.dy * u"m^2"
    return NamedTuple{names}(ntuple(k -> counts[k] * cell, length(names)))
end

# Deterministic near-uniform directions on the unit sphere (Fibonacci spiral).
# Deterministic (no RNG) so view partitions are reproducible.
function _fibonacci_sphere(n::Integer)
    golden = π * (3 - sqrt(5.0))
    dirs = Vector{NTuple{3,Float64}}(undef, n)
    for i in 0:(n - 1)
        z = 1 - 2 * (i + 0.5) / n
        r = sqrt(max(0.0, 1 - z * z))
        θ = golden * i
        dirs[i + 1] = (r * cos(θ), r * sin(θ), z)
    end
    return dirs
end

"""
    Sky(size)
    Sky(size, tilt)

The sky as a region of the direction sphere: a spherical cap covering `size` (a
fraction `0`–`1` of the whole sphere), centred straight up. `Sky(0.5)` is a
flat-ground hemisphere, `Sky(0.7)` a mountaintop (sky bulging past horizontal),
`Sky(0.0)` a sealed burrow. Pass a `tilt` direction (a 3-tuple) to point the cap off
vertical — a slope whose sky leans downhill. `Ground` is its complement; hand either
to [`silhouette_factors`](@ref) to split every direction into sky and ground.
"""
_check_size(name, size) = 0 <= size <= 1 ||
    throw(ArgumentError("$name size must be a fraction of the sphere in [0, 1], got $size"))

struct Sky{S,A}
    size::S
    axis::A
    Sky(size::S, tilt) where {S} =
        (_check_size("Sky", size); a = _normalize3(tilt); new{S,typeof(a)}(size, a))
end
Sky(size) = Sky(size, (0.0, 0.0, 1.0))

"""
    Ground(size)
    Ground(size, tilt)

The ground as a region of the direction sphere: a spherical cap covering `size` (a
fraction `0`–`1` of the whole sphere), centred straight down. `Ground(g)` is the
complement of `Sky(1 - g)`, so `Ground(0.5)` is flat ground and `Ground(1.0)` a
sealed burrow. `tilt` points the cap off vertical for a slope.
"""
struct Ground{S,A}
    size::S
    axis::A
    Ground(size::S, tilt) where {S} =
        (_check_size("Ground", size); a = _normalize3(tilt); new{S,typeof(a)}(size, a))
end
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
    Horizon(angles::A) where {A} = isempty(angles) ?
        throw(ArgumentError("Horizon needs at least one elevation angle")) : new{A}(angles)
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
by the *other* parts.

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
"""
function silhouette_factors(body::CompositeBody, region;
                            ndirections::Integer = 256, resolution::Integer = 96)
    _check_count("ndirections", ndirections)
    _check_count("resolution", resolution)
    names = propertynames(body.parts)
    npart = length(names)
    part_tris = _part_triangles(body) # only the projection changes per direction

    sky   = zeros(npart)
    grnd  = zeros(npart)
    neigh = zeros(npart, npart)

    for d in _fibonacci_sphere(ndirections)
        g = _project(part_tris, d, resolution)
        g === nothing && continue
        depths = _part_depths(g, resolution)
        # Pixel area; the shared dω cancels in the per-part normalisation, so it is omitted.
        cell = g.dx * g.dy
        sky_fraction = _sky_share(region, d) # 0–1; split each exposed pixel sky/ground

        # At each pixel every covering part looks toward the source: the frontmost
        # sees sky/ground, and every other part sees the part *directly* in front of
        # it — the covering part with the smallest depth beyond its own. (For
        # source → A → B → C, C sees B, not A.)
        @inbounds for j in 1:resolution, i in 1:resolution
            for p in 1:npart
                dp = depths[p][i, j]
                dp == -Inf && continue
                blocker, nearest = 0, Inf
                for k in 1:npart
                    dk = depths[k][i, j]
                    k != p && dp < dk < nearest && ((blocker, nearest) = (k, dk))
                end
                if blocker == 0
                    sky[p]  += cell * sky_fraction
                    grnd[p] += cell * (1 - sky_fraction)
                else
                    neigh[p, blocker] += cell
                end
            end
        end
    end

    # Normalise each part's shares to sum to 1, and package neighbours by name.
    return NamedTuple{names}(ntuple(npart) do p
        total = sky[p] + grnd[p] + sum(@view neigh[p, :])
        inv = total > 0 ? 1 / total : 0.0
        other = filter(!=(names[p]), collect(names))
        nb = NamedTuple{Tuple(other)}(ntuple(length(other)) do j
            neigh[p, findfirst(==(other[j]), collect(names))] * inv
        end)
        (; sky = sky[p] * inv, ground = grnd[p] * inv, neighbours = nb)
    end)
end
