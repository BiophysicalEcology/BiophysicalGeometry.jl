# Rasterised silhouette area for a CompositeBody.
#
# Project every part's outer mesh triangles onto a plane perpendicular to
# the sun direction, rasterise into a 2D bitmap, and count covered pixels.
# This handles part-on-part overlap correctly (unlike the per-part summed
# fallback in composition.jl, which double-counts shadowed regions).

# Pick two orthonormal basis vectors perpendicular to a unit vector `d`.
# Aligns the screen-`v` axis with world `+z` whenever possible (so the
# image axis labels stay meaningful and the basis varies continuously as
# `d` sweeps the upper hemisphere). Falls back to a `+x`-aligned reference
# when `d` is nearly vertical, to avoid a singularity at the zenith.
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
    n > 0 || error("sun direction must be non-zero")
    (d[1]/n, d[2]/n, d[3]/n)
end

# Rasterise one 2D triangle into a z-buffer, tagging each won pixel with `part`.
# `depth` is the triangle's distance along the view direction (larger = nearer the
# source); a pixel is claimed only when this triangle is in front of whatever holds
# it. This is the per-part occlusion core: the frontmost part at each pixel wins, so
# a shadowed part contributes nothing there. Signature otherwise mirrors
# `_rasterize_triangle!`.
function _rasterize_triangle_zbuf!(part_buf, depth_buf, part, depth, p1, p2, p3, x0, y0, dx, dy, n)
    s = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    abs(s) < 1e-18 && return  # degenerate
    sgn = sign(s)

    xmin_f = min(p1[1], p2[1], p3[1])
    xmax_f = max(p1[1], p2[1], p3[1])
    ymin_f = min(p1[2], p2[2], p3[2])
    ymax_f = max(p1[2], p2[2], p3[2])
    imin = max(1, floor(Int, (xmin_f - x0) / dx) + 1)
    imax = min(n, ceil(Int, (xmax_f - x0) / dx))
    jmin = max(1, floor(Int, (ymin_f - y0) / dy) + 1)
    jmax = min(n, ceil(Int, (ymax_f - y0) / dy))

    @inbounds for j in jmin:jmax
        y = y0 + (j - 0.5) * dy
        for i in imin:imax
            depth <= depth_buf[i, j] && continue   # something nearer already holds this pixel
            x = x0 + (i - 0.5) * dx
            e1 = sgn * ((p2[1] - p1[1]) * (y - p1[2]) - (p2[2] - p1[2]) * (x - p1[1]))
            e1 < 0 && continue
            e2 = sgn * ((p3[1] - p2[1]) * (y - p2[2]) - (p3[2] - p2[2]) * (x - p2[1]))
            e2 < 0 && continue
            e3 = sgn * ((p1[1] - p3[1]) * (y - p3[2]) - (p1[2] - p3[2]) * (x - p3[1]))
            e3 < 0 && continue
            depth_buf[i, j] = depth
            part_buf[i, j] = part
        end
    end
end

# Rasterise one 2D triangle into `covered`. (x0, y0) is grid origin, (dx, dy)
# the cell size, and `n` the side length. Pixels whose centre lies inside
# the triangle are flagged true.
function _rasterize_triangle!(covered, p1, p2, p3, x0, y0, dx, dy, n)
    s = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    abs(s) < 1e-18 && return  # degenerate
    sgn = sign(s)

    xmin_f = min(p1[1], p2[1], p3[1])
    xmax_f = max(p1[1], p2[1], p3[1])
    ymin_f = min(p1[2], p2[2], p3[2])
    ymax_f = max(p1[2], p2[2], p3[2])
    imin = max(1, floor(Int, (xmin_f - x0) / dx) + 1)
    imax = min(n, ceil(Int, (xmax_f - x0) / dx))
    jmin = max(1, floor(Int, (ymin_f - y0) / dy) + 1)
    jmax = min(n, ceil(Int, (ymax_f - y0) / dy))

    @inbounds for j in jmin:jmax
        y = y0 + (j - 0.5) * dy
        for i in imin:imax
            covered[i, j] && continue
            x = x0 + (i - 0.5) * dx
            e1 = sgn * ((p2[1] - p1[1]) * (y - p1[2]) - (p2[2] - p1[2]) * (x - p1[1]))
            e1 < 0 && continue
            e2 = sgn * ((p3[1] - p2[1]) * (y - p2[2]) - (p3[2] - p2[2]) * (x - p2[1]))
            e2 < 0 && continue
            e3 = sgn * ((p1[1] - p3[1]) * (y - p3[2]) - (p1[2] - p3[2]) * (x - p3[1]))
            e3 < 0 && continue
            covered[i, j] = true
        end
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
    d = _normalize3(sun_direction)
    u, v = _ortho_basis(d)

    # Collect 2D triangles + bbox first, then rasterise.
    tris = NTuple{3, NTuple{2, Float64}}[]
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    for name in propertynames(body.parts)
        part = getfield(body.parts, name)
        pose = getfield(body.poses, name)
        for grid in _part_outer_meshes(part.shape, part, 1.0)  # sc=1 → metres
            X, Y, Z = _transform_mesh(grid..., pose, 1.0)
            for tri in _each_triangle(X, Y, Z)
                p1, p2, p3 = tri
                q1 = (p1[1]*u[1] + p1[2]*u[2] + p1[3]*u[3],
                      p1[1]*v[1] + p1[2]*v[2] + p1[3]*v[3])
                q2 = (p2[1]*u[1] + p2[2]*u[2] + p2[3]*u[3],
                      p2[1]*v[1] + p2[2]*v[2] + p2[3]*v[3])
                q3 = (p3[1]*u[1] + p3[2]*u[2] + p3[3]*u[3],
                      p3[1]*v[1] + p3[2]*v[2] + p3[3]*v[3])
                push!(tris, (q1, q2, q3))
                xmin = min(xmin, q1[1], q2[1], q3[1])
                xmax = max(xmax, q1[1], q2[1], q3[1])
                ymin = min(ymin, q1[2], q2[2], q3[2])
                ymax = max(ymax, q1[2], q2[2], q3[2])
            end
        end
    end
    isempty(tris) && return return_image ?
        SilhouetteResult(falses(resolution, resolution), (0.0, 0.0), (0.0, 0.0), 0.0u"m^2") :
        0.0u"m^2"

    # 2% margin so triangles don't clip pixel-grid boundaries.
    pad = 0.02 * max(xmax - xmin, ymax - ymin)
    x0 = xmin - pad;  x1 = xmax + pad
    y0 = ymin - pad;  y1 = ymax + pad
    dx = (x1 - x0) / resolution
    dy = (y1 - y0) / resolution

    covered = falses(resolution, resolution)
    for (q1, q2, q3) in tris
        _rasterize_triangle!(covered, q1, q2, q3, x0, y0, dx, dy, resolution)
    end
    area = count(covered) * dx * dy * u"m^2"
    return return_image ? SilhouetteResult(covered, (x0, x1), (y0, y1), area) : area
end

"""
    silhouette_area_per_part(body::CompositeBody, view_direction; resolution=256)

Per-part lit (unshadowed) silhouette area projected along `view_direction` — a
`NamedTuple` keyed like `body.parts`, each a `Quantity` (m²). Every part's posed
mesh is projected and rasterised into a shared depth buffer, so at each pixel only
the frontmost part (nearest along `view_direction`) is counted. A part occluded by
another — a ground-facing half under a sky-facing half toward an overhead sun —
therefore reports (near) zero, and the parts' areas sum to the composite silhouette
(no double counting), unlike the per-part `silhouette_area`.

Run it toward the sun for direct-beam exposure, or toward the sky / ground
hemispheres for the diffuse view fractions each part sees.
"""
function silhouette_area_per_part(body::CompositeBody, view_direction::NTuple{3,<:Real};
                                  resolution::Integer = 256)
    d = _normalize3(view_direction)
    u, v = _ortho_basis(d)
    names = propertynames(body.parts)

    # (part index, projected 2D triangle, depth along d). Depth uses the triangle
    # centroid: larger = nearer the source, so it wins the pixel.
    tris = Tuple{Int,NTuple{3,NTuple{2,Float64}},Float64}[]
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    for (pidx, name) in enumerate(names)
        part = getfield(body.parts, name)
        pose = getfield(body.poses, name)
        for grid in _part_outer_meshes(part.shape, part, 1.0)
            X, Y, Z = _transform_mesh(grid..., pose, 1.0)
            for tri in _each_triangle(X, Y, Z)
                p1, p2, p3 = tri
                q1 = (p1[1]*u[1] + p1[2]*u[2] + p1[3]*u[3], p1[1]*v[1] + p1[2]*v[2] + p1[3]*v[3])
                q2 = (p2[1]*u[1] + p2[2]*u[2] + p2[3]*u[3], p2[1]*v[1] + p2[2]*v[2] + p2[3]*v[3])
                q3 = (p3[1]*u[1] + p3[2]*u[2] + p3[3]*u[3], p3[1]*v[1] + p3[2]*v[2] + p3[3]*v[3])
                depth = ((p1[1]+p2[1]+p3[1])*d[1] + (p1[2]+p2[2]+p3[2])*d[2] + (p1[3]+p2[3]+p3[3])*d[3]) / 3
                push!(tris, (pidx, (q1, q2, q3), depth))
                xmin = min(xmin, q1[1], q2[1], q3[1]); xmax = max(xmax, q1[1], q2[1], q3[1])
                ymin = min(ymin, q1[2], q2[2], q3[2]); ymax = max(ymax, q1[2], q2[2], q3[2])
            end
        end
    end
    if isempty(tris)
        return NamedTuple{names}(ntuple(_ -> 0.0u"m^2", length(names)))
    end

    pad = 0.02 * max(xmax - xmin, ymax - ymin)
    x0 = xmin - pad; x1 = xmax + pad
    y0 = ymin - pad; y1 = ymax + pad
    dx = (x1 - x0) / resolution
    dy = (y1 - y0) / resolution

    part_buf  = zeros(Int, resolution, resolution)
    depth_buf = fill(-Inf, resolution, resolution)
    for (pidx, (q1, q2, q3), depth) in tris
        _rasterize_triangle_zbuf!(part_buf, depth_buf, pidx, depth, q1, q2, q3, x0, y0, dx, dy, resolution)
    end

    cell = dx * dy
    areas = ntuple(length(names)) do k
        count(==(k), part_buf) * cell * u"m^2"
    end
    return NamedTuple{names}(areas)
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
    view_partition(body::CompositeBody; ndirections=256, resolution=96)

Occlusion-aware partition of every part's radiative view into `sky`, `ground`, and
per-`neighbour` fractions that **sum to 1**. A `NamedTuple` keyed like `body.parts`;
each entry is `(; sky, ground, neighbours)` where `neighbours` is a `NamedTuple` keyed
by the *other* parts.

Integrates the depth-buffered per-part silhouette over a Fibonacci-sphere set of
directions: toward each direction a part's *unoccluded* projected area counts as sky
(direction above the horizon) or ground (below), while the projected area it loses to a
frontmost neighbour `k` accrues to that neighbour. So the blocked solid angle is never
lost — it becomes the part-to-part term — and the three shares exhaust the hemisphere.
General for any parts at any pose (no shape/orientation assumption); internal mated
faces score zero because the neighbour buries them in the depth buffer.

`ndirections` sets the angular quadrature, `resolution` the raster grid; both trade
accuracy for cost. Runs once per pose/solar configuration (an `init!`-time quantity).
"""
function view_partition(body::CompositeBody; ndirections::Integer = 256, resolution::Integer = 96)
    names = propertynames(body.parts)
    npart = length(names)

    # Each part's posed 3D triangles, collected once (only the projection changes per direction).
    part_tris = [NTuple{3,NTuple{3,Float64}}[] for _ in 1:npart]
    for (pidx, name) in enumerate(names)
        part = getfield(body.parts, name)
        pose = getfield(body.poses, name)
        for grid in _part_outer_meshes(part.shape, part, 1.0)
            X, Y, Z = _transform_mesh(grid..., pose, 1.0)
            for tri in _each_triangle(X, Y, Z)
                push!(part_tris[pidx], tri)
            end
        end
    end

    sky   = zeros(npart)
    grnd  = zeros(npart)
    neigh = zeros(npart, npart)

    for d in _fibonacci_sphere(ndirections)
        u, v = _ortho_basis(d)
        # Project every part's triangles onto the plane ⟂ d, tracking a shared bbox.
        proj = [Tuple{NTuple{2,Float64},NTuple{2,Float64},NTuple{2,Float64},Float64}[] for _ in 1:npart]
        xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
        for pidx in 1:npart, (p1, p2, p3) in part_tris[pidx]
            q1 = (p1[1]*u[1]+p1[2]*u[2]+p1[3]*u[3], p1[1]*v[1]+p1[2]*v[2]+p1[3]*v[3])
            q2 = (p2[1]*u[1]+p2[2]*u[2]+p2[3]*u[3], p2[1]*v[1]+p2[2]*v[2]+p2[3]*v[3])
            q3 = (p3[1]*u[1]+p3[2]*u[2]+p3[3]*u[3], p3[1]*v[1]+p3[2]*v[2]+p3[3]*v[3])
            depth = ((p1[1]+p2[1]+p3[1])*d[1] + (p1[2]+p2[2]+p3[2])*d[2] + (p1[3]+p2[3]+p3[3])*d[3]) / 3
            push!(proj[pidx], (q1, q2, q3, depth))
            xmin = min(xmin, q1[1], q2[1], q3[1]); xmax = max(xmax, q1[1], q2[1], q3[1])
            ymin = min(ymin, q1[2], q2[2], q3[2]); ymax = max(ymax, q1[2], q2[2], q3[2])
        end
        isfinite(xmin) || continue

        pad = 0.02 * max(xmax - xmin, ymax - ymin)
        x0 = xmin - pad; y0 = ymin - pad
        dx = (xmax + pad - x0) / resolution
        dy = (ymax + pad - y0) / resolution
        cell = dx * dy

        # Frontmost part per pixel (depth buffer over all triangles).
        part_buf  = zeros(Int, resolution, resolution)
        depth_buf = fill(-Inf, resolution, resolution)
        for pidx in 1:npart, (q1, q2, q3, depth) in proj[pidx]
            _rasterize_triangle_zbuf!(part_buf, depth_buf, pidx, depth, q1, q2, q3, x0, y0, dx, dy, resolution)
        end

        # For each part: pixels it wins → sky/ground exposure; pixels it covers but a
        # neighbour won → that neighbour's part-to-part share. `cell` weights each pixel
        # by area; the shared `dω` cancels in the per-part normalisation, so it is omitted.
        upper = d[3] >= 0
        cov = falses(resolution, resolution)
        for pidx in 1:npart
            fill!(cov, false)
            for (q1, q2, q3, _) in proj[pidx]
                _rasterize_triangle!(cov, q1, q2, q3, x0, y0, dx, dy, resolution)
            end
            @inbounds for j in 1:resolution, i in 1:resolution
                cov[i, j] || continue
                w = part_buf[i, j]
                if w == pidx
                    (upper ? sky : grnd)[pidx] += cell
                elseif w != 0
                    neigh[pidx, w] += cell
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
