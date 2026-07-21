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
