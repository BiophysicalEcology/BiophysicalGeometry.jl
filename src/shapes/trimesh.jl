"""
    TriMesh(vertices, faces; mass) <: AbstractShape

A free-form triangular mesh — represents any closed-surface shape, given by
a list of vertex positions and triangle indices.

- `vertices` : `Vector{NTuple{3, <:Length}}` — vertex positions (e.g. metres).
- `faces`    : `Vector{NTuple{3, Int}}` — 1-based triangle vertex indices.
- `mass`     : organism mass; `density` is auto-computed from the mesh volume.

Volume and surface area are computed directly from the mesh (signed
tetrahedra and triangle areas respectively). The mesh should be closed and
have consistent triangle winding.

This shape is `attachment_surfaces = ()` by default — it can stand on its
own as a `Body` but cannot be joined into a `CompositeBody` without the
user defining their own named attachment surfaces.

For accurate silhouette area, wrap a `TriMesh` body in a one-part
`CompositeBody` and call `silhouette_rasterized`.
"""
mutable struct TriMesh{V,F,M,D} <: AbstractShape
    vertices::V
    faces::F
    mass::M
    density::D
end

function TriMesh(vertices, faces; mass)
    vol = _trimesh_volume(vertices, faces)
    density = mass / vol
    TriMesh(vertices, faces, mass, density)
end

# ── Closed-mesh integrals ─────────────────────────────────────────────────

# Signed-tetrahedron volume sum: |Σ (v1 × v2) · v3 / 6|.
function _trimesh_volume(verts, faces)
    s = zero(verts[1][1]^3)
    for (i, j, k) in faces
        v1, v2, v3 = verts[i], verts[j], verts[k]
        cx = v1[2]*v2[3] - v1[3]*v2[2]
        cy = v1[3]*v2[1] - v1[1]*v2[3]
        cz = v1[1]*v2[2] - v1[2]*v2[1]
        s += (cx*v3[1] + cy*v3[2] + cz*v3[3]) / 6
    end
    abs(s)
end

# Sum of triangle areas via cross products.
function _trimesh_surface_area(verts, faces)
    A = zero(verts[1][1]^2)
    for (i, j, k) in faces
        v1, v2, v3 = verts[i], verts[j], verts[k]
        e1 = (v2[1]-v1[1], v2[2]-v1[2], v2[3]-v1[3])
        e2 = (v3[1]-v1[1], v3[2]-v1[2], v3[3]-v1[3])
        cx = e1[2]*e2[3] - e1[3]*e2[2]
        cy = e1[3]*e2[1] - e1[1]*e2[3]
        cz = e1[1]*e2[2] - e1[2]*e2[1]
        A += sqrt(cx^2 + cy^2 + cz^2) / 2
    end
    A
end

# ── Geometry constructors ────────────────────────────────────────────────

function geometry(shape::TriMesh, ::Naked)
    volume = _trimesh_volume(shape.vertices, shape.faces)
    total = _trimesh_surface_area(shape.vertices, shape.faces)
    characteristic_dimension = volume^(1/3)
    return Geometry(volume, characteristic_dimension, (;),
                    SurfaceAreas(; total))
end

# Insulation other than Naked is not (yet) modelled for free-form meshes —
# growing fur uniformly outward needs vertex-normal displacement which is
# out of scope for now.
geometry(shape::TriMesh, ::Fur) =
    error("TriMesh does not yet support Fur insulation; use Naked.")
geometry(shape::TriMesh, ::Fat) =
    error("TriMesh does not yet support Fat insulation; use Naked.")
geometry(shape::TriMesh, ::CompositeInsulation) =
    error("TriMesh does not yet support composite insulation; use Naked.")

# ── Surface / silhouette / radii ──────────────────────────────────────────

function surface_area(shape::TriMesh, body::AbstractBody)
    body.geometry.area.total
end

# Coarse single-shape silhouette: characteristic-dimension² as a placeholder.
# For accurate values, wrap in a CompositeBody and use silhouette_rasterized.
function silhouette_area(::TriMesh, ::Naked, body::AbstractBody)
    cd = body.geometry.characteristic_dimension
    (; normal = cd^2, parallel = cd^2)
end
function silhouette_area(::TriMesh, ::Naked, body::AbstractBody, θ)
    body.geometry.characteristic_dimension^2
end

# Radii fall back to half the characteristic dimension. Defined for Naked
# only; insulation variants would require vertex-normal expansion which
# isn't implemented (matching geometry's Naked-only support).
skin_radius(::TriMesh, ::Naked, body::AbstractBody) =
    body.geometry.characteristic_dimension / 2
insulation_radius(::TriMesh, ::Naked, body::AbstractBody) =
    body.geometry.characteristic_dimension / 2
flesh_radius(::TriMesh, ::Naked, body::AbstractBody) =
    body.geometry.characteristic_dimension / 2

skin_radius(::TriMesh, ::AbstractInsulation, ::AbstractBody) =
    error("TriMesh does not yet support insulated bodies; use Naked.")
insulation_radius(::TriMesh, ::AbstractInsulation, ::AbstractBody) =
    error("TriMesh does not yet support insulated bodies; use Naked.")
flesh_radius(::TriMesh, ::AbstractInsulation, ::AbstractBody) =
    error("TriMesh does not yet support insulated bodies; use Naked.")

# ── Composition opt-out ──────────────────────────────────────────────────

attachment_surfaces(::TriMesh) = ()
