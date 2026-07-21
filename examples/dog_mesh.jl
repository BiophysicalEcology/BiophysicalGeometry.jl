# dog_mesh.jl
#
# A "TriMesh dog" — the same shape concept as `dog.jl`, but baked into a
# single triangulated mesh instead of a CompositeBody of primitives.
# Demonstrates the `TriMesh` shape: pass any closed mesh of vertex
# positions + triangle indices and it slots in as an `AbstractShape`.
#
# Run from the repository root:
#   julia --project=examples examples/dog_mesh.jl

using BiophysicalGeometry
using Unitful
using GLMakie
import BiophysicalGeometry: TriMesh

GLMakie.activate!()

# ── Build the mesh ────────────────────────────────────────────────────────
const VertT = NTuple{3, typeof(0.0u"m")}

function build_dog_mesh()
    verts = VertT[]
    faces = NTuple{3, Int}[]

    # Add a closed axis-aligned box to the mesh. Triangle winding follows
    # the right-hand rule so outward normals point away from the centre.
    function add_box!(min_pt, max_pt)
        x0, y0, z0 = min_pt
        x1, y1, z1 = max_pt
        i = length(verts)
        push!(verts, (x0, y0, z0))   # 1
        push!(verts, (x1, y0, z0))   # 2
        push!(verts, (x1, y1, z0))   # 3
        push!(verts, (x0, y1, z0))   # 4
        push!(verts, (x0, y0, z1))   # 5
        push!(verts, (x1, y0, z1))   # 6
        push!(verts, (x1, y1, z1))   # 7
        push!(verts, (x0, y1, z1))   # 8
        # 12 triangles, winding so outward normals point out of the box.
        push!(faces, (i+1, i+3, i+2));  push!(faces, (i+1, i+4, i+3))   # bottom (-z)
        push!(faces, (i+5, i+6, i+7));  push!(faces, (i+5, i+7, i+8))   # top    (+z)
        push!(faces, (i+1, i+2, i+6));  push!(faces, (i+1, i+6, i+5))   # front  (-y)
        push!(faces, (i+4, i+7, i+3));  push!(faces, (i+4, i+8, i+7))   # back   (+y)
        push!(faces, (i+1, i+5, i+8));  push!(faces, (i+1, i+8, i+4))   # left   (-x)
        push!(faces, (i+2, i+3, i+7));  push!(faces, (i+2, i+7, i+6))   # right  (+x)
    end

    # Body: a horizontal box, x = head/tail axis.
    add_box!((0.00u"m", -0.10u"m", 0.40u"m"),
             (0.60u"m",  0.10u"m", 0.60u"m"))
    # Head: a slightly smaller box just in front of the body.
    add_box!((0.60u"m", -0.09u"m", 0.43u"m"),
             (0.82u"m",  0.09u"m", 0.65u"m"))
    # Tail: a thin stub at the back.
    add_box!((-0.20u"m", -0.025u"m", 0.50u"m"),
             ( 0.00u"m",  0.025u"m", 0.55u"m"))
    # Four legs hanging down from the body to the ground (z=0).
    for (x, y) in ((0.10, -0.07), (0.10, 0.07),
                   (0.50, -0.07), (0.50, 0.07))
        add_box!((x*u"m" - 0.025u"m", y*u"m" - 0.025u"m", 0.0u"m"),
                 (x*u"m" + 0.025u"m", y*u"m" + 0.025u"m", 0.4u"m"))
    end

    return verts, faces
end

verts, faces = build_dog_mesh()
println("vertices: ", length(verts), "  faces: ", length(faces))

mesh_shape = TriMesh(verts, faces; mass = 20.0u"kg")
body       = Body(mesh_shape, Naked())

println("volume       = ", body.geometry.volume)
println("surface area = ", body.geometry.area.total)
println("density      = ", mesh_shape.density)

# Wrap in a single-part CompositeBody so we can use the interactive
# rasterised-silhouette plot.
struct DogBody end
dog = CompositeBody(; parts = (DogBody() => body,), joins = (), root = DogBody())

fig = plot_body_silhouette(dog; resolution = 128)
screen = display(fig)
println()
println("Drag the θ / φ sliders to rotate the sun direction.")
println("Close the window to exit.")
wait(screen)
