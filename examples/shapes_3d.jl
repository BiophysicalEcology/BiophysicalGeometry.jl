# shapes_3d.jl
#
# 3D surface-mesh cutaway plot for a single BiophysicalGeometry body.
#
# The outer layers are drawn over 270° (leaving one quarter open) so the
# interior is visible.  The innermost flesh core is always drawn in full (360°).
#
# Edit shape, insulation, and interactive at the top.
#
# Run from the repository root:
#   julia --project=examples examples/shapes_3d.jl
#
# Output: examples/shapes_3d.png  (when interactive = false)

using BiophysicalGeometry
using Unitful

# ── Backend ───────────────────────────────────────────────────────────────────
# false → save PNG with CairoMakie
# true  → open a rotatable GLMakie window (close it or press Enter to exit)
interactive = false

if interactive
    using GLMakie
else
    using CairoMakie
end

# Explicitly import BiophysicalGeometry types that conflict with GeometryBasics
# (which CairoMakie / GLMakie re-exports: Sphere, Cylinder, Ellipsoid).
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid

# ── Parameters ────────────────────────────────────────────────────────────────
mass    = 1.0u"kg"
density = 1000.0u"kg/m^3"
b_ratio = 4.0   # length-to-diameter ratio (used by Cylinder and Ellipsoid)

shape = Cylinder(mass, density, b_ratio)
# shape = Sphere(mass, density)
# shape = Ellipsoid(mass, density, b_ratio, b_ratio)

insulation = CompositeInsulation(Fur(10.0u"mm", 30.0u"μm", 3000u"cm^-2"),
                                  Fat(0.1, 901.0u"kg/m^3"))
# insulation = Fur(10.0u"mm", 30.0u"μm", 3000u"cm^-2")
# insulation = Fat(0.1, 901.0u"kg/m^3")
# insulation = Naked()

body = Body(shape, insulation)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plot_body(body)

if interactive
    screen = display(fig)
    println("Rotate with left-click drag, zoom with scroll.  Close the window to exit.")
    wait(screen)
else
    outfile = joinpath(@__DIR__, "shapes_3d.png")
    save(outfile, fig; px_per_unit=2)
    println("Saved: $outfile")
end
