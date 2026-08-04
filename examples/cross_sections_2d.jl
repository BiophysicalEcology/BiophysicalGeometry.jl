# cross_sections_2d.jl
#
# 2D cross-section diagrams for a single BiophysicalGeometry body.
#
# Edit shape and insulation to change the body being visualised.
#
# Run from the repository root:
#   julia --project=examples examples/cross_sections_2d.jl
#
# Output: examples/cross_sections_2d.png

using CairoMakie
using BiophysicalGeometry
using Unitful

import BiophysicalGeometry: Cylinder, Ellipsoid

# ── Body parameters ────────────────────────────────────────────────────────────
mass    = 1.0u"kg"
density = 1000.0u"kg/m^3"

shape      = Cylinder(mass, density, 4.0)
# shape    = Ellipsoid(mass, density, 4.0, 1.0)
# shape    = Sphere(mass, density)

insulation = CompositeInsulation(FibrousLayer(10.0u"mm", 30.0u"μm", 3000u"cm^-2"),
                                  FatLayer(0.1, 901.0u"kg/m^3"))
# insulation = FibrousLayer(10.0u"mm", 30.0u"μm", 3000u"cm^-2")
# insulation = FatLayer(0.1, 901.0u"kg/m^3")
# insulation = Naked()

body = Body(shape, insulation)

# ── Plot ───────────────────────────────────────────────────────────────────────
fig = plot_cross_sections(body)

outfile = joinpath(@__DIR__, "cross_sections_2d.png")
save(outfile, fig; px_per_unit=2)
println("Saved: $outfile")
