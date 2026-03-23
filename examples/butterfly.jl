# butterfly.jl
#
# Butterfly body: cylindrical thorax/abdomen + two flat wings (Plate).
# Demonstrates the dihedral_angle field of Join for wing elevation angle.
# Uses prone posture (body horizontal).
#
# Wing angle PHI (dihedral_angle in Join):
#   PHI = π/2  → wings horizontal
#   PHI > π/2  → open V-shape (e.g. 2π/3 = 120° = basking posture)
#   PHI < π/2  → tent/closed (e.g. π/3 = 60°)
#
# Run from the repository root:
#   julia --project=examples examples/butterfly.jl
#
# Output: examples/butterfly.png  (when interactive = false)

using BiophysicalGeometry
using Unitful
using Printf

# ── Backend ───────────────────────────────────────────────────────────────────
interactive = false

if interactive
    using GLMakie
else
    using CairoMakie
end

import BiophysicalGeometry: Cylinder, Plate

# ── Body parameters ───────────────────────────────────────────────────────────
# A medium-sized butterfly (e.g. Monarch, ~0.5 g)
mass    = 0.5e-3u"kg"
density = 800.0u"kg/m^3"

# Wing geometry: b = length/width, c = length/height (thin flat wing)
wing_b = 2.5   # wing length ≈ 2.5× wing width
wing_c = 8.0   # thin (height ≈ length/8)

# Wing angle: open basking posture, 120°
phi = deg2rad(120.0)

butterfly = Organism(
    parts = (
        body  = BodyPart(Body(Cylinder(mass * 0.30, density, 4.0), Naked()), 1),
        wings = BodyPart(Body(Plate(mass * 0.35, density, wing_b, wing_c), Naked()), 2),
    ),
    joins = (
        # Wings attach laterally; dihedral_angle = PHI (wing opening angle)
        Join(:body, :wings, :lateral, 0.0, phi),
    )
)

# ── Geometry summary ──────────────────────────────────────────────────────────
geom = geometry(butterfly)
println("Butterfly geometry (0.5 g, wing angle = $(round(rad2deg(phi), digits=1))°):")
println("  Total surface area : ", round(typeof(1.0u"cm^2"), geom.total_area, digits=3))
println("  Total volume       : ", round(typeof(1.0u"cm^3"), geom.total_volume * 1e6u"cm^3/m^3", digits=3))
println("\nArea fractions by part:")
for (name, frac) in pairs(geom.area_fractions)
    @printf "  %-8s  %.1f %%\n" name (frac * 100)
end

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plot_organism(butterfly, Prone())

if interactive
    screen = display(fig)
    println("Close the window to exit.")
    wait(screen)
else
    outfile = joinpath(@__DIR__, "butterfly.png")
    save(outfile, fig; px_per_unit=2)
    println("\nSaved: $outfile")
end
