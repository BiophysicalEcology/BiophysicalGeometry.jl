# human.jl
#
# Multi-part human body matching NicheMapR's plot_human.R default parameters:
#   MASS = 70 kg, density = 1050 kg/m³
#   FATPCTs = c(5, 36, 10, 23) × 0.7  (% mass as fat)
#   Fur (hair) depths: 10 mm head, 6 mm body
#   SHAPE_Bs = c(1.6, 1.9, 12, 7)
#
# Run from the repository root:
#   julia --project=examples examples/human.jl
#
# Output: examples/human.png  (when interactive = false)

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

import BiophysicalGeometry: Cylinder, Ellipsoid

# ── Body parameters (matching NicheMapR plot_human.R defaults) ────────────────
mass    = 70.0u"kg"
density = 1050.0u"kg/m^3"     # DENSITYs = rep(1050, 4)
fat_ρ   = 900.0u"kg/m^3"

# FATPCTs = c(5, 36, 10, 23) * 0.7  → mass fractions
fat_head  = Fat(0.035, fat_ρ)   #  3.5 %
fat_torso = Fat(0.252, fat_ρ)   # 25.2 %
fat_arm   = Fat(0.070, fat_ρ)   #  7.0 %
fat_leg   = Fat(0.161, fat_ρ)   # 16.1 %

# INSDEPDs = c(1e-2, 6e-3, 6e-3, 6e-3) m  (dorsal hair depth)
# DHAIRDs  = c(7.5e-5, 1e-6, 1e-6, 1e-6) m (hair diameter)
# INSDENDs = rep(3e8, 4) m⁻²               (hair density)
fur_head  = Fur(10e-3u"m", 7.5e-5u"m", 3e8u"m^-2")
fur_body  = Fur( 6e-3u"m",   1e-6u"m", 3e8u"m^-2")

human = Organism(
    parts = (
        head  = BodyPart(Body(Ellipsoid(mass * 0.0761, density, 1.6, 1.0),
                              CompositeInsulation(fur_head, fat_head)), 1),
        torso = BodyPart(Body(Cylinder(mass * 0.501,  density, 1.9),
                              CompositeInsulation(fur_body, fat_torso)), 1),
        arms  = BodyPart(Body(Cylinder(mass * 0.049,  density, 12.0),
                              CompositeInsulation(fur_body, fat_arm)), 2),
        legs  = BodyPart(Body(Cylinder(mass * 0.162,  density, 7.0),
                              CompositeInsulation(fur_body, fat_leg)), 2),
    ),
    joins = (
        Join(:torso, :head, :top,     0.15),
        Join(:torso, :arms, :lateral, 0.0),
        Join(:torso, :legs, :bottom,  0.0),
    )
)

# ── Geometry summary ──────────────────────────────────────────────────────────
geom = geometry(human)
println("Human geometry (70 kg, NicheMapR parameters):")
println("  Total surface area : ", round(typeof(1.0u"cm^2"), geom.total_area, digits=1))
println("  Total volume       : ", round(typeof(1.0u"L"),    geom.total_volume * 1000u"L/m^3", digits=2))
println("\nArea fractions by part:")
for (name, frac) in pairs(geom.area_fractions)
    @printf "  %-8s  %.1f %%\n" name (frac * 100)
end

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plot_organism(human, Upright())

if interactive
    screen = display(fig)
    println("Close the window to exit.")
    wait(screen)
else
    outfile = joinpath(@__DIR__, "human.png")
    save(outfile, fig; px_per_unit=2)
    println("\nSaved: $outfile")
end
