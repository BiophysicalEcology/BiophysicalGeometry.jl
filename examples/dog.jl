# dog.jl
#
# Multi-part quadruped (dog) body: head + torso + 4 legs, prone posture.
# Head is an Ellipsoid; torso and legs are Cylinders.
#
# Run from the repository root:
#   julia --project=examples examples/dog.jl
#
# Output: examples/dog.png  (when interactive = false)

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

# ── Body parameters ───────────────────────────────────────────────────────────
mass    = 25.0u"kg"
density = 990.0u"kg/m^3"
fat_ρ   = 900.0u"kg/m^3"

# Typical short-haired dog fur: 20 mm depth, 50 µm diameter, 1e8 hairs/m²
fur_dog  = Fur(20e-3u"m", 5e-5u"m", 1e8u"m^-2")
fat_body = Fat(0.10, fat_ρ)   # 10% body fat
fat_head = Fat(0.05, fat_ρ)   #  5% head fat
fat_legs = Fat(0.05, fat_ρ)   #  5% leg fat

dog = Organism(
    parts = (
        head       = BodyPart(Body(Ellipsoid(mass * 0.10,  density, 2.0, 1.0), CompositeInsulation(fur_dog, fat_head)), 1),
        torso      = BodyPart(Body(Cylinder( mass * 0.52,  density, 3.5),       CompositeInsulation(fur_dog, fat_body)), 1),
        legs_fore  = BodyPart(Body(Cylinder( mass * 0.04,  density, 6.0),       CompositeInsulation(fur_dog, fat_legs)), 2),
        legs_hind  = BodyPart(Body(Cylinder( mass * 0.05,  density, 5.5),       CompositeInsulation(fur_dog, fat_legs)), 2),
        tail       = BodyPart(Body(Cylinder( mass * 0.01,  density, 8.0),       CompositeInsulation(fur_dog, Fat(0.02, fat_ρ))), 1),
        ears       = BodyPart(Body(Plate(    mass * 0.003, density, 0.6, 2.5),  CompositeInsulation(fur_dog, fat_head)), 2),
    ),
    joins = (
        Join(:torso, :head,      :forward, 0.5),           # head in front, neck truncation at 50%
        Join(:torso, :legs_fore, :lateral, 0.0, π/2, 0.7), # fore legs at shoulder (70% toward front)
        Join(:torso, :legs_hind, :lateral, 0.0, π/2, -0.7),# hind legs at hip (70% toward rear)
        Join(:torso, :tail,      :back,    0.0),            # tail extends from rear of torso
        Join(:head,  :ears,      :top,     0.0),            # ears sit on top of head
    )
)

# ── Geometry summary ──────────────────────────────────────────────────────────
geom = geometry(dog)
println("Dog geometry (25 kg):")
println("  Total surface area : ", round(typeof(1.0u"cm^2"), geom.total_area, digits=1))
println("  Total volume       : ", round(typeof(1.0u"L"),    geom.total_volume * 1000u"L/m^3", digits=2))
println("\nArea fractions by part:")
for (name, frac) in pairs(geom.area_fractions)
    @printf "  %-12s  %.1f %%\n" name (frac * 100)
end

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plot_organism(dog, Prone())

if interactive
    screen = display(fig)
    println("Close the window to exit.")
    wait(screen)
else
    outfile = joinpath(@__DIR__, "dog.png")
    save(outfile, fig; px_per_unit=2)
    println("\nSaved: $outfile")
end
