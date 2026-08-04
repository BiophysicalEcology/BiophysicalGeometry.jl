# insulation_properties.jl
#
# Visualises the properties of fur/fibre insulation as modelled in
# BiophysicalGeometry.jl.  Four panels:
#
#   1. FibrousLayer schematic — side-view cartoon showing fibre diameter d, fur depth
#      (thickness), fibre length L, and the inter-fibre spacing implied by N.
#      When L > fur depth the fibres are drawn as tilted parallelograms.
#
#   2. Coverage fraction heatmap — fraction of skin area occupied by fibre
#      cross-sections as a function of d and N.
#      Formula: f = π (d/2)² × N
#      Contour lines mark f = 0.25, 0.50, 0.75, 1.0.
#
#   3. Surface area comparison — bar chart of total, skin and evaporation areas
#      for a Cylinder under the four insulation configurations.
#
#   4. Silhouette area vs. solar zenith angle — how the projected area seen by
#      the sun changes with angle θ (θ = 0° = normal, θ = 90° = parallel).
#
# Run from the repository root:
#   julia --project=examples examples/insulation_properties.jl
#
# Output: examples/insulation_properties.png

using CairoMakie
using BiophysicalGeometry
using Unitful

# Explicitly import BiophysicalGeometry types that conflict with GeometryBasics
# (which CairoMakie re-exports: Sphere, Cylinder, Ellipsoid).
import BiophysicalGeometry: Cylinder

# ── Parameters ────────────────────────────────────────────────────────────────
mass    = 1.0u"kg"
density = 1000.0u"kg/m^3"
b_ratio = 4.0

insulation_depth = 10.0u"mm"
fibre_length     = 14.0u"mm"   # actual fibre length (may exceed fur depth → fibres tilt)
fibre_diameter   = 30.0u"μm"
fibre_density    = 3000u"cm^-2"

fur  = FibrousLayer(insulation_depth, fibre_diameter, fibre_density)
fat  = FatLayer(0.1, 901.0u"kg/m^3")
comp = CompositeInsulation(fur, fat)

cyl        = Cylinder(mass, density, b_ratio)
ins_list   = [Naked(), fat, fur, comp]
ins_labels = ["Naked", "FatLayer", "FibrousLayer", "FibrousLayer + FatLayer"]

# ── Panel 3: Surface area comparison bar chart ────────────────────────────────

function draw_area_comparison!(ax, shape, insulations, labels)
    total_vec = Float64[]
    skin_vec = Float64[]
    evap_vec = Float64[]

    for ins in insulations
        body = Body(shape, ins)
        push!(total_vec, ustrip(u"cm^2", total_area(body)))
        push!(skin_vec, ustrip(u"cm^2", skin_area(body)))
        push!(evap_vec, ustrip(u"cm^2", evaporation_area(body)))
    end

    xs = 1:length(insulations)
    width = 0.26

    barplot!(ax, xs .- width, total_vec; width, color=:steelblue,
             label="Total (outer surface)")
    barplot!(ax, xs, skin_vec;  width, color=:salmon,
             label="Skin (under insulation)")
    barplot!(ax, xs .+ width, evap_vec;  width, color=:olivedrab,
             label="Evaporation (skin − hair)")

    ax.xticks = (collect(xs), labels)
    ax.title = "Cylinder surface areas by insulation type  (mass = $(shape.mass))"
    ax.xlabel = "Insulation configuration"
    ax.ylabel = "Area (cm²)"
    axislegend(ax; position=:lt, labelsize=9, framevisible=false)
    hidespines!(ax, :t, :r)
end

# ── Panel 4: Silhouette area vs. solar zenith angle ───────────────────────────

function draw_silhouette_vs_zenith!(ax, shape, insulations, labels)
    colours = [:grey50, :goldenrod, :sienna, :firebrick]
    linestyles = [:solid, :dash, :dot, :dashdot]
    θ_deg = LinRange(0.0, 90.0, 181)

    for (ins, lbl, col, ls) in zip(insulations, labels, colours, linestyles)
        body = Body(shape, ins)
        sil = [ustrip(u"cm^2",
                       silhouette_area(body, ZenithAngleVarying(), θ * u"°"))
                for θ in θ_deg]
        lines!(ax, θ_deg, sil; color=col, linestyle=ls, linewidth=2, label=lbl)
    end

    vlines!(ax, [0.0];  color=:black, linestyle=:dot, linewidth=0.8)
    vlines!(ax, [90.0]; color=:black, linestyle=:dot, linewidth=0.8)
    text!(ax, 2.0, 0.5; text="Normal\nto sun", fontsize=8, align=(:left, :bottom))
    text!(ax, 88.0, 0.5; text="Parallel\nto sun", fontsize=8, align=(:right, :bottom))

    ax.title = "Cylinder silhouette area vs. solar zenith angle"
    ax.xlabel = "Solar zenith angle θ (°)"
    ax.ylabel = "Silhouette area (cm²)"
    axislegend(ax; position=:rt, labelsize=9, framevisible=false)
    hidespines!(ax, :t, :r)
end

# ── Assemble figure ───────────────────────────────────────────────────────────
fig = Figure(size=(1120, 920), backgroundcolor=:white)

Label(fig[0, 1:3],
      "BiophysicalGeometry.jl — Insulation Properties";
      fontsize=16, font=:bold, padding=(0,0,10,0))

ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[2, 1])
ax4 = Axis(fig[2, 2:3])

draw_insulation_schematic!(ax1, fur; fibre_length)
hm = draw_insulation_coverage!(ax2, fur)
draw_area_comparison!(ax3, cyl, ins_list, ins_labels)
draw_silhouette_vs_zenith!(ax4, cyl, [Naked(), fur, comp], ["Naked", "FibrousLayer", "FibrousLayer + FatLayer"])

Colorbar(fig[1, 3], hm; label="Coverage fraction f", width=14, labelsize=10)

rowgap!(fig.layout, 16)
colgap!(fig.layout, 12)
colsize!(fig.layout, 3, Auto(0.05))

outfile = joinpath(@__DIR__, "insulation_properties.png")
save(outfile, fig; px_per_unit=2)
println("Saved: $outfile")
