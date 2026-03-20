# insulation_properties.jl
#
# Visualises the properties of fur/fibre insulation as modelled in
# BiophysicalGeometry.jl.  Four panels:
#
#   1. Fur schematic — side-view cartoon showing fibre diameter d, fur depth
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

fur_thick    = 10.0u"mm"
fibre_length = 14.0u"mm"   # actual fibre length (may exceed fur depth → fibres tilt)
fibre_d      = 30.0u"μm"
fibre_n      = 3000u"cm^-2"

fur  = Fur(fur_thick, fibre_d, fibre_n)
fat  = Fat(0.1, 901.0u"kg/m^3")
comp = CompositeInsulation(fur, fat)

# ── Panel 1: Fur schematic ────────────────────────────────────────────────────
# A to-scale side view is impractical (fibres are 30 μm wide but ~10 mm tall),
# so the figure uses a display_d that is exaggerated for clarity while keeping
# the spacing and thickness accurate.  When fibre_length > fur_thick the fibres
# are drawn as tilted parallelograms; tilt angle = acos(fur_thick / fibre_length).

function draw_fur_schematic!(ax)
    thick_mm     = ustrip(u"mm", fur_thick)
    fibre_len_mm = ustrip(u"mm", fibre_length)
    d_μm         = ustrip(u"μm", fibre_d)
    n_cm2        = ustrip(u"cm^-2", fibre_n)

    # Square-packing inter-centre spacing: spacing = 1 / √(N_mm²)
    # N in mm⁻² = N_cm² / 100
    spacing_mm = 1.0 / sqrt(n_cm2 / 100.0)

    # Exaggerate fibre width in the plot so they are visible
    d_display = spacing_mm * 0.40    # 40 % of spacing

    # Horizontal offset of fibre top relative to base: dx = √(L² − h²)
    # When L == h the fibres are vertical (dx = 0); when L > h they lean.
    dx = sqrt(max(0.0, fibre_len_mm^2 - thick_mm^2))

    n_show = 8
    skin_h = thick_mm * 0.12          # height of skin band
    W      = n_show * spacing_mm      # horizontal extent of fibre bases

    # ── Skin surface ──────────────────────────────────────────────────────────
    poly!(ax,
          [Point2f(0, 0), Point2f(W, 0),
           Point2f(W, -skin_h), Point2f(0, -skin_h)],
          color=RGBf(0.88, 0.48, 0.42),
          strokecolor=:black, strokewidth=0.5)
    text!(ax, W / 2, -skin_h / 2;
          text="Skin surface", fontsize=9, align=(:center, :center))

    # ── Fibres (parallelograms when dx > 0) ───────────────────────────────────
    for i in 0:(n_show - 1)
        xc = (i + 0.5) * spacing_mm
        poly!(ax,
              [Point2f(xc - d_display/2,      0.0),
               Point2f(xc + d_display/2,      0.0),
               Point2f(xc + d_display/2 + dx, thick_mm),
               Point2f(xc - d_display/2 + dx, thick_mm)],
              color=RGBf(0.76, 0.62, 0.42),
              strokecolor=(:saddlebrown, 0.6), strokewidth=0.4)
    end

    # ── Fur-layer bounding box (parallelogram) ────────────────────────────────
    lines!(ax, [0.0, W, W + dx, dx, 0.0],
               [0.0, 0.0, thick_mm, thick_mm, 0.0];
           color=(:grey40, 0.5), linewidth=0.8, linestyle=:dash)

    # ── Fur depth annotation (vertical, right of bounding box) ────────────────
    x_ann = W + dx + 0.15 * W
    lines!(ax, [x_ann, x_ann], [0.0, thick_mm]; color=:black, linewidth=1)
    scatter!(ax, [x_ann, x_ann], [0.0, thick_mm];
             color=:black, markersize=6, marker=:rect)
    text!(ax, x_ann + 0.04*W, thick_mm / 2;
          text="fur depth\n$(round(Int, thick_mm)) mm",
          fontsize=9, align=(:left, :center))

    # ── Fibre length annotation (along right edge of a central fibre) ─────────
    if dx > 1e-6
        i_ann  = n_show ÷ 2
        xc_ann = (i_ann + 0.5) * spacing_mm
        xfl0   = xc_ann + d_display / 2
        xfl1   = xfl0 + dx
        tilt_deg = round(Int, atand(dx, thick_mm))
        lines!(ax, [xfl0, xfl1], [0.0, thick_mm];
               color=:purple, linewidth=1.4)
        scatter!(ax, [xfl0, xfl1], [0.0, thick_mm];
                 color=:purple, markersize=5, marker=:vline)
        text!(ax, (xfl0 + xfl1)/2, thick_mm * 1.02;
              text="L = $(round(fibre_len_mm, digits=1)) mm  ($(tilt_deg)° from vertical)",
              fontsize=8, align=(:center, :bottom), color=:purple)
    end

    # ── Diameter annotation (below first fibre) ───────────────────────────────
    y_d = -skin_h * 2.4
    x0  = 0.5 * spacing_mm - d_display / 2
    x1  = 0.5 * spacing_mm + d_display / 2
    lines!(ax, [x0, x1], [y_d, y_d]; color=:darkblue, linewidth=1.2)
    scatter!(ax, [x0, x1], [y_d, y_d]; color=:darkblue, markersize=5, marker=:vline)
    text!(ax, (x0 + x1)/2, y_d - 0.3;
          text="d = $(round(Int, d_μm)) μm\n(displayed $(round(Int, d_display*1e3)) μm for clarity)",
          fontsize=8, align=(:center, :top), color=:darkblue)

    # ── Spacing annotation (between fibre tops) ────────────────────────────────
    x_sp0 = dx + 0.5 * spacing_mm
    x_sp1 = dx + 1.5 * spacing_mm
    y_sp  = thick_mm * 1.22
    lines!(ax, [x_sp0, x_sp1], [y_sp, y_sp]; color=:darkgreen, linewidth=1.2)
    scatter!(ax, [x_sp0, x_sp1], [y_sp, y_sp]; color=:darkgreen, markersize=5, marker=:vline)
    text!(ax, (x_sp0 + x_sp1)/2, y_sp + 0.2;
          text="spacing ≈ $(round(spacing_mm, digits=2)) mm  (N = $n_cm2 cm⁻²)",
          fontsize=8, align=(:center, :bottom), color=:darkgreen)

    ax.title  = "Fur schematic  (fibre width exaggerated for clarity)"
    ax.xlabel = "position (mm)"
    ax.ylabel = "height above skin (mm)"
    ylims!(ax, -skin_h * 5, thick_mm * 1.55)
    xlims!(ax, -0.1 * W, x_ann + 0.45 * W)
    hidespines!(ax, :t, :r)
end

# ── Panel 2: Coverage fraction heatmap ───────────────────────────────────────

function draw_coverage_heatmap!(ax)
    d_range = LinRange(10.0, 120.0, 200)    # fibre diameter  (μm)
    N_range = LinRange(200.0, 9000.0, 200)  # fibre density   (cm⁻²)

    # Coverage fraction: f = π (d/2)² × N
    # d in m = d_μm × 1e-6;  N in m⁻² = N_cm2 × 1e4
    cov = [π * ((d_μm * 1e-6) / 2)^2 * (N_cm2 * 1e4)
           for d_μm in d_range, N_cm2 in N_range]  # size: (n_d, n_N)

    hm = heatmap!(ax, d_range, N_range, cov;
                  colormap=:plasma, colorrange=(0.0, 1.0), highclip=:white)

    # Contour lines at key coverage fractions
    for (level, clr) in [(0.25, :white), (0.50, :white),
                          (0.75, :white), (1.00, :cyan)]
        contour!(ax, d_range, N_range, cov;
                 levels=[level], color=clr, linewidth=1.2, linestyle=:dash)
        # Label contour at the right edge
        # Find N for which cov(d_range[end], N) ≈ level
        N_label = level / (π * ((d_range[end] * 1e-6) / 2)^2 * 1e4)
        N_label_cm2 = N_label * 1e-4
        if 200 < N_label_cm2 < 9000
            text!(ax, d_range[end] - 3, N_label_cm2;
                  text="f=$(level)",
                  fontsize=8, align=(:right, :bottom), color=clr)
        end
    end

    # Reference point for default fibre parameters
    d_ref = ustrip(u"μm", fibre_d)
    N_ref = ustrip(u"cm^-2", fibre_n)
    cov_ref = π * ((d_ref * 1e-6) / 2)^2 * (N_ref * 1e4)
    scatter!(ax, [d_ref], [N_ref]; color=:lime, markersize=11,
             strokecolor=:black, strokewidth=1)
    text!(ax, d_ref + 2, N_ref + 150;
          text="$(d_ref) μm, $(N_ref) cm⁻²\nf ≈ $(round(cov_ref, digits=3))",
          fontsize=8, color=:lime)

    ax.title  = "Fibre coverage fraction  f = π(d/2)² × N"
    ax.xlabel = "Fibre diameter d (μm)"
    ax.ylabel = "Fibre density N (cm⁻²)"

    return hm  # returned so the caller can attach a Colorbar
end

# ── Panel 3: Surface area comparison bar chart ────────────────────────────────

function draw_area_comparison!(ax)
    ins_list   = [Naked(), fat, fur, comp]
    ins_labels = ["Naked", "Fat", "Fur", "Fur + Fat"]
    cyl        = Cylinder(mass, density, b_ratio)

    total_vec = Float64[]
    skin_vec  = Float64[]
    evap_vec  = Float64[]

    for ins in ins_list
        body = Body(cyl, ins)
        push!(total_vec, ustrip(u"cm^2", total_area(body)))
        push!(skin_vec,  ustrip(u"cm^2", skin_area(body)))
        push!(evap_vec,  ustrip(u"cm^2", evaporation_area(body)))
    end

    xs    = 1:length(ins_list)
    width = 0.26

    barplot!(ax, xs .- width, total_vec; width, color=:steelblue,
             label="Total (outer surface)")
    barplot!(ax, xs,          skin_vec;  width, color=:salmon,
             label="Skin (under insulation)")
    barplot!(ax, xs .+ width, evap_vec;  width, color=:olivedrab,
             label="Evaporation (skin − hair)")

    ax.xticks = (collect(xs), ins_labels)
    ax.title  = "Cylinder surface areas by insulation type  (mass = $(mass))"
    ax.xlabel = "Insulation configuration"
    ax.ylabel = "Area (cm²)"
    axislegend(ax; position=:lt, labelsize=9, framevisible=false)
    hidespines!(ax, :t, :r)
end

# ── Panel 4: Silhouette area vs. solar zenith angle ───────────────────────────

function draw_silhouette_vs_zenith!(ax)
    ins_list   = [Naked(), fat, fur, comp]
    ins_labels = ["Naked", "Fat", "Fur", "Fur + Fat"]
    colours    = [:grey50, :goldenrod, :sienna, :firebrick]
    linestyles = [:solid, :dash, :dot, :dashdot]

    θ_deg = LinRange(0.0, 90.0, 181)
    cyl   = Cylinder(mass, density, b_ratio)

    for (ins, lbl, col, ls) in zip(ins_list, ins_labels, colours, linestyles)
        body = Body(cyl, ins)
        sil = [ustrip(u"cm^2",
                       silhouette_area(body, ZenithAngleVarying(), θ * u"°"))
               for θ in θ_deg]
        lines!(ax, θ_deg, sil;
               color=col, linestyle=ls, linewidth=2, label=lbl)
    end

    vlines!(ax, [0.0];  color=:black, linestyle=:dot, linewidth=0.8)
    vlines!(ax, [90.0]; color=:black, linestyle=:dot, linewidth=0.8)
    text!(ax,  2.0, 0.5; text="Normal\nto sun", fontsize=8,
          align=(:left, :bottom))
    text!(ax, 88.0, 0.5; text="Parallel\nto sun", fontsize=8,
          align=(:right, :bottom))

    ax.title  = "Cylinder silhouette area vs. solar zenith angle"
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

draw_fur_schematic!(ax1)
hm = draw_coverage_heatmap!(ax2)
draw_area_comparison!(ax3)
draw_silhouette_vs_zenith!(ax4)

# Colorbar for the heatmap (placed in column 3, row 1)
Colorbar(fig[1, 3], hm; label="Coverage fraction f", width=14, labelsize=10)

rowgap!(fig.layout, 16)
colgap!(fig.layout, 12)
colsize!(fig.layout, 3, Auto(0.05))  # narrow colorbar column

outfile = joinpath(@__DIR__, "insulation_properties.png")
save(outfile, fig; px_per_unit=2)
println("Saved: $outfile")
