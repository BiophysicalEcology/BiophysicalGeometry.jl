# cross_sections_2d.jl
#
# 2D cross-section diagrams for a single BiophysicalGeometry body, built up
# the same way as the R plot_human example: draw layers outer-to-inner on true
# x-y axes with asp = 1 (DataAspect).
#
# Geometry stays as Unitful quantities throughout; units are stripped only at
# the Point2f boundary where CairoMakie requires plain numbers.
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

#shape      = Cylinder(mass, density, 4.0)
shape    = Ellipsoid(mass, density, 4.0, 1.0)
# shape    = Sphere(mass, density)

insulation = CompositeInsulation(Fur(10.0u"mm", 30.0u"μm", 3000u"cm^-2"),
                                  Fat(0.1, 901.0u"kg/m^3"))
# insulation = Fur(10.0u"mm", 30.0u"μm", 3000u"cm^-2")
# insulation = Fat(0.1, 901.0u"kg/m^3")
# insulation = Naked()

body = Body(shape, insulation)

# ── Colours ────────────────────────────────────────────────────────────────────
flesh_col = RGBf(0.88, 0.48, 0.42)
fat_col   = RGBf(1.00, 0.97, 0.60)
fur_col   = RGBf(0.76, 0.62, 0.42)

# ── Coordinate helpers ─────────────────────────────────────────────────────────
# Accept Unitful lengths; strip to cm only at Point2f.
# (CairoMakie requires plain Float32 coordinates — this is the only place
#  unit information is discarded.  cm gives readable axis tick values.)
pu(x) = ustrip(u"cm", x)

circle_pts(r; n=300) =
    [Point2f(pu(r)*cos(t), pu(r)*sin(t)) for t in LinRange(0, 2π, n+1)]

ellipse_pts(xr, yr; n=300) =
    [Point2f(pu(xr)*cos(t), pu(yr)*sin(t)) for t in LinRange(0, 2π, n+1)]

rect_pts(hw, hh) =
    [Point2f(-pu(hw), -pu(hh)), Point2f(pu(hw), -pu(hh)),
     Point2f( pu(hw),  pu(hh)), Point2f(-pu(hw),  pu(hh))]

# Draw a filled polygon; strokecolor matches fill to suppress CairoMakie's halo.
layer!(ax, pts, col) = poly!(ax, pts, color=col, strokecolor=col, strokewidth=1)

# ── Build plots ────────────────────────────────────────────────────────────────
shape_name = string(nameof(typeof(body.shape)))
ins_name   = string(nameof(typeof(body.insulation)))

fig = Figure(backgroundcolor=:white)
Label(fig[0, 1:2],
      "$(shape_name) · $(ins_name)  (mass = $(mass))";
      fontsize=13, font=:bold, padding=(0,0,8,0))

ax_long = Axis(fig[1, 1];
               title="Longitudinal section",
               xlabel="x (cm)", ylabel="y (cm)",
               aspect=DataAspect())
ax_tran = Axis(fig[1, 2];
               title="Transverse section",
               xlabel="x (cm)", ylabel="y (cm)",
               aspect=DataAspect())

r_f = flesh_radius(body)       # Unitful, metres
r_s = skin_radius(body)        # Unitful, metres
r_i = insulation_radius(body)  # Unitful, metres

# ── Longitudinal section ───────────────────────────────────────────────────────
if body.shape isa Cylinder
    gl   = body.geometry.length
    hl_s = gl.length_skin / 2
    hl_i = hasproperty(gl, :length_fur) ? gl.length_fur / 2 : hl_s
    # Draw outermost layer first; each inner layer paints over it
    r_i > r_s && layer!(ax_long, rect_pts(r_i, hl_i), fur_col)
    r_s > r_f && layer!(ax_long, rect_pts(r_s, hl_s), fat_col)
                 layer!(ax_long, rect_pts(r_f, hl_s), flesh_col)
    ey = hl_i
else  # Ellipsoid or Sphere
    ratio = hasproperty(body.shape, :b) ? Float64(body.shape.b) : 1.0
    r_i > r_s && layer!(ax_long, ellipse_pts(r_i, r_i * ratio), fur_col)
    r_s > r_f && layer!(ax_long, ellipse_pts(r_s, r_s * ratio), fat_col)
                 layer!(ax_long, ellipse_pts(r_f, r_f * ratio), flesh_col)
    ey = r_i * ratio
end

# ── Transverse section ─────────────────────────────────────────────────────────
r_i > r_s && layer!(ax_tran, circle_pts(r_i), fur_col)
r_s > r_f && layer!(ax_tran, circle_pts(r_s), fat_col)
             layer!(ax_tran, circle_pts(r_f), flesh_col)

# ── Axis limits (same y-scale on both panels for true to-scale comparison) ─────
pad = 0.12
xlims!(ax_long, -pu(r_i)*(1+pad),  pu(r_i)*(1+pad))
ylims!(ax_long, -pu(ey)*(1+pad),   pu(ey)*(1+pad))
xlims!(ax_tran, -pu(r_i)*(1+pad),  pu(r_i)*(1+pad))
ylims!(ax_tran, -pu(ey)*(1+pad),   pu(ey)*(1+pad))   # same y-scale as longitudinal

# ── Legend ─────────────────────────────────────────────────────────────────────
Legend(fig[2, 1:2],
    [PolyElement(polycolor=flesh_col, strokecolor=flesh_col, strokewidth=1),
     PolyElement(polycolor=fat_col,   strokecolor=fat_col,   strokewidth=1),
     PolyElement(polycolor=fur_col,   strokecolor=fur_col,   strokewidth=1)],
    ["Flesh / muscle", "Subcutaneous fat", "Fur / fibre insulation"];
    orientation=:horizontal, framevisible=false)

rowgap!(fig.layout, 8)
colgap!(fig.layout, 30)

outfile = joinpath(@__DIR__, "cross_sections_2d.png")
save(outfile, fig; px_per_unit=2)
println("Saved: $outfile")
