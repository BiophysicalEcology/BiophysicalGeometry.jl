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
# shape = Ellipsoid(mass, density, b_ratio, 1.0)

insulation = CompositeInsulation(Fur(10.0u"mm", 30.0u"μm", 3000u"cm^-2"),
                                  Fat(0.1, 901.0u"kg/m^3"))
# insulation = Fur(10.0u"mm", 30.0u"μm", 3000u"cm^-2")
# insulation = Fat(0.1, 901.0u"kg/m^3")
# insulation = Naked()

sc = 100.0  # m → cm

# RGBA colours for each layer
flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00)
fat_col   = RGBAf(1.00, 0.97, 0.60, 0.75)  # semi-transparent
fur_col   = RGBAf(0.76, 0.62, 0.42, 0.45)  # more transparent outer layer

# ── Parametric mesh generators ────────────────────────────────────────────────

function cylinder_tube(r, L; nθ=72, nz=2)
    θ = LinRange(0.0, 2π, nθ);  z = LinRange(0.0, Float64(L), nz)
    [r*cos(θi) for θi in θ, _ in z], [r*sin(θi) for θi in θ, _ in z],
    [zi        for _  in θ, zi in z]
end

function cylinder_tube_partial(r, L; nθ=60, θ_end=3π/2, z0=0.0)
    θ = LinRange(0.0, θ_end, nθ);  z = [z0, z0 + Float64(L)]
    [r*cos(θi) for θi in θ, _ in z], [r*sin(θi) for θi in θ, _ in z],
    [zi        for _  in θ, zi in z]
end

function cylinder_cap(r, z0; nθ=72, nr=12)
    θ = LinRange(0.0, 2π, nθ);  rv = LinRange(0.0, r, nr)
    [rr*cos(θi) for θi in θ, rr in rv], [rr*sin(θi) for θi in θ, rr in rv],
    fill(Float64(z0), nθ, nr)
end

function cylinder_cap_partial(r, z0; nθ=60, nr=12, θ_end=3π/2)
    θ = LinRange(0.0, θ_end, nθ);  rv = LinRange(0.0, r, nr)
    [rr*cos(θi) for θi in θ, rr in rv], [rr*sin(θi) for θi in θ, rr in rv],
    fill(Float64(z0), nθ, nr)
end

function sphere_mesh(r; n=60)
    θ = LinRange(0.0, 2π, n);  φ = LinRange(0.0, π, n)
    [r*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [r*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [r*cos(φj)         for _  in θ, φj in φ]
end

function sphere_mesh_partial(r; n=60, θ_end=3π/2)
    θ = LinRange(0.0, θ_end, n);  φ = LinRange(0.0, π, n)
    [r*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [r*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [r*cos(φj)         for _  in θ, φj in φ]
end

function ellipsoid_mesh(a, b; n=60)
    θ = LinRange(0.0, 2π, n);  φ = LinRange(0.0, π, n)
    [a*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [b*cos(φj)         for _  in θ, φj in φ]
end

function ellipsoid_mesh_partial(a, b; n=60, θ_end=3π/2)
    θ = LinRange(0.0, θ_end, n);  φ = LinRange(0.0, π, n)
    [a*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [b*cos(φj)         for _  in θ, φj in φ]
end

# ── Surface helper ────────────────────────────────────────────────────────────
function draw_surface!(ax, mesh_tuple, color)
    X, Y, Z = mesh_tuple
    surface!(ax, X, Y, Z;
             color=fill(color, size(X)...), shading=true, backlight=0.4f0)
end

# ── Unified cutaway drawing ───────────────────────────────────────────────────
function draw_cutaway!(ax, body)
    r_f = ustrip(u"m", flesh_radius(body))      * sc
    r_s = ustrip(u"m", skin_radius(body))       * sc
    r_i = ustrip(u"m", insulation_radius(body)) * sc
    has_fur = r_i > r_s + 1e-9
    has_fat = r_s > r_f + 1e-9

    if body.shape isa Cylinder
        gl    = body.geometry.length
        L_s   = ustrip(u"m", gl.length_skin) * sc
        L_i   = hasproperty(gl, :length_fur) ? ustrip(u"m", gl.length_fur) * sc : L_s
        z0_fur = -(L_i - L_s) / 2

        draw_surface!(ax, cylinder_tube(r_f, L_s), flesh_col)
        draw_surface!(ax, cylinder_cap(r_f, 0.0),  flesh_col)
        draw_surface!(ax, cylinder_cap(r_f, L_s),  flesh_col)
        if has_fat
            draw_surface!(ax, cylinder_tube_partial(r_s, L_s), fat_col)
            draw_surface!(ax, cylinder_cap_partial(r_s, 0.0),  fat_col)
            draw_surface!(ax, cylinder_cap_partial(r_s, L_s),  fat_col)
        end
        if has_fur
            draw_surface!(ax, cylinder_tube_partial(r_i, L_i; z0=z0_fur), fur_col)
            draw_surface!(ax, cylinder_cap_partial(r_i, z0_fur),           fur_col)
            draw_surface!(ax, cylinder_cap_partial(r_i, z0_fur + L_i),     fur_col)
        end

    elseif body.shape isa Sphere
        draw_surface!(ax, sphere_mesh(r_f), flesh_col)
        has_fat && draw_surface!(ax, sphere_mesh_partial(r_s), fat_col)
        has_fur && draw_surface!(ax, sphere_mesh_partial(r_i), fur_col)

    else  # Ellipsoid: flesh_radius returns minor semi-axis; major = minor * shape.b
        ratio = Float64(body.shape.b)
        draw_surface!(ax, ellipsoid_mesh(r_f * ratio, r_f), flesh_col)
        has_fat && draw_surface!(ax, ellipsoid_mesh_partial(r_s * ratio, r_s), fat_col)
        has_fur && draw_surface!(ax, ellipsoid_mesh_partial(r_i * ratio, r_i), fur_col)
    end

    ax.xlabel = "x (cm)"; ax.ylabel = "y (cm)"; ax.zlabel = "z (cm)"
end

# ── Figure ────────────────────────────────────────────────────────────────────
body       = Body(shape, insulation)
shape_name = string(nameof(typeof(body.shape)))
ins_name   = string(nameof(typeof(body.insulation)))

fig = Figure(size=(600, 560), backgroundcolor=:white)

Label(fig[0, 1],
      "BiophysicalGeometry.jl — $(shape_name) · $(ins_name)  (quarter cutaway)";
      fontsize=13, font=:bold, padding=(0,0,10,0))

ax = Axis3(fig[1, 1];
           perspectiveness=0.5, viewmode=:fit,
           elevation=π/7, azimuth=5π/4)
draw_cutaway!(ax, body)

Legend(fig[2, 1],
    [PolyElement(polycolor=flesh_col, strokecolor=:saddlebrown, strokewidth=1),
     PolyElement(polycolor=fat_col,   strokecolor=:saddlebrown, strokewidth=1),
     PolyElement(polycolor=fur_col,   strokecolor=:black,       strokewidth=1)],
    ["Flesh / muscle (full 360°)",
     "Subcutaneous fat (270° cutaway)",
     "Fur / insulation (270° cutaway)"];
    orientation=:horizontal, framevisible=false)

if interactive
    display(fig)
    println("Rotate with left-click drag, zoom with scroll.  Press Enter to exit.")
    readline()
else
    outfile = joinpath(@__DIR__, "shapes_3d.png")
    save(outfile, fig; px_per_unit=2)
    println("Saved: $outfile")
end
