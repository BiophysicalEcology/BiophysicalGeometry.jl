module BiophysicalGeometryMakieExt

using Makie
using Unitful
using BiophysicalGeometry
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Plate

# ══════════════════════════════════════════════════════════════════════════════
# 3-D MESH GENERATORS (internal)
# ══════════════════════════════════════════════════════════════════════════════

# ── Cylinder ──────────────────────────────────────────────────────────────────

function _cylinder_tube(r, L; nθ=72, nz=2)
    θ = LinRange(0.0, 2π, nθ);  z = LinRange(0.0, Float64(L), nz)
    [r*cos(θi) for θi in θ, _ in z],
    [r*sin(θi) for θi in θ, _ in z],
    [zi        for _  in θ, zi in z]
end

function _cylinder_tube_partial(r, L; nθ=60, θ_end=3π/2, z0=0.0)
    θ = LinRange(0.0, θ_end, nθ);  z = [z0, z0 + Float64(L)]
    [r*cos(θi) for θi in θ, _ in z],
    [r*sin(θi) for θi in θ, _ in z],
    [zi        for _  in θ, zi in z]
end

function _cylinder_cap(r, z0; nθ=72, nr=12)
    θ = LinRange(0.0, 2π, nθ);  rv = LinRange(0.0, r, nr)
    [rr*cos(θi) for θi in θ, rr in rv],
    [rr*sin(θi) for θi in θ, rr in rv],
    fill(Float64(z0), nθ, nr)
end

function _cylinder_cap_partial(r, z0; nθ=60, nr=12, θ_end=3π/2)
    θ = LinRange(0.0, θ_end, nθ);  rv = LinRange(0.0, r, nr)
    [rr*cos(θi) for θi in θ, rr in rv],
    [rr*sin(θi) for θi in θ, rr in rv],
    fill(Float64(z0), nθ, nr)
end

# ── Sphere ────────────────────────────────────────────────────────────────────

function _sphere_mesh(r; n=60)
    θ = LinRange(0.0, 2π, n);  φ = LinRange(0.0, π, n)
    [r*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [r*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [r*cos(φj)         for _  in θ, φj in φ]
end

function _sphere_mesh_partial(r; n=60, θ_end=3π/2)
    θ = LinRange(0.0, θ_end, n);  φ = LinRange(0.0, π, n)
    [r*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [r*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [r*cos(φj)         for _  in θ, φj in φ]
end

# ── Ellipsoid ─────────────────────────────────────────────────────────────────

function _ellipsoid_mesh(a, b; n=60)
    θ = LinRange(0.0, 2π, n);  φ = LinRange(0.0, π, n)
    [a*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [b*cos(φj)         for _  in θ, φj in φ]
end

function _ellipsoid_mesh_partial(a, b; n=60, θ_end=3π/2)
    θ = LinRange(0.0, θ_end, n);  φ = LinRange(0.0, π, n)
    [a*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [b*cos(φj)         for _  in θ, φj in φ]
end

# ── Box / Plate ───────────────────────────────────────────────────────────────

_box_face_z(x1, x2, y1, y2, z) =
    ([xi for xi in (x1, x2), _ in (y1, y2)],
     [yi for _ in (x1, x2), yi in (y1, y2)],
     fill(Float64(z), 2, 2))

_box_face_y(x1, x2, y, z1, z2) =
    ([xi for xi in (x1, x2), _ in (z1, z2)],
     fill(Float64(y), 2, 2),
     [zi for _ in (x1, x2), zi in (z1, z2)])

_box_face_x(x, y1, y2, z1, z2) =
    (fill(Float64(x), 2, 2),
     [yi for yi in (y1, y2), _ in (z1, z2)],
     [zi for _ in (y1, y2), zi in (z1, z2)])

# ── Surface helper ────────────────────────────────────────────────────────────

function _draw_surface!(ax, mesh_tuple, color)
    X, Y, Z = mesh_tuple
    surface!(ax, X, Y, Z; color=fill(color, size(X)...), shading=true, backlight=0.4f0)
end

# Draw a box layer; full=true for flesh (no cutaway), false for outer layers.
# Cutaway removes the front face (y = -hw) and the right face (x = +hl).
function _draw_box_faces!(ax, hl, hw, hh, color; full=false)
    _draw_surface!(ax, _box_face_z(-hl,  hl, -hw, hw, -hh), color)  # bottom
    _draw_surface!(ax, _box_face_z(-hl,  hl, -hw, hw,  hh), color)  # top
    _draw_surface!(ax, _box_face_y(-hl,  hl,  hw, -hh, hh), color)  # back
    _draw_surface!(ax, _box_face_x(-hl, -hw,  hw, -hh, hh), color)  # left
    if full
        _draw_surface!(ax, _box_face_y(-hl, hl, -hw, -hh, hh), color)  # front
        _draw_surface!(ax, _box_face_x( hl, -hw,  hw, -hh, hh), color)  # right
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 2-D COORDINATE HELPERS (internal)
# ══════════════════════════════════════════════════════════════════════════════

_pu(x) = Float32(ustrip(u"cm", x))

_circle_pts(r; n=300) =
    [Point2f(_pu(r)*cos(t), _pu(r)*sin(t)) for t in LinRange(0, 2π, n+1)]

_ellipse_pts(xr, yr; n=300) =
    [Point2f(_pu(xr)*cos(t), _pu(yr)*sin(t)) for t in LinRange(0, 2π, n+1)]

_rect_pts(hw, hh) =
    Point2f[(-_pu(hw), -_pu(hh)), (_pu(hw), -_pu(hh)),
            ( _pu(hw),  _pu(hh)), (-_pu(hw),  _pu(hh))]

_layer!(ax, pts, col) = poly!(ax, pts; color=col, strokecolor=col, strokewidth=1)

# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC 3-D API
# ══════════════════════════════════════════════════════════════════════════════

"""
    draw_cutaway!(ax::Axis3, body; sc=100.0, flesh_col=…, fat_col=…, fur_col=…)

Draw a quarter-cutaway 3-D surface mesh of `body` into an existing `Axis3`.
`sc` converts metres to axis units (default 100 → cm labels).
"""
function BiophysicalGeometry.draw_cutaway!(ax, body;
        sc        = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col   = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col   = RGBAf(0.76, 0.62, 0.42, 0.45))

    r_f = ustrip(u"m", flesh_radius(body))      * sc
    r_s = ustrip(u"m", skin_radius(body))       * sc
    r_i = ustrip(u"m", insulation_radius(body)) * sc
    has_fat = r_s > r_f + 1e-9
    has_fur = r_i > r_s + 1e-9

    if body.shape isa Cylinder
        gl     = body.geometry.length
        L_s    = ustrip(u"m", gl.length_skin) * sc
        L_i    = hasproperty(gl, :length_fur) ? ustrip(u"m", gl.length_fur) * sc : L_s
        z0_fur = -(L_i - L_s) / 2

        _draw_surface!(ax, _cylinder_tube(r_f, L_s),              flesh_col)
        _draw_surface!(ax, _cylinder_cap(r_f, 0.0),               flesh_col)
        _draw_surface!(ax, _cylinder_cap(r_f, L_s),               flesh_col)
        if has_fat
            _draw_surface!(ax, _cylinder_tube_partial(r_s, L_s),  fat_col)
            _draw_surface!(ax, _cylinder_cap_partial(r_s, 0.0),   fat_col)
            _draw_surface!(ax, _cylinder_cap_partial(r_s, L_s),   fat_col)
        end
        if has_fur
            _draw_surface!(ax, _cylinder_tube_partial(r_i, L_i; z0=z0_fur), fur_col)
            _draw_surface!(ax, _cylinder_cap_partial(r_i, z0_fur),           fur_col)
            _draw_surface!(ax, _cylinder_cap_partial(r_i, z0_fur + L_i),     fur_col)
        end

    elseif body.shape isa Sphere
        _draw_surface!(ax, _sphere_mesh(r_f), flesh_col)
        has_fat && _draw_surface!(ax, _sphere_mesh_partial(r_s), fat_col)
        has_fur && _draw_surface!(ax, _sphere_mesh_partial(r_i), fur_col)

    elseif body.shape isa Ellipsoid
        ratio = Float64(body.shape.b)
        _draw_surface!(ax, _ellipsoid_mesh(r_f * ratio, r_f), flesh_col)
        has_fat && _draw_surface!(ax, _ellipsoid_mesh_partial(r_s * ratio, r_s), fat_col)
        has_fur && _draw_surface!(ax, _ellipsoid_mesh_partial(r_i * ratio, r_i), fur_col)

    elseif body.shape isa Plate
        gl   = body.geometry.length
        hw_s = ustrip(u"m", gl.width_skin)  / 2 * sc
        hl_s = ustrip(u"m", gl.length_skin) / 2 * sc
        hh_s = ustrip(u"m", gl.height_skin) / 2 * sc
        hw_i = hasproperty(gl, :width_fur)  ? ustrip(u"m", gl.width_fur)  / 2 * sc : hw_s
        hl_i = hasproperty(gl, :length_fur) ? ustrip(u"m", gl.length_fur) / 2 * sc : hl_s
        hh_i = hasproperty(gl, :height_fur) ? ustrip(u"m", gl.height_fur) / 2 * sc : hh_s
        hw_f = r_f
        hl_f = hw_f * Float64(body.shape.b)
        hh_f = hl_f / Float64(body.shape.c)

        _draw_box_faces!(ax, hl_f, hw_f, hh_f, flesh_col; full=true)
        has_fat && _draw_box_faces!(ax, hl_s, hw_s, hh_s, fat_col)
        has_fur && _draw_box_faces!(ax, hl_i, hw_i, hh_i, fur_col)
    end

    ax.xlabel = "x (cm)"; ax.ylabel = "y (cm)"; ax.zlabel = "z (cm)"
end

"""
    plot_body(body; sc=100.0, flesh_col=…, fat_col=…, fur_col=…) -> Figure

Create a labelled `Figure` with `Axis3`, draw a quarter-cutaway of `body`,
and attach a legend. Returns the `Figure`.
"""
function BiophysicalGeometry.plot_body(body;
        sc        = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col   = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col   = RGBAf(0.76, 0.62, 0.42, 0.45))

    shape_name = string(nameof(typeof(body.shape)))
    ins_name   = string(nameof(typeof(body.insulation)))

    fig = Figure(size=(600, 560), backgroundcolor=:white)
    Label(fig[0, 1],
          "BiophysicalGeometry.jl — $(shape_name) · $(ins_name)  (cutaway)";
          fontsize=13, font=:bold, padding=(0, 0, 10, 0))
    ax = Axis3(fig[1, 1];
               perspectiveness=0.5, viewmode=:fit, aspect=:data,
               elevation=π/7, azimuth=5π/4)
    draw_cutaway!(ax, body; sc, flesh_col, fat_col, fur_col)
    Legend(fig[2, 1],
        [PolyElement(polycolor=flesh_col, strokecolor=:saddlebrown, strokewidth=1),
         PolyElement(polycolor=fat_col,   strokecolor=:saddlebrown, strokewidth=1),
         PolyElement(polycolor=fur_col,   strokecolor=:black,       strokewidth=1)],
        ["Flesh / muscle", "Subcutaneous fat", "Fur / insulation"];
        orientation=:horizontal, framevisible=false)
    return fig
end

# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC 2-D API
# ══════════════════════════════════════════════════════════════════════════════

"""
    draw_cross_sections!(ax_long, ax_tran, body; flesh_col=…, fat_col=…, fur_col=…)

Draw longitudinal and transverse cross-section polygons into a pair of `Axis`
objects.  Layers are drawn outer-to-inner so each inner layer paints over the
outer one.
"""
function BiophysicalGeometry.draw_cross_sections!(ax_long, ax_tran, body;
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_col   = RGBf(1.00, 0.97, 0.60),
        fur_col   = RGBf(0.76, 0.62, 0.42))

    r_f = flesh_radius(body)
    r_s = skin_radius(body)
    r_i = insulation_radius(body)
    pad = 0.12

    if body.shape isa Cylinder
        gl   = body.geometry.length
        hl_s = gl.length_skin / 2
        hl_i = hasproperty(gl, :length_fur) ? gl.length_fur / 2 : hl_s

        r_i > r_s && _layer!(ax_long, _rect_pts(r_i, hl_i), fur_col)
        r_s > r_f && _layer!(ax_long, _rect_pts(r_s, hl_s), fat_col)
                     _layer!(ax_long, _rect_pts(r_f, hl_s), flesh_col)

        r_i > r_s && _layer!(ax_tran, _circle_pts(r_i), fur_col)
        r_s > r_f && _layer!(ax_tran, _circle_pts(r_s), fat_col)
                     _layer!(ax_tran, _circle_pts(r_f),  flesh_col)

        xlims!(ax_long, -_pu(r_i)*(1+pad),  _pu(r_i)*(1+pad))
        ylims!(ax_long, -_pu(hl_i)*(1+pad), _pu(hl_i)*(1+pad))
        xlims!(ax_tran, -_pu(r_i)*(1+pad),  _pu(r_i)*(1+pad))
        ylims!(ax_tran, -_pu(hl_i)*(1+pad), _pu(hl_i)*(1+pad))

    elseif body.shape isa Plate
        gl   = body.geometry.length
        hl_s = gl.length_skin / 2
        hh_s = gl.height_skin / 2
        hw_s = gl.width_skin  / 2
        hl_i = hasproperty(gl, :length_fur) ? gl.length_fur / 2 : hl_s
        hh_i = hasproperty(gl, :height_fur) ? gl.height_fur / 2 : hh_s
        hw_i = hasproperty(gl, :width_fur)  ? gl.width_fur  / 2 : hw_s
        hw_f = r_f
        hl_f = hw_f * body.shape.b
        hh_f = hl_f / body.shape.c

        # Longitudinal: length × height
        r_i > r_s && _layer!(ax_long, _rect_pts(hl_i, hh_i), fur_col)
        r_s > r_f && _layer!(ax_long, _rect_pts(hl_s, hh_s), fat_col)
                     _layer!(ax_long, _rect_pts(hl_f, hh_f), flesh_col)

        # Transverse: width × height
        r_i > r_s && _layer!(ax_tran, _rect_pts(hw_i, hh_i), fur_col)
        r_s > r_f && _layer!(ax_tran, _rect_pts(hw_s, hh_s), fat_col)
                     _layer!(ax_tran, _rect_pts(hw_f, hh_f), flesh_col)

        xlims!(ax_long, -_pu(hl_i)*(1+pad), _pu(hl_i)*(1+pad))
        ylims!(ax_long, -_pu(hh_i)*(1+pad), _pu(hh_i)*(1+pad))
        xlims!(ax_tran, -_pu(hw_i)*(1+pad), _pu(hw_i)*(1+pad))
        ylims!(ax_tran, -_pu(hh_i)*(1+pad), _pu(hh_i)*(1+pad))

    else  # Sphere or Ellipsoid
        ratio = hasproperty(body.shape, :b) ? Float64(body.shape.b) : 1.0

        r_i > r_s && _layer!(ax_long, _ellipse_pts(r_i, r_i * ratio), fur_col)
        r_s > r_f && _layer!(ax_long, _ellipse_pts(r_s, r_s * ratio), fat_col)
                     _layer!(ax_long, _ellipse_pts(r_f, r_f * ratio), flesh_col)

        r_i > r_s && _layer!(ax_tran, _circle_pts(r_i), fur_col)
        r_s > r_f && _layer!(ax_tran, _circle_pts(r_s), fat_col)
                     _layer!(ax_tran, _circle_pts(r_f),  flesh_col)

        ey = r_i * ratio
        xlims!(ax_long, -_pu(r_i)*(1+pad), _pu(r_i)*(1+pad))
        ylims!(ax_long, -_pu(ey)*(1+pad),  _pu(ey)*(1+pad))
        xlims!(ax_tran, -_pu(r_i)*(1+pad), _pu(r_i)*(1+pad))
        ylims!(ax_tran, -_pu(ey)*(1+pad),  _pu(ey)*(1+pad))
    end
end

"""
    plot_cross_sections(body; flesh_col=…, fat_col=…, fur_col=…) -> Figure

Create a two-panel `Figure` (longitudinal + transverse cross-sections) with
legend. Returns the `Figure`.
"""
function BiophysicalGeometry.plot_cross_sections(body;
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_col   = RGBf(1.00, 0.97, 0.60),
        fur_col   = RGBf(0.76, 0.62, 0.42))

    shape_name = string(nameof(typeof(body.shape)))
    ins_name   = string(nameof(typeof(body.insulation)))

    fig = Figure(backgroundcolor=:white)
    Label(fig[0, 1:2],
          "$(shape_name) · $(ins_name)  (mass = $(body.shape.mass))";
          fontsize=13, font=:bold, padding=(0, 0, 8, 0))
    ax_long = Axis(fig[1, 1];
                   title="Longitudinal section",
                   xlabel="x (cm)", ylabel="y (cm)", aspect=DataAspect())
    ax_tran = Axis(fig[1, 2];
                   title="Transverse section",
                   xlabel="x (cm)", ylabel="y (cm)", aspect=DataAspect())
    draw_cross_sections!(ax_long, ax_tran, body; flesh_col, fat_col, fur_col)
    Legend(fig[2, 1:2],
        [PolyElement(polycolor=flesh_col, strokecolor=flesh_col, strokewidth=1),
         PolyElement(polycolor=fat_col,   strokecolor=fat_col,   strokewidth=1),
         PolyElement(polycolor=fur_col,   strokecolor=fur_col,   strokewidth=1)],
        ["Flesh / muscle", "Subcutaneous fat", "Fur / fibre insulation"];
        orientation=:horizontal, framevisible=false)
    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 30)
    return fig
end

end  # module BiophysicalGeometryMakieExt
