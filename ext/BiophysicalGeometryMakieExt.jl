module BiophysicalGeometryMakieExt

using Makie
using Unitful
using BiophysicalGeometry
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Plate, Cone, Half
import BiophysicalGeometry: Naked
import BiophysicalGeometry: CompositeBody, Pose, apply_pose, silhouette_rasterized
# Mesh helpers now live in core (src/meshes.jl); reuse them here.
import BiophysicalGeometry: _cylinder_tube, _cylinder_cap, _ellipsoid_mesh,
    _half_cylinder_tube, _half_cylinder_cap, _half_cylinder_flat,
    _half_ellipsoid_dome_mesh, _half_ellipsoid_flat_mesh,
    _cone_tube, _box_face_x, _box_face_y, _box_face_z,
    _part_outer_meshes, _transform_mesh, outer_dims

# ══════════════════════════════════════════════════════════════════════════════
# GENERIC HELPERS
# ══════════════════════════════════════════════════════════════════════════════

# Two ways units leave this module:
#   `_pu(x)`      → Float32 cm, for Point2f coordinates in 2-D recipes.
#   `_m(x, sc)`   → Float64 m·sc, for 3-D drawing that uses the caller's
#                   `sc` scale factor.
_pu(x) = Float32(ustrip(u"cm", x))
_m(x, sc) = ustrip(u"m", x) * sc

_radii(body) = (flesh=flesh_radius(body), skin=skin_radius(body), ins=insulation_radius(body))

_scaled_radii(body, sc) = map(x -> _m(x, sc), _radii(body))

_layer_flags(r) = (fat=r.skin > r.flesh + 1e-9, fur=r.ins > r.skin + 1e-9)

_axis_ratio(::Any) = 1.0
_axis_ratio(s::Ellipsoid) = Float64(s.b)

_colors(p) = (flesh=p[:flesh_col][], fat=p[:fat_col][], fur=p[:fur_col][])

_limits(lx, ly, tx, ty) = (
    long_x=(-lx, lx), long_y=(-ly, ly),
    tran_x=(-tx, tx), tran_y=(-ty, ty)
)

function _draw_layers!(target, layers)
    for (cond, geom, col) in layers
        cond && _layer!(target, geom(), col)
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 3-D DRAW PRIMITIVES
# ══════════════════════════════════════════════════════════════════════════════

function _draw_surface!(target, (X, Y, Z), color)
    surface!(target, X, Y, Z; color=fill(color, size(X)...), shading=true, backlight=0.4f0)
end

function _draw_cylinder!(target, r, L, col; θ_end=2π, z0=0.0)
    _draw_surface!(target, _cylinder_tube(r, L; θ_end, z0), col)
    _draw_surface!(target, _cylinder_cap(r, z0; θ_end), col)
    _draw_surface!(target, _cylinder_cap(r, z0 + L; θ_end), col)
end

# Cutaway removes the front face (y = -hw) and the right face (x = +hl).
function _draw_box_faces!(target, hl, hw, hh, color; full=false)
    faces = [
        _box_face_z(-hl, hl, -hw, hw, -hh),
        _box_face_z(-hl, hl, -hw, hw, hh),
        _box_face_y(-hl, hl, hw, -hh, hh),
        _box_face_x(-hl, -hw, hw, -hh, hh),
    ]
    full && append!(faces, [
        _box_face_y(-hl, hl, -hw, -hh, hh),
        _box_face_x( hl, -hw, hw, -hh, hh),
    ])
    for f in faces; _draw_surface!(target, f, color); end
end

# ══════════════════════════════════════════════════════════════════════════════
# 3-D SHAPE DISPATCH
# ══════════════════════════════════════════════════════════════════════════════

function _draw_cutaway_shape!(p, sh::Cylinder, body, sc, cols)
    r = _scaled_radii(body, sc)
    fl = _layer_flags(r)

    L_s = _m(body.geometry.length.length_skin, sc)
    L_i = _m(outer_dims(sh, body).L, sc)
    z0_fibrous = -(L_i - L_s) / 2

    _draw_cylinder!(p, r.flesh, L_s, cols.flesh)
    fl.fat && _draw_cylinder!(p, r.skin, L_s, cols.fat; θ_end=3π/2)
    fl.fur && _draw_cylinder!(p, r.ins, L_i, cols.fur; θ_end=3π/2, z0=z0_fibrous)
end

function _draw_cutaway_shape!(p, shape::Union{Sphere,Ellipsoid}, body, sc, cols)
    r = _scaled_radii(body, sc)
    fl = _layer_flags(r)
    ratio = _axis_ratio(shape)

    _draw_surface!(p, _ellipsoid_mesh(r.flesh * ratio, r.flesh), cols.flesh)
    fl.fat && _draw_surface!(p, _ellipsoid_mesh(r.skin * ratio, r.skin; θ_end=3π/2), cols.fat)
    fl.fur && _draw_surface!(p, _ellipsoid_mesh(r.ins * ratio, r.ins;  θ_end=3π/2), cols.fur)
end

function _draw_cutaway_shape!(p, sh::Plate, body, sc, cols)
    r = _scaled_radii(body, sc)
    fl = _layer_flags(r)

    gl = body.geometry.length
    di = outer_dims(sh, body)
    hw_s = _m(gl.width_skin, sc) / 2
    hl_s = _m(gl.length_skin, sc) / 2
    hh_s = _m(gl.height_skin, sc) / 2
    hw_i = _m(di.W, sc) / 2
    hl_i = _m(di.L, sc) / 2
    hh_i = _m(di.H, sc) / 2
    hw_f = r.flesh
    hl_f = hw_f * Float64(sh.b)
    hh_f = hl_f / Float64(sh.c)

    _draw_box_faces!(p, hl_f, hw_f, hh_f, cols.flesh; full=true)
    fl.fat && _draw_box_faces!(p, hl_s, hw_s, hh_s, cols.fat)
    fl.fur && _draw_box_faces!(p, hl_i, hw_i, hh_i, cols.fur)
end

# ══════════════════════════════════════════════════════════════════════════════
# COMPOSITE BODY DRAWING
# ══════════════════════════════════════════════════════════════════════════════

# Pick a fill colour per part — outer insulation if any, else flesh.
_part_color(body, cols) = body.insulation isa Naked ? cols.flesh : cols.fur

function _draw_composite!(p, b::CompositeBody, sc, cols)
    for name in propertynames(b.parts)
        part = getfield(b.parts, name)
        pose = getfield(b.poses, name)
        col = _part_color(part, cols)
        for mesh in _part_outer_meshes(part.shape, part, sc)
            _draw_surface!(p, _transform_mesh(mesh..., pose, sc), col)
        end
    end
end

# Posed world-frame bounding box of a composite (in `sc` units): (mins, maxs).
function _composite_bbox(b::CompositeBody, sc)
    xmin, xmax = Inf, -Inf
    ymin, ymax = Inf, -Inf
    zmin, zmax = Inf, -Inf
    for name in propertynames(b.parts)
        part = getfield(b.parts, name)
        pose = getfield(b.poses, name)
        for grid in _part_outer_meshes(part.shape, part, sc)
            X, Y, Z = _transform_mesh(grid..., pose, sc)
            xmin = min(xmin, minimum(X)); xmax = max(xmax, maximum(X))
            ymin = min(ymin, minimum(Y)); ymax = max(ymax, maximum(Y))
            zmin = min(zmin, minimum(Z)); zmax = max(zmax, maximum(Z))
        end
    end
    ((xmin, ymin, zmin), (xmax, ymax, zmax))
end

_composite_extent(b::CompositeBody, sc) =
    let (mins, maxs) = _composite_bbox(b, sc)
        max(maxs[1]-mins[1], maxs[2]-mins[2], maxs[3]-mins[3])
    end

# Posed bounding box for a single part (in `sc` units).
function _part_bbox(shape, body, pose::Pose, sc)
    xmin, xmax = Inf, -Inf
    ymin, ymax = Inf, -Inf
    zmin, zmax = Inf, -Inf
    for grid in _part_outer_meshes(shape, body, sc)
        X, Y, Z = _transform_mesh(grid..., pose, sc)
        xmin = min(xmin, minimum(X)); xmax = max(xmax, maximum(X))
        ymin = min(ymin, minimum(Y)); ymax = max(ymax, maximum(Y))
        zmin = min(zmin, minimum(Z)); zmax = max(zmax, maximum(Z))
    end
    ((xmin, ymin, zmin), (xmax, ymax, zmax))
end

# ══════════════════════════════════════════════════════════════════════════════
# 2-D PRIMITIVES
# ══════════════════════════════════════════════════════════════════════════════

_circle_pts(r; n=300) =
    [Point2f(_pu(r)*cos(t), _pu(r)*sin(t)) for t in LinRange(0, 2π, n+1)]

_ellipse_pts(xr, yr; n=300) =
    [Point2f(_pu(xr)*cos(t), _pu(yr)*sin(t)) for t in LinRange(0, 2π, n+1)]

_rect_pts(hw, hh) =
    Point2f[(-_pu(hw), -_pu(hh)), (_pu(hw), -_pu(hh)),
            ( _pu(hw), _pu(hh)), (-_pu(hw), _pu(hh))]

_layer!(target, pts, col) = poly!(target, pts; color=col, strokecolor=col, strokewidth=1)

# ══════════════════════════════════════════════════════════════════════════════
# 2-D SECTION LAYER DISPATCH
# ══════════════════════════════════════════════════════════════════════════════

function _section_layers(sh::Cylinder, body, mode, r, cols)
    r_f, r_s, r_i = r.flesh, r.skin, r.ins
    hl_s = body.geometry.length.length_skin / 2
    hl_i = outer_dims(sh, body).L / 2
    if mode === :long
        [
            (r_i > r_s, () -> _rect_pts(r_i, hl_i), cols.fur),
            (r_s > r_f, () -> _rect_pts(r_s, hl_s), cols.fat),
            (true, () -> _rect_pts(r_f, hl_s), cols.flesh),
        ]
    else
        [
            (r_i > r_s, () -> _circle_pts(r_i), cols.fur),
            (r_s > r_f, () -> _circle_pts(r_s), cols.fat),
            (true, () -> _circle_pts(r_f), cols.flesh),
        ]
    end
end

function _section_layers(sh::Plate, body, mode, r, cols)
    gl = body.geometry.length
    di = outer_dims(sh, body)
    r_f, r_s, r_i = r.flesh, r.skin, r.ins
    if mode === :long
        d_s = (gl.length_skin / 2, gl.height_skin / 2)
        d_i = (di.L / 2, di.H / 2)
        d_f = (r_f * sh.b, r_f * sh.b / sh.c)
    else
        d_s = (gl.width_skin / 2, gl.height_skin / 2)
        d_i = (di.W / 2, di.H / 2)
        d_f = (r_f, (r_f * sh.b) / sh.c)
    end
    [
        (r_i > r_s, () -> _rect_pts(d_i...), cols.fur),
        (r_s > r_f, () -> _rect_pts(d_s...), cols.fat),
        (true, () -> _rect_pts(d_f...), cols.flesh),
    ]
end

function _section_layers(shape::Union{Sphere,Ellipsoid}, _, mode, r, cols)
    r_f, r_s, r_i = r.flesh, r.skin, r.ins
    ratio = _axis_ratio(shape)
    geom = mode === :long ? x -> _ellipse_pts(x, x * ratio) : _circle_pts
    [
        (r_i > r_s, () -> geom(r_i), cols.fur),
        (r_s > r_f, () -> geom(r_s), cols.fat),
        (true, () -> geom(r_f), cols.flesh),
    ]
end

# ══════════════════════════════════════════════════════════════════════════════
# 2-D AXIS LIMIT DISPATCH
# ══════════════════════════════════════════════════════════════════════════════

function _section_limits(sh::Cylinder, body, r, pad)
    hl_i = outer_dims(sh, body).L / 2
    ri = _pu(r.ins) * (1 + pad)
    li = _pu(hl_i) * (1 + pad)
    _limits(ri, li, ri, li)
end

function _section_limits(sh::Plate, body, r, pad)
    di = outer_dims(sh, body)
    hl_i = di.L / 2
    hh_i = di.H / 2
    hw_i = di.W / 2
    _limits(_pu(hl_i) * (1 + pad), _pu(hh_i) * (1 + pad),
            _pu(hw_i) * (1 + pad), _pu(hh_i) * (1 + pad))
end

function _section_limits(shape::Union{Sphere,Ellipsoid}, body, r, pad)
    ratio = _axis_ratio(shape)
    ey = r.ins * ratio
    rx = _pu(r.ins) * (1 + pad)
    ry = _pu(ey) * (1 + pad)
    _limits(rx, ry, rx, ry)
end

# ══════════════════════════════════════════════════════════════════════════════
# RECIPES
# ══════════════════════════════════════════════════════════════════════════════

@recipe(BodyCutaway, body) do scene
    Theme(
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col = RGBAf(0.76, 0.62, 0.42, 0.45),
        sc = 100.0,
    )
end

function Makie.plot!(p::BodyCutaway)
    body = p[:body][];  sc = p[:sc][]
    if body isa CompositeBody
        _draw_composite!(p, body, sc, _colors(p))
    else
        _draw_cutaway_shape!(p, body.shape, body, sc, _colors(p))
    end
    p
end

@recipe(BodyLongSection, body) do scene
    Theme(
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_col = RGBf(1.00, 0.97, 0.60),
        fur_col = RGBf(0.76, 0.62, 0.42),
    )
end

function Makie.plot!(p::BodyLongSection)
    body = p[:body][];  r = _radii(body)
    _draw_layers!(p, _section_layers(body.shape, body, :long, r, _colors(p)))
    p
end

@recipe(BodyTransSection, body) do scene
    Theme(
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_col = RGBf(1.00, 0.97, 0.60),
        fur_col = RGBf(0.76, 0.62, 0.42),
    )
end

function Makie.plot!(p::BodyTransSection)
    body = p[:body][];  r = _radii(body)
    _draw_layers!(p, _section_layers(body.shape, body, :trans, r, _colors(p)))
    p
end

# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC API — extends the stubs declared in BiophysicalGeometry
# ══════════════════════════════════════════════════════════════════════════════

"""
    draw_cutaway!(ax::Axis3, body; sc=100.0, flesh_col=…, fat_col=…, fur_col=…)

Draw a quarter-cutaway 3-D surface mesh of `body` into an existing `Axis3`.
`sc` converts metres to axis units (default 100 → cm labels).
"""
function BiophysicalGeometry.draw_cutaway!(ax, body;
        sc = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col = RGBAf(0.76, 0.62, 0.42, 0.45))
    bodycutaway!(ax, body; sc, flesh_col, fat_col, fur_col)
    ax.xlabel = "x (cm)"; ax.ylabel = "y (cm)"; ax.zlabel = "z (cm)"
end

"""
    plot_body(body; sc=100.0, flesh_col=…, fat_col=…, fur_col=…) -> Figure

Create a labelled `Figure` with `Axis3`, draw a quarter-cutaway of `body`,
and attach a legend. Returns the `Figure`.
"""
function BiophysicalGeometry.plot_body(body;
        sc = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col = RGBAf(0.76, 0.62, 0.42, 0.45))

    if body isa CompositeBody
        shape_name = "Composite($(length(body.parts)) parts)"
        ins_name = ""
    else
        shape_name = string(nameof(typeof(body.shape)))
        ins_name = string(nameof(typeof(body.insulation)))
    end

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
         PolyElement(polycolor=fat_col, strokecolor=:saddlebrown, strokewidth=1),
         PolyElement(polycolor=fur_col, strokecolor=:black, strokewidth=1)],
        ["Flesh / muscle", "Subcutaneous fat", "FibrousLayer / insulation"];
        orientation=:horizontal, framevisible=false)
    return fig
end

"""
    plot_body_silhouette(body::CompositeBody; resolution=128, sc=100.0, kwargs...) -> Figure

Interactive Makie figure: 3-D cutaway of `body` on the left, rasterised
silhouette projection on the right, with `zenith θ` and `azimuth φ`
sliders below. The silhouette area updates live as you drag the sliders.

The `resolution` kwarg controls the silhouette bitmap (≈ accuracy /
performance trade-off). `sc` scales metres to plot units (default 100 → cm).
"""
function BiophysicalGeometry.plot_body_silhouette(body::CompositeBody;
        resolution::Integer = 128,
        sc = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col = RGBAf(0.76, 0.62, 0.42, 0.45))

    fig = Figure(size=(1100, 640), backgroundcolor=:white)
    Label(fig[0, 1:2],
          "BiophysicalGeometry.jl — silhouette projection (interactive)";
          fontsize=13, font=:bold, padding=(0, 0, 8, 0))

    ax3 = Axis3(fig[1, 1];
                perspectiveness=0.5, viewmode=:fit, aspect=:data,
                elevation=π/7, azimuth=5π/4,
                xlabel="x (cm)", ylabel="y (cm)", zlabel="z (cm)",
                title="Body + sun direction")
    ax2 = Axis(fig[1, 2]; aspect=DataAspect(), title="Silhouette projection",
               xlabel="u (cm)", ylabel="v (cm)")

    cols = (flesh=flesh_col, fat=fat_col, fur=fur_col)
    _draw_composite!(ax3, body, sc, cols)

    sg = SliderGrid(fig[2, 1:2],
        (label="zenith θ (0=overhead, π/2=horizon)", range=range(0.0, π/2, length=91),
         startvalue=π/4, format="{:.2f} rad"),
        (label="azimuth φ", range=range(0.0, 2π, length=361),
         startvalue=0.0, format="{:.2f} rad"))
    θobs = sg.sliders[1].value
    φobs = sg.sliders[2].value

    sun_dir = lift(θobs, φobs) do θ, φ
        (sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ))
    end

    # Silhouette image + area (single rasteriser call per slider event).
    result = lift(sun_dir) do d
        silhouette_rasterized(body, d; resolution=resolution, return_image=true)
    end

    # Display bitmap as a heatmap in cm-coords on ax2.
    img_x = lift(result) do r
        LinRange(r.x_range[1] * sc, r.x_range[2] * sc, size(r.bitmap, 1))
    end
    img_y = lift(result) do r
        LinRange(r.y_range[1] * sc, r.y_range[2] * sc, size(r.bitmap, 2))
    end
    img = lift(r -> Float32.(r.bitmap), result)
    heatmap!(ax2, img_x, img_y, img; colormap=[:white, RGBAf(0.2, 0.2, 0.25, 1.0)],
             colorrange=(0.0f0, 1.0f0))
    # Refit ax2 limits to the projected bbox each time the sun direction
    # changes (Makie won't auto-shrink the limits otherwise).
    on(result) do r
        xlims!(ax2, r.x_range[1] * sc, r.x_range[2] * sc)
        ylims!(ax2, r.y_range[1] * sc, r.y_range[2] * sc)
    end
    notify(result)

    # Sun-direction indicator: line from a sun marker to the root body's
    # centre (so the target stays at the main body, not pulled around by
    # legs/head when the bbox centre moves).
    root_part = body.parts[body.root]
    root_pose = body.poses[body.root]
    (rb_min, rb_max) = _part_bbox(root_part.shape, root_part, root_pose, sc)
    body_centre = Point3f((rb_min[1] + rb_max[1]) / 2,
                          (rb_min[2] + rb_max[2]) / 2,
                          (rb_min[3] + rb_max[3]) / 2)
    # Whole-composite bbox just used to pick a sensible arrow length and
    # axis limits.
    (bb_min, bb_max) = _composite_bbox(body, sc)
    extent = max(bb_max[1]-bb_min[1], bb_max[2]-bb_min[2], bb_max[3]-bb_min[3])
    arrow_len = 0.7 * extent

    sun_marker_pos = lift(d -> Point3f(body_centre[1] + arrow_len * d[1],
                                        body_centre[2] + arrow_len * d[2],
                                        body_centre[3] + arrow_len * d[3]), sun_dir)
    sun_line = lift(p -> [p, body_centre], sun_marker_pos)
    lines!(ax3, sun_line; color=:goldenrod, linewidth=3)
    scatter!(ax3, lift(p -> [p], sun_marker_pos);
             color=:goldenrod, markersize=20, marker=:circle,
             strokecolor=:black, strokewidth=1)

    # Expand axis limits so the sun marker is always visible.
    pad = arrow_len * 1.05
    xlims!(ax3, body_centre[1] - pad, body_centre[1] + pad)
    ylims!(ax3, body_centre[2] - pad, body_centre[2] + pad)
    zlims!(ax3, body_centre[3] - pad, body_centre[3] + pad)

    Label(fig[3, 1:2],
          lift(r -> "silhouette area = " * string(round(ustrip(u"cm^2", r.area), digits=1)) * " cm²", result);
          fontsize=14, font=:bold)

    return fig
end

"""
    draw_cross_sections!(ax_long, ax_tran, body; flesh_col=…, fat_col=…, fur_col=…)

Draw longitudinal and transverse cross-section polygons into a pair of `Axis`
objects.  Layers are drawn outer-to-inner so each inner layer paints over the
outer one.
"""
function BiophysicalGeometry.draw_cross_sections!(ax_long, ax_tran, body;
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_col = RGBf(1.00, 0.97, 0.60),
        fur_col = RGBf(0.76, 0.62, 0.42))

    bodylongsection!(ax_long, body; flesh_col, fat_col, fur_col)
    bodytranssection!(ax_tran, body; flesh_col, fat_col, fur_col)

    r = _radii(body)
    lims = _section_limits(body.shape, body, r, 0.12)
    xlims!(ax_long, lims.long_x...); ylims!(ax_long, lims.long_y...)
    xlims!(ax_tran, lims.tran_x...); ylims!(ax_tran, lims.tran_y...)
end

"""
    plot_cross_sections(body; flesh_col=…, fat_col=…, fur_col=…) -> Figure

Create a two-panel `Figure` (longitudinal + transverse cross-sections) with
legend. Returns the `Figure`.
"""
function BiophysicalGeometry.plot_cross_sections(body;
        flesh_col = RGBf(0.88, 0.48, 0.42),
        fat_col = RGBf(1.00, 0.97, 0.60),
        fur_col = RGBf(0.76, 0.62, 0.42))

    shape_name = string(nameof(typeof(body.shape)))
    ins_name = string(nameof(typeof(body.insulation)))

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
         PolyElement(polycolor=fat_col, strokecolor=fat_col, strokewidth=1),
         PolyElement(polycolor=fur_col, strokecolor=fur_col, strokewidth=1)],
        ["Flesh / muscle", "Subcutaneous fat", "FibrousLayer / fibre insulation"];
        orientation=:horizontal, framevisible=false)
    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 30)
    return fig
end

"""
    draw_insulation_schematic!(ax, fur::FibrousLayer; fibre_length=fur.thickness)

Draw a side-view schematic of a `FibrousLayer` insulation layer into `ax`.  Fibre width
is exaggerated for clarity.  When `fibre_length > fur.thickness` the fibres are
drawn as tilted parallelograms.
"""
function BiophysicalGeometry.draw_insulation_schematic!(ax, fur::FibrousLayer;
        fibre_length = fur.thickness)

    thick_mm = ustrip(u"mm", fur.thickness)
    fibre_len_mm = ustrip(u"mm", fibre_length)
    d_μm = ustrip(u"μm", fur.fibre_diameter)
    n_cm2 = ustrip(u"cm^-2", fur.fibre_density)

    spacing_mm = 1.0 / sqrt(n_cm2 / 100.0)
    d_display = spacing_mm * 0.40
    dx = sqrt(max(0.0, fibre_len_mm^2 - thick_mm^2))
    n_show = 8
    skin_h = thick_mm * 0.12
    W = n_show * spacing_mm

    poly!(ax,
          [Point2f(0, 0), Point2f(W, 0),
           Point2f(W, -skin_h), Point2f(0, -skin_h)],
          color=RGBf(0.88, 0.48, 0.42),
          strokecolor=:black, strokewidth=0.5)
    text!(ax, W / 2, -skin_h / 2;
          text="Skin surface", fontsize=9, align=(:center, :center))

    for i in 0:(n_show - 1)
        xc = (i + 0.5) * spacing_mm
        poly!(ax,
              [Point2f(xc - d_display/2, 0.0),
               Point2f(xc + d_display/2, 0.0),
               Point2f(xc + d_display/2 + dx, thick_mm),
               Point2f(xc - d_display/2 + dx, thick_mm)],
              color=RGBf(0.76, 0.62, 0.42),
              strokecolor=(:saddlebrown, 0.6), strokewidth=0.4)
    end

    lines!(ax, [0.0, W, W + dx, dx, 0.0],
               [0.0, 0.0, thick_mm, thick_mm, 0.0];
           color=(:grey40, 0.5), linewidth=0.8, linestyle=:dash)

    x_ann = W + dx + 0.15 * W
    lines!(ax, [x_ann, x_ann], [0.0, thick_mm]; color=:black, linewidth=1)
    scatter!(ax, [x_ann, x_ann], [0.0, thick_mm];
             color=:black, markersize=6, marker=:rect)
    text!(ax, x_ann + 0.04*W, thick_mm / 2;
          text="fur depth\n$(round(Int, thick_mm)) mm",
          fontsize=9, align=(:left, :center))

    if dx > 1e-6
        i_ann = n_show ÷ 2
        xc_ann = (i_ann + 0.5) * spacing_mm
        xfl0 = xc_ann + d_display / 2
        xfl1 = xfl0 + dx
        tilt_deg = round(Int, atand(dx, thick_mm))
        lines!(ax, [xfl0, xfl1], [0.0, thick_mm]; color=:purple, linewidth=1.4)
        scatter!(ax, [xfl0, xfl1], [0.0, thick_mm];
                 color=:purple, markersize=5, marker=:vline)
        text!(ax, (xfl0 + xfl1)/2, thick_mm * 1.02;
              text="L = $(round(fibre_len_mm, digits=1)) mm  ($(tilt_deg)° from vertical)",
              fontsize=8, align=(:center, :bottom), color=:purple)
    end

    y_d = -skin_h * 2.4
    x0 = 0.5 * spacing_mm - d_display / 2
    x1 = 0.5 * spacing_mm + d_display / 2
    lines!(ax, [x0, x1], [y_d, y_d]; color=:darkblue, linewidth=1.2)
    scatter!(ax, [x0, x1], [y_d, y_d]; color=:darkblue, markersize=5, marker=:vline)
    text!(ax, (x0 + x1)/2, y_d - 0.3;
          text="d = $(round(Int, d_μm)) μm\n(displayed $(round(Int, d_display*1e3)) μm for clarity)",
          fontsize=8, align=(:center, :top), color=:darkblue)

    x_sp0 = dx + 0.5 * spacing_mm
    x_sp1 = dx + 1.5 * spacing_mm
    y_sp = thick_mm * 1.22
    lines!(ax, [x_sp0, x_sp1], [y_sp, y_sp]; color=:darkgreen, linewidth=1.2)
    scatter!(ax, [x_sp0, x_sp1], [y_sp, y_sp]; color=:darkgreen, markersize=5, marker=:vline)
    text!(ax, (x_sp0 + x_sp1)/2, y_sp + 0.2;
          text="spacing ≈ $(round(spacing_mm, digits=2)) mm  (N = $n_cm2 cm⁻²)",
          fontsize=8, align=(:center, :bottom), color=:darkgreen)

    ax.title = "FibrousLayer schematic  (fibre width exaggerated for clarity)"
    ax.xlabel = "position (mm)"
    ax.ylabel = "height above skin (mm)"
    ylims!(ax, -skin_h * 5, thick_mm * 1.55)
    xlims!(ax, -0.1 * W, x_ann + 0.45 * W)
    hidespines!(ax, :t, :r)
end

"""
    draw_insulation_coverage!(ax, fur::FibrousLayer; d_range=LinRange(10,120,200), N_range=LinRange(200,9000,200))

Draw a coverage-fraction heatmap (plasma colormap) with contour lines at f = 0.25,
0.50, 0.75, 1.0, and mark the reference point for `fur`.  Returns the `Heatmap`
object so the caller can attach a `Colorbar`.
"""
function BiophysicalGeometry.draw_insulation_coverage!(ax, fur::FibrousLayer;
        d_range = LinRange(10.0, 120.0, 200),
        N_range = LinRange(200.0, 9000.0, 200))

    cov = [π * ((d_μm * 1e-6) / 2)^2 * (N_cm2 * 1e4)
           for d_μm in d_range, N_cm2 in N_range]

    hm = heatmap!(ax, d_range, N_range, cov;
                  colormap=:plasma, colorrange=(0.0, 1.0), highclip=:white)

    for (level, clr) in [(0.25, :white), (0.50, :white),
                          (0.75, :white), (1.00, :cyan)]
        contour!(ax, d_range, N_range, cov;
                 levels=[level], color=clr, linewidth=1.2, linestyle=:dash)
        N_label_cm2 = level / (π * ((d_range[end] * 1e-6) / 2)^2 * 1e4) * 1e-4
        if 200 < N_label_cm2 < 9000
            text!(ax, d_range[end] - 3, N_label_cm2;
                  text="f=$(level)", fontsize=8, align=(:right, :bottom), color=clr)
        end
    end

    d_ref = ustrip(u"μm", fur.fibre_diameter)
    N_ref = ustrip(u"cm^-2", fur.fibre_density)
    cov_ref = π * ((d_ref * 1e-6) / 2)^2 * (N_ref * 1e4)
    scatter!(ax, [d_ref], [N_ref]; color=:lime, markersize=11,
             strokecolor=:black, strokewidth=1)
    text!(ax, d_ref + 2, N_ref + 150;
          text="$(d_ref) μm, $(N_ref) cm⁻²\nf ≈ $(round(cov_ref, digits=3))",
          fontsize=8, color=:lime)

    ax.title = "Fibre coverage fraction  f = π(d/2)² × N"
    ax.xlabel = "Fibre diameter d (μm)"
    ax.ylabel = "Fibre density N (cm⁻²)"
    return hm
end

"""
    plot_insulation_properties(fur::FibrousLayer; fibre_length=fur.thickness, kwargs...) → Figure

Create a two-panel `Figure` showing a fur schematic and coverage heatmap for
the supplied `FibrousLayer` object.  `fibre_length` may exceed `fur.thickness` to show
tilted fibres.  `d_range` and `N_range` control the heatmap axes (μm / cm⁻²).
"""
function BiophysicalGeometry.plot_insulation_properties(fur::FibrousLayer;
        fibre_length = fur.thickness,
        d_range = LinRange(10.0, 120.0, 200),
        N_range = LinRange(200.0, 9000.0, 200))

    fig = Figure(size=(900, 420), backgroundcolor=:white)
    Label(fig[0, 1:3],
          "BiophysicalGeometry.jl — FibrousLayer Properties";
          fontsize=14, font=:bold, padding=(0, 0, 8, 0))
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    draw_insulation_schematic!(ax1, fur; fibre_length)
    hm = draw_insulation_coverage!(ax2, fur; d_range, N_range)
    Colorbar(fig[1, 3], hm; label="Coverage fraction f", width=14, labelsize=10)
    colgap!(fig.layout, 12)
    colsize!(fig.layout, 3, Auto(0.05))
    return fig
end

end  # module BiophysicalGeometryMakieExt
