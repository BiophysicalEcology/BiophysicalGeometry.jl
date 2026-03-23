module BiophysicalGeometryMakieExt

using Makie
using Unitful
using BiophysicalGeometry
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Plate, Organism, BodyPart, Join,
                             OrganismPosture, Upright, Prone, PartPosition, part_positions,
                             insulation_radius, flesh_radius, skin_radius, _half_axis_m

# ══════════════════════════════════════════════════════════════════════════════
# 3-D MESH GENERATORS (internal)
# ══════════════════════════════════════════════════════════════════════════════
#
# Each generator accepts θ_end (default 2π = full surface; 3π/2 = cutaway).
# Cylinder mesh is centred at z=0 with caps at ±L/2.
# Ellipsoid long axis is along Z; φ_max<π clips the south pole (neck join).

function _cylinder_tube(r, L; nθ=72, nz=2, θ_end=2π)
    L_f = Float64(L)
    θ = LinRange(0.0, θ_end, nθ);  z = LinRange(-L_f/2, L_f/2, nz)
    [r*cos(θi) for θi in θ, _ in z],
    [r*sin(θi) for θi in θ, _ in z],
    [zi        for _  in θ, zi in z]
end

function _cylinder_cap(r, z0; nθ=72, nr=12, θ_end=2π)
    θ = LinRange(0.0, θ_end, nθ);  rv = LinRange(0.0, r, nr)
    [rr*cos(θi) for θi in θ, rr in rv],
    [rr*sin(θi) for θi in θ, rr in rv],
    fill(Float64(z0), nθ, nr)
end

function _sphere_mesh(r; n=60, θ_end=2π)
    θ = LinRange(0.0, θ_end, n);  φ = LinRange(0.0, π, n)
    [r*sin(φj)*cos(θi) for θi in θ, φj in φ],
    [r*sin(φj)*sin(θi) for θi in θ, φj in φ],
    [r*cos(φj)         for _  in θ, φj in φ]
end

function _ellipsoid_mesh(a, b; n=60, θ_end=2π, φ_max=π)
    θ = LinRange(0.0, θ_end, n);  φ = LinRange(0.0, Float64(φ_max), n)
    [b*sin(φj)*cos(θi) for θi in θ, φj in φ],   # X: semi-minor
    [b*sin(φj)*sin(θi) for θi in θ, φj in φ],   # Y: semi-minor
    [a*cos(φj)         for _  in θ, φj in φ]     # Z: semi-major (long axis)
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

# Draw box faces via ds!(mesh, color).  Cutaway omits front (y=−hw) and right (x=+hl).
function _draw_box_faces!(hl, hw, hh, color, ds!; full=false)
    ds!(_box_face_z(-hl,  hl, -hw, hw, -hh), color)  # bottom
    ds!(_box_face_z(-hl,  hl, -hw, hw,  hh), color)  # top
    ds!(_box_face_y(-hl,  hl,  hw, -hh, hh), color)  # back
    ds!(_box_face_x(-hl, -hw,  hw, -hh, hh), color)  # left
    if full
        ds!(_box_face_y(-hl, hl, -hw, -hh, hh), color)  # front
        ds!(_box_face_x( hl, -hw,  hw, -hh, hh), color)  # right
    end
end

# ── 3-D transform helpers ──────────────────────────────────────────────────────

# Rodrigues rotation: rotate the default axis (0,0,1) to `target`.
function _rotation_matrix(target)
    tx, ty, tz = Float64(target[1]), Float64(target[2]), Float64(target[3])
    tnorm = sqrt(tx^2 + ty^2 + tz^2)
    tx /= tnorm; ty /= tnorm; tz /= tnorm

    kx, ky = -ty, tx            # cross product (0,0,1)×target; kz=0
    sinθ = sqrt(kx^2 + ky^2)
    cosθ = tz

    if sinθ < 1e-10
        return cosθ > 0 ?
            [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] :
            [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
    end
    kx /= sinθ; ky /= sinθ

    c, s, ic = cosθ, sinθ, 1.0 - cosθ
    [c + kx*kx*ic    kx*ky*ic        ky*s;
     ky*kx*ic        c + ky*ky*ic   -kx*s;
    -ky*s            kx*s            c   ]
end

function _apply_transform(mesh, R, offset)
    X, Y, Z = mesh
    Xr = R[1,1].*X .+ R[1,2].*Y .+ R[1,3].*Z .+ offset[1]
    Yr = R[2,1].*X .+ R[2,2].*Y .+ R[2,3].*Z .+ offset[2]
    Zr = R[3,1].*X .+ R[3,2].*Y .+ R[3,3].*Z .+ offset[3]
    (Xr, Yr, Zr)
end

# ══════════════════════════════════════════════════════════════════════════════
# RECIPE: BodyCutaway
# Draws a quarter-cutaway 3-D surface of one body part into an Axis3.
# ══════════════════════════════════════════════════════════════════════════════

@recipe(BodyCutaway, body) do scene
    Attributes(
        sc        = 100.0,
        offset    = Vec3f(0, 0, 0),
        axis_dir  = Vec3f(0, 0, 1),   # direction of the part's long axis
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col   = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col   = RGBAf(0.76, 0.62, 0.42, 0.45),
    )
end

Makie.preferred_axis_type(::BodyCutaway) = Axis3

function Makie.plot!(p::BodyCutaway)
    body      = p[:body][]
    sc        = p.sc[]
    off       = Float64.(p.offset[])
    axis_dir  = p.axis_dir[]
    flesh_col = p.flesh_col[]
    fat_col   = p.fat_col[]
    fur_col   = p.fur_col[]

    r_f = ustrip(u"m", flesh_radius(body))      * sc
    r_s = ustrip(u"m", skin_radius(body))       * sc
    r_i = ustrip(u"m", insulation_radius(body)) * sc
    has_fat = r_s > r_f + 1e-9
    has_fur = r_i > r_s + 1e-9

    R = _rotation_matrix(axis_dir)
    ds!(mesh, color) = let (X, Y, Z) = _apply_transform(mesh, R, off)
        surface!(p, X, Y, Z; color=fill(color, size(X)...), shading=true, backlight=0.4f0)
    end

    PARTIAL = 3π/2  # cutaway angle (3/4 of the circumference)

    if body.shape isa Cylinder
        gl  = body.geometry.length
        L_s = ustrip(u"m", gl.length_skin) * sc
        L_i = hasproperty(gl, :length_fur) ? ustrip(u"m", gl.length_fur) * sc : L_s

        ds!(_cylinder_tube(r_f, L_s),                               flesh_col)
        ds!(_cylinder_cap(r_f, -L_s/2),                             flesh_col)
        ds!(_cylinder_cap(r_f,  L_s/2),                             flesh_col)
        if has_fat
            ds!(_cylinder_tube(r_s, L_s; nθ=60, θ_end=PARTIAL),    fat_col)
            ds!(_cylinder_cap( r_s, -L_s/2; nθ=60, θ_end=PARTIAL), fat_col)
            ds!(_cylinder_cap( r_s,  L_s/2; nθ=60, θ_end=PARTIAL), fat_col)
        end
        if has_fur
            ds!(_cylinder_tube(r_i, L_i; nθ=60, θ_end=PARTIAL),    fur_col)
            ds!(_cylinder_cap( r_i, -L_i/2; nθ=60, θ_end=PARTIAL), fur_col)
            ds!(_cylinder_cap( r_i,  L_i/2; nθ=60, θ_end=PARTIAL), fur_col)
        end

    elseif body.shape isa Sphere
        ds!(_sphere_mesh(r_f),                           flesh_col)
        has_fat && ds!(_sphere_mesh(r_s; θ_end=PARTIAL), fat_col)
        has_fur && ds!(_sphere_mesh(r_i; θ_end=PARTIAL), fur_col)

    elseif body.shape isa Ellipsoid
        ratio = Float64(body.shape.b)
        ds!(_ellipsoid_mesh(r_f * ratio, r_f),                           flesh_col)
        has_fat && ds!(_ellipsoid_mesh(r_s * ratio, r_s; θ_end=PARTIAL), fat_col)
        has_fur && ds!(_ellipsoid_mesh(r_i * ratio, r_i; θ_end=PARTIAL), fur_col)

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

        _draw_box_faces!(hl_f, hw_f, hh_f, flesh_col, ds!; full=true)
        has_fat && _draw_box_faces!(hl_s, hw_s, hh_s, fat_col, ds!)
        has_fur && _draw_box_faces!(hl_i, hw_i, hh_i, fur_col, ds!)
    end

    return p
end

# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC 3-D API  (wrappers around the BodyCutaway recipe)
# ══════════════════════════════════════════════════════════════════════════════

"""
    draw_cutaway!(ax::Axis3, body; sc=100.0, offset=(0,0,0), axis=(0,0,1), flesh_col=…, fat_col=…, fur_col=…)

Draw a quarter-cutaway 3-D surface mesh of `body` into an existing `Axis3`.
`sc` converts metres to axis units (default 100 → cm labels).
`offset` translates the mesh (in axis units, i.e. after `sc` scaling).
`axis` rotates the mesh so the part's long axis aligns with this direction.
"""
function BiophysicalGeometry.draw_cutaway!(ax, body;
        sc        = 100.0,
        offset    = (0.0, 0.0, 0.0),
        axis      = (0.0, 0.0, 1.0),
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col   = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col   = RGBAf(0.76, 0.62, 0.42, 0.45))
    bodycutaway!(ax, body;
        sc,
        offset   = Vec3f(Float64.(offset)...),
        axis_dir = Vec3f(Float64.(axis)...),
        flesh_col, fat_col, fur_col)
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
# 3-D ORGANISM ASSEMBLY
# ══════════════════════════════════════════════════════════════════════════════

# Return the 3-D centre(s) to draw for a part, mirroring n=2 pairs.
# Upright: paired parts are symmetric about x=0 (left/right).
# Prone:   paired parts are symmetric about y=0 (left/right).
function _part_centers(pos, n::Int, ::Upright)
    c = pos.center
    n == 2 ? [c, (-c[1], c[2], c[3])] : [c]
end
function _part_centers(pos, n::Int, ::Prone)
    c = pos.center
    n == 2 ? [c, (c[1], -c[2], c[3])] : [c]
end

"""
    draw_organism!(ax, organism, posture=Upright(); sc=100.0, flesh_col=…, fat_col=…, fur_col=…)

Draw all body parts of a multi-part `Organism` as 3-D cutaway surfaces into an
existing `Axis3`. Parts with `n=2` are mirrored symmetrically.

`Plate` shapes (wings, fins) are drawn without axis rotation so their length
lies along the spine and width extends laterally.
"""
function BiophysicalGeometry.draw_organism!(ax, organism, posture=Upright();
        sc        = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col   = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col   = RGBAf(0.76, 0.62, 0.42, 0.45))

    positions = part_positions(organism, posture)

    for (name, bp) in pairs(organism.parts)
        pos  = positions[name]
        body = bp.body
        # Plate boxes are already in world orientation — skip axis rotation.
        ax_dir = body.shape isa Plate ? (0.0, 0.0, 1.0) : pos.axis
        for c3d in _part_centers(pos, bp.n, posture)
            draw_cutaway!(ax, body; sc,
                          offset = c3d .* sc,
                          axis   = ax_dir,
                          flesh_col, fat_col, fur_col)
        end
    end
    ax.xlabel = "x (cm)"; ax.ylabel = "y (cm)"; ax.zlabel = "z (cm)"
end

"""
    plot_body(organism, posture=Upright(); sc=100.0, flesh_col=…, fat_col=…, fur_col=…) -> Figure

Create a labelled `Figure` with `Axis3`, draw all body parts of a multi-part
`Organism` as 3-D cutaway surfaces, and return the `Figure`.
"""
function BiophysicalGeometry.plot_body(organism::Organism, posture::OrganismPosture=Upright();
        sc        = 100.0,
        flesh_col = RGBAf(0.88, 0.48, 0.42, 1.00),
        fat_col   = RGBAf(1.00, 0.97, 0.60, 0.75),
        fur_col   = RGBAf(0.76, 0.62, 0.42, 0.45))

    posture_name = string(nameof(typeof(posture)))
    # Upright: frontal view (camera at −y, arms spread left/right, head up, legs down)
    # Prone:   side view   (camera at +y, spine runs left/right, legs hang below)
    elevation = posture isa Prone ? π/8  : π/7
    azimuth   = posture isa Prone ? π/2  : 3π/2
    fig = Figure(size=(700, 620), backgroundcolor=:white)
    Label(fig[0, 1],
          "BiophysicalGeometry.jl — Organism  ($(posture_name), cutaway)";
          fontsize=13, font=:bold, padding=(0, 0, 10, 0))
    ax = Axis3(fig[1, 1];
               perspectiveness=0.5, viewmode=:fit, aspect=:data,
               elevation, azimuth)
    draw_organism!(ax, organism, posture; sc, flesh_col, fat_col, fur_col)
    Legend(fig[2, 1],
        [PolyElement(polycolor=flesh_col, strokecolor=:saddlebrown, strokewidth=1),
         PolyElement(polycolor=fat_col,   strokecolor=:saddlebrown, strokewidth=1),
         PolyElement(polycolor=fur_col,   strokecolor=:black,       strokewidth=1)],
        ["Flesh / muscle", "Subcutaneous fat", "Fur / insulation"];
        orientation=:horizontal, framevisible=false)
    return fig
end

# ══════════════════════════════════════════════════════════════════════════════
# 2-D ORGANISM SECTIONS (plot_organism — sagittal / coronal / transverse)
# Ported from NicheMapR/R/plot_human.R, generalised to any Organism topology.
# ══════════════════════════════════════════════════════════════════════════════

# Single unit-stripping boundary: called only inside drawing helpers below.
_pm(x::Unitful.Length) = Float32(ustrip(u"m", x))
_pm(x::Real)           = Float32(x)                  # already dimensionless metres

# Drawing helpers — accept Unitful lengths, produce Point2f for Makie.
_rect_m(cx, cy, hw, hh) =
    Point2f[(_pm(cx-hw), _pm(cy-hh)), (_pm(cx+hw), _pm(cy-hh)),
            (_pm(cx+hw), _pm(cy+hh)), (_pm(cx-hw), _pm(cy+hh))]

_ellipse_m(cx, cy, rx, ry; n=120) =
    [Point2f(_pm(cx + rx*cos(t)), _pm(cy + ry*sin(t))) for t in LinRange(0, 2π, n+1)]

# Upper half-ellipse (dorsal, θ ∈ [0, π]) — closed along the base.
# Matches R's draw.ellipse(..., segment = c(0, 180)).
_upper_ellipse_m(cx, cy, rx, ry; n=120) =
    [Point2f(_pm(cx + rx*cos(t)), _pm(cy + ry*sin(t))) for t in LinRange(0, π, n+1)]

# Lower half-ellipse (ventral, θ ∈ [π, 2π]) — matches R's segment = c(180, 360).
_lower_ellipse_m(cx, cy, rx, ry; n=120) =
    [Point2f(_pm(cx + rx*cos(t)), _pm(cy + ry*sin(t))) for t in LinRange(π, 2π, n+1)]

# Raw-bounds rectangle (xleft, ybottom, xright, ytop) — avoids centre+half arithmetic.
_rect_bounds(xl, yb, xr, yt) =
    Point2f[(_pm(xl),_pm(yb)), (_pm(xr),_pm(yb)), (_pm(xr),_pm(yt)), (_pm(xl),_pm(yt))]

# ── body-geometry helpers (Unitful) ─────────────────────────────────────────
_body_L_s(body) = body.geometry.length.length_skin
_body_L_i(body) = let gl = body.geometry.length
    hasproperty(gl, :length_fur) ? gl.length_fur :
        _body_L_s(body) + 2*(insulation_radius(body) - skin_radius(body))
end
_body_a_s(body) = body.geometry.length.a_semi_major_skin
_body_a_i(body) = let gl = body.geometry.length
    hasproperty(gl, :a_semi_major_fur) ? gl.a_semi_major_fur :
        _body_a_s(body) + (insulation_radius(body) - skin_radius(body))
end
_body_c_s(body) = body.geometry.length.c_semi_minor_skin
_body_c_i(body) = let gl = body.geometry.length
    hasproperty(gl, :c_semi_minor_fur) ? gl.c_semi_minor_fur :
        _body_c_s(body) + (insulation_radius(body) - skin_radius(body))
end

function BiophysicalGeometry.draw_organism_prone!(ax::Axis, organism::Organism;
        flesh_col    = RGBf(1.000, 0.753, 0.796),
        skin_col     = RGBf(1.000, 0.980, 0.600),   # yellow fat for prone
        fur_col      = RGBf(0.647, 0.165, 0.165),   # brown for prone
        arm_skin_col = RGBf(0.502, 0.502, 0.502))

    parts = organism.parts
    child_set  = Set(j.child for j in organism.joins)
    root_name  = first(k for k in keys(parts) if k ∉ child_set)
    root_body  = parts[root_name].body

    # Spine horizontal (x). Ground at y=0. Lateral joins = legs hanging down.
    # Forward/top join = head extending in +x from torso front.
    lat_joins  = [j for j in organism.joins if j.attachment == :lateral]
    fwd_joins  = [j for j in organism.joins if j.attachment ∈ (:top, :forward)]
    back_joins = [j for j in organism.joins if j.attachment == :back]
    # Ear joins: :top attachments whose parent is itself a child of the torso (chained)
    fwd_child_names = Set(j.child for j in fwd_joins)
    ear_joins  = [j for j in organism.joins if j.attachment == :top &&
                                                j.parent ∈ fwd_child_names]

    r_s_root  = skin_radius(root_body)
    r_i_root  = insulation_radius(root_body)
    r_f_root  = flesh_radius(root_body)
    L_s_root  = _body_L_s(root_body)
    trunk_ins = r_i_root - r_s_root
    trunk_fat = r_s_root - r_f_root

    # Head geometry
    has_hd = !isempty(fwd_joins)
    hd_body = has_hd ? parts[fwd_joins[1].child].body : nothing
    a_i_hd  = has_hd ? _body_a_i(hd_body) : 0u"m"
    a_s_hd  = has_hd ? _body_a_s(hd_body) : 0u"m"
    b_i_hd  = has_hd ? insulation_radius(hd_body) : 0u"m"
    b_s_hd  = has_hd ? skin_radius(hd_body) : 0u"m"
    fat_hd  = has_hd ? (skin_radius(hd_body) - flesh_radius(hd_body)) : 0u"m"

    # Tail geometry (back join = Cylinder extending rearward)
    has_tail  = !isempty(back_joins)
    tail_body = has_tail ? parts[back_joins[1].child].body : nothing
    r_s_tail  = has_tail ? skin_radius(tail_body)       : 0u"m"
    r_i_tail  = has_tail ? insulation_radius(tail_body) : 0u"m"
    r_f_tail  = has_tail ? flesh_radius(tail_body)      : 0u"m"
    L_s_tail  = has_tail ? _body_L_s(tail_body)         : 0u"m"
    tail_ins  = has_tail ? (r_i_tail - r_s_tail)        : 0u"m"
    tail_fat  = has_tail ? (r_s_tail - r_f_tail)        : 0u"m"

    # Ear/plate geometry (Plate children attached :top to head)
    has_ears  = !isempty(ear_joins)
    ear_body  = has_ears ? parts[ear_joins[1].child].body : nothing
    if has_ears
        gl_ear    = ear_body.geometry.length
        L_ear_s   = gl_ear.length_skin
        W_ear_s   = gl_ear.width_skin
        H_ear_s   = gl_ear.height_skin
        L_ear_i   = hasproperty(gl_ear, :length_fur) ? gl_ear.length_fur : L_ear_s + 2*(insulation_radius(ear_body) - skin_radius(ear_body))
        W_ear_i   = hasproperty(gl_ear, :width_fur)  ? gl_ear.width_fur  : W_ear_s + 2*(insulation_radius(ear_body) - skin_radius(ear_body))
        H_ear_i   = hasproperty(gl_ear, :height_fur) ? gl_ear.height_fur : H_ear_s + 2*(insulation_radius(ear_body) - skin_radius(ear_body))
        fat_ear   = skin_radius(ear_body) - flesh_radius(ear_body)
    else
        L_ear_s = W_ear_s = H_ear_s = L_ear_i = W_ear_i = H_ear_i = fat_ear = 0u"m"
    end

    # Leg ground-clearance: tallest leg determines torso height
    has_legs  = !isempty(lat_joins)
    max_leg_h = has_legs ?
        maximum(j -> _body_L_s(parts[j.child].body) + (insulation_radius(parts[j.child].body) - skin_radius(parts[j.child].body)), lat_joins) :
        0u"m"
    max_leg_r = has_legs ? maximum(j -> insulation_radius(parts[j.child].body), lat_joins) : 0u"m"

    # y-coords: torso sits on top of legs
    y_tor_sk_bot  = max_leg_h
    y_tor_center  = y_tor_sk_bot + r_s_root
    y_tor_sk_top  = y_tor_sk_bot + 2*r_s_root
    y_tor_fur_top = y_tor_sk_top + trunk_ins

    # Panel spacing
    head_horiz  = has_hd ? a_i_hd : trunk_ins
    min_spacing = L_s_root + 2*(trunk_ins + head_horiz)
    spacing     = max(min_spacing * 1.1, 2*(r_i_root + 2*max_leg_r) + r_i_root)

    SAG_OFF = -spacing/2
    COR_OFF = +spacing/2

    vlines!(ax, [0.0];
            color=:grey80, linestyle=:dash, linewidth=0.8)

    pfill!(pts, col) = poly!(ax, pts; color=col, strokecolor=col, strokewidth=0.3)
    prect(xl, yb, xr, yt, col) = pfill!(_rect_bounds(xl, yb, xr, yt), col)
    pell(cx, cy, rx, ry, col)  = pfill!(_ellipse_m(cx, cy, rx, ry),   col)

    # ── SIDE VIEW (sagittal, centred on SAG_OFF) ──────────────────────────
    prect(SAG_OFF - L_s_root/2 - trunk_ins, y_tor_sk_bot - trunk_ins,
          SAG_OFF + L_s_root/2 + trunk_ins, y_tor_fur_top, fur_col)
    prect(SAG_OFF - L_s_root/2, y_tor_sk_bot,
          SAG_OFF + L_s_root/2, y_tor_sk_top, skin_col)
    prect(SAG_OFF - L_s_root/2 + trunk_fat, y_tor_sk_bot + trunk_fat,
          SAG_OFF + L_s_root/2 - trunk_fat, y_tor_sk_top - trunk_fat, flesh_col)
    if has_hd
        cx_hd = SAG_OFF + L_s_root/2
        cy_hd = y_tor_sk_top + b_i_hd
        pell(cx_hd, cy_hd, a_i_hd, b_i_hd, fur_col)
        pell(cx_hd, cy_hd, a_s_hd, b_s_hd, arm_skin_col)
        pell(cx_hd, cy_hd, a_s_hd - fat_hd, b_s_hd - fat_hd, flesh_col)
    end
    for j in lat_joins
        b = parts[j.child].body
        r_s_l = skin_radius(b);  r_i_l = insulation_radius(b);  r_f_l = flesh_radius(b)
        l_ins = r_i_l - r_s_l;  l_fat = r_s_l - r_f_l;  L_l = _body_L_s(b)
        xc = SAG_OFF + j.spine_fraction * L_s_root/2
        prect(xc - r_i_l, 0u"m" - l_ins, xc + r_i_l, L_l + l_ins, fur_col)
        prect(xc - r_s_l, 0u"m",          xc + r_s_l, L_l,          skin_col)
        prect(xc - r_f_l, l_fat,           xc + r_f_l, L_l - l_fat,  flesh_col)
    end
    if has_tail
        x_tail_attach = SAG_OFF - L_s_root/2
        prect(x_tail_attach - tail_ins - L_s_tail, y_tor_center - r_i_tail,
              x_tail_attach + tail_ins,             y_tor_center + r_i_tail, fur_col)
        prect(x_tail_attach - L_s_tail, y_tor_center - r_s_tail,
              x_tail_attach,            y_tor_center + r_s_tail, skin_col)
        prect(x_tail_attach - L_s_tail + tail_fat, y_tor_center - r_f_tail,
              x_tail_attach - tail_fat,             y_tor_center + r_f_tail, flesh_col)
    end
    if has_ears && has_hd
        cx_hd_s = SAG_OFF + L_s_root/2
        cy_hd_s = y_tor_sk_top + b_i_hd
        prect(cx_hd_s - L_ear_i/2, cy_hd_s - W_ear_i,
              cx_hd_s + L_ear_i/2, cy_hd_s,           fur_col)
        prect(cx_hd_s - L_ear_s/2, cy_hd_s - W_ear_s,
              cx_hd_s + L_ear_s/2, cy_hd_s,           arm_skin_col)
        prect(cx_hd_s - L_ear_s/2 + fat_ear, cy_hd_s - W_ear_s + fat_ear,
              cx_hd_s + L_ear_s/2 - fat_ear, cy_hd_s - fat_ear, flesh_col)
    end

    # ── FRONT VIEW (coronal, centred on COR_OFF) ──────────────────────────
    if has_tail
        pell(COR_OFF, y_tor_center, r_i_tail, r_i_tail, fur_col)
    end
    pell(COR_OFF, y_tor_center, r_i_root, r_i_root, fur_col)
    pell(COR_OFF, y_tor_center, r_s_root, r_s_root, skin_col)
    pell(COR_OFF, y_tor_center, r_f_root, r_f_root, flesh_col)
    if has_hd
        cy_hd_cor = y_tor_sk_top + b_i_hd
        pell(COR_OFF, cy_hd_cor, b_i_hd, b_i_hd, fur_col)
        pell(COR_OFF, cy_hd_cor, b_s_hd, b_s_hd, arm_skin_col)
        pell(COR_OFF, cy_hd_cor, b_s_hd - fat_hd, b_s_hd - fat_hd, flesh_col)
    end
    fwd_lats = filter(j -> j.spine_fraction > 0, lat_joins)
    if !isempty(fwd_lats)
        jf = first(fwd_lats)
        b  = parts[jf.child].body
        r_s_l = skin_radius(b);  r_i_l = insulation_radius(b);  r_f_l = flesh_radius(b)
        l_ins = r_i_l - r_s_l;  l_fat = r_s_l - r_f_l;  L_l = _body_L_s(b)
        cx_off = r_s_root/10 + r_i_l
        for sg in (1, -1)
            cx = COR_OFF + sg * cx_off
            prect(cx - r_i_l, 0u"m" - l_ins, cx + r_i_l, L_l + l_ins, fur_col)
            prect(cx - r_s_l, 0u"m",          cx + r_s_l, L_l,          skin_col)
            prect(cx - r_f_l, l_fat,           cx + r_f_l, L_l - l_fat,  flesh_col)
        end
    end
    if has_ears && has_hd
        cy_hd_cor = y_tor_sk_top + b_i_hd
        for sg in (1, -1)
            x_ear_inner = COR_OFF + sg * b_i_hd
            x_ear_outer = x_ear_inner + sg * H_ear_i
            prect(min(x_ear_inner, x_ear_outer), cy_hd_cor - W_ear_i,
                  max(x_ear_inner, x_ear_outer), cy_hd_cor,           fur_col)
            x_ear_inner_s = COR_OFF + sg * b_s_hd
            x_ear_outer_s = x_ear_inner_s + sg * H_ear_s
            prect(min(x_ear_inner_s, x_ear_outer_s), cy_hd_cor - W_ear_s,
                  max(x_ear_inner_s, x_ear_outer_s), cy_hd_cor,           arm_skin_col)
            prect(min(x_ear_inner_s, x_ear_outer_s) + fat_ear, cy_hd_cor - W_ear_s + fat_ear,
                  max(x_ear_inner_s, x_ear_outer_s) - fat_ear, cy_hd_cor - fat_ear, flesh_col)
        end
    end

    return (; SAG_OFF, COR_OFF, spacing, L_s_root, trunk_ins, L_s_tail, tail_ins,
              b_i_hd, H_ear_i, max_leg_h, r_s_root)
end

function BiophysicalGeometry.draw_organism_upright!(ax::Axis, organism::Organism;
        flesh_col      = RGBf(1.000, 0.753, 0.796),
        skin_col       = RGBf(0.827, 0.827, 0.827),
        fur_col        = RGBf(0.678, 0.847, 0.902),
        head_fur_col   = RGBf(0.647, 0.165, 0.165),
        beard_fur_col  = RGBf(0.647, 0.165, 0.165),
        leg_fur_col    = RGBf(0.941, 0.902, 0.549),
        arm_skin_col   = RGBf(0.502, 0.502, 0.502),
        beard_depth    = 0.0u"m")

    parts = organism.parts
    child_set  = Set(j.child for j in organism.joins)
    root_name  = first(k for k in keys(parts) if k ∉ child_set)
    root_body  = parts[root_name].body

    top_names     = Symbol[j.child for j in organism.joins if j.attachment ∈ (:top, :forward)]
    bottom_names  = Symbol[j.child for j in organism.joins if j.attachment == :bottom]
    lateral_names = Symbol[j.child for j in organism.joins if j.attachment == :lateral]

    r_s_root = skin_radius(root_body)
    r_i_root = insulation_radius(root_body)
    r_f_root = flesh_radius(root_body)
    L_s_root = _body_L_s(root_body)

    has_bottom = !isempty(bottom_names)
    leg_body   = has_bottom ? parts[bottom_names[1]].body : nothing
    r_s_leg    = has_bottom ? skin_radius(leg_body)        : 0u"m"
    r_i_leg    = has_bottom ? insulation_radius(leg_body)  : 0u"m"
    r_f_leg    = has_bottom ? flesh_radius(leg_body)       : 0u"m"
    L_s_leg    = has_bottom ? _body_L_s(leg_body)          : 0u"m"
    z_leg_top  = L_s_leg

    has_top    = !isempty(top_names)
    head_body  = has_top ? parts[top_names[1]].body : nothing
    a_s_head   = has_top ? _body_a_s(head_body)           : 0u"m"
    a_i_head   = has_top ? _body_a_i(head_body)           : 0u"m"
    b_s_head   = has_top ? skin_radius(head_body)          : 0u"m"
    b_i_head   = has_top ? insulation_radius(head_body)    : 0u"m"
    c_s_head   = has_top ? _body_c_s(head_body)           : 0u"m"
    c_i_head   = has_top ? _body_c_i(head_body)           : 0u"m"
    fat_d_head = has_top ? (skin_radius(head_body) - flesh_radius(head_body)) : 0u"m"

    has_lat    = !isempty(lateral_names)
    arm_body   = has_lat ? parts[lateral_names[1]].body : nothing
    r_s_arm    = has_lat ? skin_radius(arm_body)         : 0u"m"
    r_i_arm    = has_lat ? insulation_radius(arm_body)   : 0u"m"
    r_f_arm    = has_lat ? flesh_radius(arm_body)        : 0u"m"
    L_s_arm    = has_lat ? _body_L_s(arm_body)           : 0u"m"

    z_root_bottom = z_leg_top
    z_root_top    = z_root_bottom + L_s_root
    z_head_center = z_root_top + a_s_head

    total_height   = z_leg_top + L_s_root + 2*a_s_head
    arm_half_width = has_lat ? (r_i_root + r_s_arm + r_s_arm/4 + r_i_arm) : r_i_root
    spacing        = max(total_height * 0.6, arm_half_width * 2 + r_i_root)

    mid_height = total_height / 2
    top_height = mid_height * 1.5
    bot_height = mid_height / 2

    SAG_OFF = -spacing
    COR_OFF = 0u"m"
    TRA_OFF = +spacing

    vlines!(ax, [_pm(-spacing/2), _pm(spacing/2)];
            color=:grey80, linestyle=:dash, linewidth=0.8)

    fill!(pts, col) = poly!(ax, pts; color=col, strokecolor=col, strokewidth=0.3)
    frect(xl, yb, xr, yt, col) = fill!(_rect_bounds(xl, yb, xr, yt), col)
    fell(cx, cy, rx, ry, col)  = fill!(_ellipse_m(cx, cy, rx, ry),   col)

    trunk_ins = r_i_root - r_s_root
    trunk_fat = r_s_root - r_f_root

    # ── CORONAL PANEL ─────────────────────────────────────────────────────────
    frect(-r_i_root + COR_OFF, z_root_bottom - trunk_ins,
           r_i_root + COR_OFF, z_root_top    + trunk_ins, fur_col)
    frect(-r_s_root + trunk_ins + COR_OFF, z_root_bottom,
           r_s_root - trunk_ins + COR_OFF, z_root_top,    skin_col)
    frect(-(r_s_root - trunk_ins - trunk_fat) + COR_OFF, z_root_bottom + trunk_fat,
           (r_s_root - trunk_ins - trunk_fat) + COR_OFF, z_root_top    - trunk_fat, flesh_col)
    if has_top
        fill!(_upper_ellipse_m(COR_OFF, z_head_center, b_i_head, a_i_head), head_fur_col)
        if beard_depth > 0u"m"
            b_i_beard = b_s_head + beard_depth
            a_i_beard = a_s_head + beard_depth
            fill!(_lower_ellipse_m(COR_OFF, z_head_center, b_i_beard, a_i_beard), beard_fur_col)
        end
        fell(COR_OFF, z_head_center, b_s_head, a_s_head, arm_skin_col)
        fell(COR_OFF, z_head_center, b_s_head - fat_d_head, a_s_head - fat_d_head, flesh_col)
    end
    if has_bottom
        leg_ins = r_i_leg - r_s_leg
        leg_fat = r_s_leg - r_f_leg
        xl_leg = -r_s_root - r_s_root/5 - leg_ins
        xr_leg = xl_leg + 2*r_s_leg + 2*leg_ins
        for (xl_f, xr_f) in ((xl_leg, xr_leg), (-xr_leg, -xl_leg))
            frect(xl_f,            -leg_ins, xr_f,            z_leg_top + leg_ins, leg_fur_col)
            frect(xl_f + leg_ins,  0u"m",    xr_f - leg_ins,  z_leg_top,           skin_col)
            frect(xl_f + leg_ins + leg_fat, leg_fat,
                  xr_f - leg_ins - leg_fat, z_leg_top - leg_fat, flesh_col)
        end
    end
    if has_lat
        arm_ins = r_i_arm - r_s_arm
        arm_fat = r_s_arm - r_f_arm
        ytop_arm_fur    = z_root_top + arm_ins
        ybottom_arm_fur = a_s_head + z_leg_top + (r_i_leg - r_s_leg) +
                          L_s_root - L_s_arm - arm_ins
        xr_arm = -(r_i_root + r_s_arm/2)
        xl_arm = -(r_s_root + 1.5*r_s_arm + r_i_arm)
        for (xl_f, xr_f) in ((xl_arm, xr_arm), (-xr_arm, -xl_arm))
            frect(xl_f,            ybottom_arm_fur, xr_f,            ytop_arm_fur, fur_col)
            frect(xl_f + arm_ins,  ybottom_arm_fur + arm_ins,
                  xr_f - arm_ins,  ytop_arm_fur    - arm_ins, arm_skin_col)
            frect(xl_f + arm_ins + arm_fat, ybottom_arm_fur + arm_ins,
                  xr_f - arm_ins - arm_fat, ytop_arm_fur    - arm_ins, flesh_col)
        end
    end

    # ── SAGITTAL PANEL ────────────────────────────────────────────────────────
    frect(SAG_OFF - r_i_root, z_root_bottom + (r_i_leg - r_s_leg),
          SAG_OFF + r_i_root, z_root_top    + (r_i_leg - r_s_leg), fur_col)
    frect(SAG_OFF - r_s_root + trunk_ins, z_root_bottom + (r_i_leg - r_s_leg),
          SAG_OFF + r_s_root - trunk_ins, z_root_top    + (r_i_leg - r_s_leg), skin_col)
    frect(SAG_OFF - r_s_root + trunk_ins + trunk_fat, z_root_bottom + (r_i_leg - r_s_leg) + trunk_fat,
          SAG_OFF + r_s_root - trunk_ins - trunk_fat, z_root_top    + (r_i_leg - r_s_leg) - trunk_fat, flesh_col)
    if has_top
        fill!(_upper_ellipse_m(SAG_OFF, z_head_center, c_i_head, a_i_head), head_fur_col)
        if beard_depth > 0u"m"
            c_i_beard = c_s_head + beard_depth
            a_i_beard_s = a_s_head + beard_depth
            fill!(_lower_ellipse_m(SAG_OFF, z_head_center, c_i_beard, a_i_beard_s), beard_fur_col)
        end
        fell(SAG_OFF, z_head_center, c_s_head, a_s_head, arm_skin_col)
        fell(SAG_OFF, z_head_center, c_s_head - fat_d_head, a_s_head - fat_d_head, flesh_col)
    end
    if has_bottom
        leg_ins = r_i_leg - r_s_leg
        leg_fat = r_s_leg - r_f_leg
        frect(SAG_OFF - r_i_leg,           0u"m",         SAG_OFF + r_i_leg, z_leg_top + leg_ins, leg_fur_col)
        frect(SAG_OFF - r_s_leg,           leg_ins,        SAG_OFF + r_s_leg, z_leg_top,           skin_col)
        frect(SAG_OFF - r_s_leg + leg_fat, leg_ins + leg_fat,
              SAG_OFF + r_s_leg - leg_fat, z_leg_top - leg_fat, flesh_col)
    end
    if has_lat
        arm_ins = r_i_arm - r_s_arm
        arm_fat = r_s_arm - r_f_arm
        ytop_arm_fur    = z_root_top + arm_ins
        ybottom_arm_fur = a_s_head + z_leg_top + (r_i_leg - r_s_leg) +
                          L_s_root - L_s_arm - arm_ins
        frect(SAG_OFF - r_i_arm, ybottom_arm_fur, SAG_OFF + r_i_arm, ytop_arm_fur, fur_col)
        frect(SAG_OFF - r_s_arm, ybottom_arm_fur + arm_ins,
              SAG_OFF + r_s_arm, ytop_arm_fur    - arm_ins, arm_skin_col)
        frect(SAG_OFF - r_s_arm + arm_fat, ybottom_arm_fur + arm_ins,
              SAG_OFF + r_s_arm - arm_fat, ytop_arm_fur    - arm_ins, flesh_col)
    end

    # ── TRANSVERSE PANEL ──────────────────────────────────────────────────────
    hlines!(ax, [_pm(mid_height), _pm(bot_height)];
            color=:grey80, linestyle=:dash, linewidth=0.8)
    fell(TRA_OFF, top_height, r_i_root, r_i_root, fur_col)
    if has_lat
        arm_ins_tra = r_i_arm - r_s_arm
        for sg in (1, -1)
            cx = TRA_OFF + sg*(r_s_root + 1.5*r_s_arm + trunk_ins + arm_ins_tra)
            fell(cx, top_height, r_i_arm, r_i_arm, fur_col)
        end
    end
    if has_top
        fell(TRA_OFF, top_height, b_i_head, a_i_head, head_fur_col)
    end
    fell(TRA_OFF, mid_height, r_i_root, r_i_root, fur_col)
    fell(TRA_OFF, mid_height, r_s_root, r_s_root, skin_col)
    fell(TRA_OFF, mid_height, r_f_root, r_f_root, flesh_col)
    if has_lat
        arm_ins = r_i_arm - r_s_arm
        arm_fat = r_s_arm - r_f_arm
        for sg in (1, -1)
            cx = TRA_OFF + sg*(r_s_root + 1.5*r_s_arm + trunk_ins + arm_ins)
            fell(cx, mid_height, r_i_arm, r_i_arm, fur_col)
            fell(cx, mid_height, r_s_arm, r_s_arm, arm_skin_col)
            fell(cx, mid_height, r_f_arm, r_f_arm, flesh_col)
        end
    end
    if has_bottom
        leg_ins = r_i_leg - r_s_leg
        leg_fat = r_s_leg - r_f_leg
        for sg in (1, -1)
            cx = TRA_OFF + sg*(r_s_root/10 + r_i_leg)
            fell(cx, bot_height, r_i_leg, r_i_leg, leg_fur_col)
            fell(cx, bot_height, r_s_leg, r_s_leg, skin_col)
            fell(cx, bot_height, r_f_leg, r_f_leg, flesh_col)
        end
    end

    return (; spacing, total_height)
end

"""
    plot_organism(organism, posture=Upright(); flesh_col=…, skin_col=…, fur_col=…) -> Figure

2-D anatomical schematic with **sagittal**, **coronal**, and **transverse** panels,
matching the layout of NicheMapR's `plot_human.R`.

Parts are drawn outer-to-inner (fur → skin/fat ring → flesh) as solid fills so
each inner layer paints over the outer, revealing tissue rings — exactly as R's
`rect()` / `draw.ellipse()` layering. Ground level (`y = 0`) is at the feet.

Default colours mirror R's palette (lightblue trunk/arms, khaki legs, brown head,
lightgrey/grey skin, pink flesh). Pass keyword arguments to override.
"""
function BiophysicalGeometry.plot_organism(organism::Organism,
        posture::OrganismPosture = Upright(); kwargs...)
    if posture isa Prone
        fig = Figure(backgroundcolor=:white)
        posture_name = string(nameof(typeof(posture)))
        Label(fig[0, 1], "sagittal section / coronal section  ($posture_name)";
              fontsize=13, font=:bold, padding=(0, 0, 6, 0))
        ax = Axis(fig[1, 1]; xlabel="metres", ylabel="metres", aspect=DataAspect())
        info = draw_organism_prone!(ax, organism; kwargs...)
        Legend(fig[2, 1],
            [PolyElement(polycolor=RGBf(1,0.753,0.796),     strokecolor=:black, strokewidth=0.5),
             PolyElement(polycolor=RGBf(1,0.98,0.6),        strokecolor=:black, strokewidth=0.5),
             PolyElement(polycolor=RGBf(0.647,0.165,0.165), strokecolor=:black, strokewidth=0.5)],
            ["Flesh", "Fat", "Fur"]; orientation=:horizontal, framevisible=false)
        x_content   = _pm(-info.spacing/2 - info.L_s_root/2 - info.trunk_ins - info.L_s_tail - info.tail_ins)
        x_content_r = _pm( info.spacing/2 + info.b_i_hd + info.H_ear_i)
        y_content   = _pm(info.max_leg_h + 2*info.r_s_root + 2*info.b_i_hd)
        cell_ratio  = (x_content_r - x_content) * 1.1f0 / (y_content * 1.1f0)
        colsize!(fig.layout, 1, Aspect(1, cell_ratio))
        resize_to_layout!(fig)
        return fig
    else  # Upright
        fig = Figure(backgroundcolor=:white)
        posture_name = string(nameof(typeof(posture)))
        Label(fig[0, 1], "sagittal section / coronal section / transverse section  ($posture_name)";
              fontsize=13, font=:bold, padding=(0, 0, 6, 0))
        ax = Axis(fig[1, 1]; xlabel="metres", ylabel="metres", aspect=DataAspect())
        info = draw_organism_upright!(ax, organism; kwargs...)
        Legend(fig[2, 1],
            [PolyElement(polycolor=RGBf(1.000,0.753,0.796), strokecolor=:black, strokewidth=0.5),
             PolyElement(polycolor=RGBf(0.827,0.827,0.827), strokecolor=:black, strokewidth=0.5),
             PolyElement(polycolor=RGBf(0.678,0.847,0.902), strokecolor=:black, strokewidth=0.5),
             PolyElement(polycolor=RGBf(0.647,0.165,0.165), strokecolor=:black, strokewidth=0.5),
             PolyElement(polycolor=RGBf(0.941,0.902,0.549), strokecolor=:black, strokewidth=0.5)],
            ["Flesh", "Fat/skin", "Fur (trunk/arms)", "Fur (head)", "Fur (legs)"];
            orientation=:horizontal, framevisible=false)
        rowgap!(fig.layout, 8)
        x_content = _pm(3 * info.spacing)
        y_content = _pm(info.total_height)
        colsize!(fig.layout, 1, Aspect(1, x_content / y_content * 1.1f0))
        resize_to_layout!(fig)
        return fig
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

        r_i > r_s && _layer!(ax_long, _rect_pts(hl_i, hh_i), fur_col)
        r_s > r_f && _layer!(ax_long, _rect_pts(hl_s, hh_s), fat_col)
                     _layer!(ax_long, _rect_pts(hl_f, hh_f), flesh_col)

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
