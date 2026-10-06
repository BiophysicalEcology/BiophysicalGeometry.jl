"""
    Ellipsoid(; mass, density, volume, length, width, height, axis_ratio_b, axis_ratio_c,
              pole_a_truncation=0.0) <: AbstractShape

A general (triaxial) ellipsoid centred on the origin: `length` along `x`, `width`
along `y`, `height` along `z`, all full extents (twice the semi-axes) at skin
level. `axis_ratio_b` is length / width and `axis_ratio_c` length / height, as
for [`Plate`](@ref). Give any sufficient set of keywords and the rest is solved
for.

With `pole_a_truncation > 0`, the `+x` end is sliced flat: the cut plane sits at
`x = (1 - pole_a_truncation) * length / 2`. `pole_a_truncation = 0` is a full
ellipsoid; `1` cuts through the centre. `length`, `width` and `height` are those of
the ellipsoid before the cut; volume, mass and surface area are of the cut body.
"""
struct Ellipsoid{M,D,B,C,T} <: AbstractEllipsoidal
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
    pole_a_truncation::T
    Ellipsoid(::Resolved, mass::M, density::D, axis_ratio_b::B, axis_ratio_c::C,
              pole_a_truncation::T) where {M,D,B,C,T} =
        new{M,D,B,C,T}(mass, density, axis_ratio_b, axis_ratio_c, pole_a_truncation)
end

# volume = π·k·a·b·c = (π·k/8)·length·width·height, with k = 4/3 for a full
# ellipsoid (see `_truncated_volume_factor`); ratios as for the box.
_ellipsoid_spec(truncation) =
    ShapeSpec((:length, :width, :height), (1, 1, 1), log(π * _truncated_volume_factor(truncation) / 8),
               (:axis_ratio_b => (1, 2, 1.0), :axis_ratio_c => (1, 3, 1.0)))

function Ellipsoid(; pole_a_truncation = 0.0, kw...)
    0 <= pole_a_truncation <= 1 || throw(ArgumentError(
        "Ellipsoid `pole_a_truncation` must be in [0, 1], got $pole_a_truncation"))
    s = _resolve_shape("Ellipsoid", _ellipsoid_spec(pole_a_truncation), NamedTuple(kw))
    Ellipsoid(RESOLVED, s.mass, s.density, s.axis_ratio_b, s.axis_ratio_c, pole_a_truncation)
end

# x-position of the truncated pole_a (as a fraction of a). 1.0 = full ellipsoid.
_pole_a_x_ratio(s::Ellipsoid) = 1 - s.pole_a_truncation
# Radial scale at the truncated pole (so cut disc has y/z extents = scale * b/c).
_pole_a_radial_scale(s::Ellipsoid) = sqrt(max(0.0, 1 - _pole_a_x_ratio(s)^2))

# Volume of the cut ellipsoid as π·k·a·b·c. The cap beyond x = x_r·a holds
# π·b·c·∫(1 - x²/a²)dx = π·a·b·c·[(1 - x_r) - (1 - x_r³)/3], so
#     k = 4/3 - (1 - x_r) + (1 - x_r³)/3,
# from 4/3 uncut down to 2/3 for a cut through the centre.
function _truncated_volume_factor(truncation)
    x = 1 - truncation
    4 / 3 - (1 - x) + (1 - x^3) / 3
end
_truncated_volume_factor(s::Ellipsoid) = _truncated_volume_factor(s.pole_a_truncation)

# The geometry stores full extents (length, width, height); the formulas below
# work in semi-axes a ≥ b, c along x, y, z.
_semiaxes(d) = (d.length / 2, d.width / 2, d.height / 2)
_skin_semiaxes(l) = (l.length_skin / 2, l.width_skin / 2, l.height_skin / 2)
_extents(level, (a, b, c)) = level === :skin ?
    (; length_skin = 2a, width_skin = 2b, height_skin = 2c) :
    (; length_fibrous = 2a, width_fibrous = 2b, height_fibrous = 2c)

# ── Surface area ───────────────────────────────────────────────────────────
#
# The surface area of a general ellipsoid is not elementary. With Carlson's
# symmetric elliptic integral of the second kind it is exactly
#     S = 4π·a·b·c·R_G(a⁻², b⁻², c⁻²),
#     2R_G(x, y, z) = z·R_F(x, y, z) - (x - z)(y - z)·R_D(x, y, z)/3 + √(xy/z),
# which needs no axis ordering and reduces to the familiar prolate (asin) and
# oblate (log) spheroid formulas when two axes are equal. R_F and R_D use
# Carlson's duplication algorithm (Numer. Algorithms 10, 1995) with a fixed
# number of steps — no data-dependent branching, which keeps AD well-behaved.
# Each step quarters the spread of the arguments, so 12 steps leave the fifth-
# order series accurate to machine precision. The formula is unitless, so
# lengths are stripped to metres once here.

function _carlson_rf(x, y, z)
    for _ in 1:12
        λ = sqrt(x * y) + sqrt(y * z) + sqrt(z * x)
        x, y, z = (x + λ) / 4, (y + λ) / 4, (z + λ) / 4
    end
    A = (x + y + z) / 3
    X, Y = 1 - x / A, 1 - y / A
    Z = -(X + Y)
    E2 = X * Y - Z^2
    E3 = X * Y * Z
    (1 - E2 / 10 + E3 / 14 + E2^2 / 24 - 3 * E2 * E3 / 44) / sqrt(A)
end

function _carlson_rd(x, y, z)
    sum = zero(x)
    scale = one(x)
    for _ in 1:12
        λ = sqrt(x * y) + sqrt(y * z) + sqrt(z * x)
        sum += scale / (sqrt(z) * (z + λ))
        scale /= 4
        x, y, z = (x + λ) / 4, (y + λ) / 4, (z + λ) / 4
    end
    A = (x + y + 3 * z) / 5
    X, Y = (A - x) / A, (A - y) / A
    Z = -(X + Y) / 3
    E2 = X * Y - 6 * Z^2
    E3 = (3 * X * Y - 8 * Z^2) * Z
    E4 = 3 * (X * Y - Z^2) * Z^2
    E5 = X * Y * Z^3
    3 * sum + scale / A^(3 / 2) *
        (1 - 3 * E2 / 14 + E3 / 6 + 9 * E2^2 / 88 - 3 * E4 / 22 - 9 * E2 * E3 / 52 + 3 * E5 / 26)
end

_carlson_rg(x, y, z) =
    (z * _carlson_rf(x, y, z) - (x - z) * (y - z) * _carlson_rd(x, y, z) / 3 + sqrt(x * y / z)) / 2

function _ellipsoid_area(a, b, c)
    am, bm, cm = ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c)
    4π * am * bm * cm * _carlson_rg(1 / am^2, 1 / bm^2, 1 / cm^2) * u"m^2"
end

# Area of a cut ellipsoid: the full surface, less the cap beyond the cut plane
# x = x_r·a, plus the flat elliptical disc it leaves (semi-axes b·s, c·s with
# s = √(1 - x_r²)). The cap of a triaxial ellipsoid has no closed form, so it is
# integrated over x = a·cos α, y = b·sin α·cos β, z = c·sin α·sin β for
# α ∈ [0, acos x_r]: Gauss–Legendre in α (smooth, non-periodic) and the
# trapezoid rule in β (periodic, so spectrally accurate). The surface element is
#     |∂r/∂α × ∂r/∂β| = sin α·√(b²c²·cos²α + a²·sin²α·(c²·cos²β + b²·sin²β)).
function _ellipsoid_area(a, b, c, truncation)
    full = _ellipsoid_area(a, b, c)
    truncation == 0 && return full
    am, bm, cm = ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c)
    x = 1 - truncation
    α0 = acos(x)
    nβ = length(CAP_ANGLES)
    cap = 0.0
    for (t, w) in GAUSS_LEGENDRE
        α = α0 * (t + 1) / 2
        sα, cα = sincos(α)
        ring = 0.0
        for β in CAP_ANGLES
            sβ, cβ = sincos(β)
            ring += sqrt((bm * cm * cα)^2 + (am * sα)^2 * ((cm * cβ)^2 + (bm * sβ)^2))
        end
        cap += w * sα * ring * 2π / nβ
    end
    cap *= α0 / 2
    disc = π * bm * cm * (1 - x^2)
    full + (disc - cap) * u"m^2"
end

# Gauss–Legendre nodes and weights on [-1, 1] by Newton's method on P_n,
# computed once at load time.
function _gauss_legendre(n)
    map(1:n) do i
        t = cos(π * (i - 0.25) / (n + 0.5))
        dp = 0.0
        for _ in 1:20
            p0, p1 = 1.0, t
            for k in 2:n
                p0, p1 = p1, ((2k - 1) * t * p1 - (k - 1) * p0) / k
            end
            dp = n * (t * p1 - p0) / (t^2 - 1)
            t -= p1 / dp
        end
        (t, 2 / ((1 - t^2) * dp^2))
    end
end
const GAUSS_LEGENDRE = _gauss_legendre(24)
const CAP_ANGLES = [2π * (j - 0.5) / 64 for j in 1:64]

# ── Geometry ───────────────────────────────────────────────────────────────

# Canonicalise to m^3 to avoid weird unit ratios from mass/density (e.g. g/(kg/m^3))
# that propagate into cbrt() and trip up Enzyme's typeunstablerules path.
_ellipsoid_volume(shape::Ellipsoid) = uconvert(u"m^3", body_volume(shape))
_ellipsoid_fat_volume(shape::Ellipsoid, fat_layer::FatLayer) =
    uconvert(u"m^3", fat_volume(shape, fat_layer))

# Semi-axes enclosing `volume` at the shape's proportions and truncation: with
# a = r_b·b and a = r_c·c, volume = π·k·a·b·c = π·k·a³/(r_b·r_c).
function _ellipsoid_semiaxes(shape::Ellipsoid, volume)
    a = cbrt(volume * shape.axis_ratio_b * shape.axis_ratio_c / (π * _truncated_volume_factor(shape)))
    (a, a / shape.axis_ratio_b, a / shape.axis_ratio_c)
end

# Area of the shape's (cut) surface at semi-axes `axes`.
_shape_area(shape::Ellipsoid, axes) = _ellipsoid_area(axes..., shape.pole_a_truncation)

function _skin_level(shape::Ellipsoid, volume)
    axes = _ellipsoid_semiaxes(shape, volume)
    (; dims = _extents(:skin, axes), area = _shape_area(shape, axes))
end
function _fibrous_level(shape::Ellipsoid, skin, thickness)
    axes = _skin_semiaxes(skin) .+ thickness
    (; dims = _extents(:fibrous, axes), area = _shape_area(shape, axes))
end

# Smooth Heaviside on a length: ≈ 1 for fat ≫ ε, ≈ 0 for fat ≪ -ε, smooth
# everywhere. Used to blend the "no fat" and "with fat" geometry branches.
@inline _smooth_step_meters(fat) = 0.5 * (1 + fat / sqrt(fat*fat + (1e-9u"m")^2))

# Ellipsoid fat is a *uniform* shell: each skin semi-axis is the flesh semi-axis
# plus one fat thickness (so the skin is not a scaled flesh ellipsoid). Skin dims
# are therefore coupled to the fat solve, so the fat paths override the
# generic orchestrator rather than compute skin from volume independently.
#
# The "not enough fat" fallback (full-volume axes, fat = 0) and the flesh + fat
# axes are smooth-blended rather than switched with an if/else, which would be
# a value discontinuity at raw_fat = 0 — bad for reverse-mode AD. A single
# smooth Heaviside drives both the axes blend and the effective fat value, so
# they stay consistent at every raw_fat.
function _ellipsoid_fat_skin(shape::Ellipsoid, fat_layer::FatLayer)
    volume = _ellipsoid_volume(shape)
    fat_v = _ellipsoid_fat_volume(shape, fat_layer)
    flesh = _ellipsoid_semiaxes(shape, volume - fat_v)
    raw_fat = _uniform_shell_thickness(flesh, volume, _truncated_volume_factor(shape))
    full = _ellipsoid_semiaxes(shape, volume)
    w = _smooth_step_meters(raw_fat)
    fat = w * raw_fat # smooth-clamped to ≥ 0
    axes = w .* (flesh .+ raw_fat) .+ (1 - w) .* full
    (; volume, dims = _extents(:skin, axes), fat, area = _shape_area(shape, axes))
end

function geometry(shape::Ellipsoid, fat_layer::FatLayer)
    s = _ellipsoid_fat_skin(shape, fat_layer)
    Geometry(s.volume, merge(s.dims, (; fat = s.fat)), SurfaceAreas(; total = s.area))
end
function geometry(shape::Ellipsoid, fibrous_layer::FibrousLayer, fat_layer::FatLayer)
    s = _ellipsoid_fat_skin(shape, fat_layer)
    fibrous = _fibrous_level(shape, s.dims, fibrous_layer.thickness)
    Geometry(s.volume, merge(s.dims, fibrous.dims, (; fat = s.fat)),
             SurfaceAreas(; total = fibrous.area, skin = s.area,
                          convection = _convective_area(fibrous_layer, s.area)))
end

# Fat thickness: the uniform X with π·k·(a+X)(b+X)(c+X) = volume, given the
# flesh semi-axes (k as in `_truncated_volume_factor`; a cut body keeps its cut
# fraction). The left side is strictly increasing for X > -min(a, b, c), so
# there is exactly one root; Newton from X = 0 with a fixed 10 steps reaches
# machine precision without branching (AD-friendly). The solve is unitless:
# lengths are stripped to metres once here and the result re-wrapped.
function _uniform_shell_thickness(flesh, volume, k)
    a, b, c = ustrip.(u"m", flesh)
    target = ustrip(u"m^3", volume) / (π * k)
    X = 0.0
    for _ in 1:10
        aX, bX, cX = a + X, b + X, c + X
        X -= (aX * bX * cX - target) / (bX * cX + aX * cX + aX * bX)
    end
    return max(0.0, X) * u"m"
end

# ── Silhouette ─────────────────────────────────────────────────────────────

# For an ellipsoid (a, b, c) the silhouette projected along a unit direction d
# is an ellipse of area π·√(b²c²·d_x² + a²c²·d_y² + a²b²·d_z²). As for the
# other shapes θ runs from the long axis up towards the vertical — the sun in
# the x–z plane, d = (cos θ, 0, sin θ) — so θ = π/2 looks down on the
# length × width footprint (`normal`) and θ = 0 looks end-on (`parallel`).
function silhouette(shape::Ellipsoid, a, b, c, θ)
    shape.pole_a_truncation == 0 && return π * b * sqrt(c^2 * cos(θ)^2 + a^2 * sin(θ)^2)
    _truncated_silhouette(a, b, c, 1 - shape.pole_a_truncation, θ)
end
function silhouette(sh::Ellipsoid, ins::AbstractInsulationLayer, body::AbstractBody)
    (; normal = silhouette(sh, ins, body, π / 2), parallel = silhouette(sh, ins, body, 0.0))
end

# A cut ellipsoid is convex, so its shadow is the convex hull of the shadows of
# its boundary's extreme points: the rim where the sun grazes the curved surface
# (kept only where it isn't cut away, x ≤ x_r·a) and the edge of the cut disc.
# Through x ↦ (x/a, y/b, z/c) the rim is the great circle ⟂ (d_x/a, d_y/b, d_z/c)
# on the unit sphere, since the surface normal at a·q is ∝ q ./ (a, b, c). Both curves
# are sampled densely, projected onto the plane ⟂ d and hulled; the hull area
# converges as 1/n² in the number of samples. Unitless: metres in, m² out.
function _truncated_silhouette(a, b, c, x_ratio, θ)
    am, bm, cm = ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c)
    d = (cos(θ), 0.0, sin(θ))
    u, v = (0.0, 1.0, 0.0), (-sin(θ), 0.0, cos(θ))      # basis of the plane ⟂ d
    e = (d[1] / am, d[2] / bm, d[3] / cm)
    ne = sqrt(e[1]^2 + e[2]^2 + e[3]^2)
    e = e ./ ne
    # Orthonormal pair ⟂ e, from whichever axis is least aligned with it.
    ref = abs(e[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
    w1 = ref .- (ref[1] * e[1] + ref[2] * e[2] + ref[3] * e[3]) .* e
    w1 = w1 ./ sqrt(w1[1]^2 + w1[2]^2 + w1[3]^2)
    w2 = (e[2] * w1[3] - e[3] * w1[2], e[3] * w1[1] - e[1] * w1[3], e[1] * w1[2] - e[2] * w1[1])
    onto(p) = (p[1] * u[1] + p[2] * u[2] + p[3] * u[3], p[1] * v[1] + p[2] * v[2] + p[3] * v[3])
    n = 720
    points = NTuple{2,Float64}[]
    xcut = x_ratio * am
    for t in range(0, 2π; length = n + 1)[1:n]
        q = cos(t) .* w1 .+ sin(t) .* w2
        p = (am * q[1], bm * q[2], cm * q[3])
        p[1] <= xcut && push!(points, onto(p))
    end
    scale = sqrt(max(0.0, 1 - x_ratio^2))
    for β in range(0, 2π; length = n + 1)[1:n]
        push!(points, onto((xcut, scale * bm * cos(β), scale * cm * sin(β))))
    end
    _hull_area(points) * u"m^2"
end

# Area of the convex hull of 2D points (Andrew's monotone chain + shoelace).
function _hull_area(points)
    pts = sort(unique(points))
    length(pts) < 3 && return 0.0
    cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
    hull = NTuple{2,Float64}[]
    for pass in (pts, reverse(pts))
        start = length(hull)
        for p in pass
            while length(hull) >= start + 2 && cross(hull[end - 1], hull[end], p) <= 0
                pop!(hull)
            end
            push!(hull, p)
        end
        pop!(hull)
    end
    area = 0.0
    for i in eachindex(hull)
        p, q = hull[i], hull[mod1(i + 1, length(hull))]
        area += p[1] * q[2] - q[1] * p[2]
    end
    abs(area) / 2
end
function silhouette(sh::Ellipsoid, ::AbstractInsulationLayer, body::AbstractBody, θ)
    silhouette(sh, _semiaxes(outer_dims(sh, body))..., θ)
end

# Radius accessors — shared by every ellipsoidal shape (`Ellipsoid`,
# `HalfEllipsoid`): the equivalent radius is the y semi-axis, width / 2.

_skin_radius(::AbstractEllipsoidal, length) = length.width_skin / 2
_fibrous_radius(::AbstractEllipsoidal, length) = length.width_fibrous / 2

# Composition

attachment_surfaces(::Ellipsoid) = (PoleA, PoleB, Equator)

# Outer (insulation-aware) full extents, as for `Plate`. Insulation-dispatched.
outer_dims(sh::Ellipsoid, body::AbstractBody) =
    outer_dims(sh, outer_insulation(insulation(body)), body)
outer_dims(::Ellipsoid, ::Union{Naked,FatLayer}, body::AbstractBody) =
    (length = body.geometry.length.length_skin,
     width = body.geometry.length.width_skin,
     height = body.geometry.length.height_skin)
outer_dims(::Ellipsoid, ::FibrousLayer, body::AbstractBody) =
    (length = body.geometry.length.length_fibrous,
     width = body.geometry.length.width_fibrous,
     height = body.geometry.length.height_fibrous)

# Skin-level semi-axes — used for flesh-anchored attachment positions.
_ellipsoid_skin(body::AbstractBody) = _skin_semiaxes(body.geometry.length)

# Notional area for pole attachments. For a truncated pole_a this is the
# actual flat-disc area exposed by the cut; for an untruncated pole it's
# the cross-sectional disc through the pole (used as a loose upper bound).
function surface_area(sh::Ellipsoid, body::AbstractBody, ::PoleA)
    (_, b, c) = _semiaxes(outer_dims(sh, body))
    s = _pole_a_radial_scale(sh)
    s == 0 ? π * b * c : π * (s * b) * (s * c)
end
function surface_area(sh::Ellipsoid, body::AbstractBody, ::PoleB)
    (_, b, c) = _semiaxes(outer_dims(sh, body))
    π * b * c
end

# Loose bound for equator joins: full ellipsoid surface area.
surface_area(::Ellipsoid, body::AbstractBody, ::Equator) =
    body.geometry.area.total

# For untruncated pole_a, coordinate fields are irrelevant (point pole).
# For truncated pole_a, PoleA may be bare (centre) or PoleA(r, φ) on the disc.
validate_range(::Ellipsoid, ::AbstractBody, ::PoleA) = nothing
validate_range(::Ellipsoid, ::AbstractBody, ::PoleB) = nothing
validate_range(::Ellipsoid, ::AbstractBody, ::Equator) = nothing

function surface_point(sh::Ellipsoid, body::AbstractBody, ::PoleA)
    a, _, _ = _ellipsoid_skin(body)
    (a * _pole_a_x_ratio(sh), zero(a), zero(a))
end
function surface_point(::Ellipsoid, body::AbstractBody, ::PoleB)
    a, _, _ = _ellipsoid_skin(body)
    (-a, zero(a), zero(a))
end
function surface_point(::Ellipsoid, body::AbstractBody, loc::Equator)
    _, b, c = _ellipsoid_skin(body)
    (zero(b), b * cos(loc.angle), c * sin(loc.angle))
end

surface_normal(::Ellipsoid, ::AbstractBody, ::PoleA) = ( 1.0, 0.0, 0.0)
surface_normal(::Ellipsoid, ::AbstractBody, ::PoleB) = (-1.0, 0.0, 0.0)
surface_normal(::Ellipsoid, ::AbstractBody, loc::Equator) =
    (0.0, cos(loc.angle), sin(loc.angle))

function surface_centroid(sh::Ellipsoid, body::AbstractBody, ::PoleA)
    a, _, _ = _ellipsoid_skin(body); (a * _pole_a_x_ratio(sh), zero(a), zero(a))
end
function surface_centroid(::Ellipsoid, body::AbstractBody, ::PoleB)
    a, _, _ = _ellipsoid_skin(body); (-a, zero(a), zero(a))
end
function surface_centroid(::Ellipsoid, body::AbstractBody, ::Equator)
    _, b, _ = _ellipsoid_skin(body); (zero(b), b, zero(b))  # arbitrary point on equator
end
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::PoleA) = ( 1.0, 0.0, 0.0)
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::PoleB) = (-1.0, 0.0, 0.0)
surface_centroid_normal(::Ellipsoid, ::AbstractBody, ::Equator) = (0.0, 1.0, 0.0)
