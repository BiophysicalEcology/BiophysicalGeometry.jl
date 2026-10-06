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
ellipsoid; `1` cuts through the centre. Surface area and volume are reported for
the full (untruncated) ellipsoid.
"""
struct Ellipsoid{M,D,B,C,T} <: AbstractEllipsoidal
    mass::M
    density::D
    axis_ratio_b::B
    axis_ratio_c::C
    pole_a_truncation::T
    Ellipsoid(::_Resolved, mass::M, density::D, axis_ratio_b::B, axis_ratio_c::C,
              pole_a_truncation::T) where {M,D,B,C,T} =
        new{M,D,B,C,T}(mass, density, axis_ratio_b, axis_ratio_c, pole_a_truncation)
end

# volume = (π/6)·length·width·height; ratios as for the box.
const _ELLIPSOID_SPEC = _ShapeSpec((:length, :width, :height), (1, 1, 1), log(π / 6),
                                   (:axis_ratio_b => (1, 2, 1.0), :axis_ratio_c => (1, 3, 1.0)))

function Ellipsoid(; pole_a_truncation = 0.0, kw...)
    0 <= pole_a_truncation <= 1 || throw(ArgumentError(
        "Ellipsoid `pole_a_truncation` must be in [0, 1], got $pole_a_truncation"))
    s = _resolve_shape("Ellipsoid", _ELLIPSOID_SPEC, NamedTuple(kw))
    Ellipsoid(_RESOLVED, s.mass, s.density, s.ratios..., pole_a_truncation)
end

# x-position of the truncated pole_a (as a fraction of a). 1.0 = full ellipsoid.
_pole_a_x_ratio(s::Ellipsoid) = 1 - s.pole_a_truncation
# Radial scale at the truncated pole (so cut disc has y/z extents = scale * b/c).
_pole_a_radial_scale(s::Ellipsoid) = sqrt(max(0.0, 1 - _pole_a_x_ratio(s)^2))

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

# ── Geometry ───────────────────────────────────────────────────────────────

# Canonicalise to m^3 to avoid weird unit ratios from mass/density (e.g. g/(kg/m^3))
# that propagate into cbrt() and trip up Enzyme's typeunstablerules path.
_ellipsoid_volume(shape::Ellipsoid) = uconvert(u"m^3", body_volume(shape))
_ellipsoid_fat_volume(shape::Ellipsoid, fat_layer::FatLayer) =
    uconvert(u"m^3", fat_volume(shape, fat_layer))

# Semi-axes enclosing `volume` at the shape's proportions: with a = r_b·b and
# a = r_c·c, volume = (4π/3)·a·b·c = (4π/3)·a³/(r_b·r_c).
function _ellipsoid_semiaxes(shape::Ellipsoid, volume)
    a = cbrt(3 * volume * shape.axis_ratio_b * shape.axis_ratio_c / (4π))
    (a, a / shape.axis_ratio_b, a / shape.axis_ratio_c)
end

function _skin_level(shape::Ellipsoid, volume)
    axes = _ellipsoid_semiaxes(shape, volume)
    (; dims = _extents(:skin, axes), area = _ellipsoid_area(axes...))
end
function _fibrous_level(shape::Ellipsoid, skin, thickness)
    axes = _skin_semiaxes(skin) .+ thickness
    (; dims = _extents(:fibrous, axes), area = _ellipsoid_area(axes...))
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
    raw_fat = _uniform_shell_thickness(flesh, volume)
    full = _ellipsoid_semiaxes(shape, volume)
    w = _smooth_step_meters(raw_fat)
    fat = w * raw_fat # smooth-clamped to ≥ 0
    axes = w .* (flesh .+ raw_fat) .+ (1 - w) .* full
    (; volume, dims = _extents(:skin, axes), fat, area = _ellipsoid_area(axes...))
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

# Fat thickness: the uniform X with (a+X)(b+X)(c+X) = 3·volume/(4π), given the
# flesh semi-axes. The left side is strictly increasing for X > -min(a, b, c),
# so there is exactly one root; Newton from X = 0 with a fixed 10 steps reaches
# machine precision without branching (AD-friendly). The solve is unitless:
# lengths are stripped to metres once here and the result re-wrapped.
function _uniform_shell_thickness(flesh, volume)
    a, b, c = ustrip.(u"m", flesh)
    target = 3 * ustrip(u"m^3", volume) / (4π)
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
silhouette(shape::Ellipsoid, a, b, c, θ) = π * b * sqrt(c^2 * cos(θ)^2 + a^2 * sin(θ)^2)
function silhouette(sh::Ellipsoid, ::AbstractInsulationLayer, body::AbstractBody)
    (a, b, c) = _semiaxes(outer_dims(sh, body))
    (; normal = π * a * b, parallel = π * b * c)
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
