"""
    AnimalBuilder

The animals of the "Build an animal" page, built from numbers alone: [`build`](@ref AnimalBuilder.build) makes one
of `ANIMALS` from its `settings`, with the package's machine-level `Unchecked` constructors, so that the page can run
it compiled to wasm.
"""
module AnimalBuilder

using Unitful
using ..BiophysicalGeometry
using ..BiophysicalGeometry: Pose, child_pose, rotation_align, @SMatrix

# ── The animals' settings ──────────────────────────────────────────────────────────────────────────────────────────
#
# Masses are fractions of the animal's, which is in kg; densities are in kg/m³, fur depths in m and angles in
# degrees. An animal's settings for parts it doesn't have are not used.

const DEFAULTS = (;
    mass = 20.0, density = 1000.0, fatDensity = 901.0, pitch = 0.0,
    torsoRatio = 3.0, fat = 0.1, backFur = 0.015, bellyFur = 0.005, limbFur = 0.008,
    headRatio = 1.5, headFraction = 0.08,
    neckFraction = 0.04, neckRatio = 1.0, neckAngle = 0.0,
    noseFraction = 0.005,
    earFraction = 0.002, earRatio = 1.5, earFlatness = 10.0, earAngle = 0.0,
    legFraction = 0.03, legRatio = 5.0, legTop = 0.5, legAngle = 0.0, legSpread = 0.0, hindFraction = 0.06, hindRatio = 6.0,
    armFraction = 0.05, armRatio = 12.0,
    wingFraction = 0.04, wingFold = 0.0,
    tailFraction = 0.01, tailRatio = 6.0, tailAngle = 0.0,
)

# The settings of the page's controls, in the order the page passes them.
const SETTING_NAMES = keys(DEFAULTS)

const PRESETS = (;
    Dog = (;),
    Mouse = (; mass = 0.02, torsoRatio = 2.0, fat = 0.05, backFur = 0.004, bellyFur = 0.003, headFraction = 0.12,
        legFraction = 0.02, legRatio = 4.0, legTop = 1.0, limbFur = 0.002, earFraction = 0.004, tailFraction = 0.01,
        tailRatio = 12.0),
    Elephant = (; mass = 4000.0, torsoRatio = 1.8, fat = 0.0, backFur = 0.0, bellyFur = 0.0, limbFur = 0.0,
        headRatio = 1.3, headFraction = 0.08, noseFraction = 0.02, earFraction = 0.01, earRatio = 1.2,
        earFlatness = 30.0, legFraction = 0.04, legRatio = 3.5, legTop = 0.8, tailFraction = 0.002, tailRatio = 12.0),
    Human = (; mass = 70.0, density = 1050.0, fatDensity = 1050.0, torsoRatio = 1.9, fat = 0.252, backFur = 0.006,
        bellyFur = 0.006, limbFur = 0.006, headRatio = 1.6, headFraction = 0.0761, legFraction = 0.1623,
        legRatio = 7.0, legTop = 1.0, armFraction = 0.0493, armRatio = 12.0),
    Kangaroo = (; mass = 50.0, pitch = 40.0, torsoRatio = 2.2, fat = 0.05, backFur = 0.01, bellyFur = 0.006,
        limbFur = 0.005, headRatio = 1.8, headFraction = 0.04, neckFraction = 0.03, neckRatio = 1.5,
        earFraction = 0.001, earRatio = 2.5, earFlatness = 12.0, legFraction = 0.01, legRatio = 6.0, legTop = 0.5,
        hindFraction = 0.1, hindRatio = 4.5, tailFraction = 0.08, tailRatio = 7.0),
    Tyrannosaur = (; mass = 7000.0, pitch = 5.0, torsoRatio = 2.2, fat = 0.0, backFur = 0.0, bellyFur = 0.0,
        limbFur = 0.0, headRatio = 1.8, headFraction = 0.07, neckFraction = 0.04, neckRatio = 1.0,
        legFraction = 0.002, legRatio = 5.0, legTop = 0.5, hindFraction = 0.12, hindRatio = 4.0, tailFraction = 0.12,
        tailRatio = 5.0),
    Giraffe = (; mass = 800.0, torsoRatio = 1.8, fat = 0.02, backFur = 0.003, bellyFur = 0.003, limbFur = 0.003,
        headRatio = 2.0, headFraction = 0.02, neckFraction = 0.1, neckRatio = 6.0, neckAngle = 70.0,
        earFraction = 0.0005, earRatio = 2.0, earFlatness = 12.0, legFraction = 0.05, legRatio = 10.0, legTop = 0.5,
        tailFraction = 0.002, tailRatio = 15.0),
    Cow = (; mass = 682.0, torsoRatio = 2.4, fat = 0.0, backFur = 0.0034, bellyFur = 0.0034, headRatio = 1.8,
        headFraction = 0.04, legFraction = 0.02, legRatio = 4.3, legTop = 0.5, limbFur = 0.0034, neckFraction = 0.05,
        earFraction = 0.001, tailFraction = 0.003, tailRatio = 12.0),
    Bird = (; mass = 0.05, torsoRatio = 1.6, fat = 0.05, backFur = 0.006, bellyFur = 0.006, headFraction = 0.1,
        legFraction = 0.02, legRatio = 8.0, legTop = 1.0, limbFur = 0.0, neckFraction = 0.03, noseFraction = 0.01,
        tailFraction = 0.02, tailRatio = 3.0, wingFraction = 0.06, wingFold = 80.0),
    Seal = (; mass = 100.0, torsoRatio = 4.0, fat = 0.35, backFur = 0.003, bellyFur = 0.003, headFraction = 0.05,
        limbFur = 0.003, tailFraction = 0.02, tailRatio = 2.0),
)

const ANIMAL_NAMES = map(string, keys(PRESETS))

"""
    settings(name) -> NamedTuple

The settings of the animal called `name`, one of `ANIMAL_NAMES`: the `DEFAULTS`, changed by its preset.
"""
settings(name::AbstractString) = map(Float64, merge(DEFAULTS, getfield(PRESETS, Symbol(name))))

# ── Vectors and dimensions ─────────────────────────────────────────────────────────────────────────────────────────

dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
unit3(v) = (n = sqrt(dot3(v, v)); (v[1] / n, v[2] / n, v[3] / n))
rotate3(R, v) = BiophysicalGeometry.apply_rotation(R, Tuple(v))

# Skin dimensions of a body (m), as the package holds them: semi-axes and a radius for ellipsoids and spheres,
# radius and length for cylinders and cones, and length, width and height for plates.
dims(b) = dims(shape(b), b.geometry.length)
function dims(::Union{BiophysicalGeometry.AbstractEllipsoidal,Half{<:BiophysicalGeometry.AbstractEllipsoidal}}, g)
    m(x) = ustrip(u"m", x)
    (a = m(g.length_skin) / 2, b = m(g.width_skin) / 2, r = m(g.width_skin) / 2)
end
function dims(::Union{BiophysicalGeometry.AbstractSpherical,Half{<:BiophysicalGeometry.AbstractSpherical}}, g)
    r = ustrip(u"m", g.radius_skin)
    (r = r, a = r, b = r)
end
dims(::Union{BiophysicalGeometry.AbstractCylindrical,Half{<:BiophysicalGeometry.AbstractCylindrical}}, g) =
    (r = ustrip(u"m", g.radius_skin), L = ustrip(u"m", g.length_skin))
dims(::BiophysicalGeometry.AbstractSlab, g) =
    (L = ustrip(u"m", g.length_skin), W = ustrip(u"m", g.width_skin), H = ustrip(u"m", g.height_skin))

normal_of(b, att) = BiophysicalGeometry.attach_normal(shape(b), b, att)

# How far either side of straight down legs are set on a cylindrical torso, in radians.
const LEG_SPLAY = 0.35

include("animals.jl")

end
