"""
    AnimalBuilder

The recipe of the "Build an animal" page and app: builder settings in, Julia code out.

The code is the model. Every body, attachment and join is written as text and evaluated with the package, and
the numbers the recipe needs to place the next part (radii, lengths, poses, twists) are read back from the
evaluated objects. An animal is built by evaluating the finished text in a fresh module with
[`build_animal`](@ref), so the code shown with an animal is exactly what made it.
"""
module AnimalBuilder

using Unitful
using ..BiophysicalGeometry
using ..BiophysicalGeometry: Pose

const DEFAULTS = (;
    mass = 20.0, density = 1000.0, fatDensity = 901.0,
    posture = "Horizontal", pitch = 0.0,
    torsoShape = "Cylinder", torsoRatio = 3.0, fat = 0.1, backFur = 0.015, bellyFur = 0.005, limbFur = 0.008,
    headShape = "Ellipsoid", headRatio = 1.5, headFraction = 0.08,
    neck = false, neckFraction = 0.04, neckRatio = 1.0, neckPosture = "Forward",
    nose = "None", noseFraction = 0.005,
    ears = "None", earFraction = 0.002, earPosture = "Up", earRatio = 1.5, earFlatness = 10.0,
    legs = 4, legScaling = "Manual", legFraction = 0.03, legRatio = 5.0, legTop = 0.5,
    hindLegs = "Same", hindFraction = 0.06, hindRatio = 6.0,
    arms = false, armFraction = 0.05, armRatio = 12.0,
    wings = "None", wingFraction = 0.04,
    tail = false, tailFraction = 0.01, tailRatio = 6.0,
)

# A preset: the defaults, changed.
preset(; changes...) = merge(DEFAULTS, NamedTuple(changes))

const PRESETS = (
    "Dog" => preset(; neck = true, nose = "Nose", ears = "Cone", tail = true),
    "Mouse" => preset(; mass = 0.02, torsoShape = "Ellipsoid", torsoRatio = 2.0, fat = 0.05, backFur = 0.004,
        bellyFur = 0.003, headFraction = 0.12, legFraction = 0.02, legRatio = 4.0, legTop = 1.0,
        limbFur = 0.002, nose = "Nose", ears = "Cone", earFraction = 0.004, tail = true,
        tailFraction = 0.01, tailRatio = 12.0),
    "Elephant" => preset(; mass = 4000.0, torsoRatio = 1.8, fat = 0.0, backFur = 0.0, bellyFur = 0.0,
        limbFur = 0.0, headRatio = 1.3, headFraction = 0.08, nose = "Nose", noseFraction = 0.02,
        ears = "Plate", earFraction = 0.01, earPosture = "Up", earRatio = 1.2, earFlatness = 30.0,
        legFraction = 0.04, legRatio = 3.5, legTop = 0.8, tail = true, tailFraction = 0.002,
        tailRatio = 12.0),
    "Human" => preset(; mass = 70.0, density = 1050.0, fatDensity = 1050.0, posture = "Upright",
        torsoRatio = 1.9, fat = 0.252, backFur = 0.006, bellyFur = 0.006, limbFur = 0.006,
        headRatio = 1.6, headFraction = 0.0761, legs = 2, legFraction = 0.1623, legRatio = 7.0,
        legTop = 1.0, arms = true, armFraction = 0.0493, armRatio = 12.0),
    "Kangaroo" => preset(; mass = 50.0, pitch = 40.0, torsoRatio = 2.2, fat = 0.05, backFur = 0.01,
        bellyFur = 0.006, limbFur = 0.005, headRatio = 1.8, headFraction = 0.04, neck = true,
        neckFraction = 0.03, neckRatio = 1.5, ears = "Plate", earFraction = 0.001, earRatio = 2.5,
        earFlatness = 12.0, legFraction = 0.01, legRatio = 6.0, legTop = 0.5, hindLegs = "Different",
        hindFraction = 0.1, hindRatio = 4.5, tail = true, tailFraction = 0.08, tailRatio = 7.0),
    "Tyrannosaur" => preset(; mass = 7000.0, pitch = 5.0, torsoRatio = 2.2, fat = 0.0, backFur = 0.0,
        bellyFur = 0.0, limbFur = 0.0, headRatio = 1.8, headFraction = 0.07, neck = true,
        neckFraction = 0.04, neckRatio = 1.0, legFraction = 0.002, legRatio = 5.0, legTop = 0.5,
        hindLegs = "Different", hindFraction = 0.12, hindRatio = 4.0, tail = true, tailFraction = 0.12,
        tailRatio = 5.0),
    "Giraffe" => preset(; mass = 800.0, torsoRatio = 1.8, fat = 0.02, backFur = 0.003, bellyFur = 0.003,
        limbFur = 0.003, headRatio = 2.0, headFraction = 0.02, neck = true, neckFraction = 0.1,
        neckRatio = 6.0, neckPosture = "Up", ears = "Plate", earFraction = 0.0005, earRatio = 2.0,
        earFlatness = 12.0, legFraction = 0.05, legRatio = 10.0, legTop = 0.5, tail = true,
        tailFraction = 0.002, tailRatio = 15.0),
    "Cow" => preset(; mass = 682.0, torsoRatio = 2.4, fat = 0.0, backFur = 0.0034, bellyFur = 0.0034,
        headRatio = 1.8, headFraction = 0.04, legFraction = 0.02, legRatio = 4.3, legTop = 0.5,
        limbFur = 0.0034, neck = true, neckFraction = 0.05, ears = "Cone", earFraction = 0.001,
        tail = true, tailFraction = 0.003, tailRatio = 12.0),
    "Bird" => preset(; mass = 0.05, torsoShape = "Ellipsoid", torsoRatio = 1.6, fat = 0.05, backFur = 0.006,
        bellyFur = 0.006, headShape = "Sphere", headFraction = 0.1, legs = 2, legFraction = 0.02,
        legRatio = 8.0, legTop = 1.0, limbFur = 0.0, neck = true, neckFraction = 0.03, nose = "Beak",
        noseFraction = 0.01, tail = true, tailFraction = 0.02, tailRatio = 3.0, wings = "Folded",
        wingFraction = 0.06),
    "Seal" => preset(; nose = "Nose", tail = true, tailFraction = 0.02, tailRatio = 2.0, mass = 100.0,
        torsoShape = "Ellipsoid", torsoRatio = 4.0, fat = 0.35, backFur = 0.003, bellyFur = 0.003,
        headShape = "Sphere", headFraction = 0.05, legs = 0, limbFur = 0.003),
)

# ── Checking settings ──────────────────────────────────────────────────────────────────────────────────────────
#
# Every setting is checked before any code is written: a choice must be one of its options, which are the only
# strings written into the code, and a number must be finite and in range.

const CHOICES = (;
    posture = ("Horizontal", "Upright"), torsoShape = ("Cylinder", "Ellipsoid"),
    headShape = ("None", "Sphere", "Ellipsoid"), neckPosture = ("Forward", "Up"), nose = ("None", "Nose", "Beak"),
    ears = ("None", "Cone", "Plate"), earPosture = ("Up", "Flat"), legScaling = ("Manual", "Elastic", "Geometric"),
    hindLegs = ("Same", "Different"), wings = ("None", "Folded", "Spread"), legs = (0, 2, 4),
    neck = (false, true), arms = (false, true), tail = (false, true),
)
const POSITIVE = (:mass, :density, :fatDensity, :torsoRatio, :headRatio, :headFraction, :neckFraction, :neckRatio,
    :noseFraction, :earFraction, :earRatio, :earFlatness, :legFraction, :legRatio, :legTop, :hindFraction, :hindRatio,
    :armFraction, :armRatio, :wingFraction, :tailFraction, :tailRatio)
const NONNEGATIVE = (:fat, :backFur, :bellyFur, :limbFur)
const ANGLE = (:pitch,)

function check(p::NamedTuple)
    for key in keys(p)
        x = p[key]
        if haskey(CHOICES, key)
            x in CHOICES[key] || throw(ArgumentError("$key must be one of $(CHOICES[key]); got $(repr(x))"))
        elseif key in POSITIVE
            x isa Real && isfinite(x) && x > 0 || throw(ArgumentError("$key must be a finite positive number; got $(repr(x))"))
        elseif key in NONNEGATIVE
            x isa Real && isfinite(x) && x >= 0 || throw(ArgumentError("$key must be a finite number, zero or more; got $(repr(x))"))
        elseif key in ANGLE
            x isa Real && -90 <= x <= 90 || throw(ArgumentError("$key must be an angle from -90 to 90 degrees; got $(repr(x))"))
        else
            throw(ArgumentError("unknown setting $key"))
        end
    end
    p.fat < 1 || throw(ArgumentError("fat must be a fraction less than 1; got $(p.fat)"))
    p.legTop <= 1 || throw(ArgumentError("legTop must be at most 1; got $(p.legTop)"))
    return p
end

# ── Numbers as Julia text ──────────────────────────────────────────────────────────────────────────────────────

num(x, digits = 5) = string(round(Float64(x); sigdigits = digits))
metres(x) = "$(num(x))u\"m\""
kilograms(x) = "$(num(x, 6))u\"kg\""
# An angle, written exactly where it is a simple fraction of π, so that rounding cannot put it out of range.
function angle(x)
    for (value, text) in ((π, "π"), (-π, "-π"), (π / 2, "π / 2"), (-π / 2, "-π / 2"), (0.0, "0.0"))
        abs(x - value) < 1e-9 && return text
    end
    return num(x, 7)
end
fibres(depth) = depth > 0 ? "FibrousLayer($(num(depth * 1000))u\"mm\", 30.0u\"μm\", 3000u\"cm^-2\")" : "Naked()"
# A patch radius, rounded down to the digits written, so that the patch never comes out larger than its surface.
function disc(radius)
    isfinite(radius) && radius > 0 || throw(ArgumentError("a patch radius must be finite and positive; got $radius"))
    scale = 10.0^(4 - floor(Int, log10(radius)))
    "Disc($(metres(floor(radius * scale) / scale)))"
end

# ── Vectors ────────────────────────────────────────────────────────────────────────────────────────────────────

dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
unit3(v) = (n = sqrt(dot3(v, v)); (v[1] / n, v[2] / n, v[3] / n))
rotate3(R, v) = BiophysicalGeometry.apply_rotation(R, Tuple(v))
unrotate3(R, v) = BiophysicalGeometry.apply_rotation(transpose(R), Tuple(v))

# ── The recipe's working state ─────────────────────────────────────────────────────────────────────────────────

mutable struct Recipe
    mod::Module                              # where every line is evaluated
    definitions::Vector{String}              # body lines after the torso halves
    parts::Vector{Pair{Symbol,Symbol}}       # part name => the body it uses
    joins::Vector{String}
    poses::NamedTuple                        # part name => its pose
end

function Recipe()
    mod = Module(:AnimalBuilder)
    Core.eval(mod, :(using BiophysicalGeometry, Unitful))
    Recipe(mod, String[], Pair{Symbol,Symbol}[], String[], (;))
end

evaluate(r::Recipe, code::AbstractString) = Core.eval(r.mod, Meta.parseall(code))
# A value defined by the code. Read through `eval`, as the binding is newer than the running recipe.
binding(r::Recipe, name::Symbol) = Core.eval(r.mod, name)
function body(r::Recipe, name::Symbol)
    i = findfirst(p -> first(p) == name, r.parts)
    i === nothing && throw(ArgumentError("the recipe has no part $name"))
    return binding(r, last(r.parts[i]))
end

# Define a body under `name`; `show` adds the line to the definitions written after the torso.
function define!(r::Recipe, name::Symbol, code::AbstractString; show::Bool = true)
    line = "$name = $code"
    evaluate(r, line)
    show && push!(r.definitions, line)
    return binding(r, name)
end

# Skin dimensions of a body (m), as the package holds them: semi-axes and a radius for
# ellipsoids and spheres, radius and length for cylinders and cones, and length,
# width and height for plates.
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

attachment(r::Recipe, location, patch) = (text = "Attachment($location, $patch)"; (text, evaluate(r, text)))
normal_of(b, att) = BiophysicalGeometry._attach_normal(shape(b), b, att)

# The twist that turns the child's own direction `axis` as near as it can to the world direction `towards`; none if
# `towards` lies along the join, where any twist serves.
function twist_to(r::Recipe, parent::Symbol, on, child_body, at, axis, towards)
    target = .-rotate3(r.poses[parent].rotation, normal_of(body(r, parent), on))
    now = rotate3(BiophysicalGeometry.rotation_align(normal_of(child_body, at), target), axis)
    along = dot3(towards, target)
    across = Tuple(towards) .- along .* target
    dot3(across, across) > 1e-24 || return 0.0
    wanted = unit3(across)
    return atan(dot3(cross3(now, wanted), target), dot3(now, wanted))
end

# Join `child` (a new part using the body `alias`) to `parent`. `on` and `at` are surface texts; `patch` is the
# patch on the parent, `child_patch` on the child. `twist` may be a number, or a function of the two attachments.
function join!(r::Recipe, parent::Symbol, on, child::Symbol, alias::Symbol, at, patch;
               twist = 0.0, child_patch = patch)
    pb = body(r, parent)
    push!(r.parts, child => alias)
    cb = body(r, child)
    on_text, on_att = attachment(r, on, patch)
    at_text, at_att = attachment(r, at, child_patch)
    t = twist isa Function ? twist(on_att, cb, at_att) : twist
    twist_text = abs(t) > 1e-9 ? "; twist = $(num(t, 7))" : ""
    text = "Join($parent = $on_text, $child = $at_text$twist_text),"
    j = evaluate(r, text[1:end-1])
    pose = BiophysicalGeometry._child_pose(pb, r.poses[parent], cb, j)
    r.poses = merge(r.poses, NamedTuple{(child,)}((pose,)))
    push!(r.joins, text)
    return cb
end

# ── The animal ─────────────────────────────────────────────────────────────────────────────────────────────────

const IDENTITY = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
const LEG_SPLAY = 0.35
leg_names(n) = n == 4 ? [:leg_fl, :leg_fr, :leg_bl, :leg_br] : n == 2 ? [:leg_l, :leg_r] : Symbol[]
# (position along the torso, side) of each leg
leg_places(n) = n == 4 ? [(0.85, 1), (0.85, -1), (0.15, 1), (0.15, -1)] : n == 2 ? [(0.5, 1), (0.5, -1)] : Tuple{Float64,Int}[]

# A cylinder (`top` of 1) or a cone or frustum (`top` less than 1).
axial(mass, ratio, top) = top >= 1 ?
    "Body(Cylinder(; mass = $(kilograms(mass)), density, axis_ratio_b = $(num(ratio))), coat)" :
    "Body(Cone(; mass = $(kilograms(mass)), density, axis_ratio_b = $(num(ratio)), top_ratio = $(num(top))), coat)"
plate(mass, ratio, flatness) =
    "Body(Plate(; mass = $(kilograms(mass)), density, axis_ratio_b = $(num(ratio)), axis_ratio_c = $(num(flatness))), coat)"

"""
    animal_code(settings) -> String

The Julia code that builds the animal described by `settings`, a `NamedTuple` of the builder's controls; missing
ones are taken from `DEFAULTS`. Evaluating it defines `animal`, a `CompositeBody`.
"""
animal_code(settings) = first(recipe(settings))

function recipe(settings)
    p = check(merge(DEFAULTS, settings))
    upright = p.posture == "Upright"
    if upright   # a trunk standing on two legs, with arms: no tail or wings, and a round trunk
        p = merge(p, (; torsoShape = "Cylinder", legs = min(p.legs, 2), tail = false, wings = "None",
            neckPosture = "Forward", pitch = 0.0))
    end
    M = p.mass
    has_head = p.headShape != "None"
    mass = (;
        head = has_head ? p.headFraction * M : 0.0,
        neck = has_head && p.neck ? p.neckFraction * M : 0.0,
        nose = has_head && p.nose != "None" ? p.noseFraction * M : 0.0,
        ear = has_head && p.ears != "None" ? p.earFraction * M : 0.0,
        leg = p.legs > 0 ? p.legFraction * M : 0.0,
        tail = p.tail ? p.tailFraction * M : 0.0,
        wing = p.wings != "None" ? p.wingFraction * M : 0.0,
        arm = upright && p.arms ? p.armFraction * M : 0.0,
    )
    r = Recipe()
    evaluate(r, "density = $(num(p.density))u\"kg/m^3\"")
    evaluate(r, "coat = $(fibres(p.limbFur))")

    # Legs come first, as their mass is taken from the torso: by hand, or from the mass of the body.
    scaled = p.legs > 0 && p.legScaling != "Manual"
    leg_lines = String[]
    if scaled
        cone = p.legTop < 1
        append!(leg_lines, [
            "# legs from $(lowercase(p.legScaling)) similarity, with BiologicalScaling.jl",
            "import BiologicalScaling",
            "similarity = BiologicalScaling.$(p.legScaling)Similarity()",
            "length = BiologicalScaling.limb_length(similarity, $(kilograms(M)))",
            "radius = BiologicalScaling.limb_diameter(similarity, $(kilograms(M))) / 2",
            cone ? "leg = Body(Cone(; length, radius, density, top_ratio = $(num(p.legTop))), coat)" :
                   "leg = Body(Cylinder(; length, radius, density), coat)"])
        foreach(line -> evaluate(r, line), leg_lines)
        mass = merge(mass, (; leg = ustrip(u"kg", BiophysicalGeometry.mass(shape(binding(r, :leg))))))
    end
    # hind legs of their own size, for an animal with four legs sized by hand
    own_hind = p.legs == 4 && p.hindLegs == "Different" && p.legScaling == "Manual"
    mass = merge(mass, (; hind = own_hind ? p.hindFraction * M : mass.leg))
    leg_mass = p.legs == 4 ? 2 * mass.leg + 2 * mass.hind : p.legs * mass.leg
    torso_mass = M - mass.head - mass.neck - mass.nose - 2 * mass.ear - leg_mass - mass.tail -
        2 * mass.wing - 2 * mass.arm
    torso_mass > 0 || error("the parts weigh more than the animal: lower their fractions")
    cylinder = p.torsoShape == "Cylinder"

    # The torso: a dorsal and a ventral half, each with its own coat.
    function torso_half(fur)
        fat_layer = "FatLayer($(num(p.fat)), $(num(p.fatDensity))u\"kg/m^3\")"
        layers = p.fat > 0 ? (fur > 0 ? "CompositeInsulation($(fibres(fur)), $fat_layer)" : fat_layer) : fibres(fur)
        # a half ellipsoid's own height is half the full one's, so its length / height is twice
        cylinder ? "Body(HalfCylinder(; mass = $(kilograms(torso_mass / 2)), density, axis_ratio_b = $(num(p.torsoRatio))), $layers)" :
                   "Body(HalfEllipsoid(; mass = $(kilograms(torso_mass / 2)), density, axis_ratio_b = $(num(p.torsoRatio)), axis_ratio_c = $(num(2 * p.torsoRatio))), $layers)"
    end
    torso_lines = ["dorsal = $(torso_half(p.backFur))", "ventral = $(torso_half(p.bellyFur))"]
    foreach(line -> evaluate(r, line), torso_lines)
    push!(r.parts, :dorsal => :dorsal)
    dorsal = body(r, :dorsal)
    D = dims(dorsal)

    # Every torso half lies along x with its back (dome) up, tipped head-up by the pitch; or stands upright, its
    # length up and its back to +y. Standing up sends local y to world x, so a lateral angle 0 stays to the side.
    tilt = deg2rad(p.pitch)
    pitched = [cos(tilt) 0 -sin(tilt); 0 1 0; sin(tilt) 0 cos(tilt)]
    root = upright ? [0 1 0; 0 0 1; 1 0 0] : pitched
    r.poses = (; dorsal = Pose((0.0u"m", 0.0u"m", 0.0u"m"), root))

    # Where a part sits on a torso half: at its head end (`toward_head` true, or a position along it), its tail
    # end, or to one side.
    function on_torso(half, toward_head, side = 0)
        d = dims(half)
        upright && return "EndB(0.0u\"m\", π / 2)"
        if cylinder
            side == 0 || return "Lateral($(metres(toward_head * d.L)), $(angle(π / 2 + side * LEG_SPLAY)))"
            return toward_head == true ? "EndB($(metres(d.r / 2)), π / 2)" : "EndA($(metres(d.r / 2)), π / 2)"
        end
        side == 0 && return "Dome($(angle(toward_head == true ? 0.3 : π - 0.3)), π / 2)"
        α = toward_head > 0.6 ? 0.9 : toward_head < 0.4 ? π - 0.9 : π / 2
        return "Dome($(angle(α)), $(angle(π / 2 + side * 0.5)))"
    end
    # The place on a head that faces the world direction `towards`.
    function facing(towards)
        l = unrotate3(r.poses[:head].rotation, towards)
        p.headShape == "Sphere" ?
            "Radial($(angle(acos(clamp(l[3] / sqrt(dot3(l, l)), -1, 1)))), $(angle(atan(l[2], l[1]))))" :
            "Equator($(angle(atan(l[3], l[2]))))"
    end

    # The ventral half is turned about the join until its long axis lies along that of the dorsal half.
    long_axis = (1.0, 0.0, 0.0)
    join!(r, :dorsal, "Flat()", :ventral, :ventral, "Flat()", "FullCover()";
        twist = (on, cb, at) -> twist_to(r, :dorsal, on, cb, at, long_axis, rotate3(root, long_axis)))

    if has_head
        head_code = p.headShape == "Sphere" ? "Body(Sphere(; mass = $(kilograms(mass.head)), density), coat)" :
            "Body(Ellipsoid(; mass = $(kilograms(mass.head)), density, axis_ratio_b = $(num(p.headRatio)), axis_ratio_c = $(num(p.headRatio))), coat)"
        if p.neck
            define!(r, :neck, axial(mass.neck, p.neckRatio, 0.7))
        end
        head = define!(r, :head, head_code)
        H = dims(head)
        back = p.headShape == "Sphere" ? "Radial(π / 2, π)" : "PoleB()"
        if p.neck
            neck = binding(r, :neck)
            N = dims(neck)
            up = p.neckPosture == "Up"
            base = !up ? on_torso(dorsal, true) :
                cylinder ? "Lateral($(metres(0.9 * D.L)), π / 2)" : "Dome(0.6, π / 2)"
            join!(r, :dorsal, base, :neck, :neck, "EndA(0.0u\"m\", 0.0)", disc(min(N.r, 0.45 * D.r)))
            patch = disc(min(0.7 * N.r, 0.45 * H.r))
            if up && p.headShape == "Ellipsoid"   # the head sits across the top of the neck, facing forward
                join!(r, :neck, "EndB(0.0u\"m\", 0.0)", :head, :head, "Equator(-π / 2)", patch;
                    twist = (on, cb, at) -> twist_to(r, :neck, on, cb, at, (1.0, 0.0, 0.0), (1.0, 0.0, 0.0)))
            else
                join!(r, :neck, "EndB(0.0u\"m\", 0.0)", :head, :head, back, patch)
            end
        else
            join!(r, :dorsal, on_torso(dorsal, true), :head, :head, back, disc(0.45 * min(H.r, D.r)))
        end
        if p.nose != "None"
            name = p.nose == "Beak" ? :beak : :nose
            nose = define!(r, name, p.nose == "Beak" ? axial(mass.nose, 2.0, 0.0) : axial(mass.nose, 1.0, 0.5))
            front = upright ? facing((0.0, -1.0, 0.0)) : p.headShape == "Sphere" ? "Radial(π / 2, 0.0)" : "PoleA()"
            join!(r, :head, front, name, name, "EndA(0.0u\"m\", 0.0)", disc(min(dims(nose).r, 0.45 * H.r)))
        end
        if p.ears != "None"
            flat_ear = p.ears == "Plate"
            ear = define!(r, :ear, flat_ear ? plate(mass.ear, p.earRatio, p.earFlatness) : axial(mass.ear, 1.5, 0.3))
            E = dims(ear)
            forward = upright ? (0.0, -1.0, 0.0) : rotate3(r.poses[:head].rotation, (1.0, 0.0, 0.0))
            for (name, side) in ((:ear_l, 1), (:ear_r, -1))
                on = upright ? facing((side * 0.9, 0.0, 0.45)) :
                    p.headShape == "Sphere" ? "Radial(0.6, $(angle(side * π / 2)))" : "Equator($(angle(π / 2 - side * 0.6)))"
                if !flat_ear
                    join!(r, :head, on, name, :ear, "EndA(0.0u\"m\", 0.0)", disc(min(E.r, 0.45 * H.r)))
                elseif p.earPosture == "Up"   # standing on its edge, its flat side to the front
                    join!(r, :head, on, name, :ear, "SideB(0.0u\"m\", 0.0u\"m\")", disc(0.9 * sqrt(E.W * E.H / π));
                        twist = (o, cb, a) -> twist_to(r, :head, o, cb, a, (0.0, 0.0, 1.0), forward))
                else                            # laid back against the head, which hides one whole face
                    join!(r, :head, on, name, :ear, "Bottom()", "Disc(sqrt(surface_area(ear.shape, ear, Bottom()) / π))";
                        child_patch = "FullCover()",
                        twist = (o, cb, a) -> twist_to(r, :head, o, cb, a, (1.0, 0.0, 0.0), .-forward))
                end
            end
        end
    end
    if p.legs > 0
        scaled ? append!(r.definitions, leg_lines) : define!(r, :leg, axial(mass.leg, p.legRatio, p.legTop))
        own_hind && define!(r, :hind_leg, axial(mass.hind, p.hindRatio, p.legTop))
        for ((along, side), name) in zip(leg_places(p.legs), leg_names(p.legs))
            alias = own_hind && along < 0.4 ? :hind_leg : :leg
            radius = dims(binding(r, alias)).r
            if upright   # under the trunk, side by side
                join!(r, :dorsal, "EndA($(metres(1.1 * radius)), $(angle(side > 0 ? 0.0 : π)))", name, alias,
                    "EndA(0.0u\"m\", 0.0)", disc(radius))
            else
                join!(r, :ventral, on_torso(body(r, :ventral), along, side), name, alias, "EndA(0.0u\"m\", 0.0)", disc(radius))
            end
        end
    end
    if mass.arm > 0   # hung from the shoulders by their sides, and turned to point down
        A = dims(define!(r, :arm, axial(mass.arm, p.armRatio, 1.0)))
        for (name, side) in ((:arm_l, 1), (:arm_r, -1))
            join!(r, :dorsal, "Lateral($(metres(D.L - A.r)), $(angle(side > 0 ? 0.0 : π)))", name, :arm,
                "Lateral($(metres(A.r)), $(angle(side > 0 ? π : 0.0)))", disc(A.r);
                twist = (o, cb, a) -> twist_to(r, :dorsal, o, cb, a, (1.0, 0.0, 0.0), (0.0, 0.0, -1.0)))
        end
    end
    if p.wings != "None"
        W = dims(define!(r, :wing, plate(mass.wing, 2.5, 15.0)))
        for (name, side) in ((:wing_l, 1), (:wing_r, -1))
            on = cylinder ? "Lateral($(metres(0.6 * D.L)), $(angle(π / 2 - side * 1.1)))" :
                "Dome(1.2, $(angle(π / 2 - side * 1.1)))"
            if p.wings == "Folded"   # flat against the body, lying along it
                join!(r, :dorsal, on, name, :wing, "Bottom(0.0u\"m\", 0.0u\"m\")", disc(0.3 * W.W);
                    twist = (o, cb, a) -> twist_to(r, :dorsal, o, cb, a, (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)))
            else                       # out to the side, flat side up
                join!(r, :dorsal, on, name, :wing, "SideB(0.0u\"m\", 0.0u\"m\")", disc(0.9 * sqrt(W.W * W.H / π));
                    twist = (o, cb, a) -> twist_to(r, :dorsal, o, cb, a, (0.0, 0.0, 1.0), (0.0, 0.0, 1.0)))
            end
        end
    end
    if p.tail
        T = dims(define!(r, :tail, axial(mass.tail, p.tailRatio, 0.3)))
        join!(r, :dorsal, on_torso(dorsal, false), :tail, :tail, "EndA(0.0u\"m\", 0.0)", disc(min(T.r, 0.45 * D.r)))
    end

    # ── The code ──
    lines = ["using BiophysicalGeometry, Unitful", "", "density = $(num(p.density))u\"kg/m^3\""]
    isempty(r.definitions) || push!(lines, "coat = $(fibres(p.limbFur))")
    append!(lines, torso_lines)
    append!(lines, r.definitions)
    names = [name == alias ? string(name) : "$name = $alias" for (name, alias) in r.parts]
    push!(lines, "", "animal = CompositeBody(;", "    parts = (; $(join(names, ", "))),", "    joins = (")
    append!(lines, "        " .* r.joins)
    push!(lines, "    ),")
    if !(root ≈ IDENTITY)
        matrix = join((join((num(abs(x) < 1e-12 ? 0.0 : x, 9) for x in row), " ") for row in eachrow(root)), "; ")
        push!(lines, "    root_pose = Pose((0.0u\"m\", 0.0u\"m\", 0.0u\"m\"), [$matrix]),")
    end
    push!(lines, ")", "", "total_area(animal), skin_area(animal)")
    return join(lines, "\n"), p
end

"""
    build_animal(code) -> CompositeBody

Evaluate the builder's `code` in a fresh module and return its `animal`.
"""
function build_animal(code::AbstractString)
    mod = Module(:Animal)
    Core.eval(mod, Meta.parseall(code))
    return Core.eval(mod, :animal)
end

end
