# ── The animals, at machine level ──────────────────────────────────────────────────────────────────────────────────
#
# The animals the builder page starts from, built straight from numbers with the package's `Unchecked` constructors:
# no strings, no `eval` and nothing that throws, so that the page can run them compiled to wasm. An animal's type
# says what it is made of, through the methods below; its settings are all numbers, the page's sliders. Joints that
# rise, fold or lie back are `bend`s of their joins.

abstract type Animal end
struct Dog <: Animal end
struct Mouse <: Animal end
struct Elephant <: Animal end
struct Human <: Animal end
struct Kangaroo <: Animal end
struct Tyrannosaur <: Animal end
struct Giraffe <: Animal end
struct Cow <: Animal end
struct Bird <: Animal end
struct Seal <: Animal end

const ANIMALS = (Dog(), Mouse(), Elephant(), Human(), Kangaroo(), Tyrannosaur(), Giraffe(), Cow(), Bird(), Seal())

# What each animal is made of. Parts an animal lacks are `nothing`.
struct CylinderTorso end
struct EllipsoidTorso end
torso_kind(::Animal) = CylinderTorso()
torso_kind(::Union{Mouse,Bird,Seal}) = EllipsoidTorso()

struct EllipsoidHead end
struct SphereHead end
head_kind(::Animal) = EllipsoidHead()
head_kind(::Union{Bird,Seal}) = SphereHead()

struct Neck end
neck_kind(::Animal) = nothing
neck_kind(::Union{Dog,Kangaroo,Tyrannosaur,Giraffe,Cow,Bird}) = Neck()

struct Snout end
struct Beak end
nose_kind(::Animal) = nothing
nose_kind(::Union{Dog,Mouse,Elephant,Seal}) = Snout()
nose_kind(::Bird) = Beak()

struct ConeEars end
struct PlateEars end
ear_kind(::Animal) = nothing
ear_kind(::Union{Dog,Mouse,Cow}) = ConeEars()
ear_kind(::Union{Elephant,Kangaroo,Giraffe}) = PlateEars()

struct FourLegs end
struct TwoLegs end
struct UprightLegs end
leg_kind(::Animal) = FourLegs()
leg_kind(::Bird) = TwoLegs()
leg_kind(::Human) = UprightLegs()
leg_kind(::Seal) = nothing

# Hind legs of their own size.
struct OwnHind end
hind_kind(::Animal) = nothing
hind_kind(::Union{Kangaroo,Tyrannosaur}) = OwnHind()

struct Arms end
arm_kind(::Animal) = nothing
arm_kind(::Human) = Arms()

struct Wings end
wing_kind(::Animal) = nothing
wing_kind(::Bird) = Wings()

struct Tail end
tail_kind(::Animal) = Tail()
tail_kind(::Human) = nothing

struct Level end
struct Upright end
stance(::Animal) = Level()
stance(::Human) = Upright()

"""
    animal(name) -> Animal

The animal called `name`, one of `ANIMAL_NAMES`.
"""
animal(name::AbstractString) = ANIMALS[findfirst(==(name), ANIMAL_NAMES)]

# ── Building ───────────────────────────────────────────────────────────────────────────────────────────────────────

const FIBRE_DIAMETER = 30.0u"μm"
const FIBRE_DENSITY = 3000.0u"cm^-2"
coat(depth) = FibrousLayer(depth * u"m", FIBRE_DIAMETER, FIBRE_DENSITY)

# A cylinder is a cone with a top as wide as its base, so every axial part is a cone.
axial_body(s, mass, ratio, top) =
    Body(Cone(Unchecked(); mass = mass * u"kg", density = s.density * u"kg/m^3", axis_ratio_b = ratio, top_ratio = top),
         coat(s.limbFur))
plate_body(s, mass, ratio, flatness) =
    Body(Plate(Unchecked(); mass = mass * u"kg", density = s.density * u"kg/m^3", axis_ratio_b = ratio,
               axis_ratio_c = flatness), coat(s.limbFur))

torso_half(s, ::CylinderTorso, mass, fur) =
    Body(HalfCylinder(Unchecked(); mass = mass * u"kg", density = s.density * u"kg/m^3", axis_ratio_b = s.torsoRatio),
         torso_layers(s, fur))
# A half ellipsoid's own height is half the full one's, so its length / height is twice.
torso_half(s, ::EllipsoidTorso, mass, fur) =
    Body(HalfEllipsoid(Unchecked(); mass = mass * u"kg", density = s.density * u"kg/m^3", axis_ratio_b = s.torsoRatio,
                       axis_ratio_c = 2 * s.torsoRatio), torso_layers(s, fur))
torso_layers(s, fur) = CompositeInsulation(coat(fur), FatLayer(s.fat, s.fatDensity * u"kg/m^3"))

head_body(s, ::EllipsoidHead, mass) =
    Body(Ellipsoid(Unchecked(); mass = mass * u"kg", density = s.density * u"kg/m^3", axis_ratio_b = s.headRatio,
                   axis_ratio_c = s.headRatio), coat(s.limbFur))
head_body(s, ::SphereHead, mass) =
    Body(Sphere(Unchecked(); mass = mass * u"kg", density = s.density * u"kg/m^3"), coat(s.limbFur))

# The masses of the parts, from the animal's: what its parts don't take is the torso's.
function part_masses(a::Animal, s)
    M = s.mass
    head = s.headFraction * M
    neck = _has(neck_kind(a)) * s.neckFraction * M
    nose = _has(nose_kind(a)) * s.noseFraction * M
    ear = _has(ear_kind(a)) * s.earFraction * M
    leg = _has(leg_kind(a)) * s.legFraction * M
    hind = _has(hind_kind(a)) ? s.hindFraction * M : leg
    arm = _has(arm_kind(a)) * s.armFraction * M
    wing = _has(wing_kind(a)) * s.wingFraction * M
    tail = _has(tail_kind(a)) * s.tailFraction * M
    legs = _leg_count(leg_kind(a)) == 4 ? 2 * leg + 2 * hind : _leg_count(leg_kind(a)) * leg
    torso = M - head - neck - nose - 2 * ear - legs - tail - 2 * wing - 2 * arm
    (; torso, head, neck, nose, ear, leg, hind, arm, wing, tail)
end
_has(::Nothing) = false
_has(_) = true
_leg_count(::FourLegs) = 4
_leg_count(::Union{TwoLegs,UprightLegs}) = 2
_leg_count(::Nothing) = 0

# The animal so far: its parts, joins and the poses of its parts.
struct Assembly{P<:NamedTuple,J<:Tuple,Q<:NamedTuple}
    parts::P
    joins::J
    poses::Q
end

# Join the new part `C`, `child`, to the part `P` already in place.
function attach(a::Assembly, ::Val{P}, ::Val{C}, child, on, at;
                twist = 0.0, bend = 0.0, hinge = (0.0, 0.0, 0.0)) where {P,C}
    parent = getfield(a.parts, P)
    j = Join{P,C,typeof(on),typeof(at),Float64}(on, at, twist, bend, hinge)
    pose = child_pose(parent, getfield(a.poses, P), child, j)
    Assembly(merge(a.parts, NamedTuple{(C,)}((child,))), (a.joins..., j), merge(a.poses, NamedTuple{(C,)}((pose,))))
end

# The twist that turns the child's own direction `axis` as near as it can to the world direction `towards`, on a
# join of `on` on the part `P` to `at` on `child`; none if `towards` lies along the join, where any twist serves.
function twist_to(a::Assembly, ::Val{P}, on, child, at, axis, towards) where {P}
    parent = getfield(a.parts, P)
    target = .-rotate3(getfield(a.poses, P).rotation, normal_of(parent, on))
    now = rotate3(rotation_align(normal_of(child, at), target), axis)
    along = dot3(towards, target)
    across = towards .- along .* target
    wanted = unit3(across)
    ifelse(dot3(across, across) > 1e-24, atan(dot3(cross3(now, wanted), target), dot3(now, wanted)), 0.0)
end

m(x) = x * u"m"
patch(radius) = Disc(m(radius))
const ZERO = 0.0u"m"

# Where parts go on a torso half: its head end, its tail end, or under it to one side.
front(::CylinderTorso, d) = EndB(m(d.r / 2), π / 2)
front(::EllipsoidTorso, d) = Dome(0.3, π / 2)
back(::CylinderTorso, d) = EndA(m(d.r / 2), π / 2)
back(::EllipsoidTorso, d) = Dome(π - 0.3, π / 2)
side(::CylinderTorso, d, along, side) = Lateral(m(along * d.L), π / 2 + side * LEG_SPLAY)
side(::EllipsoidTorso, d, along, side) =
    Dome(along > 0.6 ? 0.9 : along < 0.4 ? π - 0.9 : π / 2, π / 2 + side * 0.5)

# A head's front and back, and the places of its ears.
head_front(::EllipsoidHead) = PoleA()
head_front(::SphereHead) = Radial(π / 2, 0.0)
head_back(::EllipsoidHead) = PoleB()
head_back(::SphereHead) = Radial(π / 2, π)
ear_place(::EllipsoidHead, side) = Equator(π / 2 - side * 0.6)
ear_place(::SphereHead, side) = Radial(0.6, side * π / 2)

root_rotation(::Level, s) = (t = deg2rad(s.pitch); @SMatrix [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)])
root_rotation(::Upright, s) = @SMatrix [0.0 1 0; 0 0 1; 1 0 0]

"""
    build(animal, settings) -> CompositeBody

Build `animal`, one of `ANIMALS`, from `settings`, a `NamedTuple` with the `SETTING_NAMES`, with machine-level
constructors only: nothing is checked, and nothing can throw.
"""
function build(an::Animal, s)
    masses = part_masses(an, s)
    kind = torso_kind(an)
    dorsal = torso_half(s, kind, masses.torso / 2, s.backFur)
    ventral = torso_half(s, kind, masses.torso / 2, s.bellyFur)
    root = root_rotation(stance(an), s)
    a = Assembly((; dorsal), (), (; dorsal = Pose((ZERO, ZERO, ZERO), root)))
    # The ventral half is turned about the join until its long axis lies along that of the dorsal half.
    x = (1.0, 0.0, 0.0)
    on, at = Attachment(Flat(), FullCover()), Attachment(Flat(), FullCover())
    a = attach(a, Val(:dorsal), Val(:ventral), ventral, on, at;
               twist = twist_to(a, Val(:dorsal), on, ventral, at, x, rotate3(root, x)))
    a = add_head(a, an, s, masses, neck_kind(an))
    a = add_nose(a, an, s, masses, nose_kind(an))
    a = add_ears(a, an, s, masses, ear_kind(an))
    a = add_legs(a, an, s, masses, leg_kind(an))
    a = add_arms(a, an, s, masses, arm_kind(an))
    a = add_wings(a, an, s, masses, wing_kind(an))
    a = add_tail(a, an, s, masses, tail_kind(an))
    CompositeBody(Unchecked(); parts = a.parts, joins = a.joins, root_pose = a.poses.dorsal)
end

# The head, straight on the torso or on a neck, which rises by `neckAngle`; the head bends back by as much, so it
# stays level.
function add_head(a, an, s, masses, ::Nothing)
    head = head_body(s, head_kind(an), masses.head)
    D, H = dims(a.parts.dorsal), dims(head)
    on = stance(an) isa Upright ? EndB(ZERO, π / 2) : front(torso_kind(an), D)
    p = patch(0.45 * min(H.r, D.r))
    attach(a, Val(:dorsal), Val(:head), head, Attachment(on, p), Attachment(head_back(head_kind(an)), p))
end
function add_head(a, an, s, masses, ::Neck)
    neck = axial_body(s, masses.neck, s.neckRatio, 0.7)
    head = head_body(s, head_kind(an), masses.head)
    D, N, H = dims(a.parts.dorsal), dims(neck), dims(head)
    angle = deg2rad(s.neckAngle)
    p = patch(min(N.r, 0.45 * D.r))
    a = attach(a, Val(:dorsal), Val(:neck), neck, Attachment(front(torso_kind(an), D), p),
               Attachment(EndA(ZERO, 0.0), p); bend = -angle, hinge = (0.0, 1.0, 0.0))
    p = patch(min(0.7 * N.r, 0.45 * H.r))
    attach(a, Val(:neck), Val(:head), head, Attachment(EndB(ZERO, 0.0), p),
           Attachment(head_back(head_kind(an)), p); bend = angle, hinge = (0.0, 1.0, 0.0))
end

add_nose(a, an, s, masses, ::Nothing) = a
add_nose(a, an, s, masses, ::Snout) = _add_nose(a, an, Val(:nose), axial_body(s, masses.nose, 1.0, 0.5))
add_nose(a, an, s, masses, ::Beak) = _add_nose(a, an, Val(:beak), axial_body(s, masses.nose, 2.0, 0.0))
function _add_nose(a, an, name, nose)
    p = patch(min(dims(nose).r, 0.45 * dims(a.parts.head).r))
    attach(a, Val(:head), name, nose, Attachment(head_front(head_kind(an)), p), Attachment(EndA(ZERO, 0.0), p))
end

# Ears stand up from the head, and lie back by `earAngle`. A plate ear stands on its edge, flat side to the front.
add_ears(a, an, s, masses, ::Nothing) = a
function add_ears(a, an, s, masses, ::ConeEars)
    ear = axial_body(s, masses.ear, 1.5, 0.3)
    p = patch(min(dims(ear).r, 0.45 * dims(a.parts.head).r))
    lie = -deg2rad(s.earAngle)
    a = attach(a, Val(:head), Val(:ear_l), ear, Attachment(ear_place(head_kind(an), 1), p),
               Attachment(EndA(ZERO, 0.0), p); bend = lie, hinge = (0.0, 1.0, 0.0))
    attach(a, Val(:head), Val(:ear_r), ear, Attachment(ear_place(head_kind(an), -1), p),
           Attachment(EndA(ZERO, 0.0), p); bend = lie, hinge = (0.0, 1.0, 0.0))
end
function add_ears(a, an, s, masses, ::PlateEars)
    ear = plate_body(s, masses.ear, s.earRatio, s.earFlatness)
    E = dims(ear)
    p = patch(0.9 * sqrt(E.W * E.H / π))
    forward = rotate3(a.poses.head.rotation, (1.0, 0.0, 0.0))
    lie = -deg2rad(s.earAngle)
    a = _add_plate_ear(a, an, Val(:ear_l), ear, 1, p, forward, lie)
    _add_plate_ear(a, an, Val(:ear_r), ear, -1, p, forward, lie)
end
function _add_plate_ear(a, an, name, ear, side, p, forward, lie)
    on, at = Attachment(ear_place(head_kind(an), side), p), Attachment(SideB(ZERO, ZERO), p)
    attach(a, Val(:head), name, ear, on, at; twist = twist_to(a, Val(:head), on, ear, at, (0.0, 0.0, 1.0), forward),
           bend = lie, hinge = (0.0, 1.0, 0.0))
end

add_legs(a, an, s, masses, ::Nothing) = a
function add_legs(a, an, s, masses, ::FourLegs)
    leg = axial_body(s, masses.leg, s.legRatio, s.legTop)
    hind = _hind_leg(s, masses, leg, hind_kind(an))
    a = _add_leg(a, an, Val(:leg_fl), leg, 0.85, 1)
    a = _add_leg(a, an, Val(:leg_fr), leg, 0.85, -1)
    a = _add_leg(a, an, Val(:leg_bl), hind, 0.15, 1)
    _add_leg(a, an, Val(:leg_br), hind, 0.15, -1)
end
function add_legs(a, an, s, masses, ::TwoLegs)
    leg = axial_body(s, masses.leg, s.legRatio, s.legTop)
    a = _add_leg(a, an, Val(:leg_l), leg, 0.5, 1)
    _add_leg(a, an, Val(:leg_r), leg, 0.5, -1)
end
# Standing: under the trunk, side by side.
function add_legs(a, an, s, masses, ::UprightLegs)
    leg = axial_body(s, masses.leg, s.legRatio, s.legTop)
    r = dims(leg).r
    p = patch(r)
    a = attach(a, Val(:dorsal), Val(:leg_l), leg, Attachment(EndA(m(1.1 * r), 0.0), p), Attachment(EndA(ZERO, 0.0), p))
    attach(a, Val(:dorsal), Val(:leg_r), leg, Attachment(EndA(m(1.1 * r), π), p), Attachment(EndA(ZERO, 0.0), p))
end
_hind_leg(s, masses, leg, ::Nothing) = leg
_hind_leg(s, masses, leg, ::OwnHind) = axial_body(s, masses.hind, s.hindRatio, s.legTop)
function _add_leg(a, an, name, leg, along, side_)
    p = patch(dims(leg).r)
    on = side(torso_kind(an), dims(a.parts.ventral), along, side_)
    attach(a, Val(:ventral), name, leg, Attachment(on, p), Attachment(EndA(ZERO, 0.0), p))
end

# Arms hang from the shoulders by their sides, turned to point down.
add_arms(a, an, s, masses, ::Nothing) = a
function add_arms(a, an, s, masses, ::Arms)
    arm = axial_body(s, masses.arm, s.armRatio, 1.0)
    D, A = dims(a.parts.dorsal), dims(arm)
    p = patch(A.r)
    a = _add_arm(a, Val(:arm_l), arm, Attachment(Lateral(m(D.L - A.r), 0.0), p), Attachment(Lateral(m(A.r), π), p))
    _add_arm(a, Val(:arm_r), arm, Attachment(Lateral(m(D.L - A.r), π), p), Attachment(Lateral(m(A.r), 0.0), p))
end
_add_arm(a, name, arm, on, at) =
    attach(a, Val(:dorsal), name, arm, on, at; twist = twist_to(a, Val(:dorsal), on, arm, at, (1.0, 0.0, 0.0), (0.0, 0.0, -1.0)))

# Wings stand out to the sides, flat side up, and fold back along the body by `wingFold`.
add_wings(a, an, s, masses, ::Nothing) = a
function add_wings(a, an, s, masses, ::Wings)
    wing = plate_body(s, masses.wing, 2.5, 15.0)
    W = dims(wing)
    p = patch(0.9 * sqrt(W.W * W.H / π))
    fold = deg2rad(s.wingFold)
    a = _add_wing(a, an, Val(:wing_l), wing, 1, p, fold)
    _add_wing(a, an, Val(:wing_r), wing, -1, p, fold)
end
function _add_wing(a, an, name, wing, side_, p, fold)
    D = dims(a.parts.dorsal)
    on = Attachment(_wing_place(torso_kind(an), D, side_), p)
    at = Attachment(SideB(ZERO, ZERO), p)
    attach(a, Val(:dorsal), name, wing, on, at;
           twist = twist_to(a, Val(:dorsal), on, wing, at, (0.0, 0.0, 1.0), (0.0, 0.0, 1.0)),
           bend = side_ * fold, hinge = (0.0, 0.0, 1.0))
end
_wing_place(::CylinderTorso, D, side) = Lateral(m(0.6 * D.L), π / 2 - side * 1.1)
_wing_place(::EllipsoidTorso, D, side) = Dome(1.2, π / 2 - side * 1.1)

add_tail(a, an, s, masses, ::Nothing) = a
function add_tail(a, an, s, masses, ::Tail)
    tail = axial_body(s, masses.tail, s.tailRatio, 0.3)
    D = dims(a.parts.dorsal)
    p = patch(min(dims(tail).r, 0.45 * D.r))
    attach(a, Val(:dorsal), Val(:tail), tail, Attachment(back(torso_kind(an), D), p), Attachment(EndA(ZERO, 0.0), p))
end
