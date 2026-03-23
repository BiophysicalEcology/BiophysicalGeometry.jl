using BiophysicalGeometry
using Unitful
using Test

import BiophysicalGeometry: Cylinder, Ellipsoid, Plate, Sphere

density = 1000.0u"kg/m^3"

# ── Fixture bodies ────────────────────────────────────────────────────────────
c_large = Body(Cylinder(10.0u"kg", density, 4.0), Naked())
c_small = Body(Cylinder( 2.0u"kg", density, 8.0), Naked())
ell     = Body(Ellipsoid(5.0u"kg", density, 1.6, 1.0), Naked())
plt     = Body(Plate(1.0u"kg", density, 3.0, 5.0), Naked())

# ── junction_radius ───────────────────────────────────────────────────────────
@testset "junction_radius" begin
    # Cylinder always returns cap (insulation) radius regardless of neck_fraction
    r = insulation_radius(c_large)
    @test junction_radius(c_large)       ≈ r
    @test junction_radius(c_large, 0.5)  ≈ r

    # Ellipsoid: neck_fraction=0 → full equatorial radius
    r_ell = insulation_radius(ell)
    @test junction_radius(ell, 0.0) ≈ r_ell
    # neck_fraction=0.5 → half of equatorial radius
    @test junction_radius(ell, 0.5) ≈ 0.5 * r_ell
end

# ── junction_area ─────────────────────────────────────────────────────────────
@testset "junction_area" begin
    j = Join(:a, :b)

    # Cylinder–Cylinder: π × min(r1, r2)²
    r_min = min(insulation_radius(c_large), insulation_radius(c_small))
    @test junction_area(c_large, c_small, j) ≈ π * r_min^2

    # Cylinder child attached to any parent: π × r_child²
    r_s = insulation_radius(c_small)
    @test junction_area(ell, c_small, j) ≈ π * r_s^2

    # Ellipsoid child with neck_fraction=0: π × min(r_parent, r_child)²
    j0 = Join(:a, :b, :top, 0.0)
    r_c = insulation_radius(c_large)
    r_e = insulation_radius(ell)
    @test junction_area(c_large, ell, j0) ≈ π * min(r_c, r_e)^2

    # Ellipsoid child with neck_fraction=0.5: π × (0.5 × r_ell)²
    j5 = Join(:a, :b, :top, 0.5)
    @test junction_area(c_large, ell, j5) ≈ π * (0.5 * r_e)^2

    # Plate child lateral (ear/wing): width × height
    j_lat = Join(:a, :b, :lateral)
    gl = plt.geometry.length
    @test junction_area(c_large, plt, j_lat) ≈ gl.width_skin * gl.height_skin

    # Plate child top: length × width
    j_top = Join(:a, :b, :top)
    @test junction_area(c_large, plt, j_top) ≈ gl.length_skin * gl.width_skin
end

# ── arm/leg fraction verification (NicheMapR PJOINs check) ───────────────────
@testset "arm and leg fraction vs NicheMapR" begin
    # Arm: b=12 cylinder.  Cap fraction = 1/(2×12+2) = 1/50 = 2.0 %
    arm  = Body(Cylinder(1.0u"kg", density, 12.0), Naked())
    sa   = arm.geometry.area.total
    r    = insulation_radius(arm)
    @test π * r^2 / sa ≈ 1.0/50.0  atol=1e-6

    # Leg: b=7 cylinder.  SA = 2πr²(2b+1), cap fraction = 1/(2(2b+1)) = 1/(4b+2) = 1/30 ≈ 3.333 %
    leg  = Body(Cylinder(1.0u"kg", density, 7.0), Naked())
    sa   = leg.geometry.area.total
    r    = insulation_radius(leg)
    @test π * r^2 / sa ≈ 1.0/30.0  atol=1e-6
end

# ── Human organism (area fractions sum to 1, total area > 0) ─────────────────
@testset "Human organism" begin
    mass = 70.0u"kg"
    human = Organism(
        parts = (
            head  = BodyPart(Body(Ellipsoid(mass * 0.076, density, 1.6, 1.0), Naked()), 1),
            torso = BodyPart(Body(Cylinder(mass * 0.501, density, 1.9),        Naked()), 1),
            arms  = BodyPart(Body(Cylinder(mass * 0.049, density, 12.0),       Naked()), 2),
            legs  = BodyPart(Body(Cylinder(mass * 0.162, density, 7.0),        Naked()), 2),
        ),
        joins = (
            Join(:torso, :head, :top,     0.5),
            Join(:torso, :arms, :lateral, 0.0),
            Join(:torso, :legs, :bottom,  0.0),
        )
    )

    geom = geometry(human)

    @test sum(values(geom.area_fractions)) ≈ 1.0  atol=1e-10
    @test geom.total_area > 0u"m^2"
    @test geom.total_volume > 0u"m^3"

    # All individual net areas must be positive
    for (name, a) in pairs(geom.area_fractions)
        @test a > 0
    end
end

# ── Validation: unknown part symbol raises error ──────────────────────────────
@testset "Organism validation" begin
    @test_throws ErrorException Organism(
        parts = (a = BodyPart(c_large, 1),),
        joins = (Join(:a, :b),)
    )
    @test_throws ErrorException Organism(
        parts = (a = BodyPart(c_large, 1),),
        joins = (Join(:b, :a),)
    )
    @test_throws ErrorException Organism(
        parts = (a = BodyPart(c_large, 1),),
        joins = (Join(:a, :a),)
    )
end

# ── total_volume and total_area accessors ─────────────────────────────────────
@testset "Organism accessors" begin
    org = Organism(
        parts = (p = BodyPart(c_large, 1), q = BodyPart(c_small, 2)),
        joins = (Join(:p, :q),)
    )
    geom = geometry(org)
    @test total_area(org)   ≈ geom.total_area
    @test total_volume(org) ≈ geom.total_volume
    # total_volume = sum of (volume × n) without junction correction
    expected_vol = c_large.geometry.volume * 1 + c_small.geometry.volume * 2
    @test total_volume(org) ≈ expected_vol
end

# ── part_positions: root at origin, Upright spine = z ────────────────────────
@testset "part_positions Upright" begin
    org = Organism(
        parts = (torso = BodyPart(c_large, 1), head = BodyPart(ell, 1)),
        joins = (Join(:torso, :head, :top, 0.5),)
    )
    pos = part_positions(org, Upright())

    # root (torso) at origin
    @test pos.torso.center == (0.0, 0.0, 0.0)
    @test pos.torso.axis   == (0.0, 0.0, 1.0)  # z-axis for Upright

    # head center above torso center (z > 0)
    @test pos.head.center[3] > 0.0
    @test pos.head.axis      == (0.0, 0.0, 1.0)
end

# ── part_positions: paired parts have positive side coordinate ────────────────
@testset "part_positions paired parts" begin
    org = Organism(
        parts = (torso = BodyPart(c_large, 1), arms = BodyPart(c_small, 2)),
        joins = (Join(:torso, :arms, :lateral),)
    )
    pos = part_positions(org, Upright())
    # Arms should be offset in the +x direction (Upright lateral = x)
    @test pos.arms.center[1] > 0.0
    @test pos.arms.center[2] ≈ 0.0  atol=1e-12
end
