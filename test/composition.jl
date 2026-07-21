using BiophysicalGeometry
using Unitful
using Test

const ρ = 1000.0u"kg/m^3"

# Per-test part tags. Users define singleton types for each part they want
# to reference in a CompositeBody. There is no shared registry — each test
# defines what it needs.
struct Dorsal end
struct Ventral end
struct Torso end
struct Head end
struct LegFL end
struct LegFR end
struct LegBL end
struct LegBR end
struct OnlyPart end
struct SnowmanTop end
struct SnowmanBottom end
struct PartA end
struct PartB end
struct Big end
struct Small end
struct Frog end

@testset "HalfCylinder geometry" begin
    h = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    full = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    # Two halves of mass m sized the same way as one full of mass 2m.
    @test h.geometry.length.radius_skin ≈ full.geometry.length.radius_skin
    @test h.geometry.length.length_skin ≈ full.geometry.length.length_skin
    @test flesh_volume(h) ≈ flesh_volume(full) / 2

    # Total area = π r L (lateral) + π r² (two semicircular ends) + 2 r L (flat).
    r = h.geometry.length.radius_skin
    L = h.geometry.length.length_skin
    @test total_area(h) ≈ π*r*L + π*r^2 + 2*r*L
end

@testset "HalfEllipsoid geometry" begin
    h = Body(HalfEllipsoid(5u"kg", ρ, 1.5, 1.0), Naked())
    full = Body(Ellipsoid(10u"kg", ρ, 1.5, 1.0), Naked())
    @test h.geometry.length.b_semi_minor_skin ≈ full.geometry.length.b_semi_minor_skin
    @test h.geometry.length.a_semi_major_skin ≈ full.geometry.length.a_semi_major_skin
    @test flesh_volume(h) ≈ flesh_volume(full) / 2
    # Half dome + flat ≈ total
    flat = π * h.geometry.length.b_semi_minor_skin * h.geometry.length.c_semi_minor_skin
    @test total_area(h) ≈ total_area(full)/2 + flat
end

@testset "Dorsal/ventral torso == full cylinder" begin
    dorsal_body  = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    ventral_body = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    torso = CompositeBody(;
        parts = (Dorsal() => dorsal_body, Ventral() => ventral_body),
        joins = (
            Join(Dorsal(),  Attachment(Flat(), FullCover()),
                 Ventral(), Attachment(Flat(), FullCover())),
        ),
        root = Dorsal(),
    )
    full = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    @test total_area(torso) ≈ total_area(full)
    @test flesh_volume(torso) ≈ flesh_volume(full)
end

@testset "Dorsal/ventral torso == full ellipsoid" begin
    dorsal_body  = Body(HalfEllipsoid(5u"kg", ρ, 1.5, 1.0), Naked())
    ventral_body = Body(HalfEllipsoid(5u"kg", ρ, 1.5, 1.0), Naked())
    torso = CompositeBody(;
        parts = (Dorsal() => dorsal_body, Ventral() => ventral_body),
        joins = (
            Join(Dorsal(),  Attachment(Flat(), FullCover()),
                 Ventral(), Attachment(Flat(), FullCover())),
        ),
        root = Dorsal(),
    )
    full = Body(Ellipsoid(10u"kg", ρ, 1.5, 1.0), Naked())
    @test total_area(torso) ≈ total_area(full)
    @test flesh_volume(torso) ≈ flesh_volume(full)
end

@testset "Dorsal furred, ventral naked" begin
    fur = Fur(0.01u"m", 30u"μm", 3000u"cm^-2")
    dorsal_body  = Body(HalfCylinder(10u"kg", ρ, 3.0), fur)
    ventral_body = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    torso = CompositeBody(;
        parts = (Dorsal() => dorsal_body, Ventral() => ventral_body),
        joins = (
            Join(Dorsal(),  Attachment(Flat(), FullCover()),
                 Ventral(), Attachment(Flat(), FullCover())),
        ),
        root = Dorsal(),
    )
    # Sum of unjoined parts minus 2× the (matching) flat area.
    flat = 2 * dorsal_body.geometry.length.radius_skin * dorsal_body.geometry.length.length_skin
    @test total_area(torso) ≈ total_area(dorsal_body) + total_area(ventral_body) - 2*flat
end

@testset "Dog example" begin
    body = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    leg  = Body(Cylinder(1u"kg",  ρ, 5.0), Naked())
    head = Body(Ellipsoid(2u"kg", ρ, 1.5, 1.0), Naked())
    L_body = body.geometry.length.length_skin
    r_leg  = skin_radius(leg)
    r_head_attach = 0.04u"m"

    dog = CompositeBody(;
        parts = (
            Torso() => body,
            Head()  => head,
            LegFL() => leg,
            LegFR() => leg,
            LegBL() => leg,
            LegBR() => leg,
        ),
        joins = (
            Join(Torso(), Attachment(EndA(0.0u"m", 0.0),                Disc(r_head_attach)),
                 Head(),  Attachment(PoleA(),                            Disc(r_head_attach))),
            Join(Torso(), Attachment(Lateral(0.2*L_body, π/2 - π/6),    Disc(r_leg)),
                 LegFL(), Attachment(EndA(0.0u"m", 0.0),                 Disc(r_leg))),
            Join(Torso(), Attachment(Lateral(0.2*L_body, -π/2 + π/6),   Disc(r_leg)),
                 LegFR(), Attachment(EndA(0.0u"m", 0.0),                 Disc(r_leg))),
            Join(Torso(), Attachment(Lateral(0.8*L_body, π/2 - π/6),    Disc(r_leg)),
                 LegBL(), Attachment(EndA(0.0u"m", 0.0),                 Disc(r_leg))),
            Join(Torso(), Attachment(Lateral(0.8*L_body, -π/2 + π/6),   Disc(r_leg)),
                 LegBR(), Attachment(EndA(0.0u"m", 0.0),                 Disc(r_leg))),
        ),
        root = Torso(),
    )

    sum_parts = total_area(body) + total_area(head) + 4 * total_area(leg)
    sub = 2 * (π * r_head_attach^2 + 4 * π * r_leg^2)
    @test total_area(dog) ≈ sum_parts - sub
    @test flesh_volume(dog) ≈ flesh_volume(body) + flesh_volume(head) + 4 * flesh_volume(leg)

    # Pose: legs sit on the body lateral surface at the expected (z, φ).
    R_body = insulation_radius(body)
    expected_fl_x = R_body * cos(π/2 - π/6)
    expected_fl_y = R_body * sin(π/2 - π/6)
    expected_fl_z = 0.2 * L_body
    leg_fl_pose = BiophysicalGeometry._lookup(dog.poses, LegFL())
    @test leg_fl_pose.translation[1] ≈ expected_fl_x atol=1e-9*u"m"
    @test leg_fl_pose.translation[2] ≈ expected_fl_y atol=1e-9*u"m"

    # Skin radius / insulation_radius delegate to the root part.
    @test skin_radius(dog) == skin_radius(body)
    @test insulation_radius(dog) == insulation_radius(body)
end

@testset "Validation" begin
    a_body = Body(Cylinder(1u"kg", ρ, 2.0), Naked())
    b_body = Body(Cylinder(1u"kg", ρ, 2.0), Naked())

    # Patch radius too large for the surface.
    @test_throws ErrorException CompositeBody(;
        parts = (PartA() => a_body, PartB() => b_body),
        joins = (Join(PartA(), Attachment(EndA(0.0u"m", 0.0), Disc(10.0u"m")),
                      PartB(), Attachment(EndA(0.0u"m", 0.0), Disc(10.0u"m"))),),
        root = PartA(),
    )

    # Mismatched FullCover areas (different sizes).
    big_body = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    small_body = Body(HalfCylinder(1u"kg",  ρ, 3.0), Naked())
    @test_throws ErrorException CompositeBody(;
        parts = (Big() => big_body, Small() => small_body),
        joins = (Join(Big(),   Attachment(Flat(), FullCover()),
                      Small(), Attachment(Flat(), FullCover())),),
        root = Big(),
    )

    # Empirical shape (LeopardFrog) cannot be joined.
    frog_body = Body(LeopardFrog(0.04u"kg", ρ), Naked())
    @test_throws ErrorException CompositeBody(;
        parts = (Frog() => frog_body, PartB() => b_body),
        joins = (Join(Frog(),  Attachment(EndA(0.0u"m", 0.0), Disc(0.001u"m")),
                      PartB(), Attachment(EndA(0.0u"m", 0.0), Disc(0.001u"m"))),),
        root = Frog(),
    )
end

@testset "Single-part composite is identity" begin
    only_body = Body(Cylinder(2u"kg", ρ, 2.0), Naked())
    cb = CompositeBody(; parts = (OnlyPart() => only_body,), joins = (), root = OnlyPart())
    @test total_area(cb) ≈ total_area(only_body)
    @test flesh_volume(cb) ≈ flesh_volume(only_body)
    @test skin_radius(cb) == skin_radius(only_body)
end

@testset "silhouette_rasterized" begin
    # Sphere of radius r → silhouette is π·r² regardless of sun direction.
    sph = Body(Sphere(1u"kg", ρ), Naked())
    cb = CompositeBody(; parts = (OnlyPart() => sph,), joins = (), root = OnlyPart())
    R = skin_radius(sph)
    expected = π * R^2
    A1 = silhouette_rasterized(cb, (1.0, 0.0, 0.0); resolution=256)
    A2 = silhouette_rasterized(cb, (0.0, 1.0, 0.0); resolution=256)
    A3 = silhouette_rasterized(cb, (1.0, 1.0, 1.0); resolution=256)
    # Rasterisation accuracy ~ 1/resolution in each dimension; allow 1.5%.
    @test abs(A1 - expected) / expected < 0.015
    @test abs(A2 - expected) / expected < 0.015
    @test abs(A3 - expected) / expected < 0.015

    # Two equal spheres joined at the +z pole make a snowman along z.
    # Projected from +z (along the join axis) → discs overlap → ≈ π·R².
    # Projected from +x (across the axis) → side-by-side discs → ≈ 2·π·R².
    # Per-part summed silhouette is always 2·π·R² regardless of view, so
    # along-axis it overcounts by the overlap.
    s1 = Body(Sphere(1u"kg", ρ), Naked())
    s2 = Body(Sphere(1u"kg", ρ), Naked())
    snowman = CompositeBody(;
        parts = (SnowmanBottom() => s1, SnowmanTop() => s2),
        joins = (Join(SnowmanBottom(), Attachment(Radial(0.0, 0.0), Disc(0.001u"m")),
                      SnowmanTop(),    Attachment(Radial(0.0, 0.0), Disc(0.001u"m"))),),
        root = SnowmanBottom(),
    )
    Aalong  = silhouette_rasterized(snowman, (0.0, 0.0, 1.0); resolution=256)
    Across  = silhouette_rasterized(snowman, (1.0, 0.0, 0.0); resolution=256)
    summed  = silhouette_area(snowman).normal
    @test abs(Aalong - expected) / expected < 0.02
    @test abs(Across - 2*expected) / (2*expected) < 0.02
    @test summed ≈ 2 * expected            # summed counts both spheres
    @test Aalong < summed                  # overlap → rasteriser < summed
end
