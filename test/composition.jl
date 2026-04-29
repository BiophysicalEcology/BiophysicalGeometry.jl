using BiophysicalGeometry
using Unitful
using Test

const ρ = 1000.0u"kg/m^3"

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
    dorsal  = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    ventral = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    torso = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (
            Join(:dorsal,  Attachment(:flat, (;), FullCover()),
                 :ventral, Attachment(:flat, (;), FullCover())),
        ),
        root = :dorsal,
    )
    full = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    @test total_area(torso) ≈ total_area(full)
    @test flesh_volume(torso) ≈ flesh_volume(full)
end

@testset "Dorsal/ventral torso == full ellipsoid" begin
    dorsal  = Body(HalfEllipsoid(5u"kg", ρ, 1.5, 1.0), Naked())
    ventral = Body(HalfEllipsoid(5u"kg", ρ, 1.5, 1.0), Naked())
    torso = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (
            Join(:dorsal,  Attachment(:flat, (;), FullCover()),
                 :ventral, Attachment(:flat, (;), FullCover())),
        ),
        root = :dorsal,
    )
    full = Body(Ellipsoid(10u"kg", ρ, 1.5, 1.0), Naked())
    @test total_area(torso) ≈ total_area(full)
    @test flesh_volume(torso) ≈ flesh_volume(full)
end

@testset "Dorsal furred, ventral naked" begin
    fur = Fur(0.01u"m", 30u"μm", 3000u"cm^-2")
    dorsal  = Body(HalfCylinder(10u"kg", ρ, 3.0), fur)
    ventral = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    torso = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (
            Join(:dorsal,  Attachment(:flat, (;), FullCover()),
                 :ventral, Attachment(:flat, (;), FullCover())),
        ),
        root = :dorsal,
    )
    # Sum of unjoined parts minus 2× the (matching) flat area.
    flat = 2 * dorsal.geometry.length.radius_skin * dorsal.geometry.length.length_skin
    @test total_area(torso) ≈ total_area(dorsal) + total_area(ventral) - 2*flat
end

@testset "Dog example" begin
    body = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    leg  = Body(Cylinder(1u"kg",  ρ, 5.0), Naked())
    head = Body(Ellipsoid(2u"kg", ρ, 1.5, 1.0), Naked())
    L_body = body.geometry.length.length_skin
    r_leg  = skin_radius(leg)
    r_head_attach = 0.04u"m"

    dog = CompositeBody(;
        parts = (; body, head, leg_fl=leg, leg_fr=leg, leg_bl=leg, leg_br=leg),
        joins = (
            Join(:body, Attachment(:end_a,   (r=0.0u"m", φ=0.0),         Disc(r_head_attach)),
                 :head, Attachment(:pole_a,  (;),                         Disc(r_head_attach))),
            Join(:body, Attachment(:lateral, (z=0.2*L_body, φ=π/2 - π/6), Disc(r_leg)),
                 :leg_fl, Attachment(:end_a, (r=0.0u"m", φ=0.0),          Disc(r_leg))),
            Join(:body, Attachment(:lateral, (z=0.2*L_body, φ=-π/2 + π/6), Disc(r_leg)),
                 :leg_fr, Attachment(:end_a, (r=0.0u"m", φ=0.0),          Disc(r_leg))),
            Join(:body, Attachment(:lateral, (z=0.8*L_body, φ=π/2 - π/6), Disc(r_leg)),
                 :leg_bl, Attachment(:end_a, (r=0.0u"m", φ=0.0),          Disc(r_leg))),
            Join(:body, Attachment(:lateral, (z=0.8*L_body, φ=-π/2 + π/6), Disc(r_leg)),
                 :leg_br, Attachment(:end_a, (r=0.0u"m", φ=0.0),          Disc(r_leg))),
        ),
        root = :body,
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
    @test dog.poses.leg_fl.translation[1] ≈ expected_fl_x atol=1e-9*u"m"
    @test dog.poses.leg_fl.translation[2] ≈ expected_fl_y atol=1e-9*u"m"

    # Skin radius/insulation_radius delegate to the root part.
    @test skin_radius(dog) == skin_radius(body)
    @test insulation_radius(dog) == insulation_radius(body)
end

@testset "Validation" begin
    a = Body(Cylinder(1u"kg", ρ, 2.0), Naked())
    b = Body(Cylinder(1u"kg", ρ, 2.0), Naked())

    # Unknown surface name.
    @test_throws ErrorException CompositeBody(;
        parts = (; a, b),
        joins = (Join(:a, Attachment(:bogus, (z=0.0u"m", φ=0.0), Disc(0.001u"m")),
                      :b, Attachment(:end_a, (r=0.0u"m", φ=0.0), Disc(0.001u"m"))),),
        root = :a,
    )

    # Position with wrong field set.
    @test_throws ErrorException CompositeBody(;
        parts = (; a, b),
        joins = (Join(:a, Attachment(:lateral, (foo=0.0u"m",), Disc(0.001u"m")),
                      :b, Attachment(:end_a, (r=0.0u"m", φ=0.0), Disc(0.001u"m"))),),
        root = :a,
    )

    # Patch radius too large for the surface.
    @test_throws ErrorException CompositeBody(;
        parts = (; a, b),
        joins = (Join(:a, Attachment(:end_a, (r=0.0u"m", φ=0.0), Disc(10.0u"m")),
                      :b, Attachment(:end_a, (r=0.0u"m", φ=0.0), Disc(10.0u"m"))),),
        root = :a,
    )

    # Mismatched FullCover areas (different sizes).
    big = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    small = Body(HalfCylinder(1u"kg",  ρ, 3.0), Naked())
    @test_throws ErrorException CompositeBody(;
        parts = (; big, small),
        joins = (Join(:big,   Attachment(:flat, (;), FullCover()),
                      :small, Attachment(:flat, (;), FullCover())),),
        root = :big,
    )

    # Empirical shape (LeopardFrog) cannot be joined.
    frog = Body(LeopardFrog(0.04u"kg", ρ), Naked())
    @test_throws ErrorException CompositeBody(;
        parts = (; frog, b),
        joins = (Join(:frog, Attachment(:end_a, (r=0.0u"m", φ=0.0), Disc(0.001u"m")),
                      :b,    Attachment(:end_a, (r=0.0u"m", φ=0.0), Disc(0.001u"m"))),),
        root = :frog,
    )
end

@testset "Single-part composite is identity" begin
    only = Body(Cylinder(2u"kg", ρ, 2.0), Naked())
    cb = CompositeBody(; parts = (; only), joins = (), root = :only)
    @test total_area(cb) ≈ total_area(only)
    @test flesh_volume(cb) ≈ flesh_volume(only)
    @test skin_radius(cb) == skin_radius(only)
end
