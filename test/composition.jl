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
    dorsal = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    ventral = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    torso = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (Join(dorsal = Attachment(Flat(), FullCover()),
                      ventral = Attachment(Flat(), FullCover())),),
    )
    full = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    @test total_area(torso) ≈ total_area(full)
    @test flesh_volume(torso) ≈ flesh_volume(full)
end

@testset "Dorsal/ventral torso == full ellipsoid" begin
    dorsal = Body(HalfEllipsoid(5u"kg", ρ, 1.5, 1.0), Naked())
    ventral = Body(HalfEllipsoid(5u"kg", ρ, 1.5, 1.0), Naked())
    torso = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (Join(dorsal = Attachment(Flat(), FullCover()),
                      ventral = Attachment(Flat(), FullCover())),),
    )
    full = Body(Ellipsoid(10u"kg", ρ, 1.5, 1.0), Naked())
    @test total_area(torso) ≈ total_area(full)
    @test flesh_volume(torso) ≈ flesh_volume(full)
end

@testset "Dorsal furred, ventral naked" begin
    fur = FibrousLayer(0.01u"m", 30u"μm", 3000u"cm^-2")
    dorsal = Body(HalfCylinder(10u"kg", ρ, 3.0), fur)
    ventral = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    torso = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (Join(dorsal = Attachment(Flat(), FullCover()),
                      ventral = Attachment(Flat(), FullCover())),),
    )
    # Sum of unjoined parts minus 2× the (matching) flat area.
    flat = 2 * dorsal.geometry.length.radius_skin * dorsal.geometry.length.length_skin
    @test total_area(torso) ≈ total_area(dorsal) + total_area(ventral) - 2*flat
end

@testset "Dog example" begin
    torso = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    leg = Body(Cylinder(1u"kg", ρ, 5.0), Naked())
    head = Body(Ellipsoid(2u"kg", ρ, 1.5, 1.0), Naked())
    L_torso = torso.geometry.length.length_skin
    r_leg = skin_radius(leg)
    r_head_attach = 0.04u"m"

    dog = CompositeBody(;
        parts = (; torso, head, leg_fl=leg, leg_fr=leg, leg_bl=leg, leg_br=leg),
        joins = (
            Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_head_attach)),
                 head = Attachment(PoleA(), Disc(r_head_attach))),
            Join(torso = Attachment(Lateral(0.2*L_torso, π/2 - π/6), Disc(r_leg)),
                 leg_fl = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
            Join(torso = Attachment(Lateral(0.2*L_torso, -π/2 + π/6), Disc(r_leg)),
                 leg_fr = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
            Join(torso = Attachment(Lateral(0.8*L_torso, π/2 - π/6), Disc(r_leg)),
                 leg_bl = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
            Join(torso = Attachment(Lateral(0.8*L_torso, -π/2 + π/6), Disc(r_leg)),
                 leg_br = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
        ),
    )

    sum_parts = total_area(torso) + total_area(head) + 4 * total_area(leg)
    sub = 2 * (π * r_head_attach^2 + 4 * π * r_leg^2)
    @test total_area(dog) ≈ sum_parts - sub
    @test flesh_volume(dog) ≈ flesh_volume(torso) + flesh_volume(head) + 4 * flesh_volume(leg)

    # Pose: legs sit on the torso lateral surface at the expected (z, φ).
    R_torso = insulation_radius(torso)
    expected_fl_x = R_torso * cos(π/2 - π/6)
    expected_fl_y = R_torso * sin(π/2 - π/6)
    @test dog.poses.leg_fl.translation[1] ≈ expected_fl_x atol=1e-9*u"m"
    @test dog.poses.leg_fl.translation[2] ≈ expected_fl_y atol=1e-9*u"m"

    # Scalar accessors delegate to the root part (first in `parts`).
    @test skin_radius(dog) == skin_radius(torso)
    @test insulation_radius(dog) == insulation_radius(torso)
end

@testset "Join accessors" begin
    torso = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    leg   = Body(Cylinder(1u"kg", ρ, 5.0), Naked())
    head  = Body(Ellipsoid(2u"kg", ρ, 1.5, 1.0), Naked())
    L_torso = torso.geometry.length.length_skin
    R_torso = skin_radius(torso)
    r_leg  = skin_radius(leg)
    r_head = 0.04u"m"

    head_join = Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_head)),
                     head = Attachment(PoleA(), Disc(r_head)))
    leg_join  = Join(torso = Attachment(Lateral(0.2*L_torso, π/2), Disc(r_leg)),
                     leg_fl = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg)))
    dog = CompositeBody(; parts = (; torso, head, leg_fl = leg),
                          joins = (head_join, leg_join))

    # join_partners — names from the type parameters.
    @test join_partners(head_join) == (:torso, :head)
    @test join_partners(leg_join)  == (:torso, :leg_fl)

    # join_area — π r² for a Disc; parent/child agree by construction.
    @test join_area(head_join, dog) ≈ π * r_head^2
    @test join_area(leg_join, dog)  ≈ π * r_leg^2

    # join_position — world centre of the interface (identity root pose).
    hp = join_position(head_join, dog)             # torso EndA at r=0 → origin
    @test hp[1] ≈ 0.0u"m" atol=1e-9u"m"
    @test hp[2] ≈ 0.0u"m" atol=1e-9u"m"
    @test hp[3] ≈ 0.0u"m" atol=1e-9u"m"
    lp = join_position(leg_join, dog)              # torso Lateral at φ=π/2 → (0, R, 0.2L)
    @test lp[1] ≈ 0.0u"m" atol=1e-9u"m"
    @test lp[2] ≈ R_torso
    @test lp[3] ≈ 0.2*L_torso

    # internal_distance — closed forms per surface.
    hd = internal_distance(head_join, dog)
    @test hd.parent ≈ L_torso / 2                                  # cylinder centroid → end cap
    @test hd.child  ≈ head.geometry.length.a_semi_major_skin       # ellipsoid centre → pole
    ld = internal_distance(leg_join, dog)
    @test ld.parent ≈ R_torso                                      # cylinder centroid → lateral
    @test ld.child  ≈ leg.geometry.length.length_skin / 2

    # Half-cylinder flat face: half-disc centroid sits 4R/(3π) off the plane.
    dorsal = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    Rd = dorsal.geometry.length.radius_skin
    @test internal_distance(dorsal, Flat()) ≈ 4Rd / (3π)

    # Sphere: centroid → any surface point is the radius.
    sph = Body(Sphere(2u"kg", ρ), Naked())
    @test internal_distance(sph, Radial(0.0, 0.0)) ≈ skin_radius(sph)

    # TriMesh has no thermal family — flesh_centroid falls through to the
    # generic error, so a lumped-resistance path length fails loudly.
    verts = [(0.0u"m",0.0u"m",0.0u"m"), (1.0u"m",0.0u"m",0.0u"m"),
             (0.0u"m",1.0u"m",0.0u"m"), (0.0u"m",0.0u"m",1.0u"m")]
    faces = [(1,3,2), (1,2,4), (1,4,3), (2,3,4)]
    mesh_body = Body(TriMesh(verts, faces; mass=1u"kg"), Naked())
    @test_throws ErrorException internal_distance(mesh_body, Radial(0.0, 0.0))
end

@testset "Validation" begin
    a = Body(Cylinder(1u"kg", ρ, 2.0), Naked())
    b = Body(Cylinder(1u"kg", ρ, 2.0), Naked())

    # Patch radius too large for the surface.
    @test_throws ErrorException CompositeBody(;
        parts = (; a, b),
        joins = (Join(a = Attachment(EndA(0.0u"m", 0.0), Disc(10.0u"m")),
                      b = Attachment(EndA(0.0u"m", 0.0), Disc(10.0u"m"))),),
    )

    # Mismatched FullCover areas (different sizes).
    big = Body(HalfCylinder(10u"kg", ρ, 3.0), Naked())
    small = Body(HalfCylinder(1u"kg", ρ, 3.0), Naked())
    @test_throws ErrorException CompositeBody(;
        parts = (; big, small),
        joins = (Join(big = Attachment(Flat(), FullCover()),
                      small = Attachment(Flat(), FullCover())),),
    )

    # Empirical shape (LeopardFrog) cannot be joined.
    frog = Body(LeopardFrog(0.04u"kg", ρ), Naked())
    @test_throws ErrorException CompositeBody(;
        parts = (; frog, b),
        joins = (Join(frog = Attachment(EndA(0.0u"m", 0.0), Disc(0.001u"m")),
                      b = Attachment(EndA(0.0u"m", 0.0), Disc(0.001u"m"))),),
    )
end

@testset "Single-part composite is identity" begin
    sphere = Body(Sphere(2u"kg", ρ), Naked())
    cb = CompositeBody(; parts = (; sphere), joins = ())
    @test total_area(cb) ≈ total_area(sphere)
    @test flesh_volume(cb) ≈ flesh_volume(sphere)
    @test skin_radius(cb) == skin_radius(sphere)
end

@testset "silhouette_rasterized" begin
    # Sphere of radius r → silhouette is π·r² regardless of sun direction.
    sphere = Body(Sphere(1u"kg", ρ), Naked())
    cb = CompositeBody(; parts = (; sphere), joins = ())
    R = skin_radius(sphere)
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
    bottom = Body(Sphere(1u"kg", ρ), Naked())
    top = Body(Sphere(1u"kg", ρ), Naked())
    snowman = CompositeBody(;
        parts = (; bottom, top),
        joins = (Join(bottom = Attachment(Radial(0.0, 0.0), Disc(0.001u"m")),
                      top = Attachment(Radial(0.0, 0.0), Disc(0.001u"m"))),),
    )
    Aalong = silhouette_rasterized(snowman, (0.0, 0.0, 1.0); resolution=256)
    Across = silhouette_rasterized(snowman, (1.0, 0.0, 0.0); resolution=256)
    summed = silhouette_area(snowman).normal
    @test abs(Aalong - expected) / expected < 0.02
    @test abs(Across - 2*expected) / (2*expected) < 0.02
    @test summed ≈ 2 * expected            # summed counts both spheres
    @test Aalong < summed                  # overlap → rasteriser < summed
end
