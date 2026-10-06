using BiophysicalGeometry
using Unitful
using Test

const BG = BiophysicalGeometry
const ρ = 1000.0u"kg/m^3"
const fur = FibrousLayer(10.0u"mm", 30.0u"μm", 3000u"cm^-2")
const fat = FatLayer(0.1, 901.0u"kg/m^3")
const insulations = (Naked(), fur, fat, CompositeInsulation(fur, fat))

single(b) = CompositeBody(; parts = (; p = b), joins = ())

@testset "HalfCone geometry" begin
    for t in (0.0, 0.4, 1.0), ins in insulations
        h = Body(HalfCone(5u"kg", ρ, 2.0, t), ins)
        full = Body(Cone(10u"kg", ρ, 2.0, t), ins)
        l = h.geometry.length
        @test l == full.geometry.length                     # dims inherited from Cone(2m)
        @test flesh_volume(h) ≈ flesh_volume(full) / 2
        # Total = half the cone's surface + the trapezoidal axial cut face.
        flat = (1 + t) * l.radius_skin * l.length_skin
        @test total_area(h) ≈ total_area(full) / 2 + flat
        @test BG.surface_area(h.shape, h, Flat()) ≈ flat
    end
    # t = 1 is a half cylinder.
    hc = Body(HalfCone(5u"kg", ρ, 2.0, 1.0), fur)
    hy = Body(HalfCylinder(5u"kg", ρ, 2.0), fur)
    @test total_area(hc) ≈ total_area(hy)
    @test silhouette(hc, 0.7) ≈ silhouette(hy, 0.7)
    @test all(flesh_centroid(hc.shape, hc) .≈ flesh_centroid(hy.shape, hy))
end

@testset "Dorsal/ventral cone == full cone" begin
    dorsal = Body(HalfCone(5u"kg", ρ, 2.0, 0.3), Naked())
    ventral = Body(HalfCone(5u"kg", ρ, 2.0, 0.3), Naked())
    cone = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (Join(dorsal = Attachment(Flat(), FullCover()),
                      ventral = Attachment(Flat(), FullCover())),),
    )
    full = Body(Cone(10u"kg", ρ, 2.0, 0.3), Naked())
    @test total_area(cone) ≈ total_area(full)
    @test flesh_volume(cone) ≈ flesh_volume(full)
end

@testset "Half cylindrical surfaces follow the parent" begin
    t = 0.4
    h = Body(HalfCone(5u"kg", ρ, 2.0, t), fur)
    full = Body(Cone(10u"kg", ρ, 2.0, t), fur)
    for S in (EndA, EndB, Lateral)
        @test BG.surface_area(h.shape, h, S()) ≈ BG.surface_area(full.shape, full, S()) / 2
    end
    R = h.geometry.length.radius_skin
    L = h.geometry.length.length_skin
    # A point on the slant is the cone's point; φ is restricted to the dome side.
    loc = Lateral(L / 2, π / 3)
    @test all(BG.surface_point(h.shape, h, loc) .≈ BG.surface_point(full.shape, full, loc))
    @test all(BG.surface_normal(h.shape, h, loc) .≈ BG.surface_normal(full.shape, full, loc))
    @test_throws ErrorException BG.validate_range(h.shape, h, Lateral(L / 2, -0.1))
    # The top end and the flat face narrow to the top radius.
    @test (BG.validate_range(h.shape, h, EndB(0.9t * R, 0.5)); true)
    @test_throws ErrorException BG.validate_range(h.shape, h, EndB(1.1t * R, 0.5))
    @test_throws ErrorException BG.validate_range(h.shape, h, Flat(L, 1.1t * R))
    @test BG.surface_centroid(h.shape, h, EndB())[2] ≈ t * R / 2
end

@testset "Half-frustum flesh centroid" begin
    # Against direct integration of half-disc slices of radius ρ(z).
    for t in (0.0, 0.5, 1.0)
        h = Body(HalfCone(5u"kg", ρ, 2.0, t), Naked())
        R = ustrip(u"m", h.geometry.length.radius_skin)
        L = ustrip(u"m", h.geometry.length.length_skin)
        zs = range(0, L; length = 20001)
        r = @. R * (1 - (1 - t) * zs / L)
        w = r .^ 2                                    # slice area ∝ ρ²
        ȳ = sum(w .* (4 .* r ./ (3π))) / sum(w)
        z̄ = sum(w .* zs) / sum(w)
        c = ustrip.(u"m", flesh_centroid(h.shape, h))
        @test c[1] == 0
        @test c[2] ≈ ȳ rtol = 1e-3
        @test c[3] ≈ z̄ rtol = 1e-3
    end
end

@testset "HalfEllipsoid faces" begin
    h = Body(HalfEllipsoid(5u"kg", ρ, 3.0, 3.0), fur)
    l = h.geometry.length
    flat = π * l.a_semi_major_skin * l.b_semi_minor_skin
    @test BG.surface_area(h.shape, h, Flat()) ≈ flat
    @test BG.surface_area(h.shape, h, Dome()) ≈ total_area(h) - flat
end

@testset "Silhouette: analytic vs rasterized" begin
    # Local sun direction for angle θ from the long axis. Halves put the sun on
    # the dome side, in the plane of the axis and the flat-face normal.
    axial(θ) = (0.0, sin(θ), cos(θ))       # cylindrical: axis z, dome +y
    domed(θ) = (cos(θ), 0.0, sin(θ))       # domed: axis x, dome +z
    shapes = [
        Cylinder(10u"kg", ρ, 3.0)            => axial,
        Cone(10u"kg", ρ, 2.0, 0.0)           => axial,
        Cone(10u"kg", ρ, 2.0, 0.4)           => axial,
        Ellipsoid(10u"kg", ρ, 3.0, 3.0)      => domed,
        HalfCylinder(10u"kg", ρ, 3.0)        => axial,
        HalfCone(10u"kg", ρ, 2.0, 0.0)       => axial,
        HalfCone(10u"kg", ρ, 2.0, 0.4)       => axial,
        HalfEllipsoid(10u"kg", ρ, 3.0, 3.0)  => domed,
        HalfSphere(10u"kg", ρ)               => domed,
    ]
    for (sh, dir) in shapes, ins in (Naked(), fur)
        b = Body(sh, ins)
        for θ in (0.0, 0.3, 1.0, π / 2, 2.2, 3.0)
            @test silhouette(b, θ) ≈ silhouette_rasterized(single(b), dir(θ); resolution = 400) rtol = 0.01
        end
        s = silhouette(b)
        @test s.normal ≈ silhouette(b, π / 2)
        @test s.parallel ≈ silhouette(b, 0.0)
    end
end

@testset "Half silhouette closed forms" begin
    # Seen face-on, a dorsal half casts the full shape's shadow.
    for (half, full) in ((HalfCylinder(5u"kg", ρ, 3.0), Cylinder(10u"kg", ρ, 3.0)),
                         (HalfEllipsoid(5u"kg", ρ, 3.0, 3.0), Ellipsoid(10u"kg", ρ, 3.0, 3.0)),
                         (HalfSphere(5u"kg", ρ), Sphere(10u"kg", ρ)))
        @test silhouette(Body(half, Naked())).normal ≈ silhouette(Body(full, Naked())).normal
    end
    # Hemisphere: (π r²/2)(1 + sin θ).
    h = Body(HalfSphere(5u"kg", ρ), Naked())
    r = h.geometry.length.radius_skin
    @test silhouette(h, 0.4) ≈ π * r^2 / 2 * (1 + sin(0.4))
    # Frustum side view is a trapezoid.
    c = Body(Cone(10u"kg", ρ, 2.0, 0.4), Naked())
    R, L = c.geometry.length.radius_skin, c.geometry.length.length_skin
    @test silhouette(c).normal ≈ (1 + 0.4) * R * L
    # A cylinder's shadow is the same from either end.
    y = Body(Cylinder(10u"kg", ρ, 3.0), Naked())
    @test silhouette(y, 0.5) ≈ silhouette(y, π - 0.5)
end

@testset "Accessors across shapes" begin
    shapes = (Cylinder(10u"kg", ρ, 3.0), Sphere(10u"kg", ρ), Ellipsoid(10u"kg", ρ, 3.0, 3.0),
              Plate(10u"kg", ρ, 3.0, 4.0), Cone(10u"kg", ρ, 2.0, 0.3),
              HalfCylinder(5u"kg", ρ, 3.0), HalfCone(5u"kg", ρ, 2.0, 0.3),
              HalfEllipsoid(5u"kg", ρ, 3.0, 3.0), HalfSphere(5u"kg", ρ))
    for sh in shapes, ins in insulations
        b = Body(sh, ins)
        @test surface_area(b) == total_area(b)
        @test outer_dims(b.shape, b) isa NamedTuple
        @test flesh_radius(b) ≤ skin_radius(b) ≤ insulation_radius(b)
        @test occursin("mass", sprint(show, MIME"text/plain"(), b))
    end
    @test outer_dims(Sphere(10u"kg", ρ), Body(Sphere(10u"kg", ρ), fur)).r ≈
          insulation_radius(Body(Sphere(10u"kg", ρ), fur))
    @test BG.mass(HalfCone(5u"kg", ρ, 2.0)) == 5u"kg"
    @test occursin("Half (mass 5.0 kg)", sprint(show, MIME"text/plain"(), HalfCone(5u"kg", ρ, 2.0)))
end

@testset "Plate silhouette" begin
    # Mesh axes: length along x, width along y, height along z. The sun moves
    # from the smallest face's normal (θ = 0) to the largest face's (θ = π/2).
    unit(i) = ntuple(j -> j == i ? 1.0 : 0.0, 3)
    for (b_ratio, c_ratio) in ((3.0, 4.0), (0.5, 2.0), (2.0, 0.7)), ins in (Naked(), fur)
        b = Body(Plate(10u"kg", ρ, b_ratio, c_ratio), ins)
        dims = values(outer_dims(b.shape, b))
        long, short = argmax(dims), argmin(dims)
        for θ in (0.0, 0.4, 1.0, π / 2, 2.5)
            d = cos(θ) .* unit(long) .+ sin(θ) .* unit(short)
            @test silhouette(b, θ) ≈ silhouette_rasterized(single(b), d; resolution = 400) rtol = 0.01
        end
        s = silhouette(b)
        @test s.normal ≈ silhouette(b, π / 2)
        @test s.parallel ≈ silhouette(b, 0.0)
        @test silhouette(b, ZenithAngleVarying(), 90u"°") ≈ s.normal
    end
end

@testset "A half plate is a plate" begin
    @test_throws ErrorException Half(Plate(10u"kg", ρ, 3.0, 4.0))
    # The suggested replacement is the same box with half the height.
    full = Body(Plate(10u"kg", ρ, 3.0, 4.0), Naked()).geometry.length
    half = Body(Plate(5u"kg", ρ, 3.0, 8.0), Naked()).geometry.length
    @test half.length_skin ≈ full.length_skin
    @test half.width_skin ≈ full.width_skin
    @test half.height_skin ≈ full.height_skin / 2
end
