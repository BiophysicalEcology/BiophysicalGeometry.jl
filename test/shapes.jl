using BiophysicalGeometry
using Unitful
using Test

const BG = BiophysicalGeometry
const density = 1000.0u"kg/m^3"
const fur = FibrousLayer(10.0u"mm", 30.0u"μm", 3000u"cm^-2")
const fat = FatLayer(0.1, 901.0u"kg/m^3")
const insulations = (Naked(), fur, fat, CompositeInsulation(fur, fat))

single(b) = CompositeBody(; parts = (; p = b), joins = ())

@testset "HalfCone geometry" begin
    for t in (0.0, 0.4, 1.0), ins in insulations
        h = Body(HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0, top_ratio = t), ins)
        full = Body(Cone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = t), ins)
        l = h.geometry.length
        @test l == full.geometry.length # dims inherited from the full cone of mass 2m
        @test flesh_volume(h) ≈ flesh_volume(full) / 2
        # Total = half the cone's surface + the trapezoidal axial cut face.
        flat = (1 + t) * l.radius_skin * l.length_skin
        @test total_area(h) ≈ total_area(full) / 2 + flat
        @test BG.surface_area(h.shape, h, Flat()) ≈ flat
    end
    # t = 1 is a half cylinder.
    hc = Body(HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0, top_ratio = 1.0), fur)
    hy = Body(HalfCylinder(; mass = 5u"kg", density, axis_ratio_b = 2.0), fur)
    @test total_area(hc) ≈ total_area(hy)
    @test silhouette(hc, 0.7) ≈ silhouette(hy, 0.7)
    @test all(flesh_centroid(hc.shape, hc) .≈ flesh_centroid(hy.shape, hy))
end

@testset "Dorsal/ventral cone == full cone" begin
    dorsal = Body(HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.3), Naked())
    ventral = Body(HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.3), Naked())
    cone = CompositeBody(;
        parts = (; dorsal, ventral),
        joins = (Join(dorsal = Attachment(Flat(), FullCover()),
                      ventral = Attachment(Flat(), FullCover())),),
    )
    full = Body(Cone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.3), Naked())
    @test total_area(cone) ≈ total_area(full)
    @test flesh_volume(cone) ≈ flesh_volume(full)
end

@testset "Half cylindrical surfaces follow the parent" begin
    t = 0.4
    h = Body(HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0, top_ratio = t), fur)
    full = Body(Cone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = t), fur)
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
    @test BG.surface_centroid(h.shape, h, EndB())[3] ≈ t * R / 2
end

@testset "Half-frustum flesh centroid" begin
    # Against direct integration of half-disc slices of radius ρ(z).
    for t in (0.0, 0.5, 1.0)
        h = Body(HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0, top_ratio = t), Naked())
        R = ustrip(u"m", h.geometry.length.radius_skin)
        L = ustrip(u"m", h.geometry.length.length_skin)
        zs = range(0, L; length = 20001)
        r = @. R * (1 - (1 - t) * zs / L)
        w = r .^ 2 # slice area ∝ ρ²
        ȳ = sum(w .* (4 .* r ./ (3π))) / sum(w)
        x̄ = sum(w .* zs) / sum(w)
        c = ustrip.(u"m", flesh_centroid(h.shape, h)) # axis along x, dome along z
        @test c[1] ≈ x̄ rtol = 1e-3
        @test c[2] == 0
        @test c[3] ≈ ȳ rtol = 1e-3
    end
end

@testset "HalfEllipsoid faces" begin
    h = Body(HalfEllipsoid(; mass = 5u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 6.0), fur)
    l = h.geometry.length
    flat = π * l.length_skin / 2 * l.width_skin / 2
    @test BG.surface_area(h.shape, h, Flat()) ≈ flat
    @test BG.surface_area(h.shape, h, Dome()) ≈ total_area(h) - flat
end

@testset "Silhouette: analytic vs rasterized" begin
    # Local sun direction for angle θ from the long axis: every shape lies along
    # x with height along z, and halves put the sun on the dome side.
    axial(θ) = (cos(θ), 0.0, sin(θ))
    domed = axial
    shapes = [
        Cylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0)            => axial,
        Cone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.0)           => axial,
        Cone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.4)           => axial,
        Ellipsoid(; mass = 10u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 3.0)      => domed,
        HalfCylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0)        => axial,
        HalfCone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.0)       => axial,
        HalfCone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.4)       => axial,
        HalfEllipsoid(; mass = 10u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 6.0)  => domed,
        HalfSphere(; mass = 10u"kg", density)               => domed,
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
    for (half, full) in ((HalfCylinder(; mass = 5u"kg", density, axis_ratio_b = 3.0), Cylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0)),
                         (HalfEllipsoid(; mass = 5u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 6.0), Ellipsoid(; mass = 10u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 3.0)),
                         (HalfSphere(; mass = 5u"kg", density), Sphere(; mass = 10u"kg", density)))
        @test silhouette(Body(half, Naked())).normal ≈ silhouette(Body(full, Naked())).normal
    end
    # Hemisphere: (π r²/2)(1 + sin θ).
    h = Body(HalfSphere(; mass = 5u"kg", density), Naked())
    r = h.geometry.length.radius_skin
    @test silhouette(h, 0.4) ≈ π * r^2 / 2 * (1 + sin(0.4))
    # Frustum side view is a trapezoid.
    c = Body(Cone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.4), Naked())
    R, L = c.geometry.length.radius_skin, c.geometry.length.length_skin
    @test silhouette(c).normal ≈ (1 + 0.4) * R * L
    # A cylinder's shadow is the same from either end.
    y = Body(Cylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0), Naked())
    @test silhouette(y, 0.5) ≈ silhouette(y, π - 0.5)
end

@testset "Accessors across shapes" begin
    shapes = (Cylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0), Sphere(; mass = 10u"kg", density), Ellipsoid(; mass = 10u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 3.0),
              Plate(; mass = 10u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 4.0), Cone(; mass = 10u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.3),
              HalfCylinder(; mass = 5u"kg", density, axis_ratio_b = 3.0), HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.3),
              HalfEllipsoid(; mass = 5u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 6.0), HalfSphere(; mass = 5u"kg", density))
    for sh in shapes, ins in insulations
        b = Body(sh, ins)
        @test surface_area(b) == total_area(b)
        @test outer_dims(b.shape, b) isa NamedTuple
        @test flesh_radius(b) ≤ skin_radius(b) ≤ insulation_radius(b)
        @test occursin("mass", sprint(show, MIME"text/plain"(), b))
    end
    @test outer_dims(Sphere(; mass = 10u"kg", density), Body(Sphere(; mass = 10u"kg", density), fur)).radius ≈
          insulation_radius(Body(Sphere(; mass = 10u"kg", density), fur))
    @test BG.mass(HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0)) == 5u"kg"
    @test occursin("Half (mass 5.0 kg)", sprint(show, MIME"text/plain"(), HalfCone(; mass = 5u"kg", density, axis_ratio_b = 2.0)))
end

@testset "Plate silhouette" begin
    # Mesh axes: length along x, width along y, height along z. The sun moves
    # from the smallest face's normal (θ = 0) to the largest face's (θ = π/2).
    unit(i) = ntuple(j -> j == i ? 1.0 : 0.0, 3)
    for (axis_ratio_b, axis_ratio_c) in ((3.0, 4.0), (0.5, 2.0), (2.0, 0.7)), ins in (Naked(), fur)
        b = Body(Plate(; mass = 10u"kg", density, axis_ratio_b, axis_ratio_c), ins)
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
    @test_throws ErrorException Half(Plate(; mass = 10u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 4.0))
    # The suggested replacement is the same box with half the height.
    full = Body(Plate(; mass = 10u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 4.0), Naked()).geometry.length
    half = Body(Plate(; mass = 5u"kg", density, axis_ratio_b = 3.0, axis_ratio_c = 8.0), Naked()).geometry.length
    @test half.length_skin ≈ full.length_skin
    @test half.width_skin ≈ full.width_skin
    @test half.height_skin ≈ full.height_skin / 2
end

@testset "Keyword construction" begin
    # Every sufficient set of keywords gives the same shape.
    ref = Cylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0)
    l = Body(ref, Naked()).geometry.length
    length, radius, volume = l.length_skin, l.radius_skin, 10u"kg" / density
    for alt in (Cylinder(; length, radius, density), Cylinder(; length, radius, mass = 10u"kg"),
                Cylinder(; volume, density, axis_ratio_b = 3.0), Cylinder(; radius, mass = 10u"kg", density),
                Cylinder(; length, axis_ratio_b = 3.0, density))
        @test BG.mass(alt) ≈ BG.mass(ref)
        @test alt.density ≈ ref.density
        @test alt.axis_ratio_b ≈ ref.axis_ratio_b
    end
    # Given values pass through unchanged, units included.
    p = Plate(; mass = 500u"g", density, axis_ratio_b = 3.0, axis_ratio_c = 4.0)
    @test p.mass === 500u"g"
    # Dimensions round-trip for every shape.
    e = Ellipsoid(; length = 0.6u"m", width = 0.2u"m", height = 0.1u"m", density)
    g = Body(e, Naked()).geometry.length
    @test all((g.length_skin, g.width_skin, g.height_skin) .≈ (0.6u"m", 0.2u"m", 0.1u"m"))
    @test BG.mass(e) ≈ density * π / 6 * 0.6u"m" * 0.2u"m" * 0.1u"m"
    s = Sphere(; radius = 0.2u"m", mass = 10u"kg")
    @test Body(s, Naked()).geometry.length.radius_skin ≈ 0.2u"m"
    c = Cone(; length = 0.5u"m", radius = 0.1u"m", top_ratio = 0.5, density)
    @test Body(c, Naked()).geometry.volume ≈ π / 3 * (1 + 0.5 + 0.25) * (0.1u"m")^2 * 0.5u"m"
    # Half shapes describe the half: its own mass, and for a HalfEllipsoid its dome height.
    h = HalfEllipsoid(; length = 0.6u"m", width = 0.2u"m", height = 0.1u"m", density)
    # Its full ellipsoid is twice e's height, so the half weighs the same as e.
    @test BG.mass(h) ≈ BG.mass(e)
    @test BG._domed_semiaxes(h, Body(h, Naked()))[3] ≈ 0.1u"m"
    @test BG.mass(HalfCylinder(; mass = 5u"kg", density, axis_ratio_b = 3.0)) == 5u"kg"
    # Positional construction is gone; bad keyword sets say what's wrong.
    @test_throws MethodError Cylinder(10u"kg", density, 3.0)
    @test_throws ArgumentError Cylinder(; mass = 10u"kg") # under-determined
    @test_throws ArgumentError Cylinder(; mass = 10u"kg", density, volume = 1u"m^3", axis_ratio_b = 3.0)
    @test_throws ArgumentError Cylinder(; mass = 10u"kg", density, axis_ratio = 3.0) # unknown keyword
    @test_throws ArgumentError Sphere(; mass = -1u"kg", density)
    @test_throws ArgumentError Cone(; mass = 1u"kg", density, axis_ratio_b = 2.0, top_ratio = 1.5)
end

@testset "Triaxial ellipsoid" begin
    e = Body(Ellipsoid(; length = 0.6u"m", width = 0.3u"m", height = 0.1u"m", density), Naked())
    a, b, c = 0.3, 0.15, 0.05
    # Surface area against direct integration over (θ, φ) of |∂r/∂θ × ∂r/∂φ|.
    n = 2000
    S = 0.0
    for i in 1:n, j in 1:2n
        θ = (i - 0.5) * π / n
        φ = (j - 0.5) * π / n
        S += sin(θ) * sqrt((b * c * sin(θ) * cos(φ))^2 + (a * c * sin(θ) * sin(φ))^2 + (a * b * cos(θ))^2)
    end
    @test ustrip(u"m^2", total_area(e)) ≈ S * (π / n)^2 rtol = 1e-6
    # The uniform fat shell keeps the body's volume.
    f = Body(Ellipsoid(; length = 0.6u"m", width = 0.3u"m", height = 0.1u"m", density), fat)
    gl = f.geometry.length
    @test 4π / 3 * gl.length_skin * gl.width_skin * gl.height_skin / 8 ≈ f.geometry.volume
    flesh = (gl.length_skin, gl.width_skin, gl.height_skin) ./ 2 .- gl.fat
    @test 4π / 3 * prod(flesh) ≈ flesh_volume(f)
end

@testset "Insulation layer order" begin
    sh = Cylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0)
    @test Body(sh, CompositeInsulation(fat, fur)).geometry == Body(sh, CompositeInsulation(fur, fat)).geometry
end

@testset "Mixed mass units" begin
    # Grams with kg/m³ used to leave g^1/3 kg^-1/3 m lengths that broke poses.
    plate() = Body(Plate(; mass = 500u"g", density, axis_ratio_b = 3.0, axis_ratio_c = 4.0), Naked())
    x, y = plate(), plate()
    @test unit(x.geometry.length.length_skin) == u"m"
    cb = CompositeBody(; parts = (; x, y),
        joins = (Join(x = Attachment(Top(0.0u"m", 0.0u"m"), Disc(1u"mm")),
                      y = Attachment(Bottom(0.0u"m", 0.0u"m"), Disc(1u"mm"))),))
    @test total_area(cb) ≈ 2 * total_area(x) - 2π * (1u"mm")^2
end

@testset "View factors: immediate blocker" begin
    # Three touching spheres stacked along z. Each sees the sphere next to it, never
    # the one beyond: for source → A → B → C, C's blocked view is B, not A.
    ball() = Body(Sphere(; radius = 0.1u"m", density), Naked())
    a, b, c = ball(), ball(), ball()
    stack = CompositeBody(; parts = (; a, b, c),
        joins = (Join(a = Attachment(Radial(π, 0.0), Disc(1u"mm")), b = Attachment(Radial(0.0, 0.0), Disc(1u"mm"))),
                 Join(b = Attachment(Radial(π, 0.0), Disc(1u"mm")), c = Attachment(Radial(0.0, 0.0), Disc(1u"mm")))))
    f = silhouette_factors(stack, Sky(0.5); ndirections = 200, resolution = 64)
    @test f.a.neighbours.c == 0
    @test f.c.neighbours.a == 0
    @test f.c.neighbours.b > 0
    # Symmetric up to the direction and pixel sampling.
    @test f.a.neighbours.b ≈ f.c.neighbours.b rtol = 0.01
    @test f.b.neighbours.a ≈ f.b.neighbours.c rtol = 0.01
    for p in f
        @test p.sky + p.ground + sum(p.neighbours) ≈ 1
    end
    # Overhead beam: only the top sphere is lit.
    lit = silhouette(stack, Beam(0.0, 0.0, 1.0))
    @test lit.a > 0u"m^2"
    @test lit.b == lit.c == 0u"m^2"
end

@testset "Argument checks" begin
    cb = single(Body(Sphere(; radius = 0.1u"m", density), Naked()))
    @test_throws ArgumentError Beam(0.0, 0.0, 0.0)
    @test_throws ArgumentError Beam(NaN, 0.0, 1.0)
    @test_throws ArgumentError Sky(1.5)
    @test_throws ArgumentError Ground(-0.1)
    @test_throws ArgumentError Horizon(Float64[])
    @test_throws ArgumentError silhouette_factors(cb, Sky(0.5); ndirections = 0)
    @test_throws ArgumentError silhouette_rasterized(cb, (0.0, 0.0, 1.0); resolution = 0)
end
