using BiophysicalGeometry
using Unitful
using Test

const BG = BiophysicalGeometry

# Common test inputs

const density = 1000.0u"kg/m^3"
const mass = 65.0u"kg"
const axis_ratio_b = 5.0
const axis_ratio_c = 5.0
const fibrous_layer_thickness = 10.0u"mm"
const fibre_diameter = 30.0u"μm"
const fibre_density = 3000u"cm^-2"
const fat_fraction = 0.1
const fat_density = 901.0u"kg/m^3"

const fibrous_layer = FibrousLayer(fibrous_layer_thickness, fibre_diameter, fibre_density)
const fat_layer = FatLayer(fat_fraction, fat_density)
const composite = CompositeInsulation(fibrous_layer, fat_layer)

@testset "helpers" begin
    sphere = Sphere(mass, density)
    @test BG.body_volume(sphere) == mass / density
    @test BG.fat_volume(sphere, fat_layer) == mass * fat_fraction / fat_density

    skin_args = (0.25u"m",)
    fur_args = (0.26u"m",)
    fa = BG.fibrous_areas(sphere, fibrous_layer, skin_args, fur_args)
    expected_skin = 4 * π * (0.25u"m")^2
    expected_total = 4 * π * (0.26u"m")^2
    expected_hair = π * (fibre_diameter / 2)^2 * fibre_density * expected_skin
    @test fa.skin == expected_skin
    @test fa.total == expected_total
    @test fa.convection ≈ expected_skin - expected_hair
end

# Reference values captured from the original implementation. They guard the
# refactor against any behaviour change in any public accessor.

@testset "Plate / Naked" begin
    body = Body(Plate(mass, density, axis_ratio_b, axis_ratio_c), Naked())
    @test total_area(body) ≈ 1.2163304590093516u"m^2"
    @test skin_area(body) ≈ 1.2163304590093516u"m^2"
    @test evaporation_area(body) ≈ 1.2163304590093516u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_volume(body) ≈ 0.065u"m^3"
    @test body.geometry.volume ≈ 0.065u"m^3"
    @test silhouette_area(body, NormalToSun()) ≈ 0.27643874068394353u"m^2"
    @test silhouette_area(body, ParallelToSun()) ≈ 0.05528774813678871u"m^2"
    @test silhouette_area(body, Intermediate()) ≈ 0.16586324441036612u"m^2"
end

@testset "Plate / FibrousLayer" begin
    body = Body(Plate(mass, density, axis_ratio_b, axis_ratio_c), fibrous_layer)
    @test total_area(body) ≈ 1.3504052015217138u"m^2"
    @test skin_area(body) ≈ 1.2163304590093516u"m^2"
    @test evaporation_area(body) ≈ 1.190537258877413u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.12756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_volume(body) ≈ 0.065u"m^3"
    @test silhouette_area(body, NormalToSun()) ≈ 0.3050547569365926u"m^2"
    @test silhouette_area(body, ParallelToSun()) ≈ 0.06509308688767174u"m^2"
end

@testset "Plate / FatLayer" begin
    body = Body(Plate(mass, density, axis_ratio_b, axis_ratio_c), fat_layer)
    @test total_area(body) ≈ 1.2163304590093516u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11304560873468857u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Plate / CompositeInsulation" begin
    body = Body(Plate(mass, density, axis_ratio_b, axis_ratio_c), composite)
    @test total_area(body) ≈ 1.3504052015217138u"m^2"
    @test skin_area(body) ≈ 1.2163304590093516u"m^2"
    @test evaporation_area(body) ≈ 1.190537258877413u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.12756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11304560873468857u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Cylinder / Naked" begin
    body = Body(Cylinder(mass, density, axis_ratio_b), Naked())
    @test total_area(body) ≈ 1.1222291434482228u"m^2"
    @test skin_radius(body) ≈ 0.12742495668987025u"m"
    @test insulation_radius(body) ≈ 0.12742495668987025u"m"
    @test flesh_radius(body) ≈ 0.12742495668987025u"m"
    @test flesh_volume(body) ≈ 0.065u"m^3"
    @test silhouette_area(body, NormalToSun()) ≈ 0.32474239174830616u"m^2"
    @test silhouette_area(body, ParallelToSun()) ≈ 0.05101041561128287u"m^2"
    @test silhouette_area(body, Intermediate()) ≈ 0.18787640367979452u"m^2"
    @test silhouette_area(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.20654751165112634u"m^2"
end

@testset "Cylinder / FibrousLayer" begin
    body = Body(Cylinder(mass, density, axis_ratio_b), fibrous_layer)
    @test total_area(body) ≈ 1.2362029452302272u"m^2"
    @test skin_area(body) ≈ 1.1222291434482228u"m^2"
    @test evaporation_area(body) ≈ 1.098431432327489u"m^2"
    @test insulation_radius(body) ≈ 0.13742495668987026u"m"
    @test flesh_radius(body) ≈ 0.12742495668987025u"m"
    @test silhouette_area(body, NormalToSun()) ≈ 0.355724381353875u"m^2"
    @test silhouette_area(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.22924427552149568u"m^2"
end

@testset "Cylinder / FatLayer" begin
    body = Body(Cylinder(mass, density, axis_ratio_b), fat_layer)
    @test flesh_radius(body) ≈ 0.12252472497618691u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Cylinder / CompositeInsulation" begin
    body = Body(Cylinder(mass, density, axis_ratio_b), composite)
    @test total_area(body) ≈ 1.2362029452302272u"m^2"
    @test insulation_radius(body) ≈ 0.13742495668987026u"m"
    @test flesh_radius(body) ≈ 0.12252472497618691u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Sphere / Naked" begin
    body = Body(Sphere(mass, density), Naked())
    @test total_area(body) ≈ 0.7817952526648283u"m^2"
    @test skin_radius(body) ≈ 0.24942591981125847u"m"
    @test silhouette_area(body, NormalToSun()) ≈ 0.19544881316620707u"m^2"
    @test silhouette_area(body, ParallelToSun()) ≈ 0.19544881316620707u"m^2"
    @test silhouette_area(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.19544881316620707u"m^2"
end

@testset "Sphere / FibrousLayer" begin
    body = Body(Sphere(mass, density), fibrous_layer)
    @test total_area(body) ≈ 0.8457394607097782u"m^2"
    @test skin_area(body) ≈ 0.7817952526648283u"m^2"
    @test evaporation_area(body) ≈ 0.7652166976637417u"m^2"
    @test insulation_radius(body) ≈ 0.25942591981125845u"m"
    @test silhouette_area(body, NormalToSun()) ≈ 0.21143486517744456u"m^2"
end

@testset "Sphere / FatLayer" begin
    body = Body(Sphere(mass, density), fat_layer)
    @test flesh_radius(body) ≈ 0.2398340405261943u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Sphere / CompositeInsulation" begin
    body = Body(Sphere(mass, density), composite)
    @test total_area(body) ≈ 0.8457394607097782u"m^2"
    @test insulation_radius(body) ≈ 0.25942591981125845u"m"
    @test flesh_radius(body) ≈ 0.2398340405261943u"m"
end

@testset "Ellipsoid / Naked" begin
    body = Body(Ellipsoid(mass, density, axis_ratio_b, axis_ratio_c), Naked())
    @test total_area(body) ≈ 1.0679282565991794u"m^2"
    @test skin_radius(body) ≈ 0.14586516277963593u"m"
    @test silhouette_area(body, NormalToSun()) ≈ 0.3342127693207218u"m^2"
    @test silhouette_area(body, ParallelToSun()) ≈ 0.06684255386414435u"m^2"
    # NOTE: pinned to the value produced by `silhouette_area(::Ellipsoid, a, b, c, θ)`,
    # which uses an Euler-angle algorithm with a hardcoded body-orientation angle.
    # This value does NOT correspond to the textbook "long axis horizontal, sun
    # zenith θ" projection of an ellipsoid (which would give ≈ 0.2914 m² here).
    # See the comment block on `silhouette_area(shape::Ellipsoid, a, b, c, θ)`
    # in src/shapes/ellipsoid.jl for the open question.
    @test silhouette_area(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.07667366774262614u"m^2"
end

@testset "Ellipsoid / FibrousLayer" begin
    body = Body(Ellipsoid(mass, density, axis_ratio_b, axis_ratio_c), fibrous_layer)
    @test total_area(body) ≈ 1.1587845516526847u"m^2"
    @test skin_area(body) ≈ 1.0679282565991794u"m^2"
    @test evaporation_area(body) ≈ 1.045282036532102u"m^2"
    @test insulation_radius(body) ≈ 0.15586516277963594u"m"
    @test silhouette_area(body, NormalToSun()) ≈ 0.3620218640142718u"m^2"
    # NOTE: see the matching note in the Ellipsoid / Naked testset; pinned to the
    # current algorithm output, not the textbook horizontal-ellipsoid projection
    # (which would be ≈ 0.3158 m²).
    @test silhouette_area(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.08748304512630967u"m^2"
end

@testset "Ellipsoid / FatLayer" begin
    body = Body(Ellipsoid(mass, density, axis_ratio_b, axis_ratio_c), fat_layer)
    @test flesh_radius(body) ≈ 0.14586516277963593u"m"
    @test flesh_volume(body) ≈ 0.065u"m^3"
end

@testset "Ellipsoid / CompositeInsulation" begin
    body = Body(Ellipsoid(mass, density, axis_ratio_b, axis_ratio_c), composite)
    @test total_area(body) ≈ 1.1587845516526847u"m^2"
    @test insulation_radius(body) ≈ 0.15586516277963594u"m"
end

@testset "DesertIguana / Naked" begin
    body = Body(DesertIguana(mass, density), Naked())
    @test total_area(body) ≈ 2.144297264544543u"m^2"
    @test skin_radius(body) ≈ 0.1729422736164134u"m"
    @test silhouette_area(body, NormalToSun()) ≈ 0.7358254132114781u"m^2"
    @test silhouette_area(body, ParallelToSun()) ≈ 0.26142920031318184u"m^2"
    @test silhouette_area(body, Intermediate()) ≈ 0.49862730676232997u"m^2"
end

@testset "LeopardFrog / Naked" begin
    body = Body(LeopardFrog(mass, density), Naked())
    @test total_area(body) ≈ 1.055591874642603u"m^2"
    @test skin_radius(body) ≈ 0.1729422736164134u"m"
end
