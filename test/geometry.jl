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
    sphere = Sphere(; mass, density)
    @test BG.body_volume(sphere) == mass / density
    @test BG.fat_volume(sphere, fat_layer) == mass * fat_fraction / fat_density
end

# Reference values captured from the original implementation. They guard the
# refactor against any behaviour change in any public accessor.

@testset "Plate / Naked" begin
    body = Body(Plate(; mass, density, axis_ratio_b, axis_ratio_c), Naked())
    @test total_area(body) ≈ 1.2163304590093516u"m^2"
    @test skin_area(body) ≈ 1.2163304590093516u"m^2"
    @test evaporation_area(body) ≈ 1.2163304590093516u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_volume(body) ≈ 0.065u"m^3"
    @test body.geometry.volume ≈ 0.065u"m^3"
    @test silhouette(body, NormalToSun()) ≈ 0.27643874068394353u"m^2"
    @test silhouette(body, ParallelToSun()) ≈ 0.05528774813678871u"m^2"
    @test silhouette(body, Intermediate()) ≈ 0.16586324441036612u"m^2"
end

@testset "Plate / FibrousLayer" begin
    body = Body(Plate(; mass, density, axis_ratio_b, axis_ratio_c), fibrous_layer)
    @test total_area(body) ≈ 1.3504052015217138u"m^2"
    @test skin_area(body) ≈ 1.2163304590093516u"m^2"
    @test evaporation_area(body) ≈ 1.190537258877413u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.12756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_volume(body) ≈ 0.065u"m^3"
    @test silhouette(body, NormalToSun()) ≈ 0.3050547569365926u"m^2"
    @test silhouette(body, ParallelToSun()) ≈ 0.06509308688767174u"m^2"
end

@testset "Plate / FatLayer" begin
    body = Body(Plate(; mass, density, axis_ratio_b, axis_ratio_c), fat_layer)
    @test total_area(body) ≈ 1.2163304590093516u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.11756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11304560873468857u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Plate / CompositeInsulation" begin
    body = Body(Plate(; mass, density, axis_ratio_b, axis_ratio_c), composite)
    @test total_area(body) ≈ 1.3504052015217138u"m^2"
    @test skin_area(body) ≈ 1.2163304590093516u"m^2"
    @test evaporation_area(body) ≈ 1.190537258877413u"m^2"
    @test skin_radius(body) ≈ 0.11756673438603786u"m"
    @test insulation_radius(body) ≈ 0.12756673438603786u"m"
    @test flesh_radius(body) ≈ 0.11304560873468857u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Cylinder / Naked" begin
    body = Body(Cylinder(; mass, density, axis_ratio_b), Naked())
    @test total_area(body) ≈ 1.1222291434482228u"m^2"
    @test skin_radius(body) ≈ 0.12742495668987025u"m"
    @test insulation_radius(body) ≈ 0.12742495668987025u"m"
    @test flesh_radius(body) ≈ 0.12742495668987025u"m"
    @test flesh_volume(body) ≈ 0.065u"m^3"
    @test silhouette(body, NormalToSun()) ≈ 0.32474239174830616u"m^2"
    @test silhouette(body, ParallelToSun()) ≈ 0.05101041561128287u"m^2"
    @test silhouette(body, Intermediate()) ≈ 0.18787640367979452u"m^2"
    @test silhouette(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.20654751165112634u"m^2"
end

@testset "Cylinder / FibrousLayer" begin
    body = Body(Cylinder(; mass, density, axis_ratio_b), fibrous_layer)
    @test total_area(body) ≈ 1.2362029452302272u"m^2"
    @test skin_area(body) ≈ 1.1222291434482228u"m^2"
    @test evaporation_area(body) ≈ 1.098431432327489u"m^2"
    @test insulation_radius(body) ≈ 0.13742495668987026u"m"
    @test flesh_radius(body) ≈ 0.12742495668987025u"m"
    @test silhouette(body, NormalToSun()) ≈ 0.355724381353875u"m^2"
    @test silhouette(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.22924427552149568u"m^2"
end

@testset "Cylinder / FatLayer" begin
    body = Body(Cylinder(; mass, density, axis_ratio_b), fat_layer)
    @test flesh_radius(body) ≈ 0.12252472497618691u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Cylinder / CompositeInsulation" begin
    body = Body(Cylinder(; mass, density, axis_ratio_b), composite)
    @test total_area(body) ≈ 1.2362029452302272u"m^2"
    @test insulation_radius(body) ≈ 0.13742495668987026u"m"
    @test flesh_radius(body) ≈ 0.12252472497618691u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Sphere / Naked" begin
    body = Body(Sphere(; mass, density), Naked())
    @test total_area(body) ≈ 0.7817952526648283u"m^2"
    @test skin_radius(body) ≈ 0.24942591981125847u"m"
    @test silhouette(body, NormalToSun()) ≈ 0.19544881316620707u"m^2"
    @test silhouette(body, ParallelToSun()) ≈ 0.19544881316620707u"m^2"
    @test silhouette(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.19544881316620707u"m^2"
end

@testset "Sphere / FibrousLayer" begin
    body = Body(Sphere(; mass, density), fibrous_layer)
    @test total_area(body) ≈ 0.8457394607097782u"m^2"
    @test skin_area(body) ≈ 0.7817952526648283u"m^2"
    @test evaporation_area(body) ≈ 0.7652166976637417u"m^2"
    @test insulation_radius(body) ≈ 0.25942591981125845u"m"
    @test silhouette(body, NormalToSun()) ≈ 0.21143486517744456u"m^2"
end

@testset "Sphere / FatLayer" begin
    body = Body(Sphere(; mass, density), fat_layer)
    @test flesh_radius(body) ≈ 0.2398340405261943u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Sphere / CompositeInsulation" begin
    body = Body(Sphere(; mass, density), composite)
    @test total_area(body) ≈ 0.8457394607097782u"m^2"
    @test insulation_radius(body) ≈ 0.25942591981125845u"m"
    @test flesh_radius(body) ≈ 0.2398340405261943u"m"
end

@testset "Ellipsoid / Naked" begin
    body = Body(Ellipsoid(; mass, density, axis_ratio_b, axis_ratio_c), Naked())
    @test total_area(body) ≈ 1.0679282565991794u"m^2"
    @test skin_radius(body) ≈ 0.14586516277963593u"m"
    @test silhouette(body, NormalToSun()) ≈ 0.3342127693207218u"m^2"
    @test silhouette(body, ParallelToSun()) ≈ 0.06684255386414435u"m^2"
    # Zenith θ is measured from the long (a) axis, as for `Cylinder`:
    # π·c·√(b²cos²θ + a²sin²θ) — see test/hand_check.jl.
    @test silhouette(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.17684877452096545u"m^2"
end

@testset "Ellipsoid / FibrousLayer" begin
    body = Body(Ellipsoid(; mass, density, axis_ratio_b, axis_ratio_c), fibrous_layer)
    @test total_area(body) ≈ 1.1587845516526847u"m^2"
    @test skin_area(body) ≈ 1.0679282565991794u"m^2"
    @test evaporation_area(body) ≈ 1.045282036532102u"m^2"
    @test insulation_radius(body) ≈ 0.15586516277963594u"m"
    @test silhouette(body, NormalToSun()) ≈ 0.3620218640142718u"m^2"
    @test silhouette(body, ZenithAngleVarying(), 30.0u"°") ≈ 0.19270108448901746u"m^2"
end

@testset "Ellipsoid / FatLayer" begin
    body = Body(Ellipsoid(; mass, density, axis_ratio_b, axis_ratio_c), fat_layer)
    # Newton fat solve; the Cardano solve these were first pinned against
    # returned the wrong root and clamped the fat layer to zero.
    @test flesh_radius(body) ≈ 0.14025579774517113u"m"
    @test flesh_volume(body) ≈ 0.057785793562708104u"m^3"
end

@testset "Ellipsoid / CompositeInsulation" begin
    body = Body(Ellipsoid(; mass, density, axis_ratio_b, axis_ratio_c), composite)
    @test total_area(body) ≈ 1.143581562260096u"m^2"
    @test insulation_radius(body) ≈ 0.15794460094666143u"m"
end

