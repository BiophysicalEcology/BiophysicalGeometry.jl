using BiophysicalGeometry, Unitful, Random, Test
const AB = BiophysicalGeometry.AnimalBuilder

# Settings drawn at random over the ranges the builder page offers.
function random_settings(rng, s)
    pick(xs) = xs[rand(rng, 1:length(xs))]
    merge(s, (;
        mass = 10.0^(-2 + 6 * rand(rng)), pitch = pick(-30.0:80), torsoRatio = pick(1.1:0.1:6), fat = pick(0:0.05:0.5),
        backFur = pick(0:0.005:0.06), bellyFur = pick(0:0.005:0.06), limbFur = pick(0:0.005:0.06),
        headRatio = pick(1.0:0.1:3), headFraction = pick(0.01:0.01:0.15), neckRatio = pick(0.5:0.5:8),
        neckAngle = pick(-30.0:90), earAngle = pick(0.0:90), legRatio = pick(1.0:12), legTop = pick(0.1:0.1:1),
        wingFold = pick(0.0:90), tailRatio = pick(1.0:20)))
end

@testset "every animal" begin
    for name in AB.ANIMAL_NAMES
        s = AB.settings(name)
        animal = AB.build(AB.animal(name), s)
        @test animal isa CompositeBody
        @test isfinite(ustrip(total_area(animal))) && total_area(animal) > skin_area(animal) * 0.5
        # The parts weigh what the animal does.
        total = sum(BiophysicalGeometry.mass(part.shape) for part in values(animal.parts))
        @test ustrip(u"kg", total) ≈ s.mass rtol = 1e-9
    end
end

@testset "random settings" begin
    rng = MersenneTwister(1)
    for _ in 1:200
        name = rand(rng, AB.ANIMAL_NAMES)
        s = random_settings(rng, AB.settings(name))
        animal = AB.build(AB.animal(name), s)
        @test total_area(animal) > 0u"m^2"
        @test isfinite(ustrip(silhouette_rasterized(animal, (0.3, 0.2, 1.0); resolution = 32)))
    end
end

@testset "tiny animals" begin
    animal = AB.build(AB.Bird(), merge(AB.settings("Bird"), (; mass = 1e-4)))
    @test total_area(animal) > 0u"m^2"
end

@testset "built from numbers alone" begin
    s = AB.settings("Dog")
    build_dog(s) = AB.build(AB.Dog(), s)
    build_dog(s)
    @test @allocated(build_dog(s)) == 0
end

@testset "joint angles" begin
    # A giraffe's neck rises, and its head stays level.
    giraffe = AB.build(AB.Giraffe(), AB.settings("Giraffe"))
    along(body, part) = BiophysicalGeometry.apply_rotation(getfield(body.poses, part).rotation, (1.0, 0.0, 0.0))
    @test along(giraffe, :neck)[3] ≈ sind(70) atol = 1e-6
    @test along(giraffe, :head)[3] ≈ 0 atol = 1e-6
    # Legs swing forward, toward the head at +x, and tails rise.
    dog = AB.build(AB.Dog(), merge(AB.settings("Dog"), (; legAngle = 30.0, tailAngle = 30.0)))
    @test along(dog, :leg_fl)[1] ≈ 0.5 atol = 1e-6
    @test along(dog, :leg_bl)[1] ≈ 0.5 atol = 1e-6
    @test along(dog, :tail)[3] ≈ 0.5 atol = 1e-6
    # A human faces -y.
    human = AB.build(AB.Human(), merge(AB.settings("Human"), (; legAngle = 30.0)))
    @test along(human, :leg_l)[2] ≈ -0.5 atol = 1e-6
    # Spread legs go out to their own sides, from straight down for a human.
    wide = AB.build(AB.Human(), merge(AB.settings("Human"), (; legSpread = 30.0)))
    @test along(wide, :leg_l)[1] ≈ 0.5 atol = 1e-6
    @test along(wide, :leg_r)[1] ≈ -0.5 atol = 1e-6
    wide = AB.build(AB.Dog(), merge(AB.settings("Dog"), (; legSpread = 30.0)))
    straight = AB.build(AB.Dog(), AB.settings("Dog"))
    @test along(wide, :leg_fl)[2] > along(straight, :leg_fl)[2] > 0
    @test along(wide, :leg_fr)[2] < along(straight, :leg_fr)[2] < 0
    # Joint angles move parts, so they change the silhouette and not the areas.
    up = AB.build(AB.Giraffe(), merge(AB.settings("Giraffe"), (; neckAngle = 0.0)))
    @test total_area(up) ≈ total_area(giraffe)
    @test silhouette_rasterized(up, (0.0, 0.0, 1.0)) != silhouette_rasterized(giraffe, (0.0, 0.0, 1.0))
end
