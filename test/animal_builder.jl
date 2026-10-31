using BiophysicalGeometry, Unitful, Random, Test
const AB = BiophysicalGeometry.AnimalBuilder

# Settings drawn at random over the ranges the app offers. Legs sized by hand: similarity needs BiologicalScaling.jl.
function random_settings(rng)
    pick(xs) = xs[rand(rng, 1:length(xs))]
    (;
        mass = 10.0^(-2 + 6 * rand(rng)), posture = pick(["Horizontal", "Upright"]), pitch = pick(-30.0:80),
        torsoShape = pick(["Cylinder", "Ellipsoid"]), torsoRatio = pick(1.1:0.1:6), fat = pick(0:0.05:0.5),
        backFur = pick(0:0.005:0.06), bellyFur = pick(0:0.005:0.06), limbFur = pick(0:0.005:0.06),
        headShape = pick(["None", "Sphere", "Ellipsoid"]), headFraction = pick(0.01:0.01:0.15),
        neck = rand(rng, Bool), neckPosture = pick(["Forward", "Up"]), neckRatio = pick(0.5:0.5:8),
        nose = pick(["None", "Nose", "Beak"]), ears = pick(["None", "Cone", "Plate"]),
        earPosture = pick(["Up", "Flat"]), legs = pick([0, 2, 4]), legRatio = pick(1.0:12),
        legTop = pick(0.1:0.1:1), hindLegs = pick(["Same", "Different"]), arms = rand(rng, Bool),
        wings = pick(["None", "Folded", "Spread"]), tail = rand(rng, Bool), tailRatio = pick(1.0:20),
    )
end

@testset "presets" begin
    for (name, preset) in AB.PRESETS
        code = AB.animal_code(preset)
        animal = AB.build_animal(code)
        @test animal isa CompositeBody
        @test isfinite(ustrip(total_area(animal))) && total_area(animal) > skin_area(animal) * 0.5
        # The code is the model: running it again gives the same animal
        @test total_area(AB.build_animal(code)) == total_area(animal)
        @test AB.animal_code(preset) == code
    end
end

@testset "random settings" begin
    rng = MersenneTwister(1)
    built = 0
    for _ in 1:200
        s = random_settings(rng)
        code = try
            AB.animal_code(s)
        catch e
            # The only refusal the recipe makes itself: parts heavier than the animal
            @test occursin("weigh more than the animal", sprint(showerror, e))
            continue
        end
        animal = AB.build_animal(code)
        @test total_area(animal) > 0u"m^2"
        built += 1
    end
    @test built > 150
end

@testset "settings are checked" begin
    bad = [(; legs = 3), (; mass = -1.0), (; mass = NaN), (; backFur = -0.001), (; torsoRatio = 0.0),
           (; torsoShape = "Cube"), (; legScaling = "Elastic\nrun(`ls`)"), (; pitch = Inf), (; fat = 1.0),
           (; legTop = 1.5), (; colour = "red")]
    for s in bad
        @test_throws ArgumentError AB.animal_code(s)
    end
end

@testset "the parts weigh what the animal does" begin
    for (name, preset) in AB.PRESETS
        animal = AB.build_animal(AB.animal_code(preset))
        total = sum(BiophysicalGeometry.mass(part.shape) for part in values(animal.parts))
        @test ustrip(u"kg", total) ≈ preset.mass rtol = 1e-4
    end
end

@testset "tiny animals" begin
    animal = AB.build_animal(AB.animal_code((; mass = 1e-4, neck = true, nose = "Beak", ears = "Plate", tail = true)))
    @test total_area(animal) > 0u"m^2"
end
