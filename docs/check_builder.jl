# Check the JavaScript geometry of the "Build an animal" page against BiophysicalGeometry.jl.
#
# Each preset is built in JavaScript, which also writes the Julia code for it. That code is run here, and the
# areas and silhouettes from the two are compared.
#
#   julia --project=docs docs/check_builder.jl

using BiophysicalGeometry, Unitful
using NodeJS_20_jll: node

const MODULE = joinpath(@__DIR__, "src", "components", "animalGeometry.mjs")

const PRESETS = [
    "{}",
    """{"torsoShape": "Ellipsoid"}""",
    """{"mass": 0.02, "torsoShape": "Ellipsoid", "torsoRatio": 2, "fat": 0.05, "backFur": 0.004, "bellyFur": 0.003, "headFraction": 0.12, "legFraction": 0.02, "legRatio": 4, "legTop": 1, "limbFur": 0.002}""",
    """{"mass": 682, "torsoRatio": 2.4, "fat": 0, "backFur": 0.0034, "bellyFur": 0.0034, "headRatio": 1.8, "headFraction": 0.05, "legFraction": 0.02, "legRatio": 4.3, "legTop": 0.5, "limbFur": 0.0034}""",
    """{"mass": 0.05, "torsoShape": "Ellipsoid", "torsoRatio": 1.6, "fat": 0.05, "backFur": 0.006, "bellyFur": 0.006, "headShape": "Sphere", "headFraction": 0.1, "legs": 2, "legFraction": 0.02, "legRatio": 8, "legTop": 1, "limbFur": 0}""",
    """{"mass": 100, "torsoShape": "Ellipsoid", "torsoRatio": 4, "fat": 0.35, "backFur": 0.003, "bellyFur": 0.003, "headShape": "Sphere", "headFraction": 0.05, "legs": 0, "limbFur": 0.003}""",
    """{"mass": 5, "fat": 0.3, "backFur": 0, "bellyFur": 0, "headShape": "None", "legs": 2, "limbFur": 0}""",
    """{"neck": true, "nose": "Nose", "ears": "Cone", "tail": true}""",
    """{"wings": "Folded"}""",
    """{"legScaling": "Elastic"}""",
    """{"legScaling": "Elastic", "mass": 4000, "legTop": 1, "limbFur": 0}""",
    """{"legScaling": "Geometric", "mass": 0.02, "torsoShape": "Ellipsoid", "legs": 2}""",
    """{"mass": 800, "torsoRatio": 1.8, "headRatio": 2, "headFraction": 0.02, "neck": true, "neckFraction": 0.1, "neckRatio": 6, "neckPosture": "Up", "ears": "Plate", "earFraction": 0.0005, "legFraction": 0.05, "legRatio": 10, "tail": true}""",
    """{"torsoShape": "Ellipsoid", "neck": true, "neckRatio": 3, "neckPosture": "Up", "nose": "Nose", "ears": "Cone"}""",
    """{"neck": true, "neckRatio": 4, "neckPosture": "Up", "headShape": "Sphere", "nose": "Beak"}""",
    """{"ears": "Plate", "earPosture": "Up", "earFraction": 0.01}""",
    """{"ears": "Plate", "earPosture": "Flat", "earFraction": 0.01, "earRatio": 2.5, "earFlatness": 20}""",
    """{"mass": 4000, "torsoRatio": 1.8, "fat": 0, "backFur": 0, "bellyFur": 0, "limbFur": 0, "headShape": "Sphere", "ears": "Plate", "earPosture": "Flat", "earFraction": 0.01, "earRatio": 1.2, "earFlatness": 30}""",
    """{"torsoShape": "Ellipsoid", "neck": true, "ears": "Plate", "earPosture": "Up", "earFraction": 0.005, "limbFur": 0.004}""",
    """{"wings": "Spread", "limbFur": 0}""",
    """{"mass": 0.05, "torsoShape": "Ellipsoid", "torsoRatio": 1.6, "headShape": "Sphere", "legs": 2, "neck": true, "nose": "Beak", "tail": true, "wings": "Folded", "wingFraction": 0.06, "backFur": 0.006, "bellyFur": 0.006, "limbFur": 0.002}""",
    """{"mass": 0.05, "torsoShape": "Ellipsoid", "torsoRatio": 1.6, "headShape": "Sphere", "legs": 2, "wings": "Spread", "wingFraction": 0.1, "backFur": 0.006, "bellyFur": 0.006, "limbFur": 0.002}""",
    """{"torsoShape": "Ellipsoid", "neck": true, "nose": "Beak", "ears": "Cone", "tail": true, "headShape": "Sphere", "legs": 2}""",
    """{"torsoShape": "Ellipsoid", "nose": "Nose", "ears": "Cone", "tail": true, "tailRatio": 12, "limbFur": 0}""",
    """{"mass": 682, "headShape": "Sphere", "neck": true, "nose": "Beak", "ears": "Cone", "tail": true, "tailFraction": 0.05, "tailRatio": 2}""",
]

const DIRECTIONS = [(0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0), (0.5, 0.3, 0.6)]

# Run a preset through the JavaScript: total area, skin area, silhouettes, then the Julia code.
function javascript(preset)
    script = """
    import { build, silhouette } from $(repr("file:///" * replace(MODULE, "\\" => "/")));
    const animal = build($preset);
    const directions = $(replace(repr(collect(collect.(DIRECTIONS))), "(" => "[", ")" => "]"));
    console.log([animal.total, animal.skin, ...directions.map((d) => silhouette(animal.triangles, d, 256).area)].join(" "));
    console.log(animal.code);
    """
    lines = split(read(`$(node()) --input-type=module -e $script`, String), '\n'; limit = 2)
    return parse.(Float64, split(lines[1])), lines[2]
end

function check_builder(; area_tolerance = 1e-4, silhouette_tolerance = 0.03)
    worst_area = 0.0; worst_silhouette = 0.0
    for preset in PRESETS
        numbers, code = javascript(preset)
        animal = Core.eval(Module(), Meta.parseall(code * "\nanimal"))
        package = [ustrip(u"m^2", total_area(animal)), ustrip(u"m^2", skin_area(animal)),
                   (ustrip(u"m^2", silhouette_rasterized(animal, d)) for d in DIRECTIONS)...]
        errors = abs.(numbers .- package) ./ package
        worst_area = max(worst_area, maximum(errors[1:2]))
        worst_silhouette = max(worst_silhouette, maximum(errors[3:end]))
        if maximum(errors[1:2]) > area_tolerance || maximum(errors[3:end]) > silhouette_tolerance
            error("The animal builder disagrees with the package for $preset:\n  page    $numbers\n  package $package")
        end
    end
    @info "Animal builder agrees with the package" presets = length(PRESETS) worst_area worst_silhouette
end

check_builder()
