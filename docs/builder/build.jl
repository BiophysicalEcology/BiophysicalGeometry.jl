# Compile the builder page's entry point to wasm, and write what the page needs to know about the animals:
#     julia --project=docs docs/builder/build.jl
# writes docs/src/public/animal_builder.wasm and docs/src/components/animals.json.

using Whisk, JSON
include(joinpath(@__DIR__, "entry.jl"))

const OUT = joinpath(@__DIR__, "..", "src")

# Per animal: its settings, and its parts in the order the wasm module numbers them.
animals = map(AB.ANIMAL_NAMES) do name
    s = AB.settings(name)
    body = AB.build(AB.animal(name), s)
    (; name, settings = s, parts = collect(string.(keys(body.parts))))
end
mkpath(joinpath(OUT, "components"))
open(joinpath(OUT, "components", "animals.json"), "w") do io
    JSON.print(io, (; settings = collect(string.(AB.SETTING_NAMES)), animals, triangle_floats = TRIANGLE_FLOATS,
                      numbers = (; triangle_count = TRIANGLE_COUNT, outer_area = OUTER_AREA, skin_area = SKIN_AREA,
                                   hidden_area = HIDDEN_AREA, volume = VOLUME, meeh = MEEH, shadow_area = SHADOW_AREA,
                                   shadow_extent = SHADOW_EXTENT, parts = PARTS)), 2)
end

mkpath(joinpath(OUT, "public"))
Whisk.compile_to_wasm(animal_builder, Tuple{Int,Int,Int,Int,Int,Int,Int,Float64,Float64};
                      out_path = joinpath(OUT, "public", "animal_builder.wasm"), name = "animal_builder", optimize = true)
println("wrote animal_builder.wasm, ", filesize(joinpath(OUT, "public", "animal_builder.wasm")), " bytes")
