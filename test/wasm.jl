# Models compiled to wasm with `compile_wasm` give, run under node, exactly what the same entry point gives run in
# Julia. This checks that the machine-level paths compile: Unchecked construction, joins with bends, the meshes,
# rasterising and the areas.
using BiophysicalGeometry, Unitful, Whisk, Test
const BG = BiophysicalGeometry
const Ext = Base.get_extension(BG, :BiophysicalGeometryWhiskExt)

const density = 1000.0u"kg/m^3"

# A torso with a neck bent up by `bend` degrees, and a truncated head.
function necked(s)
    torso = Body(Cylinder(Unchecked(); mass = s.mass * u"kg", density, axis_ratio_b = s.ratio), Naked())
    neck = Body(Cone(Unchecked(); mass = 0.05 * s.mass * u"kg", density, axis_ratio_b = 2.0, top_ratio = 0.7), Naked())
    head = Body(Ellipsoid(Unchecked(); mass = 0.1 * s.mass * u"kg", density, axis_ratio_b = 1.5, axis_ratio_c = 1.5,
                          pole_a_truncation = 0.2), Naked())
    patch = Disc(0.01u"m" * cbrt(s.mass))
    CompositeBody(Unchecked(); parts = (; torso, neck, head), joins = (
        Join(torso = Attachment(EndB(0.0u"m", 0.0), patch), neck = Attachment(EndA(0.0u"m", 0.0), patch);
             bend = -deg2rad(s.bend), hinge = (0.0, 1.0, 0.0)),
        Join(neck = Attachment(EndB(0.0u"m", 0.0), patch), head = Attachment(PoleB(), patch))))
end

# A furred dorsal and a naked ventral half, and a round head.
function halves(s)
    half() = HalfEllipsoid(Unchecked(); mass = s.mass / 2 * u"kg", density, axis_ratio_b = s.ratio,
                           axis_ratio_c = 2 * s.ratio)
    dorsal = Body(half(), FibrousLayer(5.0u"mm", 30.0u"μm", 3000.0u"cm^-2"))
    ventral = Body(half(), Naked())
    head = Body(Sphere(Unchecked(); mass = 0.1 * s.mass * u"kg", density), Naked())
    patch = Disc(0.01u"m" * cbrt(s.mass))
    CompositeBody(Unchecked(); parts = (; dorsal, ventral, head), joins = (
        Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),
        Join(dorsal = Attachment(Dome(0.3, π / 2), patch), head = Attachment(Radial(π / 2, π), patch))))
end

const MODELS = (; necked = (necked, (; mass = 10.0, ratio = 3.0, bend = 30.0)), halves = (halves, (; mass = 5.0, ratio = 2.0)))
const CASES = [(1, (; mass = 10.0, ratio = 3.0, bend = 30.0), (0.0, 0.0, 1.0)),
               (1, (; mass = 2.0, ratio = 4.0, bend = 70.0), (0.6, 0.3, 0.74)),
               (2, (; mass = 5.0, ratio = 2.0), (0.0, 0.5, 0.86)),
               (2, (; mass = 80.0, ratio = 1.5), (1.0, 0.0, 0.0))]
const CAPACITY, SHADOW = 20000, 64

path = mktempdir()
compile_wasm(path, MODELS)
@test isfile(joinpath(path, "model.wasm")) && isfile(joinpath(path, "biophysical.mjs"))

# Each case under node: its numbers, the number of shadow pixels, and the sum of the triangles' floats.
js(settings) = "{" * join(("$k: $v" for (k, v) in pairs(settings)), ", ") * "}"
cases = join(("[$(i - 1), $(js(s)), [$(join(d, ", "))]]" for (i, s, d) in CASES), ", ")
script = """
import fs from 'node:fs';
import { loadBiophysicalModel } from '$(joinpath(path, "biophysical.mjs"))';
const spec = JSON.parse(fs.readFileSync('$(joinpath(path, "model.json"))'));
const model = await loadBiophysicalModel(fs.readFileSync('$(joinpath(path, "model.wasm"))'), spec, { capacity: $CAPACITY, shadow: $SHADOW });
for (const [index, settings, sun] of [$cases]) {
  const r = model.run(index, settings, sun);
  const covered = r.shadow.covered.reduce((a, b) => a + b, 0);
  const sum = r.triangles.reduce((a, b) => a + b, 0);
  console.log([r.count, r.total, r.skin, r.hidden, r.volume, r.shadow.area, ...r.parts.flatMap((p) => [p.area, p.mass]), covered, sum].join(' '));
}
"""
node = split(read(`$(Sys.which("node")) --input-type=module -e $script`, String), '\n'; keepempty = false)

entry = Ext.Entry{map(first, values(MODELS)), map(m -> map(Float64, last(m)), values(MODELS))}()
for ((index, settings, sun), line) in zip(CASES, node)
    s = collect(values(map(Float64, settings)))
    triangles = zeros(Float32, Ext.TRIANGLE_FLOATS, CAPACITY)
    numbers = zeros(Ext.NUMBERS.parts - 1 + 2 * Ext.MAX_PARTS)
    shadow = zeros(Bool, SHADOW, SHADOW)
    nparts = GC.@preserve s triangles numbers shadow entry(index, Int(pointer(s)), Int(pointer(triangles)), CAPACITY,
        Int(pointer(numbers)), Int(pointer(shadow)), SHADOW, sun...)
    at(k) = numbers[k]
    ntriangles = Int(at(Ext.NUMBERS.triangle_count))
    julia = [ntriangles, at(Ext.NUMBERS.outer_area), at(Ext.NUMBERS.skin_area), at(Ext.NUMBERS.hidden_area),
             at(Ext.NUMBERS.volume), at(Ext.NUMBERS.shadow_area),
             (at(Ext.NUMBERS.parts + i) for i in 0:2nparts - 1)..., count(shadow),
             sum(Float64, view(triangles, :, 1:ntriangles))]
    wasm = parse.(Float64, split(line))
    @test length(wasm) == length(julia)
    @test wasm[1:end-1] ≈ julia[1:end-1] rtol = 1e-12
    # The triangles are summed in Float32 in one and Float64 in the other.
    @test wasm[end] ≈ julia[end] rtol = 1e-5
end
