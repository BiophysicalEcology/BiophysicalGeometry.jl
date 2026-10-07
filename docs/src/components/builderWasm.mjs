// Loads the animal builder's wasm module: BiophysicalGeometry.jl's own machine-level builder, compiled with
// Whisk.jl by docs/builder/build.jl. The page's buffers live in the module's memory, allocated once.

import spec from './animals.json'

const CAPACITY = 80000   // triangles
const SHADOW = 160       // shadow pixels across

export { spec }

export async function loadBuilder(url) {
  const bytes = await (await fetch(url)).arrayBuffer()
  let memory
  const N = (b) => (typeof b === 'bigint' ? Number(b) : b)
  const env = {
    js_console(level, ptr, len) {
      console.error('[wasm]', new TextDecoder().decode(new Uint8Array(memory.buffer, N(ptr), N(len))))
    },
    memset(p, v, n) { new Uint8Array(memory.buffer).fill(Number(v) & 0xff, N(p), N(p) + N(n)); return p },
    memcpy(d, s, n) { new Uint8Array(memory.buffer).copyWithin(N(d), N(s), N(s) + N(n)); return d },
    memmove(d, s, n) { new Uint8Array(memory.buffer).copyWithin(N(d), N(s), N(s) + N(n)); return d },
  }
  const { instance } = await WebAssembly.instantiate(bytes, { env })
  const w = instance.exports
  memory = w.memory
  // pinned, so that the module's collector never takes them back
  const alloc = (n) => { const p = w.whisk_gc_alloc(BigInt(n)); w.whisk_pin(p); return p }
  const nnumbers = spec.numbers.parts + 64
  const buffers = {
    settings: alloc(8 * spec.settings.length),
    triangles: alloc(4 * spec.triangle_floats * CAPACITY),
    numbers: alloc(8 * nnumbers),
    shadow: alloc(SHADOW * SHADOW),
  }

  // Build animal `index` (0-based, as in spec.animals) from `settings`, keyed like spec.settings, with the sun at
  // `zenith` and `azimuth` (degrees). The memory can move when the module grows it, so every view is made afresh
  // and copied out.
  function run(index, settings, zenith, azimuth) {
    const view = (Type, address, length) => new Type(memory.buffer, Number(address), length)
    view(Float64Array, buffers.settings, spec.settings.length).set(spec.settings.map((k) => settings[k]))
    const nparts = N(w.animal_builder(BigInt(index + 1), buffers.settings, buffers.triangles, BigInt(CAPACITY),
      buffers.numbers, buffers.shadow, BigInt(SHADOW), (zenith * Math.PI) / 180, (azimuth * Math.PI) / 180))
    const x = view(Float64Array, buffers.numbers, nnumbers)
    const at = (k) => x[k - 1]
    const count = at(spec.numbers.triangle_count)
    const extent = [0, 1, 2, 3].map((i) => at(spec.numbers.shadow_extent + i))
    return {
      triangles: view(Float32Array, buffers.triangles, spec.triangle_floats * count).slice(),
      count,
      parts: spec.animals[index].parts.map((name, k) => ({
        name, area: at(spec.numbers.parts + 2 * k), mass: at(spec.numbers.parts + 2 * k + 1) })),
      total: at(spec.numbers.outer_area),
      skin: at(spec.numbers.skin_area),
      hidden: at(spec.numbers.hidden_area),
      volume: at(spec.numbers.volume),
      meeh: at(spec.numbers.meeh),
      shadow: { covered: view(Uint8Array, buffers.shadow, SHADOW * SHADOW).slice(), n: SHADOW,
                width: extent[1] - extent[0], height: extent[3] - extent[2], area: at(spec.numbers.shadow_area) },
      nparts,
    }
  }
  return { run }
}
