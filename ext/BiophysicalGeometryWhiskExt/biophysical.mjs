// Loads models compiled by BiophysicalGeometry.jl's `compile_wasm`. No framework: use it from Vue, Bonito or plain
// HTML.
//
//     import { loadBiophysicalModel } from './biophysical.mjs'
//     const model = await loadBiophysicalModel('model.wasm', spec)    // spec: the contents of model.json
//     const result = model.run(0, { mass: 20 }, [0, 0.5, 1])
//
// `loadBiophysicalModel` takes the module as a URL, an ArrayBuffer, a typed array or a fetch Response. `run(model, settings, sun)`
// builds model `model` (its number from 0, or its name) from `settings`, keyed by name, missing ones taken from its
// defaults, with the sun in the direction `sun`, [x, y, z]. It returns the body's triangles to draw (Float32Array,
// `spec.triangle_floats` each: three corners in metres, then the part's number), each part's exposed area (m²) and
// mass (kg), the totals, and the shadow: an n × n grid of 0 and 1, row j from the bottom and column i from the left
// at index j * n + i, over `width` × `height` metres.

export async function loadBiophysicalModel(source, spec, { capacity = 80000, shadow = 160 } = {}) {
  const bytes = await asBytes(source)
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
  const nsettings = Math.max(...spec.models.map((m) => m.settings.length))
  const nnumbers = spec.numbers.parts - 1 + 2 * spec.max_parts
  const buffers = {
    settings: alloc(8 * nsettings),
    triangles: alloc(4 * spec.triangle_floats * capacity),
    numbers: alloc(8 * nnumbers),
    shadow: alloc(shadow * shadow),
  }

  function run(model, settings = {}, sun = [0, 0, 1]) {
    const index = typeof model === 'number' ? model : spec.models.findIndex((m) => m.name === model)
    const m = spec.models[index]
    if (!m) throw new Error(`no model ${model}`)
    // The memory can move when the module grows it, so every view is made afresh, and results are copied out.
    const view = (Type, address, length) => new Type(memory.buffer, Number(address), length)
    view(Float64Array, buffers.settings, m.settings.length)
      .set(m.settings.map((name, i) => (name in settings ? settings[name] : m.defaults[i])))
    const nparts = N(w.biophysical_run(BigInt(index + 1), buffers.settings, buffers.triangles, BigInt(capacity),
      buffers.numbers, buffers.shadow, BigInt(shadow), sun[0], sun[1], sun[2]))
    const x = view(Float64Array, buffers.numbers, nnumbers)
    const at = (k) => x[k - 1]
    const count = at(spec.numbers.triangle_count)
    const extent = [0, 1, 2, 3].map((i) => at(spec.numbers.shadow_extent + i))
    return {
      triangles: view(Float32Array, buffers.triangles, spec.triangle_floats * count).slice(),
      count,
      parts: m.parts.slice(0, nparts).map((name, k) => ({
        name, area: at(spec.numbers.parts + 2 * k), mass: at(spec.numbers.parts + 2 * k + 1) })),
      total: at(spec.numbers.outer_area),
      skin: at(spec.numbers.skin_area),
      hidden: at(spec.numbers.hidden_area),
      volume: at(spec.numbers.volume),
      shadow: { covered: view(Uint8Array, buffers.shadow, shadow * shadow).slice(), n: shadow,
                width: extent[1] - extent[0], height: extent[3] - extent[2], area: at(spec.numbers.shadow_area) },
    }
  }
  return { run, spec }
}

async function asBytes(source) {
  if (typeof source === 'string' || source instanceof URL) source = await fetch(source)
  if (typeof Response !== 'undefined' && source instanceof Response) return source.arrayBuffer()
  return source
}
