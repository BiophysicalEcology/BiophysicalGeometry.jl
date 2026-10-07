// Run every animal of the builder's wasm module under node and print its numbers:
//     node docs/builder/test_node.mjs
import fs from 'node:fs';
const here = new URL('.', import.meta.url).pathname;
const spec = JSON.parse(fs.readFileSync(here + '../src/components/animals.json'));
const bytes = fs.readFileSync(here + '../src/public/animal_builder.wasm');
let mem;
const N = (b) => typeof b === 'bigint' ? Number(b) : b;
const env = {
  js_console(level, ptr, len) { console.error('[wasm]', new TextDecoder().decode(new Uint8Array(mem.buffer, N(ptr), N(len)))); },
  memset(p, v, n) { new Uint8Array(mem.buffer).fill(Number(v) & 0xff, N(p), N(p) + N(n)); return p; },
  memcpy(d, s, n) { new Uint8Array(mem.buffer).copyWithin(N(d), N(s), N(s) + N(n)); return d; },
  memmove(d, s, n) { new Uint8Array(mem.buffer).copyWithin(N(d), N(s), N(s) + N(n)); return d; },
};
const { instance } = await WebAssembly.instantiate(bytes, { env });
const w = instance.exports;
mem = w.memory;
const alloc = (bytes) => { const p = w.whisk_gc_alloc(BigInt(bytes)); w.whisk_pin(p); return p; };
const capacity = 60000, n = 160;
const buffers = {
  settings: alloc(8 * spec.settings.length), triangles: alloc(4 * spec.triangle_floats * capacity),
  numbers: alloc(8 * (spec.numbers.parts + 64)), shadow: alloc(n * n),
};
spec.animals.forEach((animal, i) => {
  new Float64Array(mem.buffer, Number(buffers.settings), spec.settings.length).set(spec.settings.map((k) => animal.settings[k]));
  const t = performance.now();
  const nparts = w.animal_builder(BigInt(i + 1), buffers.settings, buffers.triangles, BigInt(capacity), buffers.numbers,
                                  buffers.shadow, BigInt(n), 30 * Math.PI / 180, 90 * Math.PI / 180);
  const ms = (performance.now() - t).toFixed(1);
  const x = new Float64Array(mem.buffer, Number(buffers.numbers), spec.numbers.parts + 64);
  const at = (k) => x[k - 1];
  console.log(animal.name.padEnd(12), 'parts', N(nparts), 'triangles', at(spec.numbers.triangle_count),
              'outer', at(spec.numbers.outer_area).toFixed(6), 'shadow', at(spec.numbers.shadow_area).toFixed(6), ms, 'ms');
});
