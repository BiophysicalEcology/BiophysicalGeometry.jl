<script setup>
// Build an animal: the page's controls, and BiophysicalGeometry.jl itself, compiled to wasm by `compile_wasm` in
// docs/make.jl, building the animal and computing its areas and shadow.
import { computed, onMounted, reactive, ref, shallowRef, watch } from 'vue'
import { loadBiophysicalModel } from './builder/biophysical.mjs'
import spec from './builder/animals.json'
import wasm from './builder/animals.wasm?url'

// Each animal's settings, as an object keyed by name.
const animals = spec.models.map((m) => ({
  name: m.name, parts: m.parts, settings: Object.fromEntries(m.settings.map((name, i) => [name, m.defaults[i]])) }))
const index = ref(0)                               // the animal, in spec.models
const p = reactive({ ...animals[0].settings })     // its settings
const logMass = ref(Math.log10(p.mass))
const sun = reactive({ zenith: 30, azimuth: 90 })
const view = reactive({ azimuth: -0.9, elevation: 0.35 })
const builder = shallowRef(null)
const animal = shallowRef(null)

watch(logMass, (x) => { p.mass = Number(Math.pow(10, x).toPrecision(3)) })
function useAnimal() {
  Object.assign(p, animals[index.value].settings)
  logMass.value = Math.log10(p.mass)
}

// The parts this animal has, by kind.
const parts = computed(() => animals[index.value].parts)
const has = (prefix) => parts.value.some((name) => name.startsWith(prefix))
const upright = computed(() => animals[index.value].name === 'Human')
const hindLegs = computed(() => ['Kangaroo', 'Tyrannosaur'].includes(animals[index.value].name))
const plateEars = computed(() => ['Elephant', 'Kangaroo', 'Giraffe'].includes(animals[index.value].name))
const sphereHead = computed(() => ['Bird', 'Seal'].includes(animals[index.value].name))

function update() {
  if (!builder.value) return
  const z = (sun.zenith * Math.PI) / 180, a = (sun.azimuth * Math.PI) / 180
  const result = builder.value.run(index.value, p, [Math.sin(z) * Math.cos(a), Math.sin(z) * Math.sin(a), Math.cos(z)])
  animal.value = { ...result, meeh: result.total / Math.pow(p.mass, 2 / 3) }
}

const colours = { dorsal: [76, 140, 191], ventral: [230, 158, 51], head: [89, 173, 115], neck: [148, 115, 184],
                  nose: [140, 107, 89], beak: [217, 166, 64], tail: [128, 128, 128], ear: [217, 140, 191], wing: [100, 170, 180], arm: [160, 190, 90] }
const legColour = [204, 102, 115]
const colourOf = (part) => colours[part] || (part.startsWith('ear') ? colours.ear : part.startsWith('wing') ? colours.wing : part.startsWith('arm') ? colours.arm : legColour)

const fmt = (x, unit) => {
  const s = x >= 100 ? x.toFixed(0) : x >= 10 ? x.toFixed(1) : x >= 1 ? x.toFixed(2) : x.toPrecision(3)
  return `${s} ${unit}`
}
const area = (x) => (x >= 0.1 ? fmt(x, 'm²') : fmt(x * 1e4, 'cm²'))
const weight = (x) => (x >= 1 ? fmt(x, 'kg') : fmt(x * 1000, 'g'))
const massLabel = computed(() => weight(p.mass))
const rows = computed(() => {
  const groups = {}
  for (const part of animal.value.parts) {
    const key = part.name.startsWith('leg') ? 'legs' : part.name.startsWith('ear') ? 'ears' : part.name.startsWith('wing') ? 'wings' : part.name.startsWith('arm') ? 'arms' : part.name
    groups[key] = groups[key] || { area: 0, mass: 0 }
    groups[key].area += part.area
    groups[key].mass += part.mass
  }
  return Object.entries(groups)
})

// ── Drawing ─────────────────────────────────────────────────────────────────

const canvas = ref(null)
const shadowCanvas = ref(null)

// What drawing needs of the animal, worked out once per build rather than once per frame: each triangle's corners
// about the animal's centre, its unit normal and colour, and the size of the whole.
const mesh = computed(() => {
  if (!animal.value) return null
  const { triangles: t, count } = animal.value
  const F = spec.triangle_floats
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]
  for (let i = 0; i < count; i++) for (let c = 0; c < 3; c++) for (let k = 0; k < 3; k++) {
    const x = t[F * i + 3 * c + k]
    lo[k] = Math.min(lo[k], x); hi[k] = Math.max(hi[k], x)
  }
  const centre = lo.map((x, k) => (x + hi[k]) / 2)
  const corners = new Float32Array(9 * count), normals = new Float32Array(3 * count), keep = []
  const colour = animal.value.parts.map((part) => colourOf(part.name)), colours = []
  for (let i = 0; i < count; i++) {
    for (let c = 0; c < 3; c++) for (let k = 0; k < 3; k++) corners[9 * i + 3 * c + k] = t[F * i + 3 * c + k] - centre[k]
    const q = (c, k) => corners[9 * i + 3 * c + k]
    const e1 = [0, 1, 2].map((k) => q(1, k) - q(0, k)), e2 = [0, 1, 2].map((k) => q(2, k) - q(0, k))
    const n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
    const length = Math.hypot(...n)
    if (length === 0) continue   // the degenerate half of a triangular cell
    for (let k = 0; k < 3; k++) normals[3 * i + k] = n[k] / length
    keep.push(i)
    colours[i] = colour[t[F * i + F - 1]]
  }
  return { corners, normals, colours, keep: Uint32Array.from(keep),
           size: Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) }
})

function draw() {
  const el = canvas.value
  if (!el || !mesh.value) return
  const ctx = el.getContext('2d')
  const w = el.width, h = el.height
  ctx.clearRect(0, 0, w, h)
  const ca = Math.cos(view.azimuth), sa = Math.sin(view.azimuth), ce = Math.cos(view.elevation), se = Math.sin(view.elevation)
  // towards the viewer, to the right and up on the screen
  const toward = [ce * ca, ce * sa, se], right = [-sa, ca, 0], up = [-se * ca, -se * sa, ce]
  const { corners: q, normals: n, colours, keep, size } = mesh.value
  const scale = (0.92 * Math.min(w, h)) / size
  const light = [0, 1, 2].map((k) => toward[k] + 0.35 * right[k] + 0.6 * up[k])
  const lightLength = Math.hypot(...light)
  const depth = new Float32Array(q.length / 9)
  for (const i of keep) {
    let d = 0
    for (let c = 0; c < 3; c++) d += q[9 * i + 3 * c] * toward[0] + q[9 * i + 3 * c + 1] * toward[1] + q[9 * i + 3 * c + 2] * toward[2]
    depth[i] = d
  }
  const order = keep.slice().sort((a, b) => depth[a] - depth[b])
  const sx = (i, c) => w / 2 + scale * (q[9 * i + 3 * c] * right[0] + q[9 * i + 3 * c + 1] * right[1] + q[9 * i + 3 * c + 2] * right[2])
  const sy = (i, c) => h / 2 - scale * (q[9 * i + 3 * c] * up[0] + q[9 * i + 3 * c + 1] * up[1] + q[9 * i + 3 * c + 2] * up[2])
  for (const i of order) {
    const shade = 0.55 + 0.45 * Math.abs(n[3 * i] * light[0] + n[3 * i + 1] * light[1] + n[3 * i + 2] * light[2]) / lightLength
    const c = colours[i]
    ctx.fillStyle = ctx.strokeStyle = `rgb(${Math.round(c[0] * shade)},${Math.round(c[1] * shade)},${Math.round(c[2] * shade)})`
    ctx.beginPath()
    ctx.moveTo(sx(i, 0), sy(i, 0)); ctx.lineTo(sx(i, 1), sy(i, 1)); ctx.lineTo(sx(i, 2), sy(i, 2))
    ctx.closePath(); ctx.fill(); ctx.stroke()
  }
  // scale bar
  const nice = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5].filter((x) => x * scale < w / 3).pop() || 0.001
  ctx.strokeStyle = ctx.fillStyle = getComputedStyle(el).color
  ctx.lineWidth = 2
  ctx.beginPath(); ctx.moveTo(16, h - 16); ctx.lineTo(16 + nice * scale, h - 16); ctx.stroke()
  ctx.lineWidth = 1
  ctx.font = '17px sans-serif'
  ctx.fillText(nice >= 1 ? `${nice} m` : `${Math.round(nice * 1000) / 10} cm`, 16, h - 26)
}

function drawShadow() {
  const el = shadowCanvas.value
  if (!el || !animal.value) return
  const ctx = el.getContext('2d')
  const { covered, n, width, height } = animal.value.shadow
  ctx.clearRect(0, 0, el.width, el.height)
  const s = Math.min(el.width / width, el.height / height)
  const cw = (width / n) * s, ch = (height / n) * s
  const x0 = (el.width - width * s) / 2, y0 = (el.height - height * s) / 2
  ctx.fillStyle = getComputedStyle(el).color
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
    if (covered[j * n + i]) ctx.fillRect(x0 + i * cw, el.height - y0 - (j + 1) * ch, cw + 0.5, ch + 0.5)
  }
}

let dragging = null
function down(e) { dragging = [e.clientX, e.clientY]; e.target.setPointerCapture(e.pointerId) }
function move(e) {
  if (!dragging) return
  view.azimuth -= (e.clientX - dragging[0]) * 0.01
  view.elevation = Math.max(-1.5, Math.min(1.5, view.elevation + (e.clientY - dragging[1]) * 0.01))
  dragging = [e.clientX, e.clientY]
}
function release() { dragging = null }

onMounted(async () => {
  // a link such as builder?animal=Bird opens on that animal, and any setting can follow, as in
  // builder?animal=Kangaroo&pitch=20
  const query = new URLSearchParams(window.location.search)
  const wanted = animals.findIndex((a) => a.name === query.get('animal'))
  if (wanted >= 0) { index.value = wanted; useAnimal() }
  for (const [key, value] of query) if (key in p) p[key] = Number(value)
  logMass.value = Math.log10(p.mass)
  builder.value = await loadBiophysicalModel(wasm, spec)
  update()
})
watch([p, sun, index], update, { deep: true })
watch([mesh, view], draw, { deep: true })
// after the DOM updates: the shadow's canvas only appears with the first result
watch(animal, drawShadow, { flush: 'post' })
</script>

<template>
  <div class="builder">
    <div class="panes">
    <div class="controls">
      <label class="wide">Animal
        <select v-model.number="index" @change="useAnimal">
          <option v-for="(a, i) in animals" :key="a.name" :value="i">{{ a.name }}</option>
        </select>
      </label>

      <h4>Body</h4>
      <label>Mass <output>{{ massLabel }}</output>
        <input type="range" min="-2" max="4" step="0.01" v-model.number="logMass" /></label>
      <label v-if="!upright">Body pitch, head up <output>{{ p.pitch }}°</output>
        <input type="range" min="-30" max="80" step="1" v-model.number="p.pitch" /></label>
      <label>Torso length / width <output>{{ p.torsoRatio }}</output>
        <input type="range" min="1.1" max="6" step="0.1" v-model.number="p.torsoRatio" /></label>
      <label>Fat, fraction of torso mass <output>{{ p.fat }}</output>
        <input type="range" min="0" max="0.5" step="0.01" v-model.number="p.fat" /></label>

      <h4>Coat</h4>
      <label>Depth on the back <output>{{ (p.backFur * 1000).toFixed(1) }} mm</output>
        <input type="range" min="0" max="0.06" step="0.0005" v-model.number="p.backFur" /></label>
      <label>Depth on the belly <output>{{ (p.bellyFur * 1000).toFixed(1) }} mm</output>
        <input type="range" min="0" max="0.06" step="0.0005" v-model.number="p.bellyFur" /></label>
      <label>Depth on head and legs <output>{{ (p.limbFur * 1000).toFixed(1) }} mm</output>
        <input type="range" min="0" max="0.06" step="0.0005" v-model.number="p.limbFur" /></label>

      <h4>Head</h4>
      <label>Fraction of mass <output>{{ p.headFraction }}</output>
        <input type="range" min="0.01" max="0.25" step="0.01" v-model.number="p.headFraction" /></label>
      <label v-if="!sphereHead">Head length / width <output>{{ p.headRatio }}</output>
        <input type="range" min="1" max="3" step="0.1" v-model.number="p.headRatio" /></label>
      <template v-if="has('neck')">
        <label>Neck, fraction of mass <output>{{ p.neckFraction }}</output>
          <input type="range" min="0.01" max="0.15" step="0.005" v-model.number="p.neckFraction" /></label>
        <label>Neck length / width <output>{{ p.neckRatio }}</output>
          <input type="range" min="0.5" max="8" step="0.1" v-model.number="p.neckRatio" /></label>
        <label>Neck angle, up <output>{{ p.neckAngle }}°</output>
          <input type="range" min="-30" max="90" step="1" v-model.number="p.neckAngle" /></label>
      </template>
      <label v-if="has('nose') || has('beak')">{{ has('beak') ? 'Beak' : 'Nose' }}, fraction of mass <output>{{ p.noseFraction }}</output>
        <input type="range" min="0.001" max="0.03" step="0.001" v-model.number="p.noseFraction" /></label>
      <template v-if="has('ear')">
        <label>Ear, fraction of mass, each <output>{{ p.earFraction }}</output>
          <input type="range" min="0.0005" max="0.03" step="0.0005" v-model.number="p.earFraction" /></label>
        <label>Ears laid back <output>{{ p.earAngle }}°</output>
          <input type="range" min="0" max="90" step="1" v-model.number="p.earAngle" /></label>
        <template v-if="plateEars">
          <label>Ear length / width <output>{{ p.earRatio }}</output>
            <input type="range" min="0.5" max="4" step="0.1" v-model.number="p.earRatio" /></label>
          <label>Ear length / thickness <output>{{ p.earFlatness }}</output>
            <input type="range" min="3" max="40" step="1" v-model.number="p.earFlatness" /></label>
        </template>
      </template>

      <template v-if="has('leg')">
        <h4>Legs</h4>
        <label>Fraction of mass, each <output>{{ p.legFraction }}</output>
          <input type="range" min="0.001" max="0.2" step="0.001" v-model.number="p.legFraction" /></label>
        <label>Length / width <output>{{ p.legRatio }}</output>
          <input type="range" min="1" max="12" step="0.1" v-model.number="p.legRatio" /></label>
        <label>Taper, foot / top <output>{{ p.legTop >= 1 ? 'cylinder' : p.legTop }}</output>
          <input type="range" min="0.1" max="1" step="0.05" v-model.number="p.legTop" /></label>
        <label>Swung forward <output>{{ p.legAngle }}°</output>
          <input type="range" min="-60" max="60" step="1" v-model.number="p.legAngle" /></label>
        <template v-if="hindLegs">
          <label>Hind leg, fraction of mass, each <output>{{ p.hindFraction }}</output>
            <input type="range" min="0.005" max="0.2" step="0.005" v-model.number="p.hindFraction" /></label>
          <label>Hind leg, length / width <output>{{ p.hindRatio }}</output>
            <input type="range" min="1" max="12" step="0.1" v-model.number="p.hindRatio" /></label>
        </template>
      </template>

      <template v-if="has('arm')">
        <h4>Arms</h4>
        <label>Fraction of mass, each <output>{{ p.armFraction }}</output>
          <input type="range" min="0.005" max="0.1" step="0.005" v-model.number="p.armFraction" /></label>
        <label>Length / width <output>{{ p.armRatio }}</output>
          <input type="range" min="2" max="16" step="0.5" v-model.number="p.armRatio" /></label>
      </template>

      <template v-if="has('wing')">
        <h4>Wings</h4>
        <label>Fraction of mass, each <output>{{ p.wingFraction }}</output>
          <input type="range" min="0.005" max="0.15" step="0.005" v-model.number="p.wingFraction" /></label>
        <label>Folded back <output>{{ p.wingFold }}°</output>
          <input type="range" min="0" max="90" step="1" v-model.number="p.wingFold" /></label>
      </template>

      <template v-if="has('tail')">
        <h4>Tail</h4>
        <label>Fraction of mass <output>{{ p.tailFraction }}</output>
          <input type="range" min="0.001" max="0.25" step="0.001" v-model.number="p.tailFraction" /></label>
        <label>Length / width <output>{{ p.tailRatio }}</output>
          <input type="range" min="1" max="20" step="0.5" v-model.number="p.tailRatio" /></label>
        <label>Tail angle, up <output>{{ p.tailAngle }}°</output>
          <input type="range" min="-60" max="90" step="1" v-model.number="p.tailAngle" /></label>
      </template>

      <h4>Sun</h4>
      <label>Zenith angle <output>{{ sun.zenith }}°</output>
        <input type="range" min="0" max="90" step="1" v-model.number="sun.zenith" /></label>
      <label>Azimuth, 0° is head-on <output>{{ sun.azimuth }}°</output>
        <input type="range" min="0" max="360" step="5" v-model.number="sun.azimuth" /></label>
    </div>

    <div class="result">
      <canvas ref="canvas" class="body" width="560" height="400" @pointerdown="down" @pointermove="move"
              @pointerup="release" @pointerleave="release"></canvas>
      <p class="hint">Drag to turn the animal.</p>
      <div v-if="animal" class="numbers">
        <table>
          <tr><th>Outer area</th><td>{{ area(animal.total) }}</td></tr>
          <tr><th>Skin area</th><td>{{ area(animal.skin) }}</td></tr>
          <tr><th>Hidden by joins</th><td>{{ area(animal.hidden) }}</td></tr>
          <tr><th>Volume</th><td>{{ fmt(animal.volume * 1000, 'L') }}</td></tr>
          <tr><th>Meeh coefficient</th><td>{{ animal.meeh.toFixed(3) }}</td></tr>
        </table>
        <table>
          <tr><th>Part</th><th class="right">Exposed area</th><th class="right">Mass</th></tr>
          <tr v-for="[name, value] in rows" :key="name">
            <td class="name">{{ name }}</td><td>{{ area(value.area) }}</td><td>{{ weight(value.mass) }}</td></tr>
          <tr><th>Whole animal</th><td>{{ area(animal.total) }}</td><td>{{ weight(p.mass) }}</td></tr>
        </table>
        <div class="shadow">
          <canvas ref="shadowCanvas" width="200" height="150"></canvas>
          <p>Silhouette to the sun<br /><strong>{{ area(animal.shadow.area) }}</strong></p>
        </div>
      </div>
      <p v-else class="hint">Loading the model…</p>
    </div>
    </div>
  </div>
</template>

<style scoped>
/* The controls scroll in their own column, so that the animal stays in view while any slider is moved. */
.builder { margin: 16px 0; --pane: calc(100vh - var(--vp-nav-height, 64px) - 40px); }
.panes { display: grid; grid-template-columns: minmax(230px, 280px) 1fr; gap: 20px; }
.controls { display: flex; flex-direction: column; gap: 6px; font-size: 13px; position: sticky;
            top: calc(var(--vp-nav-height, 64px) + 20px); align-self: start; max-height: var(--pane); overflow-y: auto;
            padding-right: 10px; }
.controls h4:first-of-type { margin-top: 4px; }
.controls h4 { margin: 10px 0 0; font-size: 13px; text-transform: uppercase; letter-spacing: 0.05em;
               color: var(--vp-c-text-2); }
.controls label { display: grid; grid-template-columns: 1fr auto; align-items: center; gap: 2px 8px; }
.controls input[type='range'] { grid-column: 1 / -1; width: 100%; accent-color: var(--vp-c-brand-3); }
.controls select { border: 1px solid var(--vp-c-divider); border-radius: 6px; padding: 2px 6px;
                   background: var(--vp-c-bg-soft); color: var(--vp-c-text-1); }
.controls output { font-variant-numeric: tabular-nums; color: var(--vp-c-text-2); }
.result { min-width: 0; position: sticky; top: calc(var(--vp-nav-height, 64px) + 20px); align-self: start; }
canvas.body { width: min(100%, calc(0.5 * var(--pane) * 1.4)); height: auto; border: 1px solid var(--vp-c-divider); border-radius: 8px;
              cursor: grab; touch-action: none; color: var(--vp-c-text-2); background: var(--vp-c-bg-soft); }
.hint { margin: 2px 0 8px; font-size: 12px; color: var(--vp-c-text-3); }
.numbers { display: flex; flex-wrap: wrap; gap: 20px; align-items: flex-start; }
.numbers table { margin: 0; font-size: 13px; display: table; width: auto; }
.numbers th { text-align: left; font-weight: 500; padding: 3px 12px 3px 8px; }
.numbers th.right { text-align: right; }
.numbers td.name { text-align: left; }
.numbers td { text-align: right; font-variant-numeric: tabular-nums; padding: 3px 8px; }
.shadow { text-align: center; font-size: 13px; }
.shadow canvas { color: var(--vp-c-text-1); border: 1px solid var(--vp-c-divider); border-radius: 8px;
                 background: var(--vp-c-bg-soft); }
.shadow p { margin: 4px 0 0; line-height: 1.4; }
@media (max-width: 720px) {
  .panes { grid-template-columns: 1fr; }
  .controls { position: static; max-height: none; overflow: visible; order: 2; }
  .result { top: var(--vp-nav-height, 64px); z-index: 2; background: var(--vp-c-bg); order: 1; }
  canvas.body { width: min(100%, 56vh); }
  .numbers { display: none; }
}
</style>
