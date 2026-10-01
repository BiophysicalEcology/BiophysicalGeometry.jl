<script setup>
import { computed, onMounted, reactive, ref, watch } from 'vue'
import { build, defaults, silhouette } from './animalGeometry.mjs'

const presets = {
  Dog: { ...defaults, neck: true, nose: 'Nose', ears: 'Cone', tail: true },
  Mouse: { ...defaults, mass: 0.02, torsoShape: 'Ellipsoid', torsoRatio: 2, fat: 0.05, backFur: 0.004, bellyFur: 0.003,
           headFraction: 0.12, legFraction: 0.02, legRatio: 4, legTop: 1, limbFur: 0.002, nose: 'Nose', ears: 'Cone',
           earFraction: 0.004, tail: true, tailFraction: 0.01, tailRatio: 12 },
  Elephant: { ...defaults, mass: 4000, torsoRatio: 1.8, fat: 0, backFur: 0, bellyFur: 0, limbFur: 0, headShape: 'Ellipsoid',
              headRatio: 1.3, headFraction: 0.08, nose: 'Nose', noseFraction: 0.02, ears: 'Plate', earFraction: 0.01,
              earPosture: 'Up', earRatio: 1.2, earFlatness: 30, legFraction: 0.04, legRatio: 3.5, legTop: 0.8, tail: true,
              tailFraction: 0.002, tailRatio: 12 },
  Giraffe: { ...defaults, mass: 800, torsoRatio: 1.8, fat: 0.02, backFur: 0.003, bellyFur: 0.003, limbFur: 0.003,
             headRatio: 2, headFraction: 0.02, neck: true, neckFraction: 0.1, neckRatio: 6, neckPosture: 'Up',
             ears: 'Plate', earFraction: 0.0005, earRatio: 2, earFlatness: 12, legFraction: 0.05, legRatio: 10, legTop: 0.5,
             tail: true, tailFraction: 0.002, tailRatio: 15 },
  Cow: { ...defaults, mass: 682, torsoRatio: 2.4, fat: 0, backFur: 0.0034, bellyFur: 0.0034, headShape: 'Ellipsoid',
         headRatio: 1.8, headFraction: 0.04, legFraction: 0.02, legRatio: 4.3, legTop: 0.5, limbFur: 0.0034, neck: true,
         neckFraction: 0.05, ears: 'Cone', earFraction: 0.001, tail: true, tailFraction: 0.003, tailRatio: 12 },
  Bird: { ...defaults, mass: 0.05, torsoShape: 'Ellipsoid', torsoRatio: 1.6, fat: 0.05, backFur: 0.006, bellyFur: 0.006,
          headShape: 'Sphere', headFraction: 0.1, legs: 2, legFraction: 0.02, legRatio: 8, legTop: 1, limbFur: 0,
          neck: true, neckFraction: 0.03, nose: 'Beak', noseFraction: 0.01, tail: true, tailFraction: 0.02, tailRatio: 3,
          wings: 'Folded', wingFraction: 0.06 },
  Seal: { ...defaults, nose: 'Nose', tail: true, tailFraction: 0.02, tailRatio: 2, mass: 100, torsoShape: 'Ellipsoid', torsoRatio: 4, fat: 0.35, backFur: 0.003, bellyFur: 0.003,
          headShape: 'Sphere', headFraction: 0.05, legs: 0, limbFur: 0.003 },
}

const p = reactive({ ...defaults })
const logMass = ref(Math.log10(defaults.mass))
const sun = reactive({ zenith: 30, azimuth: 90 })
const view = reactive({ azimuth: -0.9, elevation: 0.35 })
const preset = ref('Dog')

watch(logMass, (x) => { p.mass = Number(Math.pow(10, x).toPrecision(3)) })
function usePreset() {
  Object.assign(p, defaults, presets[preset.value])
  logMass.value = Math.log10(p.mass)
}

const animal = computed(() => build(p))
const sunDirection = computed(() => {
  const z = (sun.zenith * Math.PI) / 180, a = (sun.azimuth * Math.PI) / 180
  return [Math.sin(z) * Math.cos(a), Math.sin(z) * Math.sin(a), Math.cos(z)]
})
const shadow = computed(() => silhouette(animal.value.triangles, sunDirection.value, 160))

const colours = { dorsal: [76, 140, 191], ventral: [230, 158, 51], head: [89, 173, 115], neck: [148, 115, 184],
                  nose: [140, 107, 89], beak: [217, 166, 64], tail: [128, 128, 128], ear: [217, 140, 191], wing: [100, 170, 180] }
const legColour = [204, 102, 115]
const colourOf = (part) => colours[part] || (part.startsWith('ear') ? colours.ear : part.startsWith('wing') ? colours.wing : legColour)

const fmt = (x, unit) => {
  const s = x >= 100 ? x.toFixed(0) : x >= 10 ? x.toFixed(1) : x >= 1 ? x.toFixed(2) : x.toPrecision(3)
  return `${s} ${unit}`
}
const area = (x) => (x >= 0.1 ? fmt(x, 'm²') : fmt(x * 1e4, 'cm²'))
const massLabel = computed(() => (p.mass >= 1 ? fmt(p.mass, 'kg') : fmt(p.mass * 1000, 'g')))
const rows = computed(() => {
  const groups = {}
  for (const part of animal.value.parts) {
    const key = part.name.startsWith('leg') ? 'legs' : part.name.startsWith('ear') ? 'ears' : part.name.startsWith('wing') ? 'wings' : part.name
    groups[key] = (groups[key] || 0) + part.total - part.hidden
  }
  return Object.entries(groups)
})

// ── Drawing ─────────────────────────────────────────────────────────────────

const canvas = ref(null)
const shadowCanvas = ref(null)

function draw() {
  const el = canvas.value
  if (!el) return
  const ctx = el.getContext('2d')
  const w = el.width, h = el.height
  ctx.clearRect(0, 0, w, h)
  const ca = Math.cos(view.azimuth), sa = Math.sin(view.azimuth), ce = Math.cos(view.elevation), se = Math.sin(view.elevation)
  // towards the viewer, to the right and up on the screen
  const toward = [ce * ca, ce * sa, se], right = [-sa, ca, 0], up = [-se * ca, -se * sa, ce]
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
  const tris = animal.value.triangles
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]
  for (const t of tris) for (const q of t.pts) for (let k = 0; k < 3; k++) {
    lo[k] = Math.min(lo[k], q[k]); hi[k] = Math.max(hi[k], q[k])
  }
  const centre = lo.map((x, k) => (x + hi[k]) / 2)
  const size = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])
  const scale = (0.92 * Math.min(w, h)) / size
  const light = [toward[0] + 0.35 * right[0] + 0.6 * up[0], toward[1] + 0.35 * right[1] + 0.6 * up[1], toward[2] + 0.35 * right[2] + 0.6 * up[2]]
  const lightLength = Math.hypot(...light)
  const drawn = tris.map((t) => {
    const rel = t.pts.map((q) => [q[0] - centre[0], q[1] - centre[1], q[2] - centre[2]])
    const e1 = rel[1].map((x, k) => x - rel[0][k]), e2 = rel[2].map((x, k) => x - rel[0][k])
    const n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
    const length = Math.hypot(...n) || 1
    return {
      depth: (dot(rel[0], toward) + dot(rel[1], toward) + dot(rel[2], toward)) / 3,
      xy: rel.map((q) => [w / 2 + scale * dot(q, right), h / 2 - scale * dot(q, up)]),
      shade: 0.55 + 0.45 * Math.abs(dot(n, light)) / (length * lightLength),
      colour: colourOf(t.part),
    }
  })
  drawn.sort((a, b) => a.depth - b.depth)
  for (const t of drawn) {
    const [r, g, b] = t.colour.map((c) => Math.round(c * t.shade))
    ctx.fillStyle = ctx.strokeStyle = `rgb(${r},${g},${b})`
    ctx.beginPath()
    ctx.moveTo(t.xy[0][0], t.xy[0][1]); ctx.lineTo(t.xy[1][0], t.xy[1][1]); ctx.lineTo(t.xy[2][0], t.xy[2][1])
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
  if (!el) return
  const ctx = el.getContext('2d')
  const { covered, n, width, height } = shadow.value
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

const copied = ref(false)
async function copyCode() {
  await navigator.clipboard.writeText(animal.value.code)
  copied.value = true
  setTimeout(() => (copied.value = false), 1500)
}

onMounted(() => {
  // a link such as builder?preset=Bird opens on that animal
  const wanted = new URLSearchParams(window.location.search).get('preset')
  if (wanted in presets) { preset.value = wanted; usePreset() }
  draw(); drawShadow()
})
watch([animal, view], draw, { deep: true })
watch(shadow, drawShadow)
</script>

<template>
  <div class="builder">
    <div class="controls">
      <label class="wide">Start from
        <select v-model="preset" @change="usePreset">
          <option v-for="(_, name) in presets" :key="name">{{ name }}</option>
        </select>
      </label>

      <h4>Body</h4>
      <label>Mass <output>{{ massLabel }}</output>
        <input type="range" min="-2" max="3" step="0.01" v-model.number="logMass" /></label>
      <label>Torso shape
        <select v-model="p.torsoShape"><option>Cylinder</option><option>Ellipsoid</option></select></label>
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
      <label>Shape
        <select v-model="p.headShape"><option>None</option><option>Sphere</option><option>Ellipsoid</option></select></label>
      <template v-if="p.headShape !== 'None'">
        <label>Fraction of mass <output>{{ p.headFraction }}</output>
          <input type="range" min="0.01" max="0.25" step="0.01" v-model.number="p.headFraction" /></label>
        <label class="check"><span><input type="checkbox" v-model="p.neck" /> Neck</span>
          <output v-if="p.neck">{{ p.neckFraction }}</output>
          <input v-if="p.neck" type="range" min="0.01" max="0.15" step="0.005" v-model.number="p.neckFraction" /></label>
        <template v-if="p.neck">
          <label>Neck length / width <output>{{ p.neckRatio }}</output>
            <input type="range" min="0.5" max="8" step="0.1" v-model.number="p.neckRatio" /></label>
          <label>Neck posture
            <select v-model="p.neckPosture"><option>Forward</option><option>Up</option></select></label>
        </template>
        <label>Nose or beak
          <select v-model="p.nose"><option>None</option><option>Nose</option><option>Beak</option></select></label>
        <label v-if="p.nose !== 'None'">Fraction of mass <output>{{ p.noseFraction }}</output>
          <input type="range" min="0.001" max="0.03" step="0.001" v-model.number="p.noseFraction" /></label>
        <label>Ears
          <select v-model="p.ears"><option>None</option><option>Cone</option><option>Plate</option></select></label>
        <template v-if="p.ears !== 'None'">
          <label>Fraction of mass, each <output>{{ p.earFraction }}</output>
            <input type="range" min="0.0005" max="0.03" step="0.0005" v-model.number="p.earFraction" /></label>
          <template v-if="p.ears === 'Plate'">
            <label>Ear posture
              <select v-model="p.earPosture"><option>Up</option><option>Flat</option></select></label>
            <label>Ear length / width <output>{{ p.earRatio }}</output>
              <input type="range" min="0.5" max="4" step="0.1" v-model.number="p.earRatio" /></label>
            <label>Ear length / thickness <output>{{ p.earFlatness }}</output>
              <input type="range" min="3" max="40" step="1" v-model.number="p.earFlatness" /></label>
          </template>
        </template>
      </template>

      <h4>Legs</h4>
      <label>Number
        <select v-model.number="p.legs"><option :value="0">0</option><option :value="2">2</option><option :value="4">4</option></select></label>
      <template v-if="p.legs > 0">
        <label>Fraction of mass, each <output>{{ p.legFraction }}</output>
          <input type="range" min="0.005" max="0.1" step="0.005" v-model.number="p.legFraction" /></label>
        <label>Length / width <output>{{ p.legRatio }}</output>
          <input type="range" min="1" max="12" step="0.1" v-model.number="p.legRatio" /></label>
        <label>Taper, foot / top <output>{{ p.legTop >= 1 ? 'cylinder' : p.legTop }}</output>
          <input type="range" min="0.1" max="1" step="0.05" v-model.number="p.legTop" /></label>
      </template>

      <h4>Wings</h4>
      <label>Wings
        <select v-model="p.wings"><option>None</option><option>Folded</option><option>Spread</option></select></label>
      <label v-if="p.wings !== 'None'">Fraction of mass, each <output>{{ p.wingFraction }}</output>
        <input type="range" min="0.005" max="0.15" step="0.005" v-model.number="p.wingFraction" /></label>

      <h4>Tail</h4>
      <label class="check"><span><input type="checkbox" v-model="p.tail" /> Tail</span>
        <output v-if="p.tail">{{ p.tailFraction }}</output>
        <input v-if="p.tail" type="range" min="0.001" max="0.1" step="0.001" v-model.number="p.tailFraction" /></label>
      <label v-if="p.tail">Length / width <output>{{ p.tailRatio }}</output>
        <input type="range" min="1" max="20" step="0.5" v-model.number="p.tailRatio" /></label>

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
      <div class="numbers">
        <table>
          <tr><th>Outer area</th><td>{{ area(animal.total) }}</td></tr>
          <tr><th>Skin area</th><td>{{ area(animal.skin) }}</td></tr>
          <tr><th>Hidden by joins</th><td>{{ area(2 * animal.joined) }}</td></tr>
          <tr><th>Volume</th><td>{{ fmt(animal.volume * 1000, 'L') }}</td></tr>
          <tr><th>Meeh coefficient</th><td>{{ animal.meeh.toFixed(3) }}</td></tr>
          <tr v-for="[name, value] in rows" :key="name"><th class="part">{{ name }}</th><td>{{ area(value) }}</td></tr>
        </table>
        <div class="shadow">
          <canvas ref="shadowCanvas" width="200" height="150"></canvas>
          <p>Silhouette to the sun<br /><strong>{{ area(shadow.area) }}</strong></p>
        </div>
      </div>
    </div>

    <div class="code">
      <h4>Julia code for this animal</h4>
      <button @click="copyCode">{{ copied ? 'Copied' : 'Copy' }}</button>
      <pre><code>{{ animal.code }}</code></pre>
    </div>
  </div>
</template>

<style scoped>
.builder { display: grid; grid-template-columns: minmax(230px, 280px) 1fr; gap: 20px; margin: 16px 0; }
.controls { display: flex; flex-direction: column; gap: 6px; font-size: 13px; }
.controls h4 { margin: 10px 0 0; font-size: 13px; text-transform: uppercase; letter-spacing: 0.05em;
               color: var(--vp-c-text-2); }
.controls label { display: grid; grid-template-columns: 1fr auto; align-items: center; gap: 2px 8px; }
.controls input[type='range'] { grid-column: 1 / -1; width: 100%; accent-color: var(--vp-c-brand-3); }
.controls select { border: 1px solid var(--vp-c-divider); border-radius: 6px; padding: 2px 6px;
                   background: var(--vp-c-bg-soft); color: var(--vp-c-text-1); }
.controls label.check span { display: flex; align-items: center; gap: 6px; }
.controls output { font-variant-numeric: tabular-nums; color: var(--vp-c-text-2); }
.result { min-width: 0; }
canvas.body { width: 100%; height: auto; border: 1px solid var(--vp-c-divider); border-radius: 8px;
              cursor: grab; touch-action: none; color: var(--vp-c-text-2); background: var(--vp-c-bg-soft); }
.hint { margin: 2px 0 8px; font-size: 12px; color: var(--vp-c-text-3); }
.numbers { display: flex; flex-wrap: wrap; gap: 20px; align-items: flex-start; }
.numbers table { margin: 0; font-size: 13px; display: table; width: auto; }
.numbers th { text-align: left; font-weight: 500; padding: 3px 12px 3px 8px; }
.numbers th.part { font-weight: 400; padding-left: 20px; color: var(--vp-c-text-2); }
.numbers td { text-align: right; font-variant-numeric: tabular-nums; padding: 3px 8px; }
.shadow { text-align: center; font-size: 13px; }
.shadow canvas { color: var(--vp-c-text-1); border: 1px solid var(--vp-c-divider); border-radius: 8px;
                 background: var(--vp-c-bg-soft); }
.shadow p { margin: 4px 0 0; line-height: 1.4; }
.code { grid-column: 1 / -1; position: relative; }
.code pre { margin: 0; padding: 14px 16px; border-radius: 8px; background: var(--vp-code-block-bg); overflow-x: auto;
            font-size: 12.5px; line-height: 1.5; }
.code h4 { margin: 0 0 6px; font-size: 14px; }
.code button { position: absolute; top: 36px; right: 8px; font-size: 12px; padding: 2px 10px; border-radius: 6px;
               border: 1px solid var(--vp-c-divider); background: var(--vp-c-bg); color: var(--vp-c-text-2); }
@media (max-width: 720px) { .builder { grid-template-columns: 1fr; } }
</style>
