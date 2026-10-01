// Geometry of the animal builder page.
//
// A JavaScript copy of the BiophysicalGeometry.jl calculations for the shapes the builder offers, so the page
// works without Julia. `docs/check_builder.jl` runs presets through this file and through the package, and
// compares the numbers. Lengths in m, masses in kg, areas in m².

const PI = Math.PI

// ── Shapes ──────────────────────────────────────────────────────────────────

const cylinderRadius = (volume, ratio) => Math.cbrt(volume / (2 * PI * ratio))
const cylinderArea = (r, L) => 2 * PI * r * L + 2 * PI * r * r

const coneRadius = (volume, ratio, top) => Math.cbrt((3 * volume) / (2 * PI * ratio * (1 + top + top * top)))
function coneArea(R, L, top) {
  const r = top * R
  return PI * R * R + PI * r * r + PI * (R + r) * Math.hypot(R - r, L)
}

const sphereRadius = (volume) => Math.cbrt((3 * volume) / (4 * PI))

function prolateArea(a, b) {
  if (Math.abs(a - b) < a * 1e-9) return 4 * PI * b * b
  const e = Math.sqrt(a * a - b * b) / a
  return 2 * PI * b * b + 2 * PI * ((a * b) / e) * Math.asin(e)
}

// Thickness x of an even fat layer over a prolate flesh ellipsoid: (ratio * b + x)(b + x)² = 3 V / 4π.
function ellipsoidFat(volume, ratio, bFlesh) {
  const target = (3 * volume) / (4 * PI)
  let x = 0
  for (let i = 0; i < 50; i++) {
    const f = (ratio * bFlesh + x) * (bFlesh + x) ** 2 - target
    const df = (bFlesh + x) ** 2 + 2 * (ratio * bFlesh + x) * (bFlesh + x)
    x -= f / df
  }
  return Math.max(0, x)
}

// ── Julia code for values ───────────────────────────────────────────────────

const num = (x, digits = 5) => {
  const s = Number(x.toPrecision(digits)).toString()
  return s.includes('.') || s.includes('e') ? s : s + '.0'
}
const metres = (x) => `${num(x)}u"m"`
const kilograms = (x) => `${num(x, 6)}u"kg"`
const fibres = (depth) => (depth > 0 ? `FibrousLayer(${num(depth * 1000)}u"mm", 30.0u"μm", 3000u"cm^-2")` : 'Naked()')

// ── Parts ───────────────────────────────────────────────────────────────────
//
// Each part has a `kind` that fixes its local frame, as in the package, its skin dimensions, the depth of its
// fibres, its outer and skin areas, and the Julia code of its shape.

function torsoHalf(p, fur) {
  const volume = p.torsoMass / p.density // of the whole torso, which the half is cut from
  const fleshVolume = volume - (p.torsoMass * p.fat) / p.fatDensity
  const fatLayer = `FatLayer(${num(p.fat)}, ${num(p.fatDensity)}u"kg/m^3")`
  const layers = p.fat > 0 ? (fur > 0 ? `CompositeInsulation(${fibres(fur)}, ${fatLayer})` : fatLayer) : fibres(fur)
  if (p.torsoShape === 'Cylinder') {
    const r = cylinderRadius(volume, p.torsoRatio)
    const L = 2 * p.torsoRatio * r
    const cut = 2 * r * L
    const skin = cylinderArea(r, L) / 2 + cut
    const total = fur > 0 ? cylinderArea(r + fur, L + 2 * fur) / 2 + cut : skin
    const code = `Body(HalfCylinder(${kilograms(p.torsoMass / 2)}, density, ${num(p.torsoRatio)}), ${layers})`
    return { kind: 'halfCylinder', r, L, cut, skin, total, fur, code }
  }
  let b = sphereRadius(volume / p.torsoRatio)
  let a = p.torsoRatio * b
  if (p.fat > 0) {
    const bFlesh = sphereRadius(fleshVolume / p.torsoRatio)
    const fatThickness = ellipsoidFat(volume, p.torsoRatio, bFlesh)
    if (fatThickness > 0) {
      b = bFlesh + fatThickness
      a = p.torsoRatio * bFlesh + fatThickness
    }
  }
  const cut = PI * b * b
  const skin = prolateArea(a, b) / 2 + cut
  const total = fur > 0 ? prolateArea(a + fur, b + fur) / 2 + cut : skin
  const code = `Body(HalfEllipsoid(${kilograms(p.torsoMass / 2)}, density, ${num(p.torsoRatio)}, 1.0), ${layers})`
  return { kind: 'halfEllipsoid', a, b, r: b, L: 2 * a, cut, skin, total, fur, code }
}

function sphere(mass, p) {
  const r = sphereRadius(mass / p.density), fur = p.limbFur
  return { kind: 'sphere', r, a: r, b: r, skin: 4 * PI * r * r, total: 4 * PI * (r + fur) ** 2, fur, mass,
           code: `Body(Sphere(${kilograms(mass)}, density), coat)` }
}

function ellipsoid(mass, ratio, p) {
  const b = sphereRadius(mass / p.density / ratio), a = ratio * b, fur = p.limbFur
  return { kind: 'ellipsoid', r: b, a, b, skin: prolateArea(a, b), total: prolateArea(a + fur, b + fur), fur, mass,
           code: `Body(Ellipsoid(${kilograms(mass)}, density, ${num(ratio)}, 1.0), coat)` }
}

// A cylinder (`top` of 1) or a cone or frustum (`top` less than 1).
function axial(mass, ratio, top, p) {
  const volume = mass / p.density, fur = p.limbFur
  if (top >= 1) {
    const r = cylinderRadius(volume, ratio), L = 2 * ratio * r
    return { kind: 'axial', r, L, top: 1, skin: cylinderArea(r, L), total: cylinderArea(r + fur, L + 2 * fur), fur, mass,
             code: `Body(Cylinder(${kilograms(mass)}, density, ${num(ratio)}), coat)` }
  }
  const r = coneRadius(volume, ratio, top), L = 2 * ratio * r
  return { kind: 'axial', r, L, top, skin: coneArea(r, L, top), total: coneArea(r + fur, L + 2 * fur, top), fur, mass,
           code: `Body(Cone(${kilograms(mass)}, density, ${num(ratio)}, ${num(top)}), coat)` }
}

// A flat plate: width = length / ratio, height = length / flatness.
function plate(mass, ratio, flatness, p) {
  const L = Math.cbrt((mass / p.density) * ratio * flatness), W = L / ratio, H = L / flatness, fur = p.limbFur
  const area = (l, w, h) => 2 * (l * w + l * h + w * h)
  return { kind: 'plate', L, W, H, skin: area(L, W, H), total: area(L + 2 * fur, W + 2 * fur, H + 2 * fur), fur, mass,
           code: `Body(Plate(${kilograms(mass)}, density, ${num(ratio)}, ${num(flatness)}), coat)` }
}

// ── Surfaces ────────────────────────────────────────────────────────────────
//
// A location on a part: its point and outward normal in the part's own frame, and its Julia code.

const unit = (v) => { const n = Math.hypot(...v); return v.map((x) => x / n) }

const at = {
  flat: (part) => part.kind === 'halfCylinder'
    ? { point: [0, 0, part.L / 2], normal: [0, -1, 0], code: 'Flat()' }
    : { point: [0, 0, 0], normal: [0, 0, -1], code: 'Flat()' },
  endA: (r = 0, phi = 0) => ({ point: [r * Math.cos(phi), r * Math.sin(phi), 0], normal: [0, 0, -1],
                               code: `EndA(${metres(r)}, ${num(phi)})` }),
  endB: (part, r = 0, phi = 0) => ({ point: [r * Math.cos(phi), r * Math.sin(phi), part.L], normal: [0, 0, 1],
                                     code: `EndB(${metres(r)}, ${num(phi)})` }),
  lateral: (part, z, phi) => ({ point: [part.r * Math.cos(phi), part.r * Math.sin(phi), z],
                                normal: [Math.cos(phi), Math.sin(phi), 0], code: `Lateral(${metres(z)}, ${num(phi)})` }),
  dome: (part, alpha, beta) => {
    const { a, b } = part
    return { point: [a * Math.cos(alpha), b * Math.sin(alpha) * Math.cos(beta), b * Math.sin(alpha) * Math.sin(beta)],
             normal: unit([Math.cos(alpha) * b * b, Math.sin(alpha) * Math.cos(beta) * a * b, Math.sin(alpha) * Math.sin(beta) * a * b]),
             code: `Dome(${num(alpha)}, ${num(beta)})` }
  },
  poleA: (part) => ({ point: [part.a, 0, 0], normal: [1, 0, 0], code: 'PoleA()' }),
  poleB: (part) => ({ point: [-part.a, 0, 0], normal: [-1, 0, 0], code: 'PoleB()' }),
  equator: (part, phi) => ({ point: [0, part.b * Math.cos(phi), part.b * Math.sin(phi)], normal: [0, Math.cos(phi), Math.sin(phi)],
                             code: `Equator(${num(phi)})` }),
  bottom: (part) => ({ point: [0, 0, -part.H / 2], normal: [0, 0, -1], code: 'Bottom(0.0u"m", 0.0u"m")' }),
  bottomFace: (part) => ({ point: [0, 0, -part.H / 2], normal: [0, 0, -1], code: 'Bottom()' }),
  sideB: (part) => ({ point: [-part.L / 2, 0, 0], normal: [-1, 0, 0], code: 'SideB(0.0u"m", 0.0u"m")' }),
  radial: (part, theta, phi) => {
    const n = [Math.sin(theta) * Math.cos(phi), Math.sin(theta) * Math.sin(phi), Math.cos(theta)]
    return { point: n.map((x) => part.r * x), normal: n, code: `Radial(${num(theta)}, ${num(phi)})` }
  },
}

// ── Poses ───────────────────────────────────────────────────────────────────
//
// As in the package: a child is placed so that the two points of a join meet and the two normals are opposed,
// then turned about the join by the twist.

const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
const rotate = (R, v) => R.map((row) => dot(row, v))
const multiply = (A, B) => A.map((row) => [0, 1, 2].map((j) => row[0] * B[0][j] + row[1] * B[1][j] + row[2] * B[2][j]))
const IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

function axisAngle([x, y, z], angle) {
  const c = Math.cos(angle), s = Math.sin(angle), t = 1 - c
  return [[t * x * x + c, t * x * y - s * z, t * x * z + s * y],
          [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
          [t * x * z - s * y, t * y * z + s * x, t * z * z + c]]
}

function align(a, b) {
  const d = dot(a, b)
  if (d > 1 - 1e-12) return IDENTITY
  if (d < -1 + 1e-12) {
    const axis = Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]
    const proj = dot(a, axis)
    return axisAngle(unit(axis.map((x, k) => x - proj * a[k])), PI)
  }
  return axisAngle(unit(cross(a, b)), Math.acos(d))
}

// The twist that turns the child's own direction `axis` as near as it can to the direction `towards`.
function twistTo(parent, on, child, axis, towards) {
  const target = rotate(parent.rotation, on.normal).map((x) => -x)
  const now = rotate(align(child.normal, target), axis)
  const along = dot(towards, target)
  const wanted = unit(towards.map((x, k) => x - along * target[k]))
  return Math.atan2(dot(cross(now, wanted), target), dot(now, wanted))
}

function childPose(parent, on, child, twist = 0) {
  const point = rotate(parent.rotation, on.point).map((x, k) => x + parent.translation[k])
  const target = rotate(parent.rotation, on.normal).map((x) => -x)
  const rotation = multiply(axisAngle(target, twist), align(child.normal, target))
  const moved = rotate(rotation, child.point)
  return { rotation, translation: point.map((x, k) => x - moved[k]) }
}

// ── The animal ──────────────────────────────────────────────────────────────

export const defaults = {
  mass: 20, density: 1000, fatDensity: 901,
  torsoShape: 'Cylinder', torsoRatio: 3, fat: 0.1, backFur: 0.015, bellyFur: 0.005, limbFur: 0.008,
  headShape: 'Ellipsoid', headRatio: 1.5, headFraction: 0.08,
  neck: false, neckFraction: 0.04, neckRatio: 1, neckPosture: 'Forward',
  nose: 'None', noseFraction: 0.005,
  ears: 'None', earFraction: 0.002, earPosture: 'Up', earRatio: 1.5, earFlatness: 10,
  legs: 4, legFraction: 0.03, legRatio: 5, legTop: 0.5,
  wings: 'None', wingFraction: 0.04,
  tail: false, tailFraction: 0.01, tailRatio: 6,
}

const LEG_SPLAY = 0.35
const legNames = (n) => (n === 4 ? ['leg_fl', 'leg_fr', 'leg_bl', 'leg_br'] : n === 2 ? ['leg_l', 'leg_r'] : [])
// (position along the torso, side) of each leg
const legPlaces = (n) => (n === 4 ? [[0.85, 1], [0.85, -1], [0.15, 1], [0.15, -1]] : n === 2 ? [[0.5, 1], [0.5, -1]] : [])

export function build(input) {
  const p = { ...defaults, ...input }
  const hasHead = p.headShape !== 'None'
  const mass = {
    head: hasHead ? p.headFraction * p.mass : 0,
    neck: hasHead && p.neck ? p.neckFraction * p.mass : 0,
    nose: hasHead && p.nose !== 'None' ? p.noseFraction * p.mass : 0,
    ear: hasHead && p.ears !== 'None' ? p.earFraction * p.mass : 0,
    leg: p.legs > 0 ? p.legFraction * p.mass : 0,
    tail: p.tail ? p.tailFraction * p.mass : 0,
    wing: p.wings !== 'None' ? p.wingFraction * p.mass : 0,
  }
  p.torsoMass = p.mass - mass.head - mass.neck - mass.nose - 2 * mass.ear - p.legs * mass.leg - mass.tail - 2 * mass.wing
  const cylinder = p.torsoShape === 'Cylinder'

  const dorsal = { name: 'dorsal', ...torsoHalf(p, p.backFur), mass: p.torsoMass / 2 }
  const ventral = { name: 'ventral', ...torsoHalf(p, p.bellyFur), mass: p.torsoMass / 2 }
  // the torso lies along x with the back up
  dorsal.pose = { rotation: cylinder ? [[0, 0, 1], [1, 0, 0], [0, 1, 0]] : IDENTITY, translation: [0, 0, 0] }
  const parts = [dorsal, ventral]
  const joins = []
  const definitions = [] // Julia lines for the bodies, each used under one or more names

  function join(parent, on, child, childAt, patchArea, patchCode, twist = 0, childPatchCode = patchCode) {
    child.pose = childPose(parent.pose, on, childAt, twist)
    joins.push({ parent, child, area: patchArea, twist,
                 code: `Join(${parent.name} = Attachment(${on.code}, ${patchCode}), ${child.name} = Attachment(${childAt.code}, ${childPatchCode})` +
                       (Math.abs(twist) > 1e-9 ? `; twist = ${num(twist, 7)}),` : '),') })
    parts.push(child)
  }
  const disc = (radius) => [PI * radius * radius, `Disc(${metres(radius)})`]
  const onTorso = (part, towardHead, side = 0) => cylinder
    ? (side === 0 ? at[towardHead ? 'endB' : 'endA'](...(towardHead ? [part, part.r / 2, PI / 2] : [part.r / 2, PI / 2]))
                  : at.lateral(part, towardHead * part.L, PI / 2 + side * LEG_SPLAY))
    : (side === 0 ? at.dome(part, towardHead ? 0.3 : PI - 0.3, PI / 2)
                  : at.dome(part, towardHead > 0.6 ? 0.9 : towardHead < 0.4 ? PI - 0.9 : PI / 2, PI / 2 + side * 0.5))

  parts.pop() // ventral is added by its join
  join(dorsal, at.flat(dorsal), ventral, at.flat(ventral), dorsal.cut, 'FullCover()', cylinder ? -PI / 2 : 0)

  if (hasHead) {
    const head = { name: 'head', ...(p.headShape === 'Sphere' ? sphere(mass.head, p) : ellipsoid(mass.head, p.headRatio, p)) }
    const back = p.headShape === 'Sphere' ? at.radial(head, PI / 2, PI) : at.poleB(head)
    definitions.push(`head = ${head.code}`)
    if (p.neck) {
      const neck = { name: 'neck', ...axial(mass.neck, p.neckRatio, 0.7, p) }
      definitions.unshift(`neck = ${neck.code}`)
      const up = p.neckPosture === 'Up'
      const base = !up ? onTorso(dorsal, true) : cylinder ? at.lateral(dorsal, 0.9 * dorsal.L, PI / 2) : at.dome(dorsal, 0.6, PI / 2)
      join(dorsal, base, neck, at.endA(), ...disc(Math.min(neck.r, 0.45 * dorsal.r)))
      const patch = disc(Math.min(0.7 * neck.r, 0.45 * head.r))
      if (up && p.headShape === 'Ellipsoid') {   // the head sits across the top of the neck, facing forward
        const under = at.equator(head, -PI / 2)
        join(neck, at.endB(neck), head, under, ...patch, twistTo(neck.pose, at.endB(neck), under, [1, 0, 0], [1, 0, 0]))
      } else {
        join(neck, at.endB(neck), head, back, ...patch)
      }
    } else {
      join(dorsal, onTorso(dorsal, true), head, back, ...disc(0.45 * Math.min(head.r, dorsal.r)))
    }
    if (p.nose !== 'None') {
      const nose = { name: p.nose === 'Beak' ? 'beak' : 'nose',
                     ...(p.nose === 'Beak' ? axial(mass.nose, 2.0, 0, p) : axial(mass.nose, 1.0, 0.5, p)) }
      definitions.push(`${nose.name} = ${nose.code}`)
      const front = p.headShape === 'Sphere' ? at.radial(head, PI / 2, 0) : at.poleA(head)
      join(head, front, nose, at.endA(), ...disc(Math.min(nose.r, 0.45 * head.r)))
    }
    if (p.ears !== 'None') {
      const flatEar = p.ears === 'Plate'
      const ear = flatEar ? plate(mass.ear, p.earRatio, p.earFlatness, p) : axial(mass.ear, 1.5, 0.3, p)
      definitions.push(`ear = ${ear.code}`)
      const forward = rotate(head.pose.rotation, [1, 0, 0])
      for (const [name, side] of [['ear_l', 1], ['ear_r', -1]]) {
        const on = p.headShape === 'Sphere' ? at.radial(head, 0.6, (side * PI) / 2) : at.equator(head, PI / 2 - side * 0.6)
        const child = { name, ...ear, alias: 'ear' }
        if (!flatEar) {
          join(head, on, child, at.endA(), ...disc(Math.min(ear.r, 0.45 * head.r)))
        } else if (p.earPosture === 'Up') {   // standing on its edge, its flat side to the front
          const twist = twistTo(head.pose, on, at.sideB(ear), [0, 0, 1], forward)
          join(head, on, child, at.sideB(ear), ...disc(0.9 * Math.sqrt((ear.W * ear.H) / PI)), twist)
        } else {                              // laid back against the head, which hides one whole face
          const twist = twistTo(head.pose, on, at.bottomFace(ear), [1, 0, 0], forward.map((x) => -x))
          const face = (ear.L + 2 * ear.fur) * (ear.W + 2 * ear.fur)
          join(head, on, child, at.bottomFace(ear), face, 'Disc(sqrt(surface_area(ear.shape, ear, Bottom()) / π))', twist, 'FullCover()')
        }
      }
    }
  }
  if (p.legs > 0) {
    const leg = axial(mass.leg, p.legRatio, p.legTop, p)
    definitions.push(`leg = ${leg.code}`)
    legPlaces(p.legs).forEach(([along, side], k) => {
      join(ventral, onTorso(ventral, along, side), { name: legNames(p.legs)[k], ...leg, alias: 'leg' }, at.endA(), ...disc(leg.r))
    })
  }
  if (p.wings !== 'None') {
    const wing = plate(mass.wing, 2.5, 15, p)
    definitions.push(`wing = ${wing.code}`)
    for (const [name, side] of [['wing_l', 1], ['wing_r', -1]]) {
      const child = { name, ...wing, alias: 'wing' }
      const on = cylinder ? at.lateral(dorsal, 0.6 * dorsal.L, PI / 2 - side * 1.1) : at.dome(dorsal, 1.2, PI / 2 - side * 1.1)
      if (p.wings === 'Folded') {   // flat against the body, lying along it
        const twist = twistTo(dorsal.pose, on, at.bottom(wing), [1, 0, 0], [-1, 0, 0])
        join(dorsal, on, child, at.bottom(wing), ...disc(0.3 * wing.W), twist)
      } else {                      // out to the side, flat side up
        const twist = twistTo(dorsal.pose, on, at.sideB(wing), [0, 0, 1], [0, 0, 1])
        join(dorsal, on, child, at.sideB(wing), ...disc(0.9 * Math.sqrt((wing.W * wing.H) / PI)), twist)
      }
    }
  }
  if (p.tail) {
    const tail = { name: 'tail', ...axial(mass.tail, p.tailRatio, 0.3, p) }
    definitions.push(`tail = ${tail.code}`)
    join(dorsal, onTorso(dorsal, false), tail, at.endA(), ...disc(Math.min(tail.r, 0.45 * dorsal.r)))
  }

  for (const part of parts) part.hidden = 0
  let joined = 0
  for (const j of joins) { j.parent.hidden += j.area; j.child.hidden += j.area; joined += j.area }
  const total = parts.reduce((s, part) => s + part.total, 0) - 2 * joined
  const skin = parts.reduce((s, part) => s + part.skin, 0) - 2 * joined
  return {
    params: p, parts, total, skin, joined, volume: p.mass / p.density, meeh: total / Math.cbrt(p.mass) ** 2,
    triangles: parts.flatMap(mesh), code: juliaCode(p, dorsal, ventral, definitions, parts, joins, cylinder),
  }
}

// ── Meshes ──────────────────────────────────────────────────────────────────

// Triangles of a grid of points `f(i / nu, j / nv)` in the frame of `part`, moved to its pose.
function grid(f, nu, nv, part) {
  const { rotation, translation } = part.pose
  const pts = []
  for (let i = 0; i <= nu; i++) for (let j = 0; j <= nv; j++) {
    pts.push(rotate(rotation, f(i / nu, j / nv)).map((x, k) => x + translation[k]))
  }
  const get = (i, j) => pts[i * (nv + 1) + j]
  const out = []
  for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) {
    out.push({ part: part.name, pts: [get(i, j), get(i + 1, j), get(i, j + 1)] })
    out.push({ part: part.name, pts: [get(i + 1, j), get(i + 1, j + 1), get(i, j + 1)] })
  }
  return out
}

// The outer surface of a part, over its fibres, as the package draws it.
function mesh(part) {
  const fur = part.fur
  if (part.kind === 'plate') {
    const h = [part.L / 2 + fur, part.W / 2 + fur, part.H / 2 + fur]
    const tris = []
    for (let k = 0; k < 3; k++) for (const sign of [-1, 1]) {
      const i = (k + 1) % 3, j = (k + 2) % 3
      tris.push(...grid((s, t) => { const q = [0, 0, 0]; q[k] = sign * h[k]; q[i] = (2 * s - 1) * h[i]; q[j] = (2 * t - 1) * h[j]; return q },
        4, 4, part))
    }
    return tris
  }
  if (part.kind === 'halfCylinder' || part.kind === 'axial') {
    const turn = part.kind === 'axial' ? 2 * PI : PI
    const top = part.kind === 'axial' ? part.top : 1
    const r = part.r + fur, z0 = -fur, z1 = part.L + fur
    const ring = (z, radius, t) => [radius * Math.cos(turn * t), radius * Math.sin(turn * t), z]
    const n = part.kind === 'axial' ? 24 : 12
    const tris = grid((i, j) => ring(z0 + i * (z1 - z0), r * (1 + (top - 1) * i), j), 8, n, part)
    tris.push(...grid((i, j) => ring(z0, r * i, j), 1, n, part))
    if (top > 0) tris.push(...grid((i, j) => ring(z1, top * r * i, j), 1, n, part))
    return tris
  }
  const a = part.a + fur, b = part.b + fur
  const span = part.kind === 'halfEllipsoid' ? PI / 2 : PI
  return grid((i, j) => {
    const phi = span * i, theta = 2 * PI * j
    return [a * Math.sin(phi) * Math.cos(theta), b * Math.sin(phi) * Math.sin(theta), b * Math.cos(phi)]
  }, part.kind === 'halfEllipsoid' ? 8 : 16, 24, part)
}

// ── Silhouette ──────────────────────────────────────────────────────────────

/** Area of the animal seen from `direction`, by drawing its triangles on a grid, and the grid itself. */
export function silhouette(triangles, direction, n = 240) {
  const d = unit(direction)
  const up = Math.abs(d[2]) > 0.999 ? [1, 0, 0] : [0, 0, 1]
  const proj = dot(up, d)
  const v = unit(up.map((x, k) => x - proj * d[k]))
  const u = cross(v, d)
  const flat = triangles.map((t) => t.pts.map((q) => [dot(q, u), dot(q, v)]))
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity
  for (const t of flat) for (const q of t) {
    x0 = Math.min(x0, q[0]); x1 = Math.max(x1, q[0]); y0 = Math.min(y0, q[1]); y1 = Math.max(y1, q[1])
  }
  const pad = 0.02 * Math.max(x1 - x0, y1 - y0)
  x0 -= pad; x1 += pad; y0 -= pad; y1 += pad
  const dx = (x1 - x0) / n, dy = (y1 - y0) / n
  const covered = new Uint8Array(n * n)
  for (const [a, b, c] of flat) {
    const s = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    if (Math.abs(s) < 1e-18) continue
    const sign = Math.sign(s)
    const i0 = Math.max(0, Math.floor((Math.min(a[0], b[0], c[0]) - x0) / dx))
    const i1 = Math.min(n - 1, Math.ceil((Math.max(a[0], b[0], c[0]) - x0) / dx))
    const j0 = Math.max(0, Math.floor((Math.min(a[1], b[1], c[1]) - y0) / dy))
    const j1 = Math.min(n - 1, Math.ceil((Math.max(a[1], b[1], c[1]) - y0) / dy))
    for (let j = j0; j <= j1; j++) {
      const y = y0 + (j + 0.5) * dy
      for (let i = i0; i <= i1; i++) {
        if (covered[j * n + i]) continue
        const x = x0 + (i + 0.5) * dx
        if (sign * ((b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0])) < 0) continue
        if (sign * ((c[0] - b[0]) * (y - b[1]) - (c[1] - b[1]) * (x - b[0])) < 0) continue
        if (sign * ((a[0] - c[0]) * (y - c[1]) - (a[1] - c[1]) * (x - c[0])) < 0) continue
        covered[j * n + i] = 1
      }
    }
  }
  let count = 0
  for (let k = 0; k < covered.length; k++) count += covered[k]
  return { area: count * dx * dy, covered, n, width: x1 - x0, height: y1 - y0 }
}

// ── Julia code ──────────────────────────────────────────────────────────────

function juliaCode(p, dorsal, ventral, definitions, parts, joins, cylinder) {
  const lines = ['using BiophysicalGeometry, Unitful', '', `density = ${num(p.density)}u"kg/m^3"`]
  if (definitions.length > 0) lines.push(`coat = ${fibres(p.limbFur)}`)
  lines.push(`dorsal = ${dorsal.code}`, `ventral = ${ventral.code}`, ...definitions)
  const names = parts.map((part) => (part.alias ? `${part.name} = ${part.alias}` : part.name))
  lines.push('', 'animal = CompositeBody(;', `    parts = (; ${names.join(', ')}),`, '    joins = (')
  joins.forEach((j) => lines.push('        ' + j.code))
  lines.push('    ),')
  if (cylinder) lines.push('    root_pose = Pose((0.0u"m", 0.0u"m", 0.0u"m"), [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]),')
  lines.push(')', '', 'total_area(animal), skin_area(animal)')
  return lines.join('\n')
}
