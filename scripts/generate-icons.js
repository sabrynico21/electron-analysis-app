#!/usr/bin/env node
/**
 * Generates the application icons used by electron-builder and by the
 * BrowserWindow:
 *
 *   resources/icons/icon.png    1024x1024 master (Linux + Electron window icon)
 *   resources/icons/icon.ico    multi-size Windows icon
 *   resources/icons/icon.icns   multi-size macOS icon
 *
 * Requires ImageMagick 7 (`magick`) on PATH. Re-run after changing the design:
 *
 *   node scripts/generate-icons.js
 *
 * The design below is a neutral placeholder (a microbial co-occurrence network).
 * Replace the generated files with your official artwork before publishing if
 * you have branded assets — the build only requires the three files to exist.
 */
const { execFileSync } = require('child_process')
const fs = require('fs')
const os = require('os')
const path = require('path')

const OUT_DIR = path.resolve(__dirname, '..', 'resources', 'icons')
const BASE = 1024
const RADIUS = 196

// --- Network layout (coordinates on a 1024x1024 canvas) ---------------------
// `kind` only drives the colour: seed = amber, bacteria = sky, fungi = violet.
const NODES = [
  { id: 'seed', x: 512, y: 522, r: 74, kind: 'seed' },
  { id: 'b1', x: 250, y: 330, r: 46, kind: 'bacteria' },
  { id: 'b2', x: 500, y: 248, r: 38, kind: 'bacteria' },
  { id: 'b3', x: 200, y: 700, r: 40, kind: 'bacteria' },
  { id: 'b4', x: 520, y: 806, r: 42, kind: 'bacteria' },
  { id: 'b5', x: 648, y: 566, r: 30, kind: 'bacteria' },
  { id: 'f1', x: 786, y: 322, r: 42, kind: 'fungi' },
  { id: 'f2', x: 836, y: 690, r: 48, kind: 'fungi' },
  { id: 'f3', x: 706, y: 764, r: 32, kind: 'fungi' },
]

const EDGES = [
  ['seed', 'b1'], ['seed', 'b2'], ['seed', 'b3'], ['seed', 'b4'],
  ['seed', 'b5'], ['seed', 'f1'], ['seed', 'f2'], ['seed', 'f3'],
  ['b1', 'b2'], ['b3', 'b4'], ['f1', 'b5'], ['b4', 'f2'], ['f2', 'f3'],
]

const COLOR = {
  seed: '#fbbf24',
  bacteria: '#7dd3fc',
  fungi: '#c4b5fd',
}

const byId = Object.fromEntries(NODES.map((node) => [node.id, node]))

function magick(args, label) {
  try {
    execFileSync('magick', args, { stdio: ['ignore', 'ignore', 'pipe'] })
  } catch (error) {
    const detail = error.stderr ? error.stderr.toString().trim() : error.message
    throw new Error(`ImageMagick failed during ${label}: ${detail}`)
  }
}

function draw(nodes, colorTag) {
  return nodes.flatMap(({ x, y, r }) => ['-draw', `circle ${x},${y} ${x},${y - r}`])
}

function buildMaster(file) {
  // Rounded-square backdrop with a dark navy -> teal diagonal gradient.
  magick([
    '-size', `${BASE}x${BASE}`,
    'gradient:#0b1220-#115e59',
    '(',
    '-size', `${BASE}x${BASE}`, 'xc:none',
    '-fill', 'white',
    '-draw', `roundrectangle 0,0,${BASE - 1},${BASE - 1},${RADIUS},${RADIUS}`,
    ')',
    '-compose', 'CopyOpacity', '-composite',
    file,
  ], 'background')
}

function drawNetwork(file) {
  const edgeArgs = EDGES.flatMap(([a, b]) =>
    ['-draw', `line ${byId[a].x},${byId[a].y} ${byId[b].x},${byId[b].y}`]
  )

  // Edges first so that nodes sit on top of them.
  magick([
    file,
    '-fill', 'none',
    '-stroke', 'rgba(125,211,252,0.38)',
    '-strokewidth', '7',
    ...edgeArgs,
    file,
  ], 'edges')

  for (const kind of ['bacteria', 'fungi']) {
    magick([
      file,
      '-stroke', 'none',
      '-fill', COLOR[kind],
      ...draw(NODES.filter((node) => node.kind === kind)),
      file,
    ], `${kind} nodes`)
  }

  // Seed node: amber disc with a light halo ring.
  const seed = byId.seed
  magick([
    file,
    '-stroke', 'none',
    '-fill', COLOR.seed,
    '-draw', `circle ${seed.x},${seed.y} ${seed.x},${seed.y - seed.r}`,
    '-fill', 'none',
    '-stroke', 'rgba(255,255,255,0.85)',
    '-strokewidth', '9',
    '-draw', `circle ${seed.x},${seed.y} ${seed.x},${seed.y - (seed.r + 26)}`,
    file,
  ], 'seed node')
}

function resize(source, size, target) {
  magick([source, '-resize', `${size}x${size}`, target], `resize ${size}`)
}

/**
 * Minimal ICNS writer: the modern container simply stores PNG payloads with a
 * 4-byte OSType and a big-endian length. ImageMagick cannot write ICNS here.
 */
function writeIcns(entries, target) {
  const chunks = entries.map(({ type, data }) => {
    const header = Buffer.alloc(8)
    header.write(type, 0, 4, 'ascii')
    header.writeUInt32BE(data.length + 8, 4)
    return Buffer.concat([header, data])
  })
  const body = Buffer.concat(chunks)
  const header = Buffer.alloc(8)
  header.write('icns', 0, 4, 'ascii')
  header.writeUInt32BE(body.length + 8, 4)
  fs.writeFileSync(target, Buffer.concat([header, body]))
}

function main() {
  try {
    execFileSync('magick', ['-version'], { stdio: 'ignore' })
  } catch {
    console.error('ImageMagick 7 (`magick`) is required on PATH to generate icons.')
    process.exit(1)
  }

  fs.mkdirSync(OUT_DIR, { recursive: true })
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'analysis-app-icons-'))

  try {
    const master = path.join(tmp, 'master.png')
    buildMaster(master)
    drawNetwork(master)

    const iconPng = path.join(OUT_DIR, 'icon.png')
    fs.copyFileSync(master, iconPng)

    const iconIco = path.join(OUT_DIR, 'icon.ico')
    magick([
      master,
      '-define', 'icon:auto-resize=256,128,64,48,32,16',
      iconIco,
    ], 'icon.ico')

    // macOS icon: PNG entries for each Retina-aware slot.
    const icnsSizes = [
      { size: 32, type: 'ic11' },
      { size: 64, type: 'ic12' },
      { size: 128, type: 'ic07' },
      { size: 256, type: 'ic13' },
      { size: 512, type: 'ic09' },
      { size: 1024, type: 'ic10' },
    ]
    const entries = icnsSizes.map(({ size, type }) => {
      const file = path.join(tmp, `icon_${size}.png`)
      resize(master, size, file)
      return { type, data: fs.readFileSync(file) }
    })
    writeIcns(entries, path.join(OUT_DIR, 'icon.icns'))

    console.log('Icons written to resources/icons:')
    for (const name of ['icon.png', 'icon.ico', 'icon.icns']) {
      const { size } = fs.statSync(path.join(OUT_DIR, name))
      console.log(`  ${name}  ${(size / 1024).toFixed(1)} KiB`)
    }
  } finally {
    fs.rmSync(tmp, { recursive: true, force: true })
  }
}

main()
