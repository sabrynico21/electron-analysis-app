/**
 * One-off recovery of the analysis history that was lost when the archive moved
 * out of the renderer's localStorage.
 *
 * Sources, in order of trust:
 *   1. the legacy Zustand `persist` payload recovered from the Chromium LevelDB
 *      (`<userData>/Local Storage/leveldb`) — carries analysis names and the full
 *      input file snapshots;
 *   2. the results already persisted in `<userData>/analysis-results`.
 *
 * The existing archive (if any) is backed up before anything is written.
 *
 * Usage:
 *   node scripts/restore-history.mjs <userDataDir> <legacyHistoryJson>
 */
import fs from 'node:fs'
import path from 'node:path'
import { createRequire } from 'node:module'

const require = createRequire(import.meta.url)
const { HistoryStore } = require('../src/main/analysis/historyStore.js')

const [, , userDataDir, legacyFile] = process.argv

if (!userDataDir || !legacyFile) {
  console.error('Usage: node scripts/restore-history.mjs <userDataDir> <legacyHistoryJson>')
  process.exit(1)
}

const resultsDir = path.join(userDataDir, 'analysis-results')

/** Read the history array out of a dumped legacy Zustand state. */
function readLegacyEntries(file) {
  const raw = JSON.parse(fs.readFileSync(file, 'utf-8'))
  const entries = raw?.state?.history
  return Array.isArray(entries) ? entries : []
}

/** Map `jobId -> result file mtime`, used to recover the real "last updated" time. */
function resultMtimes(dir) {
  const map = new Map()
  if (!fs.existsSync(dir)) return map
  for (const name of fs.readdirSync(dir)) {
    if (!name.endsWith('.json')) continue
    map.set(name.replace(/\.json$/, ''), fs.statSync(path.join(dir, name)).mtimeMs)
  }
  return map
}

const store = new HistoryStore({ baseDir: userDataDir, resultsDir })

if (fs.existsSync(store.historyPath)) {
  const backupPath = `${store.historyPath}.bak-${Date.now()}`
  fs.copyFileSync(store.historyPath, backupPath)
  console.log(`Existing archive backed up to ${backupPath}`)
}

const legacy = readLegacyEntries(legacyFile)
console.log(`Legacy entries recovered from localStorage: ${legacy.length}`)

// Enrich with the result file mtime so "last update" reflects the real last write
// (re-running an analysis or recomputing a cluster touches that file again).
const mtimes = resultMtimes(resultsDir)
const enriched = legacy.map((entry) => {
  const mtime = mtimes.get(entry.id)
  if (!mtime) return entry
  const created = Number(entry.createdAt ?? entry.startedAt) || mtime
  return { ...entry, createdAt: created, updatedAt: Math.max(mtime, created) }
})
console.log(`  of which still have results on disk: ${enriched.filter((e) => mtimes.has(e.id)).length}`)

const imported = await store.importLegacy(enriched)
console.log(`Imported into the archive: ${imported.imported}`)

const rebuilt = await store.rebuildFromResults()
console.log(`Recovered from the results on disk: ${rebuilt.recovered}`)

const entries = await store.list()
console.log(`Archive now holds ${entries.length} analyses -> ${store.historyPath}`)
entries.slice(0, 8).forEach((entry) => {
  console.log(
    `  - ${entry.analysisName} | created ${new Date(entry.createdAt).toLocaleString()}` +
    ` | updated ${new Date(entry.updatedAt).toLocaleString()}`
  )
})
