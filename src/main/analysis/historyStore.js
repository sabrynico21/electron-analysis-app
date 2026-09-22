/**
 * HistoryStore
 * Persistent archive of past analyses.
 *
 * The history used to live in the renderer's `localStorage` through the Zustand
 * `persist` middleware. That is a fragile place for user data: browser storage is
 * partitioned **per origin**, so the Vite dev server (`http://localhost:5173`) and
 * the packaged app (`file://…/app.asar/dist/renderer/`) keep two completely
 * separate copies — switching between the two made entire histories "disappear".
 * It is also wiped whenever the renderer profile is cleared.
 *
 * This store keeps the archive in plain JSON inside the operating system's
 * application-data directory, next to the persisted analysis results:
 *
 *   <userData>/analysis-history.json     ← the archive (this file)
 *   <userData>/analysis-results/*.json   ← the full result of every job
 *
 * Because every result is already written to disk, the history can be rebuilt
 * from scratch (`rebuildFromResults`) if the index is ever lost.
 */
const { app } = require('electron')
const path = require('path')
const fs = require('fs').promises
const os = require('os')
const log = require('electron-log')

/** Bump when the on-disk entry shape changes. */
const SCHEMA_VERSION = 2
/** Safety valve so the archive file cannot grow without bound. */
const MAX_ENTRIES = 500

const DEFAULT_SCRIPT = 'analysis.py'

function toFiniteNumber(value) {
  const numeric = Number(value)
  return Number.isFinite(numeric) ? numeric : null
}

function formatDefaultName({ createdAt, params }) {
  const conditionA = params?.conditionA || params?.condition_a
  const conditionB = params?.conditionB || params?.condition_b
  if (conditionA && conditionB) return `${conditionA} vs ${conditionB}`
  if (conditionA) return String(conditionA)
  const when = new Date(createdAt)
  return `Analysis ${when.toLocaleString()}`
}

class HistoryStore {
  /**
   * @param {{ baseDir?: string, resultsDir?: string }} [options]
   *   `baseDir` defaults to the OS application-data directory and is only
   *   resolved lazily, so the store can also be driven outside of Electron
   *   (useful for the one-off archive import/repair scripts).
   */
  constructor({ baseDir, resultsDir } = {}) {
    this.baseDir = baseDir || app.getPath('userData')
    this.historyPath = path.join(this.baseDir, 'analysis-history.json')
    this.resultsDir = resultsDir || path.join(this.baseDir, 'analysis-results')
  }

  /**
   * Normalize an arbitrarily-shaped history entry (legacy localStorage payloads
   * use `startedAt`, older snapshots use `metadataFiles`, …) into the canonical
   * shape used on disk.
   * @param {object} raw
   * @returns {object|null} `null` when the entry has no usable identifier
   */
  normalize(raw) {
    if (!raw || typeof raw !== 'object') return null
    const id = String(raw.id || raw.jobId || '').trim()
    if (!id) return null

    const createdRaw = toFiniteNumber(raw.createdAt) ?? toFiniteNumber(raw.startedAt) ?? Date.now()
    const updatedRaw = toFiniteNumber(raw.updatedAt) ?? createdRaw
    const params = raw?.payloadSnapshot?.params || raw.params || {}
    const files = raw?.payloadSnapshot?.files || raw.files || {}

    const analysisName = String(raw.analysisName || '').trim()
      || formatDefaultName({ createdAt: createdRaw, params })

    return {
      id,
      analysisName,
      scriptName: raw.scriptName || DEFAULT_SCRIPT,
      status: raw.status || 'completed',
      createdAt: createdRaw,
      updatedAt: updatedRaw,
      runCount: toFiniteNumber(raw.runCount) || 1,
      payloadSnapshot: { files, params },
      summary: raw.summary || null,
      hasResults: raw.hasResults !== false,
    }
  }

  async _ensureBaseDir() {
    await fs.mkdir(this.baseDir, { recursive: true })
  }

  /**
   * Read the archive. A missing or corrupt file is not an error: it simply means
   * "no history yet" so the app can start clean instead of failing to boot.
   * @returns {Promise<{version: number, entries: object[]}>}
   */
  async _read() {
    try {
      const raw = await fs.readFile(this.historyPath, 'utf-8')
      const parsed = JSON.parse(raw)
      const entries = Array.isArray(parsed?.entries)
        ? parsed.entries.map((entry) => this.normalize(entry)).filter(Boolean)
        : []
      return { version: SCHEMA_VERSION, entries }
    } catch (error) {
      if (error.code !== 'ENOENT') {
        log.warn(`History archive unreadable (${this.historyPath}): ${error.message}`)
      }
      return { version: SCHEMA_VERSION, entries: [] }
    }
  }

  /** Write atomically (temp file + rename) so a crash cannot truncate the archive. */
  async _write(entries) {
    await this._ensureBaseDir()
    const sorted = [...entries]
      .sort((a, b) => (b.updatedAt || 0) - (a.updatedAt || 0))
      .slice(0, MAX_ENTRIES)
    const payload = { version: SCHEMA_VERSION, updatedAt: new Date().toISOString(), entries: sorted }
    const tmpPath = `${this.historyPath}.${process.pid}.tmp`
    await fs.writeFile(tmpPath, JSON.stringify(payload, null, 2), 'utf-8')
    await fs.rename(tmpPath, this.historyPath)
    return sorted
  }

  /** @returns {Promise<object[]>} newest-updated first */
  async list() {
    const { entries } = await this._read()
    return [...entries].sort((a, b) => (b.updatedAt || 0) - (a.updatedAt || 0))
  }

  /**
   * Insert or replace one entry. An existing entry keeps its original
   * `createdAt`; unless the caller supplies an explicit `runCount` the run
   * counter is incremented, because upserts happen when a run completes.
   */
  async upsert(rawEntry) {
    const entry = this.normalize(rawEntry)
    if (!entry) return { success: false, error: 'History entry needs an id' }

    const { entries } = await this._read()
    const index = entries.findIndex((item) => item.id === entry.id)
    if (index >= 0) {
      const existing = entries[index]
      const explicitRunCount = toFiniteNumber(rawEntry?.runCount)
      entries[index] = {
        ...existing,
        ...entry,
        createdAt: existing.createdAt,
        runCount: explicitRunCount ?? (existing.runCount || 1) + 1,
      }
    } else {
      entries.push(entry)
    }
    const saved = await this._write(entries)
    return { success: true, entry: saved.find((item) => item.id === entry.id) || entry, entries: saved }
  }

  /**
   * Patch an existing entry (rename, refreshed snapshot, timestamps…).
   * @returns the updated entry, or an error when the id is unknown.
   */
  async update(id, patch = {}) {
    const { entries } = await this._read()
    const index = entries.findIndex((item) => item.id === id)
    if (index < 0) return { success: false, error: 'History entry not found' }

    const current = entries[index]
    const merged = this.normalize({
      ...current,
      ...patch,
      id: current.id,
      createdAt: current.createdAt,
      payloadSnapshot: {
        files: patch?.payloadSnapshot?.files ?? current.payloadSnapshot.files,
        params: patch?.payloadSnapshot?.params ?? current.payloadSnapshot.params,
      },
      updatedAt: patch?.updatedAt ?? Date.now(),
    })
    entries[index] = merged
    const saved = await this._write(entries)
    return { success: true, entry: saved.find((item) => item.id === id) || merged, entries: saved }
  }

  async remove(id) {
    const { entries } = await this._read()
    const saved = await this._write(entries.filter((entry) => entry.id !== id))
    return { success: true, entries: saved }
  }

  async clear() {
    const saved = await this._write([])
    return { success: true, entries: saved }
  }

  /**
   * One-time merge of entries recovered from the legacy `localStorage` archive.
   * Existing ids win, so this is safe to call repeatedly.
   * Legacy names/timestamps are preserved when present.
   */
  async importLegacy(rawEntries = []) {
    if (!Array.isArray(rawEntries) || rawEntries.length === 0) {
      return { success: true, imported: 0, entries: await this.list() }
    }

    const { entries } = await this._read()
    const known = new Set(entries.map((entry) => entry.id))
    let imported = 0

    for (const raw of rawEntries) {
      const entry = this.normalize(raw)
      if (!entry || known.has(entry.id)) continue
      // A legacy entry with a name that is just the script filename is treated as
      // unnamed so the generated label is nicer.
      if (entry.analysisName === entry.scriptName) {
        entry.analysisName = formatDefaultName({ createdAt: entry.createdAt, params: entry.payloadSnapshot.params })
      }
      entries.push(entry)
      known.add(entry.id)
      imported += 1
    }

    if (imported === 0) {
      return { success: true, imported: 0, entries: await this.list() }
    }

    const saved = await this._write(entries)
    log.info(`History: imported ${imported} legacy entries`)
    return { success: true, imported, entries: saved }
  }

  /**
   * Recover history entries straight from the persisted analysis results, and
   * backfill the cached summary of entries that are already known.
   *
   * Every job already leaves `<userData>/analysis-results/<jobId>.json` behind,
   * so the archive can always be rebuilt. The file modification time is used as
   * `updatedAt`, and the parameters are read back from the result summary (or the
   * `input` block written by newer versions of AnalysisManager).
   */
  async rebuildFromResults() {
    let fileNames = []
    try {
      fileNames = (await fs.readdir(this.resultsDir)).filter((name) => name.endsWith('.json'))
    } catch (error) {
      if (error.code !== 'ENOENT') {
        log.warn(`History rebuild: cannot read ${this.resultsDir}: ${error.message}`)
      }
      return { success: true, recovered: 0, entries: await this.list() }
    }

    const { entries } = await this._read()
    const byId = new Map(entries.map((entry) => [entry.id, entry]))
    let recovered = 0
    let enriched = 0

    for (const fileName of fileNames) {
      const jobId = fileName.replace(/\.json$/, '')
      const filePath = path.join(this.resultsDir, fileName)

      // Already archived: only fill in a summary that was never captured.
      const existing = byId.get(jobId)
      if (existing) {
        if (!existing.summary) {
          try {
            const parsed = JSON.parse(await fs.readFile(filePath, 'utf-8'))
            if (parsed?.summary) {
              existing.summary = parsed.summary
              enriched += 1
            }
          } catch (error) {
            log.warn(`History rebuild: cannot read summary of ${fileName}: ${error.message}`)
          }
        }
        continue
      }

      try {
        const [raw, stat] = await Promise.all([fs.readFile(filePath, 'utf-8'), fs.stat(filePath)])
        const parsed = JSON.parse(raw)
        const summary = parsed?.summary || null
        const input = parsed?.input || null
        const params = input?.params || {
          datasetType: summary?.dataset_type || '',
          correlationMethod: summary?.correlation_method || '',
          conditionA: summary?.condition_a || '',
          conditionB: summary?.condition_b || '',
          groupingMode: summary?.grouping_mode || '',
          pValueThreshold: summary?.p_value_threshold ?? undefined,
        }
        const createdAt = toFiniteNumber(input?.completedAt) ?? stat.mtimeMs

        const entry = {
          id: jobId,
          analysisName: String(input?.analysisName || '').trim()
            || formatDefaultName({ createdAt, params }),
          scriptName: input?.scriptName || DEFAULT_SCRIPT,
          status: summary ? 'completed' : 'unknown',
          createdAt,
          updatedAt: stat.mtimeMs,
          runCount: 1,
          payloadSnapshot: {
            files: input?.files || {},
            params,
          },
          summary,
          hasResults: true,
        }
        entries.push(entry)
        byId.set(jobId, entry)
        recovered += 1
      } catch (error) {
        log.warn(`History rebuild: skipping ${fileName}: ${error.message}`)
      }
    }

    if (recovered === 0 && enriched === 0) {
      return { success: true, recovered: 0, enriched: 0, entries: await this.list() }
    }

    const saved = await this._write(entries)
    log.info(`History: recovered ${recovered} entries from persisted results (${enriched} summaries backfilled)`)
    return { success: true, recovered, enriched, entries: saved }
  }

  /**
   * Make sure the archive exists. Called once at startup: when the index is
   * missing (fresh install, cleared profile, or the old localStorage era) it is
   * rebuilt from the results already on disk.
   */
  async ensureInitialized() {
    const { entries } = await this._read()
    if (entries.length > 0) return { rebuilt: false, entries }
    return { rebuilt: true, ...(await this.rebuildFromResults()) }
  }

  /** Diagnostic info surfaced in Settings. */
  info() {
    return {
      historyPath: this.historyPath,
      resultsDir: this.resultsDir,
      tempDir: os.tmpdir(),
    }
  }

  /** Backup copy of the raw archive, used before destructive operations. */
  async backupRaw() {
    try {
      return await fs.readFile(this.historyPath, 'utf-8')
    } catch {
      return null
    }
  }
}

module.exports = { HistoryStore, SCHEMA_VERSION, MAX_ENTRIES }
