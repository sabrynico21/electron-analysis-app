/**
 * Unit tests for the persistent history archive (the main-process side).
 *
 * These run in plain Node with no Electron instance: `HistoryStore` only needs a
 * base directory, so the tests are fast and hermetic.
 *
 *   npm test
 */
const fs = require('fs')
const os = require('os')
const path = require('path')

const { HistoryStore, MAX_ENTRIES } = require('../../src/main/analysis/historyStore.js')

const makeStore = () => {
  const baseDir = fs.mkdtempSync(path.join(os.tmpdir(), 'history-store-'))
  return { baseDir, resultsDir: path.join(baseDir, 'analysis-results'), store: new HistoryStore({ baseDir, resultsDir: path.join(baseDir, 'analysis-results') }) }
}

const entry = (overrides = {}) => ({
  id: 'job-1',
  analysisName: 'Analysis one',
  createdAt: 1_700_000_000_000,
  updatedAt: 1_700_000_000_000,
  payloadSnapshot: { files: { bacteria: ['/data/a.csv'] }, params: { conditionA: 'A', conditionB: 'B' } },
  ...overrides,
})

const writeResult = (resultsDir, jobId, result) => {
  fs.mkdirSync(resultsDir, { recursive: true })
  fs.writeFileSync(path.join(resultsDir, `${jobId}.json`), JSON.stringify(result), 'utf-8')
}

describe('HistoryStore', () => {
  let context

  beforeEach(() => {
    context = makeStore()
  })

  afterEach(() => {
    fs.rmSync(context.baseDir, { recursive: true, force: true })
  })

  describe('normalize', () => {
    it('maps a legacy entry (startedAt, no name) onto the canonical shape', () => {
      const normalized = context.store.normalize({
        id: 'legacy',
        scriptName: 'analysis.py',
        startedAt: 1_600_000_000_000,
        payloadSnapshot: { params: { conditionA: 'Ctrl', conditionB: 'Treat' } },
      })
      expect(normalized.createdAt).toBe(1_600_000_000_000)
      expect(normalized.updatedAt).toBe(1_600_000_000_000)
      expect(normalized.runCount).toBe(1)
      expect(normalized.analysisName).toBe('Ctrl vs Treat')
    })

    it('rejects entries without an identifier', () => {
      expect(context.store.normalize({ analysisName: 'no id' })).toBeNull()
      expect(context.store.normalize(null)).toBeNull()
    })
  })

  describe('upsert', () => {
    it('inserts a new entry and returns the archive sorted by last update', async () => {
      await context.store.upsert(entry({ id: 'old', updatedAt: 1_000 }))
      await context.store.upsert(entry({ id: 'new', updatedAt: 2_000 }))

      const entries = await context.store.list()
      expect(entries.map((item) => item.id)).toEqual(['new', 'old'])
    })

    it('keeps the original creation time and counts the re-runs', async () => {
      const first = await context.store.upsert(entry())
      expect(first.entry.runCount).toBe(1)

      const second = await context.store.upsert(entry({ createdAt: 1_800_000_000_000, updatedAt: 1_800_000_000_000 }))
      expect(second.entry.runCount).toBe(2)
      expect(second.entry.createdAt).toBe(1_700_000_000_000)
      expect(second.entry.updatedAt).toBe(1_800_000_000_000)
    })

    it('honours an explicit run count', async () => {
      await context.store.upsert(entry())
      const updated = await context.store.upsert(entry({ runCount: 7 }))
      expect(updated.entry.runCount).toBe(7)
    })

    it('caps the archive size', async () => {
      for (let index = 0; index < MAX_ENTRIES + 5; index += 1) {
        await context.store.upsert(entry({ id: `job-${index}`, updatedAt: 1_000 + index }))
      }
      const entries = await context.store.list()
      expect(entries).toHaveLength(MAX_ENTRIES)
      expect(entries[0].id).toBe(`job-${MAX_ENTRIES + 4}`)
    })
  })

  describe('update', () => {
    it('patches a field without touching the run count nor the creation time', async () => {
      await context.store.upsert(entry())
      const response = await context.store.update('job-1', { analysisName: 'Renamed', summary: { p_value: 0.01 } })

      expect(response.success).toBe(true)
      expect(response.entry.analysisName).toBe('Renamed')
      expect(response.entry.summary).toEqual({ p_value: 0.01 })
      expect(response.entry.runCount).toBe(1)
      expect(response.entry.createdAt).toBe(1_700_000_000_000)
    })

    it('reports an unknown id', async () => {
      const response = await context.store.update('missing', { analysisName: 'x' })
      expect(response.success).toBe(false)
      expect(response.error).toMatch(/not found/i)
    })
  })

  describe('importLegacy', () => {
    it('imports unknown entries once and is idempotent', async () => {
      const legacy = [
        { id: 'job-a', analysisName: 'prova', startedAt: 1_600_000_000_000, payloadSnapshot: { files: {}, params: {} } },
        { id: 'job-b', scriptName: 'analysis.py', startedAt: 1_600_000_001_000, payloadSnapshot: { files: {}, params: { conditionA: 'X' } } },
      ]

      const first = await context.store.importLegacy(legacy)
      expect(first.imported).toBe(2)

      const second = await context.store.importLegacy(legacy)
      expect(second.imported).toBe(0)
      expect(await context.store.list()).toHaveLength(2)

      const entries = await context.store.list()
      expect(entries.find((item) => item.id === 'job-a').analysisName).toBe('prova')
      // An entry named after the script gets a readable label instead.
      expect(entries.find((item) => item.id === 'job-b').analysisName).toBe('X')
    })

    it('keeps the existing entry when an id is already archived', async () => {
      await context.store.upsert(entry({ analysisName: 'Current name' }))
      const response = await context.store.importLegacy([entry({ analysisName: 'Legacy name' })])
      expect(response.imported).toBe(0)
      expect((await context.store.list())[0].analysisName).toBe('Current name')
    })

    it('ignores an empty payload', async () => {
      const response = await context.store.importLegacy([])
      expect(response).toEqual({ success: true, imported: 0, entries: [] })
    })
  })

  describe('rebuildFromResults', () => {
    it('recovers entries from the persisted results', async () => {
      writeResult(context.resultsDir, 'job-1', {
        summary: { condition_a: 'Ctrl', condition_b: 'Treat', correlation_method: 'spearman' },
        input: {
          analysisName: 'Recovered analysis',
          files: { bacteria: ['/data/a.csv'] },
          params: { conditionA: 'Ctrl', conditionB: 'Treat' },
          completedAt: 1_650_000_000_000,
        },
      })

      const response = await context.store.rebuildFromResults()
      expect(response.recovered).toBe(1)

      const [recovered] = await context.store.list()
      expect(recovered.id).toBe('job-1')
      expect(recovered.analysisName).toBe('Recovered analysis')
      expect(recovered.createdAt).toBe(1_650_000_000_000)
      expect(recovered.payloadSnapshot.files.bacteria).toEqual(['/data/a.csv'])
      expect(recovered.summary).toEqual({ condition_a: 'Ctrl', condition_b: 'Treat', correlation_method: 'spearman' })
    })

    it('falls back to the result summary when the input block is missing', async () => {
      writeResult(context.resultsDir, 'job-legacy', {
        summary: { condition_a: 'Ctrl', condition_b: 'Treat', correlation_method: 'sparcc' },
      })

      await context.store.rebuildFromResults()
      const [recovered] = await context.store.list()
      expect(recovered.analysisName).toBe('Ctrl vs Treat')
      expect(recovered.payloadSnapshot.params.correlationMethod).toBe('sparcc')
    })

    it('backfills the summary of an already archived entry', async () => {
      await context.store.upsert(entry({ summary: null }))
      writeResult(context.resultsDir, 'job-1', { summary: { p_value: 0.05 } })

      const response = await context.store.rebuildFromResults()
      expect(response.recovered).toBe(0)
      expect(response.enriched).toBe(1)
      expect((await context.store.list())[0].summary).toEqual({ p_value: 0.05 })
    })

    it('does nothing when there is no results directory', async () => {
      const response = await context.store.rebuildFromResults()
      expect(response.success).toBe(true)
      expect(response.recovered).toBe(0)
    })
  })

  describe('ensureInitialized', () => {
    it('rebuilds from the results only when the archive is empty', async () => {
      writeResult(context.resultsDir, 'job-1', {
        summary: { condition_a: 'A', condition_b: 'B' },
        input: { params: { conditionA: 'A', conditionB: 'B' } },
      })

      const first = await context.store.ensureInitialized()
      expect(first.rebuilt).toBe(true)
      expect(first.recovered).toBe(1)

      const second = await context.store.ensureInitialized()
      expect(second.rebuilt).toBe(false)
      expect(second.entries).toHaveLength(1)
    })
  })

  describe('robustness', () => {
    it('treats a corrupt archive as empty instead of failing', async () => {
      fs.writeFileSync(context.store.historyPath, '{ this is not json', 'utf-8')
      await expect(context.store.list()).resolves.toEqual([])
    })

    it('removes and clears entries', async () => {
      await context.store.upsert(entry())
      const removed = await context.store.remove('job-1')
      expect(removed.entries).toHaveLength(0)

      await context.store.upsert(entry())
      const cleared = await context.store.clear()
      expect(cleared.entries).toHaveLength(0)
    })
  })
})
