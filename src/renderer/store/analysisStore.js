import { create } from 'zustand'

/**
 * Analysis history archive.
 *
 * The history is **not** kept in `localStorage` any more: browser storage is
 * partitioned per origin, so the Vite dev server and the packaged build used to
 * see two different (and partially empty) histories, and clearing the renderer
 * profile wiped everything. The archive now lives in a plain JSON file inside the
 * OS application-data directory, managed by the main process
 * (`src/main/analysis/historyStore.js`).
 *
 * The legacy `localStorage` payload is still read once, so a history created by
 * an older version is imported instead of being lost.
 */

/** Key used by the previous Zustand `persist` middleware. */
const LEGACY_STORAGE_KEY = 'analysis-history'
const LEGACY_MIGRATION_KEY = 'analysis-history-migrated-to-archive'

function readLegacyHistory() {
  try {
    const raw = window.localStorage.getItem(LEGACY_STORAGE_KEY)
    if (!raw) return []
    const parsed = JSON.parse(raw)
    const history = parsed?.state?.history
    return Array.isArray(history) ? history : []
  } catch {
    return []
  }
}

function legacyAlreadyMigrated() {
  try {
    return window.localStorage.getItem(LEGACY_MIGRATION_KEY) === '1'
  } catch {
    return false
  }
}

function markLegacyMigrated() {
  try {
    window.localStorage.setItem(LEGACY_MIGRATION_KEY, '1')
  } catch {
    // Storage can be unavailable; the import is idempotent so this is harmless.
  }
}

const byRecency = (entries = []) =>
  [...entries].sort((a, b) => (b.updatedAt || b.createdAt || 0) - (a.updatedAt || a.createdAt || 0))

export const useAnalysisStore = create((set, get) => ({
  /** @type {object[]} newest first */
  history: [],
  isLoading: false,
  initialized: false,
  error: null,

  /**
   * Load the archive from disk, importing the legacy localStorage history the
   * first time the new storage is used.
   */
  loadHistory: async () => {
    set({ isLoading: true, error: null })
    try {
      let response = await window.electronAPI.listHistory()
      if (!response?.success) throw new Error(response?.error || 'Unable to load the history archive')

      if (!legacyAlreadyMigrated()) {
        const legacy = readLegacyHistory()
        if (legacy.length > 0) {
          const imported = await window.electronAPI.importHistory(legacy)
          if (imported?.success && imported.entries) {
            response = { success: true, entries: imported.entries }
          }
        }
        markLegacyMigrated()
      }

      set({
        history: byRecency(response.entries || []),
        isLoading: false,
        initialized: true,
        error: null,
      })
    } catch (error) {
      set({ isLoading: false, initialized: true, error: error?.message || 'Unable to load the history' })
    }
  },

  /** Insert or replace a history entry (used when a run completes). */
  addToHistory: async (entry) => {
    if (!entry?.id) return null
    try {
      const response = await window.electronAPI.upsertHistory(entry)
      if (response?.success && response.entries) {
        set({ history: byRecency(response.entries) })
      }
      return response
    } catch (error) {
      set({ error: error?.message || 'Unable to save the history entry' })
      return null
    }
  },

  /** Patch an entry (rename, refreshed parameters, new timestamp…). */
  updateHistoryEntry: async (id, patch) => {
    if (!id) return null
    try {
      const response = await window.electronAPI.updateHistory(id, patch)
      if (response?.success && response.entries) {
        set({ history: byRecency(response.entries) })
      }
      return response
    } catch (error) {
      set({ error: error?.message || 'Unable to update the history entry' })
      return null
    }
  },

  /** Delete an entry together with its results and cached graphs. */
  removeFromHistory: async (id) => {
    if (!id) return null
    try {
      const response = await window.electronAPI.deleteHistory(id)
      if (response?.success && response.entries) {
        set({ history: byRecency(response.entries) })
      }
      return response
    } catch (error) {
      set({ error: error?.message || 'Unable to delete the history entry' })
      return null
    }
  },

  clearHistory: async () => {
    try {
      const response = await window.electronAPI.clearHistory()
      if (response?.success) set({ history: byRecency(response.entries || []) })
      return response
    } catch (error) {
      set({ error: error?.message || 'Unable to clear the history' })
      return null
    }
  },

  /** Recover entries from the analysis results already stored on disk. */
  rebuildHistory: async () => {
    try {
      const response = await window.electronAPI.rebuildHistory()
      if (response?.success) {
        const list = await window.electronAPI.listHistory()
        set({ history: byRecency(list?.entries || []) })
      }
      return response
    } catch (error) {
      set({ error: error?.message || 'Unable to rebuild the history' })
      return null
    }
  },

  getHistoryEntry: (id) => get().history.find((entry) => entry.id === id) || null,
}))
