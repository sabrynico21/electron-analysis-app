const { ipcMain, dialog, BrowserWindow } = require('electron')
const fs = require('fs').promises
const path = require('path')
const log = require('electron-log')
const { AnalysisManager } = require('../analysis/AnalysisManager')
const { HistoryStore } = require('../analysis/historyStore')
const Store = require('electron-store')

const store = new Store()
const analysisManager = new AnalysisManager()
// The archive lives next to the persisted results so both stay in sync.
const historyStore = new HistoryStore({ resultsDir: analysisManager.resultsDir })

/**
 * Inspect a list of raw paths coming from a drag & drop or from the native
 * "open file" dialog. Returns one descriptor per path so the renderer never has
 * to guess whether an item is a real file, and never has to build a path itself.
 * @param {string[]} filePaths
 */
async function describeFiles(filePaths = []) {
  const unique = [...new Set((Array.isArray(filePaths) ? filePaths : []).filter(Boolean))]
  return Promise.all(
    unique.map(async (filePath) => {
      const baseName = path.basename(filePath)
      try {
        const stats = await fs.stat(filePath)
        return {
          path: filePath,
          name: baseName,
          size: stats.size,
          isFile: stats.isFile(),
          isDirectory: stats.isDirectory(),
          modifiedAt: stats.mtimeMs,
          exists: true,
        }
      } catch (error) {
        log.warn(`fs:describeFiles could not stat ${filePath}: ${error.message}`)
        return {
          path: filePath,
          name: baseName || filePath,
          size: 0,
          isFile: false,
          isDirectory: false,
          modifiedAt: null,
          exists: false,
          error: error.code === 'ENOENT' ? 'File not found' : error.message,
        }
      }
    })
  )
}

function setupIpcHandlers() {
  // --- File Dialogs ---
  ipcMain.handle('dialog:openFile', async (event, options = {}) => {
    const parent = BrowserWindow.fromWebContents(event.sender)
    const result = parent
      ? await dialog.showOpenDialog(parent, options)
      : await dialog.showOpenDialog(options)
    if (result.canceled) return result
    return { ...result, files: await describeFiles(result.filePaths) }
  })

  ipcMain.handle('dialog:saveFile', async (event, options = {}) => {
    const parent = BrowserWindow.fromWebContents(event.sender)
    return parent ? dialog.showSaveDialog(parent, options) : dialog.showSaveDialog(options)
  })

  // Validate/normalize a set of absolute paths (drag & drop goes through here).
  ipcMain.handle('fs:describeFiles', async (_event, filePaths) => describeFiles(filePaths))

  ipcMain.handle('fs:readFile', async (_event, filePath) => {
    try {
      if (!filePath) throw new Error('Missing file path')
      const resolvedPath = path.resolve(filePath)
      const stats = await fs.stat(resolvedPath)
      if (!stats.isFile()) {
        return { success: false, error: `${resolvedPath} is not a file` }
      }
      const content = await fs.readFile(resolvedPath, 'utf-8')
      return { success: true, content }
    } catch (error) {
      log.error('fs:readFile error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('fs:stat', async (_event, filePath) => {
    try {
      const resolvedPath = path.resolve(filePath)
      const stats = await fs.stat(resolvedPath)
      return { success: true, size: stats.size, modifiedAt: stats.mtimeMs }
    } catch (error) {
      log.error('fs:stat error:', error)
      return { success: false, error: error.message }
    }
  })

  // --- Analysis ---
  ipcMain.handle('analysis:run', async (event, payload) => {
    try {
      const settings = store.store || {}
      const withRuntimeSettings = {
        ...payload,
        params: {
          ...(payload?.params || {}),
          fastsparPath: settings.fastsparPath || '',
        },
      }
      const jobId = await analysisManager.runAnalysis(
        withRuntimeSettings,
        (progressData) => event.sender.send('analysis:progress', progressData),
        (logLine) => event.sender.send('analysis:log', logLine)
      )
      return { success: true, jobId }
    } catch (error) {
      log.error('analysis:run error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('analysis:cancel', async (_event, jobId) => {
    analysisManager.cancelJob(jobId)
    return { success: true }
  })

  // --- Results ---
  ipcMain.handle('results:get', async (_event, jobId) => {
    return await analysisManager.getResults(jobId)
  })

  ipcMain.handle('results:getConditionGraph', async (_event, cacheKey) => {
    return analysisManager.getConditionGraph(cacheKey)
  })

  ipcMain.handle('results:reclusterConditionGraph', async (_event, payload) => {
    return analysisManager.reclusterConditionGraph(payload)
  })

  ipcMain.handle('results:export', async (event, jobId, format) => {
    return analysisManager.exportResults(
      jobId,
      format,
      BrowserWindow.fromWebContents(event.sender)
    )
  })

  ipcMain.handle('results:deleteAnalysis', async (_event, jobId) => {
    return await analysisManager.deleteAnalysis(jobId)
  })

  // --- History (persistent archive in the OS application-data directory) ---
  ipcMain.handle('history:list', async () => {
    try {
      await historyStore.ensureInitialized()
      return { success: true, entries: await historyStore.list() }
    } catch (error) {
      log.error('history:list error:', error)
      return { success: false, error: error.message, entries: [] }
    }
  })

  ipcMain.handle('history:upsert', async (_event, entry) => {
    try {
      return await historyStore.upsert(entry)
    } catch (error) {
      log.error('history:upsert error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('history:update', async (_event, id, patch) => {
    try {
      return await historyStore.update(id, patch)
    } catch (error) {
      log.error('history:update error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('history:delete', async (_event, jobId) => {
    try {
      await analysisManager.deleteAnalysis(jobId)
      return await historyStore.remove(jobId)
    } catch (error) {
      log.error('history:delete error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('history:clear', async () => {
    try {
      return await historyStore.clear()
    } catch (error) {
      log.error('history:clear error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('history:import', async (_event, entries) => {
    try {
      return await historyStore.importLegacy(entries)
    } catch (error) {
      log.error('history:import error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('history:rebuild', async () => {
    try {
      return await historyStore.rebuildFromResults()
    } catch (error) {
      log.error('history:rebuild error:', error)
      return { success: false, error: error.message }
    }
  })

  ipcMain.handle('history:info', () => ({ success: true, ...historyStore.info() }))

  // --- Settings ---
  ipcMain.handle('settings:get', () => store.store)
  ipcMain.handle('settings:save', (_event, settings) => {
    store.set(settings)
    return { success: true }
  })
  ipcMain.handle('sparcc:status', async () => {
    const settings = store.store || {}
    return analysisManager.getSparccRuntimeStatus({
      customFastsparPath: settings.fastsparPath || '',
    })
  })
}

module.exports = { setupIpcHandlers }
