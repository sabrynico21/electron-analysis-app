const { contextBridge, ipcRenderer, webUtils } = require('electron')

contextBridge.exposeInMainWorld('electronAPI', {
  // File operations
  openFileDialog: (options) => ipcRenderer.invoke('dialog:openFile', options),
  saveFileDialog: (options) => ipcRenderer.invoke('dialog:saveFile', options),
  readFile: (filePath) => ipcRenderer.invoke('fs:readFile', filePath),
  getFileStats: (filePath) => ipcRenderer.invoke('fs:stat', filePath),
  describeFiles: (filePaths) => ipcRenderer.invoke('fs:describeFiles', filePaths),

  /**
   * Resolve the absolute path of a File coming from a drag & drop event.
   *
   * `File.path` was removed in Electron 32, which silently turned every dropped
   * file into an empty path and made the app fail with "file not found". The
   * supported replacement is `webUtils.getPathForFile`, which can only run in the
   * preload/renderer realm (a `File` object cannot cross the IPC boundary), so it
   * is exposed as a plain synchronous function.
   * @param {File} file
   * @returns {string} absolute path, or '' when it cannot be resolved
   */
  getPathForFile: (file) => {
    try {
      if (!file) return ''
      return webUtils.getPathForFile(file) || ''
    } catch {
      return ''
    }
  },

  // Analysis
  runAnalysis: (payload) => ipcRenderer.invoke('analysis:run', payload),
  cancelAnalysis: (jobId) => ipcRenderer.invoke('analysis:cancel', jobId),

  // Progress events (subscribe)
  onAnalysisProgress: (callback) =>
    ipcRenderer.on('analysis:progress', (_event, data) => callback(data)),
  onAnalysisLog: (callback) =>
    ipcRenderer.on('analysis:log', (_event, data) => callback(data)),

  // Results
  getResults: (jobId) => ipcRenderer.invoke('results:get', jobId),
  getConditionGraph: (cacheKey) => ipcRenderer.invoke('results:getConditionGraph', cacheKey),
  reclusterConditionGraph: (payload) => ipcRenderer.invoke('results:reclusterConditionGraph', payload),
  exportResults: (jobId, format) => ipcRenderer.invoke('results:export', jobId, format),
  deleteAnalysis: (jobId) => ipcRenderer.invoke('results:deleteAnalysis', jobId),

  // History archive (persisted in the OS application-data directory)
  listHistory: () => ipcRenderer.invoke('history:list'),
  upsertHistory: (entry) => ipcRenderer.invoke('history:upsert', entry),
  updateHistory: (id, patch) => ipcRenderer.invoke('history:update', id, patch),
  deleteHistory: (jobId) => ipcRenderer.invoke('history:delete', jobId),
  clearHistory: () => ipcRenderer.invoke('history:clear'),
  importHistory: (entries) => ipcRenderer.invoke('history:import', entries),
  rebuildHistory: () => ipcRenderer.invoke('history:rebuild'),
  getHistoryInfo: () => ipcRenderer.invoke('history:info'),

  // Settings
  getSettings: () => ipcRenderer.invoke('settings:get'),
  saveSettings: (settings) => ipcRenderer.invoke('settings:save', settings),
  getSparccStatus: () => ipcRenderer.invoke('sparcc:status'),

  // Cleanup
  removeAllListeners: (channel) => ipcRenderer.removeAllListeners(channel),
})
