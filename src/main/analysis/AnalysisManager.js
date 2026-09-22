/**
 * AnalysisManager
 * Spawns and manages the Python child process for analysis jobs.
 */
const { spawn } = require('child_process')
const path = require('path')
const fs = require('fs').promises
const fsSync = require('fs')
const crypto = require('crypto')
const { app, dialog } = require('electron')
const os = require('os')
const log = require('electron-log')
const { getPythonPath, resolveFastsparPath, resolveFastsparCompanions, getSparccRuntimeStatus } = require('./runtimeDetector')

class AnalysisManager {
  constructor() {
    /** @type {Map<string, ChildProcess>} */
    this.activeJobs = new Map()
    /** @type {Map<string, object>} */
    this.results = new Map()
    /** @type {Map<string, { cacheKey: string, path: string, version: number, jobId: string, condition: string }>} */
    this.graphCacheIndex = new Map()

    this.cacheVersion = 1
    this.pythonDepsReadyByInterpreter = new Map()
    /** @type {string|null} Interpreter resolved for this session */
    this.pythonInterpreter = null
    this.cacheDir = path.join(app.getPath('userData'), 'analysis-cache', 'graphs')
    this.cacheIndexPath = path.join(app.getPath('userData'), 'analysis-cache', 'index.json')
    this.resultsDir = path.join(app.getPath('userData'), 'analysis-results')
    this.clusterSnapshotsDir = path.join(this.resultsDir, 'clusters')

    this._loadGraphCacheIndex().catch((error) => {
      log.warn('Failed to load graph cache index:', error)
    })
  }

  async _ensureCacheDir() {
    await fs.mkdir(this.cacheDir, { recursive: true })
  }

  async _ensureResultsDir() {
    await fs.mkdir(this.resultsDir, { recursive: true })
  }

  _resultDiskPath(jobId) {
    return path.join(this.resultsDir, `${jobId}.json`)
  }

  async _persistResultToDisk(jobId, result) {
    await this._ensureResultsDir()
    const filePath = this._resultDiskPath(jobId)
    await fs.writeFile(filePath, JSON.stringify(result, null, 2), 'utf-8')
  }

  async _readResultFromDisk(jobId) {
    const filePath = this._resultDiskPath(jobId)
    const raw = await fs.readFile(filePath, 'utf-8')
    return JSON.parse(raw)
  }

  _clusterSnapshotPath(jobId) {
    return path.join(this.clusterSnapshotsDir, `${jobId}.json`)
  }

  async _readClusterSnapshot(jobId) {
    const filePath = this._clusterSnapshotPath(jobId)
    const raw = await fs.readFile(filePath, 'utf-8')
    const parsed = JSON.parse(raw)
    return parsed && typeof parsed === 'object' ? parsed : {}
  }

  async _saveClusterSnapshot(jobId, cacheKey, payload) {
    await fs.mkdir(this.clusterSnapshotsDir, { recursive: true })
    const filePath = this._clusterSnapshotPath(jobId)
    let current = {}
    try {
      current = await this._readClusterSnapshot(jobId)
    } catch {
      current = {}
    }
    current[cacheKey] = payload
    await fs.writeFile(filePath, JSON.stringify(current, null, 2), 'utf-8')
    return current
  }

  _datasetFingerprint(payload) {
    const files = payload.files || {}
    const flatFiles = [
      ...(files.bacteria || []),
      ...(files.fungi || []),
      ...(files.metadata || []),
    ]
    const basis = JSON.stringify(flatFiles.sort())
    return crypto.createHash('sha1').update(basis).digest('hex')
  }

  _graphCacheKey({ jobId, condition, correlationMethod, pValueThreshold, datasetFingerprint, datasetType }) {
    const raw = [jobId, condition, correlationMethod, pValueThreshold, datasetFingerprint, datasetType || 'unknown'].join('|')
    return crypto.createHash('sha1').update(raw).digest('hex')
  }

  async _loadGraphCacheIndex() {
    try {
      const raw = await fs.readFile(this.cacheIndexPath, 'utf-8')
      const parsed = JSON.parse(raw)
      Object.entries(parsed).forEach(([key, value]) => {
        this.graphCacheIndex.set(key, value)
      })
    } catch {
      // Missing index is expected on first run.
    }
  }

  async _persistGraphCacheIndex() {
    await fs.mkdir(path.dirname(this.cacheIndexPath), { recursive: true })
    const data = Object.fromEntries(this.graphCacheIndex.entries())
    await fs.writeFile(this.cacheIndexPath, JSON.stringify(data, null, 2), 'utf-8')
  }

  async _persistConditionGraphs(jobId, payload, result) {
    const pipeline = result?.graph_pipeline
    const conditionGraphs = pipeline?.condition_graphs || []
    if (!conditionGraphs.length) {
      return { graphRefs: [] }
    }

    await this._ensureCacheDir()
    const pValueThreshold = result?.summary?.p_value_threshold ?? payload?.params?.pValueThreshold ?? 0.05
    const correlationMethod = payload?.params?.correlationMethod || 'spearman'
    const datasetType = payload?.params?.datasetType || 'unknown'
    const datasetFingerprint = this._datasetFingerprint(payload)

    const graphRefs = []

    for (const conditionGraph of conditionGraphs) {
      const condition = conditionGraph.condition
      const cacheKey = this._graphCacheKey({
        jobId,
        condition,
        correlationMethod,
        pValueThreshold,
        datasetFingerprint,
        datasetType,
      })
      const cachePath = path.join(this.cacheDir, `${cacheKey}.json`)
      const serialized = {
        version: this.cacheVersion,
        cacheKey,
        jobId,
        condition,
        correlationMethod,
        datasetType,
        pValueThreshold,
        datasetFingerprint,
        createdAt: new Date().toISOString(),
        graph: conditionGraph.graph,
        pairFile: conditionGraph.pair_file,
      }
      await fs.writeFile(cachePath, JSON.stringify(serialized), 'utf-8')

      const ref = {
        cacheKey,
        path: cachePath,
        version: this.cacheVersion,
        jobId,
        condition,
      }
      this.graphCacheIndex.set(cacheKey, ref)
      graphRefs.push(ref)
    }

    await this._persistGraphCacheIndex()
    return { graphRefs }
  }

  async _readCachedGraph(cacheKey) {
    const ref = this.graphCacheIndex.get(cacheKey)
    if (!ref) {
      throw new Error(`Graph cache miss: ${cacheKey}`)
    }

    const raw = await fs.readFile(ref.path, 'utf-8')
    const parsed = JSON.parse(raw)
    if (parsed.version !== this.cacheVersion) {
      throw new Error(`Graph cache version mismatch for key ${cacheKey}`)
    }
    if (!parsed.graph || !Array.isArray(parsed.graph.nodes) || !Array.isArray(parsed.graph.edges)) {
      throw new Error(`Graph cache invalid payload for key ${cacheKey}`)
    }
    return parsed
  }

  /**
   * Directory containing the Python analysis scripts.
   *
   * Packaged builds ship them through `extraResources` at `<resources>/python`.
   * The copy that would live inside `app.asar` must never be used: app.asar is a
   * virtual archive that only Electron's patched `fs` can read, so handing such
   * a path to an external Python process fails with ENOTDIR.
   * @returns {string}
   */
  _pythonRoot() {
    if (app.isPackaged) {
      const packagedRoot = path.join(process.resourcesPath, 'python')
      if (fsSync.existsSync(path.join(packagedRoot, 'analysis.py'))) {
        return packagedRoot
      }
      throw new Error(
        `Python analysis scripts are missing from the installation (expected in ${packagedRoot}). ` +
        'The build is incomplete — reinstall the application.'
      )
    }
    return path.join(__dirname, '../../../python')
  }

  async _runPythonClustering({ graphData, seedNode, options = {} }) {
    const runId = crypto.randomUUID()
    const inputPath = path.join(os.tmpdir(), `cluster_input_${runId}.json`)
    const pythonScript = path.join(this._pythonRoot(), 'clustering', 'run_clustering.py')
    const interpreter = await getPythonPath()

    const payload = {
      graphData,
      seedNode,
      targetConductance: options.targetConductance ?? 0.2,
      b: Number.isFinite(options.b) ? options.b : null,
      teleportAlpha: Number.isFinite(options.teleportAlpha) ? options.teleportAlpha : null,
      weighted: options.weighted ?? false,
    }

    await fs.writeFile(inputPath, JSON.stringify(payload), 'utf-8')

    return new Promise((resolve, reject) => {
      const proc = spawn(interpreter, [pythonScript, '--input', inputPath])
      let stdout = ''
      let stderr = ''

      proc.stdout.on('data', (data) => {
        stdout += data.toString()
      })

      proc.stderr.on('data', (data) => {
        stderr += data.toString()
      })

      proc.on('close', async (code) => {
        await fs.unlink(inputPath).catch(() => {})
        if (code !== 0) {
          return reject(new Error(stderr || `Clustering process exited with code ${code}`))
        }

        const lines = stdout
          .split('\n')
          .map((line) => line.trim())
          .filter(Boolean)
        const lastLine = lines[lines.length - 1]
        if (!lastLine) {
          return reject(new Error('Clustering process returned empty output'))
        }

        try {
          resolve(JSON.parse(lastLine))
        } catch {
          reject(new Error(`Unable to parse clustering output: ${lastLine}`))
        }
      })

      proc.on('error', async (error) => {
        await fs.unlink(inputPath).catch(() => {})
        reject(error)
      })
    })
  }

  _runCommand(command, args, options = {}) {
    return new Promise((resolve, reject) => {
      const proc = spawn(command, args, options)
      let stdout = ''
      let stderr = ''

      proc.stdout.on('data', (data) => {
        stdout += data.toString()
      })

      proc.stderr.on('data', (data) => {
        stderr += data.toString()
      })

      proc.on('close', (code) => {
        resolve({ code, stdout, stderr })
      })

      proc.on('error', (error) => {
        reject(error)
      })
    })
  }

  /**
   * Fingerprint of the bundled requirements file, used to cache dependency checks.
   * @returns {Promise<string|null>}
   */
  async _requirementsFingerprint() {
    const requirementsPath = path.join(this._pythonRoot(), 'requirements.txt')
    try {
      const stat = await fs.stat(requirementsPath)
      return `${requirementsPath}:${stat.mtimeMs}:${stat.size}`
    } catch {
      return null
    }
  }

  async _hasRequiredPackages(interpreter) {
    const check = await this._runCommand(interpreter, ['-c', 'import scipy, networkx, numpy'])
    return check.code === 0
  }

  /**
   * Create (once) and populate a Python environment owned by the application.
   *
   * Distributions that mark their system interpreter as "externally managed"
   * (PEP 668 — Arch, Debian 12+, Ubuntu 23+, Fedora) reject `pip install`, so
   * installing into the system Python is not an option. A virtual environment in
   * the application's user-data directory works on every platform, needs no
   * administrator rights, and leaves the user's Python installation untouched.
   *
   * @param {string} baseInterpreter
   * @param {(line: string) => void} emitLog
   * @returns {Promise<string|null>} the environment interpreter, or null on failure
   */
  async _ensureManagedEnvironment(baseInterpreter, emitLog) {
    const envDir = path.join(app.getPath('userData'), 'python-env')
    const envPython = process.platform === 'win32'
      ? path.join(envDir, 'Scripts', 'python.exe')
      : path.join(envDir, 'bin', 'python')
    const markerPath = path.join(envDir, '.requirements-fingerprint')

    const fingerprint = await this._requirementsFingerprint()
    if (!fingerprint) return null

    // Fast path: the environment exists and was built from the same requirements
    // file, so there is nothing to install. Skipping pip keeps the first analysis
    // of every session close to instant.
    if (fsSync.existsSync(envPython)) {
      let storedFingerprint = null
      try {
        storedFingerprint = (await fs.readFile(markerPath, 'utf-8')).trim()
      } catch {
        // Missing marker: fall through and reconcile the environment.
      }
      if (storedFingerprint === fingerprint && await this._hasRequiredPackages(envPython)) {
        this.pythonDepsReadyByInterpreter.set(envPython, fingerprint)
        return envPython
      }
    }

    if (!fsSync.existsSync(envPython)) {
      emitLog(`First run: creating a private Python environment in ${envDir} …`)
      const created = await this._runCommand(baseInterpreter, ['-m', 'venv', envDir])
      if (created.code !== 0) {
        log.warn(`Could not create a virtual environment: ${created.stderr || created.stdout}`)
        return null
      }
    }

    const requirementsPath = path.join(this._pythonRoot(), 'requirements.txt')
    emitLog('Installing the required Python packages (this happens only once) …')
    const installed = await this._runCommand(envPython, [
      '-m', 'pip', 'install',
      '--disable-pip-version-check',
      '--upgrade',
      '-r', requirementsPath,
    ])
    if (installed.code !== 0) {
      log.warn(`Could not install Python dependencies: ${installed.stderr || installed.stdout}`)
      return null
    }

    if (!(await this._hasRequiredPackages(envPython))) {
      log.warn('Python dependencies are still missing after installation')
      return null
    }

    await fs.writeFile(markerPath, fingerprint, 'utf-8').catch(() => {})
    this.pythonDepsReadyByInterpreter.set(envPython, fingerprint)
    emitLog('Python environment ready.')
    return envPython
  }

  /**
   * Resolve an interpreter that can actually run the analysis, installing what
   * is missing when possible.
   * @param {(line: string) => void} emitLog
   * @returns {Promise<string>}
   */
  async _resolvePythonInterpreter(emitLog = () => {}) {
    if (this.pythonInterpreter) return this.pythonInterpreter

    const base = await getPythonPath()

    if (await this._hasRequiredPackages(base)) {
      this.pythonInterpreter = base
      return base
    }

    log.warn(`Python dependencies missing for interpreter ${base}. Setting up a managed environment.`)

    const managed = await this._ensureManagedEnvironment(base, emitLog)
    if (managed) {
      this.pythonInterpreter = managed
      return managed
    }

    throw new Error(
      'The analysis packages (scipy, numpy, networkx) are unavailable and could not be installed automatically. ' +
      `Tried interpreter: ${base}. ` +
      'Connect to the internet and retry — or create an environment with ' +
      '"python3 -m venv <folder>", run "pip install scipy numpy networkx" inside it, ' +
      'then select it in Settings → Python.'
    )
  }

  /**
   * Run an analysis job.
   *
   * `payload.jobId` may be supplied to re-run an existing history entry: the new
   * results then overwrite that entry's artifacts instead of creating an orphan
   * job, which is what makes a stored analysis fully editable.
   *
   * @param {{ scriptName: string, files: string[], params: object, jobId?: string, analysisName?: string }} payload
   * @param {(progress: object) => void} onProgress
   * @param {(line: string) => void} onLog
   * @returns {Promise<string>} jobId
   */
  async runAnalysis(payload, onProgress, onLog) {
    const requestedJobId = String(payload?.jobId || '').trim()
    if (requestedJobId && this.activeJobs.has(requestedJobId)) {
      throw new Error('This analysis is already running.')
    }
    const jobId = requestedJobId || crypto.randomUUID()
    const { scriptName, files, params } = payload

    // Write params to a temp JSON file so the script can read them
    const paramsPath = path.join(require('os').tmpdir(), `params_${jobId}.json`)
    await fs.writeFile(paramsPath, JSON.stringify({ files, params, jobId }))

    const scriptPath = path.join(this._pythonRoot(), scriptName)

    const requestedMethod = (params?.correlationMethod || 'spearman').toLowerCase()
    if (requestedMethod === 'sparcc') {
      const resolved = resolveFastsparPath({ customPath: params?.fastsparPath || '' })
      const companions = resolved.available
        ? resolveFastsparCompanions(resolved.path)
        : { fastspar_bootstrap: '', fastspar_pvalues: '' }
      const missing = ['fastspar_bootstrap', 'fastspar_pvalues'].filter((name) => !companions[name])
      if (!resolved.available || missing.length > 0) {
        const detail = missing.length > 0
          ? `missing companion binaries: ${missing.join(', ')}`
          : 'fastspar binary not found'
        throw new Error(
          `SparCC requires the FastSpar runtime (${detail}). ` +
          'Configure FastSpar in Settings or ship bundled binaries under resources/bin/<platform>/<arch>/.'
        )
      }
      params.fastsparPath = resolved.path
    }
    const interpreter = await this._resolvePythonInterpreter((line) => onLog({ jobId, line }))

    const args = [scriptPath, '--params', paramsPath]

    log.info(
      `[Job ${jobId}] Starting python analysis: ${scriptName} ` +
      `(interpreter: ${interpreter}, script: ${scriptPath})`
    )
    onProgress({ jobId, status: 'running', percent: 0 })

    return new Promise((resolve, reject) => {
      const proc = spawn(interpreter, args)
      this.activeJobs.set(jobId, proc)
      let stderr = ''

      proc.stdout.on('data', (data) => {
        const lines = data.toString().split('\n').filter(Boolean)
        lines.forEach((line) => {
          onLog({ jobId, line })
          // Scripts can emit JSON progress lines: {"progress": 42}
          try {
            const json = JSON.parse(line)
            if (json.progress !== undefined) {
              onProgress({ jobId, status: 'running', percent: json.progress })
            }
          } catch {
            // not JSON, just a log line
          }
        })
      })

      proc.stderr.on('data', (data) => {
        const chunk = data.toString()
        stderr += chunk
        onLog({ jobId, line: `[stderr] ${chunk}`, level: 'warn' })
      })

      proc.on('close', async (code) => {
        this.activeJobs.delete(jobId)
        if (code === 0) {
          try {
            const resultPath = path.join(require('os').tmpdir(), `result_${jobId}.json`)
            const raw = await fs.readFile(resultPath, 'utf-8')
            const parsed = JSON.parse(raw)
            const { graphRefs } = await this._persistConditionGraphs(jobId, payload, parsed)

            parsed.graph_pipeline = parsed.graph_pipeline || {
              version: this.cacheVersion,
              status: 'not_available',
              condition_graphs: [],
            }
            parsed.graph_pipeline.cache_version = this.cacheVersion
            parsed.graph_pipeline.graph_refs = graphRefs
            // Persist the exact inputs next to the outputs: this is what allows
            // the history archive to be rebuilt without losing any information.
            parsed.input = {
              jobId,
              scriptName,
              analysisName: payload?.analysisName || '',
              files: files || {},
              params: params || {},
              completedAt: Date.now(),
            }
            this.results.set(jobId, parsed)
            await this._persistResultToDisk(jobId, parsed)
            onProgress({ jobId, status: 'completed', percent: 100 })
            resolve(jobId)
          } catch (err) {
            reject(new Error(`Analysis completed but result file not found: ${err.message}`))
          }
        } else {
          onProgress({ jobId, status: 'failed', percent: 0 })
          const stderrTail = stderr
            .split('\n')
            .map((line) => line.trim())
            .filter(Boolean)
            .slice(-8)
            .join(' | ')
          const message = stderrTail
            ? `Process exited with code ${code}: ${stderrTail}`
            : `Process exited with code ${code}`
          reject(new Error(message))
        }
        await fs.unlink(paramsPath).catch(() => {})
      })

      proc.on('error', (err) => {
        this.activeJobs.delete(jobId)
        reject(err)
      })
    })
  }

  cancelJob(jobId) {
    const proc = this.activeJobs.get(jobId)
    if (proc) {
      proc.kill('SIGTERM')
      this.activeJobs.delete(jobId)
      log.info(`[Job ${jobId}] Cancelled`)
    }
  }

  async getResults(jobId) {
    const inMemory = this.results.get(jobId)
    if (inMemory) {
      try {
        const snapshot = await this._readClusterSnapshot(jobId)
        inMemory.graph_pipeline = inMemory.graph_pipeline || {}
        inMemory.graph_pipeline.last_clustering = {
          ...(inMemory.graph_pipeline.last_clustering || {}),
          ...snapshot,
        }
      } catch {
        // no persisted cluster snapshots available
      }
      return inMemory
    }

    try {
      const fromDisk = await this._readResultFromDisk(jobId)
      try {
        const snapshot = await this._readClusterSnapshot(jobId)
        fromDisk.graph_pipeline = fromDisk.graph_pipeline || {}
        fromDisk.graph_pipeline.last_clustering = {
          ...(fromDisk.graph_pipeline.last_clustering || {}),
          ...snapshot,
        }
      } catch {
        // no persisted cluster snapshots available
      }
      this.results.set(jobId, fromDisk)
      return fromDisk
    } catch {
      return null
    }
  }

  async getConditionGraph(cacheKey) {
    try {
      const cached = await this._readCachedGraph(cacheKey)
      return { success: true, graph: cached }
    } catch (error) {
      return {
        success: false,
        error: `${error.message}. Please rerun the analysis to rebuild graph cache.`,
      }
    }
  }

  async reclusterConditionGraph({ cacheKey, seedNode, options = {} }) {
    try {
      const cached = await this._readCachedGraph(cacheKey)
      const clusteringResult = await this._runPythonClustering({
        graphData: cached,
        seedNode,
        options,
      })
      if (!clusteringResult.success) {
        throw new Error(clusteringResult.error || 'Clustering failed')
      }

      const response = {
        success: true,
        cacheKey,
        seedNode,
        condition: cached.condition,
        algorithm: options.algorithm || 'pagerank_nibble_adaptive',
        graphVersion: cached.version,
        clusters: [clusteringResult.clusterNodes],
        clusterNodes: clusteringResult.clusterNodes,
        subgraph: clusteringResult.subgraph,
        metadata: clusteringResult.metadata,
        message: 'Cluster computed from cached graph without recomputing correlations.',
      }

      const jobId = cached?.jobId
      if (jobId) {
        const existing = await this.getResults(jobId)
        if (existing) {
          const mergedSnapshots = await this._saveClusterSnapshot(jobId, cacheKey, response)
          existing.graph_pipeline = existing.graph_pipeline || {}
          existing.graph_pipeline.last_clustering = {
            ...(existing.graph_pipeline.last_clustering || {}),
            ...mergedSnapshots,
          }
          this.results.set(jobId, existing)
          await this._persistResultToDisk(jobId, existing)
        }
      }

      return response
    } catch (error) {
      return {
        success: false,
        error: `${error.message}. Please rerun the analysis to rebuild graph cache.`,
      }
    }
  }

  async deleteAnalysis(jobId) {
    if (!jobId) {
      return { success: false, error: 'Missing jobId' }
    }

    this.results.delete(jobId)

    const resultPath = this._resultDiskPath(jobId)
    await fs.unlink(resultPath).catch(() => {})
    const clusterSnapshotPath = this._clusterSnapshotPath(jobId)
    await fs.unlink(clusterSnapshotPath).catch(() => {})

    const keysToDelete = []
    for (const [cacheKey, ref] of this.graphCacheIndex.entries()) {
      if (ref?.jobId === jobId) {
        keysToDelete.push(cacheKey)
      }
    }

    for (const cacheKey of keysToDelete) {
      const ref = this.graphCacheIndex.get(cacheKey)
      if (ref?.path) {
        await fs.unlink(ref.path).catch(() => {})
      }
      this.graphCacheIndex.delete(cacheKey)
    }

    await this._persistGraphCacheIndex().catch(() => {})
    return { success: true }
  }

  /**
   * Serialize the results table to CSV. Values are quoted only when needed and
   * the union of all row keys is used as the header, so rows with heterogeneous
   * keys survive the round-trip.
   * @param {Array<object>} rows
   * @returns {string}
   */
  _tableToCsv(rows) {
    if (!Array.isArray(rows) || rows.length === 0) return ''

    const columns = [...new Set(rows.flatMap((row) => Object.keys(row || {})))]
    const escape = (value) => {
      if (value === null || value === undefined) return ''
      const text = typeof value === 'object' ? JSON.stringify(value) : String(value)
      return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text
    }

    return [
      columns.map(escape).join(','),
      ...rows.map((row) => columns.map((column) => escape((row || {})[column])).join(',')),
    ].join('\r\n')
  }

  /**
   * Export a completed job's results to a user-chosen file.
   * @param {string} jobId
   * @param {'csv'|'json'} format
   * @param {import('electron').BrowserWindow=} browserWindow Parent for the save dialog
   */
  async exportResults(jobId, format, browserWindow = null) {
    const normalized = String(format || '').toLowerCase()
    if (!['csv', 'json'].includes(normalized)) {
      return { success: false, error: `Unsupported export format: ${format}` }
    }

    const results = await this.getResults(jobId)
    if (!results) {
      return { success: false, error: 'No results for jobId' }
    }

    const content = normalized === 'json'
      ? JSON.stringify(results, null, 2)
      : this._tableToCsv(results.table)

    if (!content) {
      return { success: false, error: 'This analysis produced no rows to export.' }
    }

    const { canceled, filePath } = await dialog.showSaveDialog(browserWindow || undefined, {
      title: 'Export results',
      defaultPath: `analysis-${jobId}.${normalized}`,
      filters: normalized === 'json'
        ? [{ name: 'JSON', extensions: ['json'] }]
        : [{ name: 'CSV', extensions: ['csv'] }],
    })

    if (canceled || !filePath) {
      return { success: false, error: 'cancelled' }
    }

    await fs.writeFile(filePath, content, 'utf-8')
    log.info(`[Job ${jobId}] Exported results as ${normalized} to ${filePath}`)

    return { success: true, path: filePath, format: normalized }
  }

  getSparccRuntimeStatus({ customFastsparPath = '' } = {}) {
    return getSparccRuntimeStatus({ customPath: customFastsparPath })
  }
}

module.exports = { AnalysisManager }
