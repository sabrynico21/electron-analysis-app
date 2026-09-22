/**
 * AnalysisWorkspace
 *
 * The editable analysis area: upload → parameters → run → taxa → results.
 *
 * It is deliberately shared by two entry points:
 *  - `mode="new"`      → creating a brand new analysis (`/`);
 *  - `mode="history"`  → reopening a stored analysis (`/history/:jobId`).
 *
 * A stored analysis is fully editable: files, metadata, parameters, name and the
 * clustering inputs can all be changed and re-run. Re-running reuses the original
 * job id (`AnalysisManager.runAnalysis` accepts `payload.jobId`), so the history
 * entry is updated in place instead of being duplicated.
 */
import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import FileUploadZone from '../FileUpload/FileUploadZone'
import ParametersForm from '../Parameters/ParametersForm'
import AnalysisProgress from '../shared/AnalysisProgress'
import ChartPanel from '../Results/ChartPanel'
import SeedSelector from '../Results/SeedSelector'
import { useAnalysis } from '../../hooks/useAnalysis'
import { useAnalysisStore } from '../../store/analysisStore'
import { useResultsStore } from '../../store/resultsStore'
import { dedupeDescriptors, describePaths, parseMetadataConditions, toPathList } from '../../utils/fileLoading'
import styles from '../../styles/AnalysisPage.module.css'

const DEFAULT_PARAMS = {
  correlationMethod: 'spearman',
  groupingMode: 'all',
  pValueThreshold: 0.05,
}

const TELEPORT_MIN = 0

function computeDefaultB(nodeCount, edgeCount) {
  const n = Math.max(1, Number(nodeCount) || 1)
  const m = Math.max(1, Number(edgeCount) || 1)

  const bMax = Math.max(1, Math.ceil(Math.log2(Math.max(m, 1))))
  const bMin = Math.max(1, Math.floor(Math.log2(Math.sqrt(Math.max(n, 1)))))
  const avgDegree = (2 * m) / n

  let bDefault
  if (avgDegree < 10) {
    bDefault = Math.max(bMin, Math.floor(Math.log2(Math.max(n / 10, 1))))
  } else if (avgDegree > 50) {
    bDefault = Math.max(bMin, Math.floor(Math.log2(Math.max(n / 3, 1))))
  } else {
    bDefault = Math.max(bMin, Math.floor(Math.log2(Math.max(Math.sqrt(n), 1))))
  }

  return Math.min(bDefault, bMax)
}

function computeDefaultTeleportAlpha(targetConductance, edgeCount) {
  const phi = Number(targetConductance)
  const m = Math.max(1, Number(edgeCount) || 1)
  const denominator = 225 * Math.log(100 * Math.sqrt(m))
  if (!Number.isFinite(phi) || denominator <= 0) return TELEPORT_MIN
  return normalizeTeleport((phi * phi) / denominator)
}

function normalizeTeleport(value) {
  if (value === '') return TELEPORT_MIN
  const numeric = Number(value)
  if (!Number.isFinite(numeric)) return TELEPORT_MIN
  const clamped = Math.max(TELEPORT_MIN, Math.min(1, numeric))
  return Number(clamped.toFixed(6))
}

const formatDateTime = (value) => (value ? new Date(value).toLocaleString() : '—')

export default function AnalysisWorkspace({
  mode = 'new',
  jobId = '',
  entry = null,
  onBack,
  onEntryChanged,
}) {
  const isHistory = mode === 'history' && !!jobId

  const [analysisName, setAnalysisName] = useState('')
  const [bacteriaFiles, setBacteriaFiles] = useState([])
  const [fungiFiles, setFungiFiles] = useState([])
  const [bacteriaMetadataFiles, setBacteriaMetadataFiles] = useState([])
  const [fungiMetadataFiles, setFungiMetadataFiles] = useState([])
  const [bacteriaMetadataConditions, setBacteriaMetadataConditions] = useState([])
  const [fungiMetadataConditions, setFungiMetadataConditions] = useState([])
  const [bacteriaMetadataError, setBacteriaMetadataError] = useState('')
  const [fungiMetadataError, setFungiMetadataError] = useState('')
  const [params, setParams] = useState(DEFAULT_PARAMS)
  const [analysisJobId, setAnalysisJobId] = useState(isHistory ? jobId : '')
  const [isBootstrapping, setIsBootstrapping] = useState(isHistory)
  const [saveState, setSaveState] = useState('')
  const [isCheckingSparcc, setIsCheckingSparcc] = useState(false)
  const [sparccStatus, setSparccStatus] = useState(null)
  const [selectedTaxa, setSelectedTaxa] = useState('')
  const [taxaSearch, setTaxaSearch] = useState('')
  const [runStatusByGraph, setRunStatusByGraph] = useState({})
  const [isClusteringRunning, setIsClusteringRunning] = useState(false)
  const [clusterParamsByGraph, setClusterParamsByGraph] = useState({})

  const { runAnalysis, resetRun, status, progress, logs, error } = useAnalysis()
  const updateHistoryEntry = useAnalysisStore((state) => state.updateHistoryEntry)
  const {
    results,
    conditionGraphs,
    reclusterByGraph,
    fetchResults,
    fetchConditionGraph,
    reclusterConditionGraph,
    hydrateReclusterByGraph,
    clearResults,
    isLoading: isResultsLoading,
  } = useResultsStore()

  // The entry object is recreated by the archive refresh, so it is read through a
  // ref: the bootstrap effect must depend on `jobId` alone, never on identity.
  const entryRef = useRef(entry)
  entryRef.current = entry
  const hydratedRef = useRef(false)

  /** Reset every piece of form state. */
  const resetWorkspace = useCallback(() => {
    setAnalysisName('')
    setBacteriaFiles([])
    setFungiFiles([])
    setBacteriaMetadataFiles([])
    setFungiMetadataFiles([])
    setBacteriaMetadataConditions([])
    setFungiMetadataConditions([])
    setBacteriaMetadataError('')
    setFungiMetadataError('')
    setParams(DEFAULT_PARAMS)
    setSelectedTaxa('')
    setTaxaSearch('')
    setRunStatusByGraph({})
    setClusterParamsByGraph({})
    setSaveState('')
  }, [])

  /**
   * Load the workspace.
   *
   * `mode="new"` always starts from a clean form — this is what the "New Analysis"
   * entry point is for. `mode="history"` restores the stored files, parameters and
   * results of the selected analysis.
   */
  useEffect(() => {
    let cancelled = false
    hydratedRef.current = false
    resetWorkspace()
    resetRun()

    if (!isHistory) {
      setAnalysisJobId('')
      clearResults()
      hydratedRef.current = true
      setIsBootstrapping(false)
      return undefined
    }

    const bootstrap = async () => {
      setIsBootstrapping(true)
      setAnalysisJobId(jobId)
      clearResults()

      const snapshot = entryRef.current?.payloadSnapshot || {}
      const storedFiles = snapshot.files || {}

      const [bacteria, bacteriaMeta, fungi, fungiMeta] = await Promise.all([
        describePaths(storedFiles.bacteria || []),
        describePaths(storedFiles.bacteria_metadata || []),
        describePaths(storedFiles.fungi || []),
        describePaths(storedFiles.fungi_metadata || []),
      ])

      if (cancelled) return

      setAnalysisName(entryRef.current?.analysisName || '')
      setBacteriaFiles(bacteria.accepted)
      setBacteriaMetadataFiles(bacteriaMeta.accepted)
      setFungiFiles(fungi.accepted)
      setFungiMetadataFiles(fungiMeta.accepted)
      setParams({ ...DEFAULT_PARAMS, ...(snapshot.params || {}) })
      hydratedRef.current = true
      setIsBootstrapping(false)

      await fetchResults(jobId)
    }

    bootstrap()
    return () => {
      cancelled = true
    }
    // Intentionally keyed on the job only: `resetRun`, `clearResults` and
    // `fetchResults` are stable store actions.
  }, [jobId, isHistory])

  /**
   * Single source of truth for the metadata conditions.
   *
   * Everything (browse, drag & drop, restored history, removal) ends up changing
   * the file lists, so deriving the conditions from them here keeps all paths
   * identical and removes the duplicated per-handler logic that used to make the
   * condition selects behave differently depending on how a file was added.
   */
  const metadataPathsKey = useMemo(
    () => [
      ...bacteriaMetadataFiles.map((file) => file.path),
      ...fungiMetadataFiles.map((file) => file.path),
    ].join('|'),
    [bacteriaMetadataFiles, fungiMetadataFiles]
  )

  useEffect(() => {
    let cancelled = false

    const parseAll = async (fileList) => {
      const conditions = new Set()
      let firstError = ''
      for (const file of fileList) {
        try {
          const parsed = await parseMetadataConditions(file.path)
          parsed.forEach((condition) => conditions.add(condition))
        } catch (err) {
          if (!firstError) firstError = `${file.name}: ${err.message}`
        }
      }
      return { conditions: [...conditions], error: firstError }
    }

    const sync = async () => {
      const bacteria = await parseAll(bacteriaMetadataFiles)
      const fungi = await parseAll(fungiMetadataFiles)
      if (cancelled) return

      setBacteriaMetadataConditions(bacteria.conditions)
      setBacteriaMetadataError(bacteria.error)
      setFungiMetadataConditions(fungi.conditions)
      setFungiMetadataError(fungi.error)

      const allConditions = [...new Set([...bacteria.conditions, ...fungi.conditions])]

      setParams((prev) => {
        // Never invalidate a restored selection: only prune once the metadata has
        // actually been read, or when there is no metadata file left at all.
        const canPrune = allConditions.length > 0
          || (metadataPathsKey === '' && hydratedRef.current)
        if (!canPrune) return prev

        const nextA = allConditions.includes(prev.conditionA) ? prev.conditionA : ''
        const nextB = allConditions.includes(prev.conditionB) ? prev.conditionB : ''
        if (nextA === (prev.conditionA ?? '') && nextB === (prev.conditionB ?? '')) return prev
        return { ...prev, conditionA: nextA, conditionB: nextB }
      })
    }

    sync()
    return () => {
      cancelled = true
    }
    // Re-derive only when the set of metadata files actually changes.
  }, [metadataPathsKey])

  useEffect(() => {
    const check = async () => {
      if (params?.correlationMethod !== 'sparcc') {
        setSparccStatus(null)
        return
      }
      setIsCheckingSparcc(true)
      try {
        const result = await window.electronAPI.getSparccStatus()
        setSparccStatus(result)
      } finally {
        setIsCheckingSparcc(false)
      }
    }
    check()
  }, [params?.correlationMethod])

  const availableDatasets = {
    bacteria: bacteriaFiles.length > 0,
    fungi: fungiFiles.length > 0,
  }

  useEffect(() => {
    setParams((prev) => {
      const { datasetType } = prev
      if (datasetType === 'bacteria' && !availableDatasets.bacteria) return { ...prev, datasetType: '' }
      if (datasetType === 'fungi' && !availableDatasets.fungi) return { ...prev, datasetType: '' }
      if (datasetType === 'both' && !(availableDatasets.bacteria && availableDatasets.fungi)) {
        return { ...prev, datasetType: '' }
      }
      return prev
    })
  }, [availableDatasets.bacteria, availableDatasets.fungi])

  const allMetadataConditions = useMemo(
    () => [...new Set([...bacteriaMetadataConditions, ...fungiMetadataConditions])],
    [bacteriaMetadataConditions, fungiMetadataConditions]
  )

  const metadataRequirementMet = (
    (availableDatasets.bacteria && bacteriaMetadataFiles.length > 0) ||
    (availableDatasets.fungi && fungiMetadataFiles.length > 0)
  )

  const handleAddDatasetFiles = useCallback((kind, descriptors) => {
    const setter = kind === 'bacteria' ? setBacteriaFiles : setFungiFiles
    setter((prev) => dedupeDescriptors(prev, descriptors).merged)
  }, [])

  // The conditions themselves are re-derived by the effect above from the file
  // lists, so adding a file only has to register it and clear a stale error.
  const handleAddMetadataFiles = useCallback((kind, descriptors) => {
    if (kind === 'bacteria') {
      setBacteriaMetadataFiles((prev) => dedupeDescriptors(prev, descriptors).merged)
      setBacteriaMetadataError('')
    } else {
      setFungiMetadataFiles((prev) => dedupeDescriptors(prev, descriptors).merged)
      setFungiMetadataError('')
    }
  }, [])

  const handleRemoveMetadataFile = useCallback((kind, index) => {
    const setter = kind === 'bacteria' ? setBacteriaMetadataFiles : setFungiMetadataFiles
    setter((prev) => prev.filter((_, idx) => idx !== index))
  }, [])

  const buildPayload = useCallback(({ overrideJobId } = {}) => ({
    scriptName: 'analysis.py',
    analysisName: analysisName.trim(),
    ...(overrideJobId ? { jobId: overrideJobId } : {}),
    files: {
      bacteria: toPathList(bacteriaFiles),
      bacteria_metadata: toPathList(bacteriaMetadataFiles),
      fungi: toPathList(fungiFiles),
      fungi_metadata: toPathList(fungiMetadataFiles),
    },
    params,
  }), [analysisName, bacteriaFiles, bacteriaMetadataFiles, fungiFiles, fungiMetadataFiles, params])

  const handleSaveDetails = useCallback(async () => {
    if (!isHistory) return
    setSaveState('Saving…')
    const payload = buildPayload()
    const response = await updateHistoryEntry(jobId, {
      analysisName: analysisName.trim(),
      updatedAt: Date.now(),
      payloadSnapshot: { files: payload.files, params },
    })
    setSaveState(response?.success ? 'Saved ✓' : 'Save failed')
    if (response?.success) onEntryChanged?.()
    setTimeout(() => setSaveState(''), 2500)
  }, [analysisName, buildPayload, isHistory, jobId, onEntryChanged, params, updateHistoryEntry])

  const handleRun = async () => {
    setSelectedTaxa('')
    setTaxaSearch('')
    setRunStatusByGraph({})
    setClusterParamsByGraph({})

    if (!isHistory) {
      // A brand new analysis must not keep showing the previous results while it runs.
      setAnalysisJobId('')
      clearResults()
    }

    // Re-running a stored analysis keeps its job id so the archive entry is
    // updated in place rather than duplicated.
    const reuseJobId = isHistory ? jobId : ''
    const result = await runAnalysis(buildPayload({ overrideJobId: reuseJobId }))

    if (result?.success) {
      setAnalysisJobId(result.jobId)
      const data = await fetchResults(result.jobId)
      // Keep the archive self-describing: the summary is what the history list
      // and the comparison view show without loading the (large) full results.
      if (data?.summary) {
        await updateHistoryEntry(result.jobId, { summary: data.summary })
      }
      if (isHistory) onEntryChanged?.()
    }
  }

  const hasInputData = availableDatasets.bacteria || availableDatasets.fungi
  const hasRunCompleted = !!analysisJobId && status !== 'running'
  const hasInlineResults = Boolean(analysisJobId && results)

  useEffect(() => {
    const refs = results?.graph_pipeline?.graph_refs || []
    refs.slice(0, 2).forEach((ref) => {
      fetchConditionGraph(ref.cacheKey)
    })
  }, [results, fetchConditionGraph])

  useEffect(() => {
    const savedClusters = results?.graph_pipeline?.last_clustering || {}
    const entries = Object.entries(savedClusters)
    if (entries.length === 0) return

    hydrateReclusterByGraph(savedClusters)
    setRunStatusByGraph((prev) => {
      const next = { ...prev }
      entries.forEach(([cacheKey, payload]) => {
        next[cacheKey] = payload?.success ? 'done' : 'error'
      })
      return next
    })

    const firstSeed = entries.find(([, payload]) => typeof payload?.seedNode === 'string' && payload.seedNode.trim())
    if (firstSeed?.[1]?.seedNode) {
      setSelectedTaxa(firstSeed[1].seedNode)
      setTaxaSearch(firstSeed[1].seedNode)
    }
  }, [results, hydrateReclusterByGraph])

  // Memoised on purpose: a fresh array on every render would make the effects
  // below re-run (and re-set their state) on every render, freezing the page.
  const graphRefs = useMemo(
    () => (results?.graph_pipeline?.graph_refs || []).slice(0, 2),
    [results]
  )
  const taxaSectionReady = hasRunCompleted && graphRefs.length > 0
  const isAnyRecomputeRunning = Object.values(runStatusByGraph).some((state) => state === 'running')
  const isPrimaryRunRunning = status === 'running' || isClusteringRunning || isResultsLoading
  const isInteractionLocked = isPrimaryRunRunning || isAnyRecomputeRunning
  const canUseRecomputeButtons = !isPrimaryRunRunning

  const canRun = hasInputData
    && metadataRequirementMet
    && allMetadataConditions.length >= 2
    && !!params.datasetType
    && !!params.correlationMethod
    && !!params.groupingMode
    && Number.isFinite(Number(params.pValueThreshold))
    && Number(params.pValueThreshold) >= 0
    && Number(params.pValueThreshold) <= 1
    && !!params.conditionA
    && !!params.conditionB
    && params.conditionA !== params.conditionB
    && status !== 'running'
    && !isInteractionLocked

  useEffect(() => {
    if (graphRefs.length === 0) return
    setClusterParamsByGraph((prev) => {
      // Only allocate a new state object when something really changes,
      // otherwise return `prev` so React bails out of the re-render.
      let next = null
      graphRefs.forEach((ref) => {
        const payload = conditionGraphs[ref.cacheKey]
        const stats = payload?.graph?.stats
        const nodeCount = stats?.node_count ?? payload?.graph?.nodes?.length ?? 1
        const edgeCount = stats?.edge_count ?? payload?.graph?.edges?.length ?? 1
        const defaultConductance = 0.2
        const prevParams = (next || prev)[ref.cacheKey]
        const hasGraphShapeChanged =
          !prevParams
          || prevParams._defaultNodeCount !== nodeCount
          || prevParams._defaultEdgeCount !== edgeCount

        if (!prevParams || (prevParams._defaultsAuto && hasGraphShapeChanged)) {
          if (!next) next = { ...prev }
          next[ref.cacheKey] = {
            conductance: defaultConductance,
            b: computeDefaultB(nodeCount, edgeCount),
            teleport: computeDefaultTeleportAlpha(defaultConductance, edgeCount),
            _defaultsAuto: true,
            _defaultNodeCount: nodeCount,
            _defaultEdgeCount: edgeCount,
          }
        }
      })
      return next || prev
    })
  }, [graphRefs, conditionGraphs])

  const selectorNodes = useMemo(() => {
    const seen = new Set()
    const merged = []
    graphRefs.forEach((ref) => {
      const payload = conditionGraphs[ref.cacheKey]
      const nodes = payload?.graph?.nodes || []
      nodes.forEach((node) => {
        if (!node?.id || seen.has(node.id)) return
        seen.add(node.id)
        merged.push(node)
      })
    })
    return merged
  }, [graphRefs, conditionGraphs])

  const runSubgraphsForConditions = async () => {
    const taxaNode = selectedTaxa.trim()
    if (!taxaNode || isClusteringRunning || graphRefs.length === 0) return

    setIsClusteringRunning(true)

    const loadingStatus = {}
    graphRefs.forEach((ref) => {
      loadingStatus[ref.cacheKey] = 'running'
    })
    setRunStatusByGraph(loadingStatus)

    const nextStatus = { ...loadingStatus }
    for (const ref of graphRefs) {
      const payload = conditionGraphs[ref.cacheKey]
      const graphNodes = payload?.graph?.nodes || []
      const hasTaxa = graphNodes.some((node) => node.id === taxaNode)
      if (!hasTaxa) {
        nextStatus[ref.cacheKey] = 'missing-taxa'
        setRunStatusByGraph({ ...nextStatus })
        continue
      }

      const response = await reclusterConditionGraph({
        cacheKey: ref.cacheKey,
        seedNode: taxaNode,
        options: {
          weighted: true,
          targetConductance: Number.isFinite(Number(clusterParamsByGraph[ref.cacheKey]?.conductance))
            ? Number(clusterParamsByGraph[ref.cacheKey]?.conductance)
            : 0.2,
          b: Number.isFinite(Number(clusterParamsByGraph[ref.cacheKey]?.b))
            ? Number(clusterParamsByGraph[ref.cacheKey]?.b)
            : 1,
          teleportAlpha: Number.isFinite(Number(clusterParamsByGraph[ref.cacheKey]?.teleport))
            ? normalizeTeleport(clusterParamsByGraph[ref.cacheKey]?.teleport)
            : normalizeTeleport(0),
        },
      })

      nextStatus[ref.cacheKey] = response?.success ? 'done' : 'error'
      setRunStatusByGraph({ ...nextStatus })
    }

    setIsClusteringRunning(false)
    if (isHistory) onEntryChanged?.()
  }

  const updateClusterParams = (cacheKey, nextValues) => {
    const payload = conditionGraphs[cacheKey]
    const stats = payload?.graph?.stats
    const nodeCount = stats?.node_count ?? payload?.graph?.nodes?.length ?? 1
    const edgeCount = stats?.edge_count ?? payload?.graph?.edges?.length ?? 1
    const bMax = Math.max(1, Math.ceil(Math.log2(Math.max(1, edgeCount))))
    const sanitized = { ...nextValues }

    if (Object.prototype.hasOwnProperty.call(sanitized, 'conductance')) {
      if (sanitized.conductance === '') {
        sanitized.conductance = ''
      } else {
        const n = Number(sanitized.conductance)
        sanitized.conductance = Number.isFinite(n) ? Math.max(0, Math.min(1, n)) : 0.2
      }
    }
    if (Object.prototype.hasOwnProperty.call(sanitized, 'b')) {
      if (sanitized.b === '') {
        sanitized.b = ''
      } else {
        const n = Number(sanitized.b)
        sanitized.b = Number.isFinite(n) ? Math.max(1, Math.min(bMax, Math.floor(n))) : 1
      }
    }
    if (Object.prototype.hasOwnProperty.call(sanitized, 'teleport')) {
      if (sanitized.teleport === '') {
        sanitized.teleport = ''
      } else if (typeof sanitized.teleport === 'string') {
        // Preserve in-progress typing in the input; clamp only when blur sends a number.
      } else {
        const n = Number(sanitized.teleport)
        sanitized.teleport = Number.isFinite(n) ? Number(Math.max(0, Math.min(1, n)).toFixed(6)) : 0
      }
    }

    setClusterParamsByGraph((prev) => ({
      ...prev,
      [cacheKey]: {
        ...(prev[cacheKey] || {
          conductance: 0.2,
          b: computeDefaultB(nodeCount, edgeCount),
          teleport: computeDefaultTeleportAlpha(0.2, edgeCount),
        }),
        ...sanitized,
        _defaultsAuto: false,
        _defaultNodeCount: nodeCount,
        _defaultEdgeCount: edgeCount,
      },
    }))
  }

  const rerunSingleCluster = async (cacheKey) => {
    const taxaNode = selectedTaxa.trim()
    if (!taxaNode || !cacheKey || isPrimaryRunRunning) return

    const ref = graphRefs.find((item) => item.cacheKey === cacheKey)
    if (!ref) return

    const payload = conditionGraphs[cacheKey]
    const graphNodes = payload?.graph?.nodes || []
    const hasTaxa = graphNodes.some((node) => node.id === taxaNode)
    if (!hasTaxa) {
      setRunStatusByGraph((prev) => ({ ...prev, [cacheKey]: 'missing-taxa' }))
      return
    }

    setRunStatusByGraph((prev) => ({ ...prev, [cacheKey]: 'running' }))

    const graphParams = clusterParamsByGraph[cacheKey] || { conductance: 0.2, b: 1, teleport: 0 }
    const response = await reclusterConditionGraph({
      cacheKey,
      seedNode: taxaNode,
      options: {
        weighted: true,
        targetConductance: Number.isFinite(Number(graphParams.conductance))
          ? Number(graphParams.conductance)
          : 0.2,
        b: Number.isFinite(Number(graphParams.b)) ? Number(graphParams.b) : 1,
        teleportAlpha: Number.isFinite(Number(graphParams.teleport))
          ? normalizeTeleport(graphParams.teleport)
          : normalizeTeleport(0),
      },
    })

    setRunStatusByGraph((prev) => ({
      ...prev,
      [cacheKey]: response?.success ? 'done' : 'error',
    }))
    if (isHistory) onEntryChanged?.()
  }

  const missingFiles = isHistory
    ? (entryRef.current?.payloadSnapshot?.files?.bacteria?.length || 0)
      + (entryRef.current?.payloadSnapshot?.files?.fungi?.length || 0)
      + (entryRef.current?.payloadSnapshot?.files?.bacteria_metadata?.length || 0)
      + (entryRef.current?.payloadSnapshot?.files?.fungi_metadata?.length || 0)
      - (bacteriaFiles.length + fungiFiles.length + bacteriaMetadataFiles.length + fungiMetadataFiles.length)
    : 0

  return (
    <div className={styles.page}>
      <header className={styles.workspaceHeader}>
        <div className={styles.headerMain}>
          {typeof onBack === 'function' && (
            <div className={styles.headerTop}>
              <button
                type='button'
                className={styles.backButton}
                data-testid='back-to-history'
                onClick={onBack}
              >
                ← Back to history
              </button>
            </div>
          )}
          <h1 className={styles.title} data-testid='workspace-title'>
            {analysisName?.trim() || (isHistory ? 'Untitled analysis' : 'Create new analysis')}
          </h1>
          {isHistory && (
            <p className={styles.headerMeta}>
              Created {formatDateTime(entry?.createdAt)} · Last update {formatDateTime(entry?.updatedAt)}
              {entry?.runCount ? ` · ${entry.runCount} run${entry.runCount > 1 ? 's' : ''}` : ''}
            </p>
          )}
        </div>
        {isHistory && (
          <div className={styles.headerActions}>
            <button type='button' className={styles.secondaryButton} onClick={handleSaveDetails} disabled={isInteractionLocked}>
              {saveState || 'Save changes'}
            </button>
          </div>
        )}
      </header>

      {isBootstrapping && (
        <p className={styles.resultsHint}>Loading the stored analysis…</p>
      )}

      {isHistory && missingFiles > 0 && (
        <p className={styles.warning}>
          {missingFiles} file{missingFiles > 1 ? 's' : ''} from this analysis could not be found on
          disk any more (they were moved or deleted). Re-add them before running the analysis again.
        </p>
      )}

      <div className={styles.workspaceBody}>
        <section className={`${styles.section} ${styles.uploadSection}`}>
          <h2>1 · Upload Files</h2>
          <div className={styles.uploadGrid}>
            <div className={styles.uploadColumn}>
              <h3>Bacteria abundance file</h3>
              <FileUploadZone
                testId='upload-zone-bacteria'
                onFilesAdded={(files) => handleAddDatasetFiles('bacteria', files)}
                files={bacteriaFiles}
                onRemove={(idx) => setBacteriaFiles((prev) => prev.filter((_, i) => i !== idx))}
                disabled={isInteractionLocked}
              />
            </div>

            <div className={styles.uploadColumn}>
              <h3>Fungi abundance file</h3>
              <FileUploadZone
                testId='upload-zone-fungi'
                onFilesAdded={(files) => handleAddDatasetFiles('fungi', files)}
                files={fungiFiles}
                onRemove={(idx) => setFungiFiles((prev) => prev.filter((_, i) => i !== idx))}
                disabled={isInteractionLocked}
              />
            </div>

            <div className={styles.uploadColumn}>
              <h3>Bacteria metadata file</h3>
              <FileUploadZone
                testId='upload-zone-bacteria-metadata'
                onFilesAdded={(files) => handleAddMetadataFiles('bacteria', files)}
                files={bacteriaMetadataFiles}
                onRemove={(idx) => handleRemoveMetadataFile('bacteria', idx)}
                disabled={isInteractionLocked}
              />
              {bacteriaMetadataError && <p className={styles.error}>{bacteriaMetadataError}</p>}
            </div>

            <div className={styles.uploadColumn}>
              <h3>Fungi metadata file</h3>
              <FileUploadZone
                testId='upload-zone-fungi-metadata'
                onFilesAdded={(files) => handleAddMetadataFiles('fungi', files)}
                files={fungiMetadataFiles}
                onRemove={(idx) => handleRemoveMetadataFile('fungi', idx)}
                disabled={isInteractionLocked}
              />
              {fungiMetadataError && <p className={styles.error}>{fungiMetadataError}</p>}
            </div>
          </div>
        </section>

        <section className={`${styles.section} ${styles.parametersSection}`}>
          <h2>2 · Configure Parameters</h2>
          <div className={styles.inputRow}>
            <label htmlFor='analysis-name'>Analysis name</label>
            <input
              id='analysis-name'
              value={analysisName}
              onChange={(e) => setAnalysisName(e.target.value)}
              onBlur={() => { if (isHistory) handleSaveDetails() }}
              placeholder='Optional name for history'
              disabled={isInteractionLocked}
            />
          </div>
          <ParametersForm
            values={params}
            onChange={setParams}
            availableDatasets={availableDatasets}
            conditionOptions={allMetadataConditions}
            disabled={isInteractionLocked}
          />
          <div className={styles.runInline}>
            <h3>3 · Run Analysis</h3>
            {status === 'running' && <AnalysisProgress progress={progress} logs={logs} />}
            {error && <p className={styles.error}>{error}</p>}
            {params.correlationMethod === 'sparcc' && isCheckingSparcc && (
              <p className={styles.resultsHint}>Checking SparCC runtime...</p>
            )}
            {params.correlationMethod === 'sparcc' && !isCheckingSparcc && !sparccStatus?.available && (
              <p className={styles.warning}>
                FastSpar is not available for {sparccStatus?.platform || 'this platform'}
                {sparccStatus?.arch ? `/${sparccStatus.arch}` : ''}. The analysis will run with{' '}
                <strong>Spearman correlation</strong> instead of SparCC, and the substitution is
                recorded in the log. Configure FastSpar in Settings to use SparCC.
              </p>
            )}
            <button className={styles.runButton} onClick={handleRun} disabled={!canRun}>
              {status === 'running' ? 'Running…' : isHistory ? 'Re-run analysis' : 'Start Analysis'}
            </button>
            {!hasInputData && (
              <p className={styles.resultsHint}>Add at least one abundance file to enable the analysis.</p>
            )}
          </div>
        </section>

        <section className={`${styles.section} ${styles.clusterControlsSection}`}>
          <h2>4 · Taxa Selection</h2>
          {taxaSectionReady ? (
            <>
              <SeedSelector
                nodes={selectorNodes}
                searchValue={taxaSearch}
                selectedTaxa={selectedTaxa}
                onSearchChange={setTaxaSearch}
                onTaxaSelect={setSelectedTaxa}
                inputId='cluster-taxa-selector'
                disabled={isInteractionLocked}
              />
              <button
                className={styles.runButton}
                onClick={runSubgraphsForConditions}
                disabled={!selectedTaxa.trim() || isInteractionLocked || graphRefs.length === 0}
              >
                {isClusteringRunning ? 'Clustering in progress...' : 'Start clustering algorithm'}
              </button>
            </>
          ) : (
            <p className={styles.resultsHint}>Taxa are available only after the analysis completes.</p>
          )}
        </section>

        {analysisJobId && isResultsLoading && (
          <section className={`${styles.section} ${styles.resultsSection}`}>
            <h2>5 · Cluster Results</h2>
            <p className={styles.resultsHint}>Loading results…</p>
          </section>
        )}

        {hasInlineResults && !isResultsLoading && (
          <section className={`${styles.section} ${styles.resultsSection}`}>
            <h2>5 · Cluster Results</h2>
            <ChartPanel
              graphRefs={graphRefs}
              conditionGraphs={conditionGraphs}
              reclusterByGraph={reclusterByGraph}
              selectedTaxa={selectedTaxa}
              runStatusByGraph={runStatusByGraph}
              clusterParamsByGraph={clusterParamsByGraph}
              onClusterParamsChange={updateClusterParams}
              onRerunCluster={rerunSingleCluster}
              isClusteringRunning={isClusteringRunning}
              interactionLocked={isInteractionLocked}
              canUseRecomputeButtons={canUseRecomputeButtons}
            />
          </section>
        )}
      </div>
    </div>
  )
}
