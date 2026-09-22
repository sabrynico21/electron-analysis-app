/**
 * ComparisonView
 *
 * Two stored analyses side by side, each with its own colour tint (blue for A,
 * amber for B) and both stretched to the same height.
 *
 * Shape of a column:
 *   - header: side badge, name, timestamps, and the run status next to the
 *     "Open analysis" button (top right);
 *   - input files, parameters and result summary;
 *   - **two separate graph areas of fixed size**: one for the condition graphs,
 *     one for the clusters (including the intersection / exclusive / complement
 *     views of the two clusters). Each area has a single `select` switcher that
 *     always shows the full selected label, a node limit, a button that opens the
 *     graph enlarged in an overlay, and a collapsible line of details.
 *
 * Fixed areas keep an analysis with fewer conditions exactly as tall as one with
 * more: only the options inside the select change.
 *
 * There is deliberately no "how many differences" counter and no correlation-value
 * table: the two states are meant to be read against each other visually.
 *
 * The two analyses are loaded independently here (own state per job) instead of
 * going through the shared results store, which by design holds a single job.
 */
import { useEffect, useMemo, useState } from 'react'
import SubgraphPlot from '../Results/SubgraphPlot'
import { buildDerivedSubgraphs } from '../../utils/subgraphComparison'
import styles from './ComparisonView.module.css'

const PARAM_LABELS = {
  datasetType: 'Data type',
  correlationMethod: 'Correlation method',
  groupingMode: 'Taxa grouping',
  pValueThreshold: 'P-value threshold',
  conditionA: 'Condition A',
  conditionB: 'Condition B',
}

const FILE_SLOTS = [
  ['bacteria', 'Bacteria abundance'],
  ['bacteria_metadata', 'Bacteria metadata'],
  ['fungi', 'Fungi abundance'],
  ['fungi_metadata', 'Fungi metadata'],
]

/** Clustering metadata entries that are huge arrays rather than scalars. */
const SKIPPED_METADATA_KEYS = new Set(['ranked_nodes', 'rankedNodes'])

const formatValue = (value) => {
  if (value === null || value === undefined || value === '') return '—'
  if (typeof value === 'number') return Number.isInteger(value) ? String(value) : value.toFixed(4)
  if (typeof value === 'boolean') return value ? 'yes' : 'no'
  if (Array.isArray(value)) return value.length === 0 ? '—' : value.join(', ')
  if (typeof value === 'object') return JSON.stringify(value)
  return String(value)
}

const formatDateTime = (value) => (value ? new Date(value).toLocaleString() : '—')

const fileName = (filePath) => String(filePath || '').split(/[/\\]/).pop() || ''

/** Scalar clustering metadata, flattened into readable `key value` pairs. */
function metadataEntries(metadata) {
  return Object.entries(metadata || {}).filter(([key, value]) => (
    !SKIPPED_METADATA_KEYS.has(key) && (value === null || typeof value !== 'object')
  ))
}

/** Load one job's results together with its cached condition graphs. */
function useJobBundle(jobId) {
  const [state, setState] = useState({ loading: true, results: null, graphs: {}, error: null })

  useEffect(() => {
    let cancelled = false

    const load = async () => {
      setState({ loading: true, results: null, graphs: {}, error: null })
      if (!jobId) {
        setState({ loading: false, results: null, graphs: {}, error: 'Missing analysis id' })
        return
      }

      const results = await window.electronAPI.getResults(jobId)
      if (cancelled) return
      if (!results) {
        setState({ loading: false, results: null, graphs: {}, error: 'No stored results for this analysis' })
        return
      }

      const refs = results?.graph_pipeline?.graph_refs || []
      const graphs = {}
      await Promise.all(
        refs.map(async (ref) => {
          const response = await window.electronAPI.getConditionGraph(ref.cacheKey)
          if (response?.success) graphs[ref.cacheKey] = response.graph
        })
      )
      if (cancelled) return
      setState({ loading: false, results, graphs, error: null })
    }

    load()
    return () => {
      cancelled = true
    }
  }, [jobId])

  return state
}

/** One option of a graph switcher. */
function makeView({ id, label, title, nodes = [], edges = [], seedNode = '', details = [] }) {
  return {
    id,
    label,
    title,
    nodes,
    edges,
    seedNode,
    details: details.filter(Boolean),
    subtitle: `${nodes.length} nodes · ${edges.length} edges`,
  }
}

/** One condition graph per condition. */
function buildConditionGraphViews(refs, graphs) {
  return refs.slice(0, 2).flatMap((ref) => {
    const graph = graphs[ref.cacheKey]?.graph || null
    const nodes = graph?.nodes || []
    const edges = graph?.edges || []
    if (nodes.length === 0) return []
    return [makeView({
      id: `${ref.cacheKey}::graph`,
      label: ref.condition,
      title: `Condition graph · ${ref.condition}`,
      nodes,
      edges,
      details: [
        graph?.stats?.density !== undefined ? `density ${formatValue(graph.stats.density)}` : '',
      ],
    })]
  })
}

/** One cluster per condition, plus the combined views when both exist. */
function buildClusterViews(refs, lastClustering, derived) {
  const views = []

  refs.slice(0, 2).forEach((ref) => {
    const cluster = lastClustering[ref.cacheKey] || null
    if (!cluster) return
    views.push(makeView({
      id: `${ref.cacheKey}::cluster`,
      label: `${ref.condition} · cluster`,
      title: `Cluster · ${ref.condition}`,
      nodes: cluster.subgraph?.nodes || [],
      edges: cluster.subgraph?.edges || [],
      seedNode: cluster.seedNode || '',
      details: [
        `seed ${cluster.seedNode || '—'}`,
        cluster.algorithm || '',
        ...metadataEntries(cluster.metadata).map(([key, value]) => `${key} ${formatValue(value)}`),
      ],
    }))
  })

  const firstRef = refs[0] || null
  const secondRef = refs[1] || null
  const firstCluster = firstRef ? lastClustering[firstRef.cacheKey] || null : null
  const secondCluster = secondRef ? lastClustering[secondRef.cacheKey] || null : null

  if (firstCluster && secondCluster) {
    const first = firstRef?.condition || 'condition 1'
    const second = secondRef?.condition || 'condition 2'
    views.push(
      makeView({
        id: 'derived::intersection',
        label: `Intersection · ${first} ∩ ${second}`,
        title: 'Intersection of the two clusters',
        nodes: derived.intersection.nodes,
        edges: derived.intersection.edges,
        seedNode: firstCluster?.seedNode || '',
      }),
      makeView({
        id: 'derived::exclusiveA',
        label: `Only in ${first}`,
        title: `Only in ${first}`,
        nodes: derived.exclusiveA.nodes,
        edges: derived.exclusiveA.edges,
        seedNode: firstCluster?.seedNode || '',
      }),
      makeView({
        id: 'derived::exclusiveB',
        label: `Only in ${second}`,
        title: `Only in ${second}`,
        nodes: derived.exclusiveB.nodes,
        edges: derived.exclusiveB.edges,
        seedNode: secondCluster?.seedNode || '',
      }),
      makeView({
        id: 'derived::complement',
        label: 'Complement (union without intersection)',
        title: 'Complement (union without intersection)',
        nodes: derived.complement.nodes,
        edges: derived.complement.edges,
        seedNode: firstCluster?.seedNode || '',
      })
    )
  }

  return views
}

function KeyValueList({ object, labels = {}, emptyLabel = 'No data.' }) {
  const entries = Object.entries(object || {})
  if (entries.length === 0) return <p className={styles.muted}>{emptyLabel}</p>
  return (
    <dl className={styles.kv}>
      {entries.map(([key, value]) => (
        <div key={key} className={styles.kvRow}>
          <dt>{labels[key] || key}</dt>
          <dd>{formatValue(value)}</dd>
        </div>
      ))}
    </dl>
  )
}

/**
 * A fixed-size graph area with a single `select` switcher, an enlarging overlay
 * and a collapsible details line.
 *
 * The viewport height does not depend on the number of options, so two analyses
 * with different numbers of conditions still line up.
 */
function GraphPanel({ side, kind, title, views, emptyMessage }) {
  const [selectedId, setSelectedId] = useState('')
  const [maxNodes, setMaxNodes] = useState(40)
  const [expanded, setExpanded] = useState(false)

  useEffect(() => {
    if (views.length === 0) {
      if (selectedId !== '') setSelectedId('')
      return
    }
    if (!views.some((view) => view.id === selectedId)) {
      setSelectedId(views[0].id)
    }
  }, [views, selectedId])

  useEffect(() => {
    if (!expanded) return undefined
    const onKey = (event) => {
      if (event.key === 'Escape') setExpanded(false)
    }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [expanded])

  const activeView = views.find((view) => view.id === selectedId) || views[0] || null
  const hasNodes = Boolean(activeView && activeView.nodes.length > 0)

  return (
    <div className={styles.graphPanel} data-testid={`graph-panel-${side}-${kind}`}>
      <h4 className={styles.graphPanelTitle}>{title}</h4>

      <div className={styles.graphControls}>
        <select
          className={styles.graphSelect}
          value={activeView?.id || ''}
          disabled={views.length === 0}
          data-testid={`graph-select-${side}-${kind}`}
          aria-label={`${title}: graph to display`}
          onChange={(event) => setSelectedId(event.target.value)}
        >
          {views.length === 0 && <option value=''>{emptyMessage}</option>}
          {views.map((view) => (
            <option key={view.id} value={view.id}>{view.label}</option>
          ))}
        </select>

        <label className={styles.nodeLimit}>
          Nodes
          <input
            type='number'
            min='1'
            max='200'
            value={maxNodes}
            data-testid={`graph-nodes-${side}-${kind}`}
            onChange={(event) => {
              const parsed = Number(event.target.value)
              setMaxNodes(Number.isFinite(parsed) ? Math.max(1, Math.min(200, parsed)) : 40)
            }}
          />
        </label>

        <button
          type='button'
          className={styles.expandButton}
          data-testid={`graph-expand-${side}-${kind}`}
          onClick={() => setExpanded(true)}
          disabled={!hasNodes}
          title='Show the graph enlarged'
          aria-label='Show the graph enlarged'
        >
          ⤢
        </button>
      </div>

      <div className={styles.plotFill} data-testid={`graph-viewport-${side}-${kind}`}>
        {hasNodes ? (
          <SubgraphPlot
            nodes={activeView.nodes}
            edges={activeView.edges}
            maxNodes={maxNodes}
            seedNode={activeView.seedNode}
          />
        ) : (
          <div className={styles.emptyViewport}>{emptyMessage}</div>
        )}
      </div>

      {activeView && (
        <details className={styles.viewDetails} data-testid={`graph-details-${side}-${kind}`}>
          <summary>{activeView.title}</summary>
          <div className={styles.viewDetailsBody}>
            <div>{activeView.subtitle}</div>
            {activeView.details.length > 0 && (
              <ul className={styles.detailList}>
                {activeView.details.map((detail) => <li key={detail}>{detail}</li>)}
              </ul>
            )}
          </div>
        </details>
      )}

      {expanded && activeView && (
        <div
          className={styles.overlay}
          data-testid={`graph-overlay-${side}-${kind}`}
          onClick={(event) => {
            if (event.target === event.currentTarget) setExpanded(false)
          }}
        >
          <div className={styles.overlayModal} role='dialog' aria-modal='true' aria-label={activeView.title}>
            <div className={styles.overlayHeader}>
              <span className={styles.overlayTitle}>{activeView.title}</span>
              <button
                type='button'
                className={styles.overlayClose}
                onClick={() => setExpanded(false)}
                title='Close (Esc)'
                aria-label='Close the enlarged graph'
              >
                ✕
              </button>
            </div>
            <div className={styles.plotFillLarge}>
              <SubgraphPlot
                nodes={activeView.nodes}
                edges={activeView.edges}
                maxNodes={maxNodes}
                seedNode={activeView.seedNode}
              />
            </div>
            <p className={styles.overlayFooter}>{activeView.subtitle}</p>
          </div>
        </div>
      )}
    </div>
  )
}

function ComparisonColumn({ side, entry, bundle, onOpen }) {
  const tintClass = side === 'A' ? styles.tintA : styles.tintB
  const params = entry?.payloadSnapshot?.params || {}
  const summary = bundle.results?.summary || entry?.summary || {}
  const refs = bundle.results?.graph_pipeline?.graph_refs || []
  const lastClustering = bundle.results?.graph_pipeline?.last_clustering || {}
  const graphs = bundle.graphs

  const firstRef = refs[0] || null
  const secondRef = refs[1] || null
  const firstCluster = firstRef ? lastClustering[firstRef.cacheKey] || null : null
  const secondCluster = secondRef ? lastClustering[secondRef.cacheKey] || null : null

  // The two clusters are stable references from the results snapshot, so they are
  // safe (and sufficient) dependencies for the derived-subgraph computation.
  const derived = useMemo(
    () => buildDerivedSubgraphs(
      { nodes: firstCluster?.subgraph?.nodes || [], edges: firstCluster?.subgraph?.edges || [] },
      { nodes: secondCluster?.subgraph?.nodes || [], edges: secondCluster?.subgraph?.edges || [] }
    ),
    [firstCluster, secondCluster]
  )

  // `refs`, `lastClustering` and `graphs` only change when another analysis is
  // loaded, and `derived` is memoised on the two clusters.
  const conditionViews = useMemo(() => buildConditionGraphViews(refs, graphs), [refs, graphs])
  const clusterViews = useMemo(
    () => buildClusterViews(refs, lastClustering, derived),
    [refs, lastClustering, derived]
  )

  return (
    <div className={`${styles.column} ${tintClass}`} data-testid={`compare-column-${side.toLowerCase()}`}>
      <header className={styles.columnHeader}>
        <div className={styles.columnHeaderTop}>
          <span className={styles.sideBadge} data-testid={`compare-badge-${side.toLowerCase()}`}>{side}</span>
          <div className={styles.columnHeaderActions}>
            <span className={styles.statusPill} data-testid={`compare-status-${side.toLowerCase()}`}>
              {entry?.status || 'unknown'}
            </span>
            <button
              type='button'
              className={styles.openButton}
              data-testid={`compare-open-${side.toLowerCase()}`}
              onClick={() => onOpen(entry.id)}
            >
              Open analysis
            </button>
          </div>
        </div>
        <h2 className={styles.columnTitle}>{entry?.analysisName || entry?.scriptName || 'Untitled'}</h2>
        <p className={styles.columnMeta}>
          Created {formatDateTime(entry?.createdAt)}<br />
          Updated {formatDateTime(entry?.updatedAt)}
          {entry?.runCount ? ` · ${entry.runCount} run${entry.runCount > 1 ? 's' : ''}` : ''}
        </p>
      </header>

      {bundle.loading && <p className={styles.muted}>Loading stored results…</p>}
      {!bundle.loading && bundle.error && <p className={styles.errorText}>{bundle.error}</p>}

      <section className={styles.block}>
        <h3 className={styles.blockTitle}>Input files</h3>
        <dl className={styles.kv}>
          {FILE_SLOTS.map(([slot, label]) => {
            const files = entry?.payloadSnapshot?.files?.[slot] || []
            return (
              <div key={slot} className={styles.kvRow}>
                <dt>{label}</dt>
                <dd title={files.join('\n')}>{files.length === 0 ? '—' : files.map(fileName).join(', ')}</dd>
              </div>
            )
          })}
        </dl>
      </section>

      <section className={styles.block}>
        <h3 className={styles.blockTitle}>Parameters</h3>
        <KeyValueList object={params} labels={PARAM_LABELS} />
      </section>

      <section className={styles.block}>
        <h3 className={styles.blockTitle}>Result summary</h3>
        <KeyValueList object={summary} emptyLabel='No stored summary.' />
      </section>

      {!bundle.loading && !bundle.error && (
        <section className={styles.graphs} data-testid={`compare-graphs-${side.toLowerCase()}`}>
          <GraphPanel
            side={side.toLowerCase()}
            kind='conditions'
            title='Condition graphs'
            views={conditionViews}
            emptyMessage={bundle.loading ? 'Loading…' : 'No condition graph stored'}
          />
          <GraphPanel
            side={side.toLowerCase()}
            kind='clusters'
            title='Clusters'
            views={clusterViews}
            emptyMessage={bundle.loading ? 'Loading…' : 'No cluster stored'}
          />
        </section>
      )}
    </div>
  )
}

export default function ComparisonView({ entryA, entryB, onClose, onOpen }) {
  const bundleA = useJobBundle(entryA?.id)
  const bundleB = useJobBundle(entryB?.id)

  return (
    <div className={styles.view}>
      <header className={styles.viewHeader}>
        <div>
          <h1 className={styles.viewTitle}>Comparison</h1>
          <p className={styles.viewSubtitle}>
            The graphs produced by the two analyses, side by side. Condition graphs and clusters have
            their own area; use the drop-down to switch, and ⤢ to enlarge a graph. Each analysis keeps
            its own colour tint: <span className={styles.legendA}>blue</span> for A,{' '}
            <span className={styles.legendB}>amber</span> for B.
          </p>
        </div>
        <button type='button' className={styles.closeButton} onClick={onClose}>
          ✕ Close comparison
        </button>
      </header>

      <div className={styles.columns} data-testid='compare-columns'>
        <ComparisonColumn side='A' entry={entryA} bundle={bundleA} onOpen={onOpen} />
        <ComparisonColumn side='B' entry={entryB} bundle={bundleB} onOpen={onOpen} />
      </div>
    </div>
  )
}
