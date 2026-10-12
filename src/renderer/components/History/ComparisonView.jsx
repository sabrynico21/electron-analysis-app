/**
 * ComparisonView
 *
 * Two stored analyses side by side, each with its own colour tint (blue for A,
 * amber for B) and both stretched to the same height.
 *
 * Shape of a column:
 *   - header: side badge, name, timestamps, and the run status next to the
 *     "Open analysis" button (top right);
 *   - input files, parameters and result summary. The summary only reports the
 *     statistics that are not already shown as parameters (the analysed
 *     conditions, the correlation method, the thresholds…), so the same value is
 *     never displayed twice;
 *   - **one graph area of fixed size** for the clusters (including the
 *     intersection / exclusive / complement views of the two clusters), with a
 *     single `select` switcher that always shows the full selected label, a node
 *     limit, a button that opens the graph enlarged in an overlay, and a
 *     collapsible line of details.
 *
 * The fixed area keeps an analysis with fewer options exactly as tall as one with
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
import ExportableTable from '../Results/ExportableTable'
import { buildDerivedSubgraphs } from '../../utils/subgraphComparison'
import {
  conditionEdgeRows,
  clusterEdgeRows,
  derivedEdgeRows,
  intersectionEdgeRows,
  nodePresenceRows,
} from '../../utils/graphExports'
import styles from './ComparisonView.module.css'

/** Turn a condition name into a filesystem-friendly prefix. */
const safeFileName = (value) => String(value || 'condition').replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '') || 'condition'

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

/**
 * Result-summary keys that merely repeat a parameter already listed in the
 * "Parameters" block (the analysed conditions above all), so the "Result
 * summary" must not report them a second time.
 */
const SUMMARY_DUPLICATE_KEYS = new Set([
  'dataset_type',        // datasetType
  'correlation_method',  // correlationMethod
  'grouping_mode',       // groupingMode
  'p_value_threshold',   // pValueThreshold
  'condition_a',         // conditionA
  'condition_b',         // conditionB
  'selected_conditions', // the very same two conditions, listed once more
])

/** The stored summary without the entries already shown as parameters. */
function summaryWithoutDuplicates(summary) {
  return Object.fromEntries(
    Object.entries(summary || {}).filter(([key]) => !SUMMARY_DUPLICATE_KEYS.has(key))
  )
}

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

/** Load one job's stored results (its clustering snapshot included). */
function useJobBundle(jobId) {
  const [state, setState] = useState({ loading: true, results: null, error: null })

  useEffect(() => {
    let cancelled = false

    const load = async () => {
      setState({ loading: true, results: null, error: null })
      if (!jobId) {
        setState({ loading: false, results: null, error: 'Missing analysis id' })
        return
      }

      const results = await window.electronAPI.getResults(jobId)
      if (cancelled) return
      if (!results) {
        setState({ loading: false, results: null, error: 'No stored results for this analysis' })
        return
      }

      setState({ loading: false, results, error: null })
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
function GraphPanel({ side, kind, title, views, emptyMessage, conflictLabels = null }) {
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
            conflictLabels={conflictLabels}
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
                conflictLabels={conflictLabels}
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
  // The summary repeats the analysed conditions, the method and the thresholds
  // already shown under "Parameters": drop those copies here.
  const summary = summaryWithoutDuplicates(bundle.results?.summary || entry?.summary || {})
  const refs = bundle.results?.graph_pipeline?.graph_refs || []
  const lastClustering = bundle.results?.graph_pipeline?.last_clustering || {}

  const firstRef = refs[0] || null
  const secondRef = refs[1] || null
  const firstCluster = firstRef ? lastClustering[firstRef.cacheKey] || null : null
  const secondCluster = secondRef ? lastClustering[secondRef.cacheKey] || null : null

  // The full condition graphs are part of the stored results (used here for the
  // exports: they tell whether a cluster node is also in the condition graph).
  const conditionGraphs = bundle.results?.graph_pipeline?.condition_graphs || []
  const graphByCondition = useMemo(
    () => new Map(conditionGraphs.map((cg) => [cg.condition, cg.graph])),
    [conditionGraphs]
  )
  const firstGraph = firstRef ? graphByCondition.get(firstRef.condition) || null : null
  const secondGraph = secondRef ? graphByCondition.get(secondRef.condition) || null : null
  const labelA = firstRef?.condition || 'condition 1'
  const labelB = secondRef?.condition || 'condition 2'
  const prefix = side.toLowerCase()

  // The two clusters are stable references from the results snapshot, so they are
  // safe (and sufficient) dependencies for the derived-subgraph computation.
  const derived = useMemo(
    () => buildDerivedSubgraphs(
      { nodes: firstCluster?.subgraph?.nodes || [], edges: firstCluster?.subgraph?.edges || [] },
      { nodes: secondCluster?.subgraph?.nodes || [], edges: secondCluster?.subgraph?.edges || [] }
    ),
    [firstCluster, secondCluster]
  )

  // `refs` and `lastClustering` only change when another analysis is loaded,
  // and `derived` is memoised on the two clusters.
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
            kind='clusters'
            title='Clusters'
            views={clusterViews}
            emptyMessage={bundle.loading ? 'Loading…' : 'No cluster stored'}
            conflictLabels={{ first: labelA, second: labelB }}
          />
        </section>
      )}

      {!bundle.loading && !bundle.error && firstRef && secondRef && (
        <section className={styles.exports} data-testid={`compare-exports-${prefix}`}>
          <h3 className={styles.blockTitle}>Downloadable tables</h3>
          <div className={styles.exportTables}>
            <ExportableTable
              title={`Full condition graph · ${labelA}`}
              description='Every edge of the condition graph with its correlation value.'
              rows={conditionEdgeRows(firstGraph, labelA)}
              fileName={`${safeFileName(labelA)}_full_graph_edges.csv`}
              testId={`compare-export-full-graph-${prefix}`}
            />
            <ExportableTable
              title={`Full condition graph · ${labelB}`}
              description='Every edge of the condition graph with its correlation value.'
              rows={conditionEdgeRows(secondGraph, labelB)}
              fileName={`${safeFileName(labelB)}_full_graph_edges.csv`}
              testId={`compare-export-full-graph-${prefix}-b`}
            />
            <ExportableTable
              title={`Cluster · ${labelA}`}
              description='Edges of the cluster subgraph, with the correlation value of each edge.'
              rows={clusterEdgeRows(firstCluster?.subgraph, labelA)}
              fileName={`${safeFileName(labelA)}_cluster_edges.csv`}
              testId={`compare-export-cluster-${prefix}`}
              emptyLabel='No cluster computed yet.'
            />
            <ExportableTable
              title={`Only in ${labelA} · edges`}
              description='Edges of the first cluster that are not part of the intersection (same as the "Only in" view).'
              rows={derivedEdgeRows(derived.exclusiveA.edges, labelA)}
              fileName={`${safeFileName(labelA)}_only_edges.csv`}
              testId={`compare-export-exclusive-${prefix}`}
              emptyLabel='No exclusive edge.'
            />
            <ExportableTable
              title={`Cluster · ${labelB}`}
              description='Edges of the cluster subgraph, with the correlation value of each edge.'
              rows={clusterEdgeRows(secondCluster?.subgraph, labelB)}
              fileName={`${safeFileName(labelB)}_cluster_edges.csv`}
              testId={`compare-export-cluster-${prefix}-b`}
              emptyLabel='No cluster computed yet.'
            />
            <ExportableTable
              title={`Only in ${labelB} · edges`}
              description='Edges of the second cluster that are not part of the intersection (same as the "Only in" view).'
              rows={derivedEdgeRows(derived.exclusiveB.edges, labelB)}
              fileName={`${safeFileName(labelB)}_only_edges.csv`}
              testId={`compare-export-exclusive-${prefix}-b`}
              emptyLabel='No exclusive edge.'
            />
            <ExportableTable
              title={`Intersection · ${labelA} ∩ ${labelB}`}
              description='Edges joining the same two nodes in both clusters. Edges with opposite signs keep both correlation values.'
              rows={intersectionEdgeRows(
                derived.intersection.edges,
                safeFileName(labelA),
                safeFileName(labelB),
              )}
              fileName={`compare_${prefix}_intersection_edges.csv`}
              testId={`compare-export-intersection-${prefix}`}
              emptyLabel='No shared edge between the two clusters.'
            />
            <ExportableTable
              title='Complement · union without intersection'
              description='Edges present in exactly one cluster (nodes shared by both clusters are excluded).'
              rows={derivedEdgeRows(derived.complement.edges, 'complement')}
              fileName={`compare_${prefix}_complement_edges.csv`}
              testId={`compare-export-complement-${prefix}`}
              emptyLabel='The complement is empty.'
            />
            <ExportableTable
              title='Nodes present in only one condition'
              description='For every node in exactly one cluster: whether it is missing only from the other cluster or from the other condition graph entirely.'
              rows={nodePresenceRows(
                firstCluster?.subgraph,
                secondCluster?.subgraph,
                firstGraph,
                secondGraph,
                labelA,
                labelB,
              )}
              fileName={`compare_${prefix}_nodes_only_in_one_condition.csv`}
              testId={`compare-export-node-presence-${prefix}`}
              emptyLabel='Every cluster node is shared by both conditions.'
            />
          </div>
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
            The stored results of the two analyses, side by side. Only the clusters are drawn — use the
            drop-down to switch between them and ⤢ to enlarge one. Values already shown as parameters
            are not repeated in the result summary. Each analysis keeps its own colour tint:{' '}
            <span className={styles.legendA}>blue</span> for A,{' '}
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
