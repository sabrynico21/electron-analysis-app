import { useMemo, useState } from 'react'
import styles from './ChartPanel.module.css'
import SubgraphPlot from './SubgraphPlot'
import ExportableTable from './ExportableTable'
import { buildDerivedSubgraphs } from '../../utils/subgraphComparison'
import {
  conditionEdgeRows,
  clusterEdgeRows,
  derivedEdgeRows,
  intersectionEdgeRows,
  nodePresenceRows,
} from '../../utils/graphExports'

/** Turn a condition name into a filesystem-friendly prefix. */
const safeFileName = (value) => String(value || 'condition').replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '') || 'condition'

/** Filesystem-friendly single token; returns '' when there is nothing usable. */
const slugToken = (value) => String(value || '').replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '')

const countBySign = (edges = []) => {
  let negative = 0
  let positive = 0
  edges.forEach((edge) => {
    const weight = Number(edge?.weight) || 0
    if (weight > 0) positive += 1
    if (weight < 0) negative += 1
  })
  return { negative, positive }
}

const formatNegativeToPositiveRatio = (edges = []) => {
  const { negative, positive } = countBySign(edges)
  if (positive === 0) {
    if (negative === 0) return '0.00'
    return 'inf'
  }
  return (negative / positive).toFixed(2)
}

const formatMeanWeight = (edges = []) => {
  if (!edges.length) return '0.0000'
  const mean = edges.reduce((total, edge) => total + (Number(edge?.weight) || 0), 0) / edges.length
  return mean.toFixed(4)
}

/**
 * Conductance φ(S) measured on the extracted cluster (clustering metadata).
 * The backend serialises the undefined-conductance fallback (seed singleton)
 * as null, so null / non-finite values are shown as ∞; while no cluster is
 * available yet the box stays empty.
 */
const formatClusterConductance = (runStatus, conductance) => {
  if (runStatus !== 'done') return '-'
  if (conductance == null) return '∞'
  const value = Number(conductance)
  if (!Number.isFinite(value)) return '∞'
  return value.toFixed(4)
}

export default function ChartPanel({
  graphRefs = [],
  conditionGraphs = {},
  reclusterByGraph = {},
  selectedTaxa = '',
  runStatusByGraph = {},
  clusterParamsByGraph = {},
  onClusterParamsChange = () => {},
  onRerunCluster = () => {},
  interactionLocked = false,
  canUseRecomputeButtons = false,
}) {
  const [maxNodesByGraph, setMaxNodesByGraph] = useState({})
  const [viewModeByGraph, setViewModeByGraph] = useState({})
  const [maxIntersectionNodes, setMaxIntersectionNodes] = useState(60)
  const [maxComplementNodes, setMaxComplementNodes] = useState(60)

  const twoConditionRefs = graphRefs.slice(0, 2)
  const firstRef = twoConditionRefs[0]
  const secondRef = twoConditionRefs[1]
  // Seed taxon used to grow the cluster, prepended to the cluster export name so
  // clusters computed from different seeds do not overwrite each other on disk.
  const seedTaxaSlug = slugToken(selectedTaxa)
  const resolveStatus = (cacheKey) => {
    if (!cacheKey) return 'idle'
    const explicit = runStatusByGraph[cacheKey]
    if (explicit) return explicit
    const recluster = reclusterByGraph[cacheKey]
    if (recluster?.success && recluster?.subgraph) return 'done'
    if (recluster?.error) return 'error'
    return 'idle'
  }
  const firstStatus = resolveStatus(firstRef?.cacheKey)
  const secondStatus = resolveStatus(secondRef?.cacheKey)

  const firstSubgraph = firstRef ? (reclusterByGraph[firstRef.cacheKey]?.subgraph || {}) : {}
  const secondSubgraph = secondRef ? (reclusterByGraph[secondRef.cacheKey]?.subgraph || {}) : {}

  const derivedSubgraphs = useMemo(
    () => buildDerivedSubgraphs(firstSubgraph, secondSubgraph),
    [firstSubgraph, secondSubgraph],
  )

  const canRenderCombined = firstStatus === 'done' && secondStatus === 'done'

  if (twoConditionRefs.length === 0) return null

  return (
    <div className={styles.panel}>
      <div className={styles.grid}>
      {twoConditionRefs.map((ref, graphIndex) => {
        const payload = conditionGraphs[ref.cacheKey]
        const graph = payload?.graph
        const stats = graph?.stats
        const weightType = stats?.weight_type || graph?.edges?.[0]?.weight_type || 'correlation'
        const isStabilityWeight = weightType === 'selection_stability'
        const methodDiagnostics = stats?.method_diagnostics || {}
        const sampleCount = stats?.condition_samples_after_dataset_filter
          ?? stats?.condition_samples_in_metadata
          ?? null
        const recluster = reclusterByGraph[ref.cacheKey]
        const subgraph = recluster?.subgraph
        const runStatus = resolveStatus(ref.cacheKey)
        const hasTaxaInGraph = (graph?.nodes || []).some((node) => node.id === selectedTaxa)
        const baseSubgraphNodes = subgraph?.nodes || []
        const baseSubgraphEdges = subgraph?.edges || []
        const cardMode = viewModeByGraph[ref.cacheKey] || 'normal'
        const useExclusive = cardMode === 'exclusive' && canRenderCombined
        const conditionDerived = ref.cacheKey === firstRef?.cacheKey ? derivedSubgraphs.exclusiveA : derivedSubgraphs.exclusiveB
        const subgraphNodes = useExclusive ? (conditionDerived?.nodes || []) : baseSubgraphNodes
        const subgraphEdges = useExclusive ? (conditionDerived?.edges || []) : baseSubgraphEdges
        const defaultMaxNodes = Math.min(40, Math.max(1, subgraphNodes.length || 1))
        const maxNodes = Object.prototype.hasOwnProperty.call(maxNodesByGraph, ref.cacheKey)
          ? maxNodesByGraph[ref.cacheKey]
          : defaultMaxNodes
        // Transient display fallback until AnalysisWorkspace seeds the real
        // per-graph defaults (reference values: phi=0.4, b=10, alpha=0.05).
        const clusterParams = clusterParamsByGraph[ref.cacheKey] || { conductance: 0.4, b: 10, teleport: 0.05 }

        return (
          <div key={ref.cacheKey} className={styles.regionCard}>
            <div className={styles.regionHeader}>
              <h3>
                {ref.condition}
                {sampleCount !== null ? ` (${sampleCount} samples)` : ''}
              </h3>
              <span className={styles.regionBadge}>Condition</span>
            </div>
            {!graph && <p>Loading condition graph...</p>}
            {graph && (
              <>
                <section className={styles.metricSection}>
                  <h4 className={styles.metricTitle}>Graph Statistics</h4>
                  {stats?.evidence_type === 'stars' && (
                    <p className={styles.infoText}>
                      Edges are selected by StARS model stability; pairwise p-values are not produced.
                    </p>
                  )}
                  {(methodDiagnostics.selected_lambda != null || methodDiagnostics.lambda != null) && (
                    <p className={styles.infoText}>
                      Selected λ: {Number(methodDiagnostics.selected_lambda ?? methodDiagnostics.lambda).toPrecision(4)}
                      {methodDiagnostics.opt_index != null ? ` · StARS index ${methodDiagnostics.opt_index}` : ''}
                      {methodDiagnostics.admm_iteration_cap_count > 0
                        ? ` · ${methodDiagnostics.admm_iteration_cap_count} CCLasso ADMM calls hit the iteration cap`
                        : ''}
                    </p>
                  )}
                  <div className={`${styles.propsGrid} ${styles.propsGridRatio}`}>
                    <div className={styles.propBox}>
                      <span>Graph nodes</span>
                      <strong>{stats?.node_count ?? graph.nodes.length}</strong>
                    </div>
                    <div className={styles.propBox}>
                      <span>Graph edges</span>
                      <strong>{stats?.edge_count ?? graph.edges.length}</strong>
                    </div>
                    <div className={styles.propBox}>
                      <span>{isStabilityWeight ? 'Mean StARS stability' : 'N/P ratio'}</span>
                      <strong>{isStabilityWeight
                        ? formatMeanWeight(graph.edges)
                        : formatNegativeToPositiveRatio(graph.edges)}</strong>
                    </div>
                    <div className={styles.propBox}>
                      <span>Taxa in graph</span>
                      <strong>{selectedTaxa ? (hasTaxaInGraph ? 'YES' : '-') : '-'}</strong>
                    </div>
                  </div>
                </section>

                <section className={styles.metricSection}>
                  <h4 className={styles.metricTitle}>Cluster Results</h4>

                  <div className={styles.statusRow}>
                    <span className={`${styles.statusDot} ${styles[`status_${runStatus}`]}`} />
                    <p className={styles.statusText}>
                      {runStatus === 'idle' && 'Waiting to start clustering.'}
                      {runStatus === 'running' && 'Algorithm running from selected taxa...'}
                      {runStatus === 'missing-taxa' && 'Selected taxa is not present in this condition.'}
                      {runStatus === 'done' && 'Cluster computed successfully.'}
                      {runStatus === 'error' && (recluster?.error || 'Error while clustering.')}
                    </p>
                  </div>

                  <div className={styles.subgraphSection}>
                    <h4>{useExclusive ? 'Cluster subgraph without intersection' : 'Cluster subgraph'}</h4>
                    <div className={`${styles.propsGrid} ${styles.propsGridRatio}`}>
                      <div className={styles.propBox}>
                        <span>Cluster nodes</span>
                        <strong>{subgraphNodes.length}</strong>
                      </div>
                      <div className={styles.propBox}>
                        <span>Cluster edges</span>
                        <strong>{subgraphEdges.length}</strong>
                      </div>
                      <div className={styles.propBox}>
                        <span>N/P ratio</span>
                        <strong>{formatNegativeToPositiveRatio(subgraphEdges)}</strong>
                      </div>
                      <div className={styles.propBox}>
                        <span>Conductance</span>
                        <strong>{formatClusterConductance(runStatus, recluster?.metadata?.conductance)}</strong>
                      </div>
                    </div>
                    {useExclusive && runStatus === 'done' && subgraphNodes.length === 0 && (
                      <p className={styles.infoText}>No exclusive nodes: this condition cluster matches the intersection.</p>
                    )}
                  </div>

                  <div className={styles.subgraphSection}>
                    {runStatus === 'missing-taxa' && (
                      <div className={styles.warningBox}>
                        The plot cannot be generated.
                      </div>
                    )}
                    <div className={styles.controlsRow}>
                      <div className={styles.nodeLimitControl}>
                        <label htmlFor={`max-nodes-${ref.cacheKey}`}>
                          Nodes to display <span className={styles.range}>(1 - {subgraphNodes.length || 0})</span>
                        </label>
                        <input
                          id={`max-nodes-${ref.cacheKey}`}
                          className={styles.nodeLimitInput}
                          type='number'
                          min='1'
                          max={String(Math.max(1, subgraphNodes.length || 1))}
                          value={maxNodes}
                          disabled={interactionLocked || runStatus === 'missing-taxa'}
                          onChange={(event) => {
                            const raw = event.target.value
                            if (raw === '') {
                              setMaxNodesByGraph((prev) => ({
                                ...prev,
                                [ref.cacheKey]: '',
                              }))
                              return
                            }
                            const value = Number(raw)
                            if (Number.isFinite(value)) {
                              const maxAvailable = Math.max(1, subgraphNodes.length || 1)
                              const clamped = Math.min(maxAvailable, value)
                              setMaxNodesByGraph((prev) => ({
                                ...prev,
                                [ref.cacheKey]: clamped,
                              }))
                            }
                          }}
                          onBlur={(event) => {
                            const parsed = Number(event.target.value)
                            const maxAvailable = Math.max(1, subgraphNodes.length || 1)
                            const normalized = Number.isFinite(parsed) ? Math.max(1, Math.min(maxAvailable, parsed)) : maxAvailable
                            setMaxNodesByGraph((prev) => ({
                              ...prev,
                              [ref.cacheKey]: normalized,
                            }))
                          }}
                        />
                        <small>Available max: {subgraphNodes.length || 0}</small>
                      </div>
                      <div className={styles.viewModeRow}>
                        <span>Visualization</span>
                        <div className={styles.viewModeButtons}>
                          <button
                            type='button'
                            className={cardMode === 'normal' ? `${styles.viewModeButton} ${styles.viewModeButtonActive}` : styles.viewModeButton}
                            onClick={() => setViewModeByGraph((prev) => ({ ...prev, [ref.cacheKey]: 'normal' }))}
                            disabled={interactionLocked}
                          >
                            Full cluster
                          </button>
                          <button
                            type='button'
                            className={cardMode === 'exclusive' ? `${styles.viewModeButton} ${styles.viewModeButtonActive}` : styles.viewModeButton}
                            onClick={() => setViewModeByGraph((prev) => ({ ...prev, [ref.cacheKey]: 'exclusive' }))}
                            disabled={interactionLocked || !canRenderCombined}
                          >
                            Without intersection
                          </button>
                        </div>
                      </div>
                    </div>
                    {runStatus !== 'missing-taxa' && (
                      <SubgraphPlot
                        nodes={subgraphNodes}
                        edges={subgraphEdges}
                        maxNodes={maxNodes}
                        seedNode={selectedTaxa}
                        controlsDisabled={interactionLocked}
                      />
                    )}
                    <div className={styles.clusterParamsBox}>
                      <h4>Clustering parameters</h4>
                      <div className={styles.clusterParamsGrid}>
                        <div className={styles.clusterParamField}>
                          <label htmlFor={`conductance-${ref.cacheKey}`}>
                            Conductance <span className={styles.range}>(0.0 - 1.0)</span>
                          </label>
                          <input
                            id={`conductance-${ref.cacheKey}`}
                            className={styles.clusterParamInput}
                            type='number'
                            min='0'
                            max='1'
                            step='0.01'
                            value={clusterParams.conductance}
                            disabled={interactionLocked}
                            onChange={(event) => {
                              const raw = event.target.value
                              if (raw === '') {
                                onClusterParamsChange(ref.cacheKey, { conductance: '' })
                                return
                              }
                              const nextValue = Number(raw)
                              if (Number.isFinite(nextValue)) {
                                onClusterParamsChange(ref.cacheKey, { conductance: nextValue })
                              }
                            }}
                            onBlur={(event) => {
                              const parsed = Number(event.target.value)
                              // Invalid input resets to the reference default.
                              const normalized = Number.isFinite(parsed) ? Math.max(0, Math.min(1, parsed)) : 0.4
                              onClusterParamsChange(ref.cacheKey, { conductance: normalized })
                            }}
                          />
                        </div>
                        <div className={styles.clusterParamField}>
                          <label htmlFor={`b-${ref.cacheKey}`}>
                            b <span className={styles.range}>(integer ≥ 1)</span>
                          </label>
                          <input
                            id={`b-${ref.cacheKey}`}
                            className={styles.clusterParamInput}
                            type='number'
                            min='1'
                            step='1'
                            value={clusterParams.b}
                            disabled={interactionLocked}
                            onChange={(event) => {
                              const raw = event.target.value
                              if (raw === '') {
                                onClusterParamsChange(ref.cacheKey, { b: '' })
                                return
                              }
                              const nextValue = Number(raw)
                              if (Number.isFinite(nextValue)) {
                                onClusterParamsChange(ref.cacheKey, { b: nextValue })
                              }
                            }}
                            onBlur={(event) => {
                              const parsed = Number(event.target.value)
                              const normalized = Math.max(1, Math.floor(Number.isFinite(parsed) ? parsed : 1))
                              onClusterParamsChange(ref.cacheKey, { b: normalized })
                            }}
                          />
                        </div>
                        <div className={styles.clusterParamField}>
                          <label htmlFor={`teleport-${ref.cacheKey}`}>
                            Teleport (alpha) <span className={styles.range}>(0 - 1.0)</span>
                          </label>
                          <input
                            id={`teleport-${ref.cacheKey}`}
                            className={styles.clusterParamInput}
                            type='number'
                            min='0'
                            max='1'
                            step='0.000001'
                            value={clusterParams.teleport}
                            disabled={interactionLocked}
                            onChange={(event) => {
                              const raw = event.target.value
                              if (raw === '') {
                                onClusterParamsChange(ref.cacheKey, { teleport: '' })
                                return
                              }
                              const parsed = Number(raw)
                              if (!Number.isFinite(parsed)) return
                              const clamped = Math.max(0, Math.min(1, parsed))
                              onClusterParamsChange(ref.cacheKey, { teleport: clamped })
                            }}
                            onBlur={() => {
                              const parsed = Number(clusterParams.teleport)
                              if (!Number.isFinite(parsed)) {
                                // Invalid input resets to the reference default.
                                onClusterParamsChange(ref.cacheKey, { teleport: 0.05 })
                                return
                              }
                              const normalized = Number(Math.max(0, Math.min(1, parsed)).toFixed(6))
                              onClusterParamsChange(ref.cacheKey, { teleport: normalized })
                            }}
                          />
                        </div>
                      </div>
                      <button
                        type='button'
                        className={styles.rerunButton}
                        onClick={() => onRerunCluster(ref.cacheKey)}
                        disabled={!selectedTaxa.trim() || !canUseRecomputeButtons || runStatus === 'running'}
                      >
                        Recompute this cluster
                      </button>
                    </div>
                  </div>
                </section>

                <section className={styles.metricSection}>
                  <h4 className={styles.metricTitle}>Downloadable tables</h4>
                  <div className={styles.exportTables}>
                    <ExportableTable
                      title={`Full condition graph · ${ref.condition}`}
                      description='Every edge with its method-specific weight, weight type, and available p-value evidence.'
                      rows={conditionEdgeRows(graph, ref.condition)}
                      fileName={`${safeFileName(ref.condition)}_full_graph_edges.csv`}
                      testId={`export-full-graph-${graphIndex}`}
                    />
                    <ExportableTable
                      title={`Cluster · ${ref.condition}`}
                      description='Edges of the cluster subgraph, with method-specific edge weights and available p-value evidence.'
                      rows={clusterEdgeRows(subgraph, ref.condition)}
                      fileName={`${seedTaxaSlug ? `${seedTaxaSlug}_` : ''}${safeFileName(ref.condition)}_cluster_edges.csv`}
                      testId={`export-cluster-${graphIndex}`}
                      emptyLabel='No cluster computed yet.'
                    />
                    <ExportableTable
                      title={`Only in ${ref.condition} · edges`}
                      description='Edges of this cluster that are not part of the intersection with the other condition (same as the "Without intersection" view).'
                      rows={derivedEdgeRows(conditionDerived?.edges || [], ref.condition)}
                      fileName={`${safeFileName(ref.condition)}_only_edges.csv`}
                      testId={`export-exclusive-${graphIndex}`}
                      emptyLabel='No exclusive edge.'
                    />
                  </div>
                </section>
              </>
            )}
          </div>
        )
      })}
      </div>

      <div className={styles.comparisonSection}>
        <h3>Comparison of the two clusters</h3>
        {!canRenderCombined && (
          <p className={styles.infoText}>Intersection and complement graphs will be available when both conditions have a valid cluster.</p>
        )}
        {canRenderCombined && (
          <>
          <div className={styles.comparisonGrid}>
            <div className={styles.regionCard}>
              <div className={styles.regionHeader}>
                <h3>Subgraph intersection</h3>
                <span className={styles.regionBadge}>Shared</span>
              </div>
              <div className={styles.propsGrid}>
                <div className={styles.propBox}>
                  <span>Nodes</span>
                  <strong>{derivedSubgraphs.intersection.nodes.length}</strong>
                </div>
                <div className={styles.propBox}>
                  <span>Edges</span>
                  <strong>{derivedSubgraphs.intersection.edges.length}</strong>
                </div>
              </div>
              <div className={styles.nodeLimitControl}>
                <label htmlFor='intersection-max-nodes'>
                  Nodes to display <span className={styles.range}>(1 - {derivedSubgraphs.intersection.nodes.length || 0})</span>
                </label>
                <input
                  id='intersection-max-nodes'
                  className={styles.nodeLimitInput}
                  type='number'
                  min='1'
                  max={String(Math.max(1, derivedSubgraphs.intersection.nodes.length || 1))}
                  value={maxIntersectionNodes}
                  disabled={interactionLocked}
                  onChange={(event) => {
                    const raw = event.target.value
                    if (raw === '') {
                      setMaxIntersectionNodes('')
                      return
                    }
                    const value = Number(raw)
                    if (Number.isFinite(value)) {
                      const maxAvailable = Math.max(1, derivedSubgraphs.intersection.nodes.length || 1)
                      setMaxIntersectionNodes(Math.min(maxAvailable, value))
                    }
                  }}
                  onBlur={(event) => {
                    const parsed = Number(event.target.value)
                    const maxAvailable = Math.max(1, derivedSubgraphs.intersection.nodes.length || 1)
                    const normalized = Number.isFinite(parsed) ? Math.max(1, Math.min(maxAvailable, parsed)) : maxAvailable
                    setMaxIntersectionNodes(normalized)
                  }}
                />
              </div>
              <SubgraphPlot
                nodes={derivedSubgraphs.intersection.nodes}
                edges={derivedSubgraphs.intersection.edges}
                maxNodes={maxIntersectionNodes}
                seedNode={selectedTaxa}
                controlsDisabled={interactionLocked}
                conflictLabels={{ first: firstRef?.condition, second: secondRef?.condition }}
              />
            </div>

            <div className={styles.regionCard}>
              <div className={styles.regionHeader}>
                <h3>Intersection complement</h3>
                <span className={styles.regionBadge}>Complement</span>
              </div>
              <div className={styles.propsGrid}>
                <div className={styles.propBox}>
                  <span>Nodes</span>
                  <strong>{derivedSubgraphs.complement.nodes.length}</strong>
                </div>
                <div className={styles.propBox}>
                  <span>Edges</span>
                  <strong>{derivedSubgraphs.complement.edges.length}</strong>
                </div>
              </div>
              <div className={styles.nodeLimitControl}>
                <label htmlFor='complement-max-nodes'>
                  Nodes to display <span className={styles.range}>(1 - {derivedSubgraphs.complement.nodes.length || 0})</span>
                </label>
                <input
                  id='complement-max-nodes'
                  className={styles.nodeLimitInput}
                  type='number'
                  min='1'
                  max={String(Math.max(1, derivedSubgraphs.complement.nodes.length || 1))}
                  value={maxComplementNodes}
                  disabled={interactionLocked}
                  onChange={(event) => {
                    const raw = event.target.value
                    if (raw === '') {
                      setMaxComplementNodes('')
                      return
                    }
                    const value = Number(raw)
                    if (Number.isFinite(value)) {
                      const maxAvailable = Math.max(1, derivedSubgraphs.complement.nodes.length || 1)
                      setMaxComplementNodes(Math.min(maxAvailable, value))
                    }
                  }}
                  onBlur={(event) => {
                    const parsed = Number(event.target.value)
                    const maxAvailable = Math.max(1, derivedSubgraphs.complement.nodes.length || 1)
                    const normalized = Number.isFinite(parsed) ? Math.max(1, Math.min(maxAvailable, parsed)) : maxAvailable
                    setMaxComplementNodes(normalized)
                  }}
                />
              </div>
              <SubgraphPlot
                nodes={derivedSubgraphs.complement.nodes}
                edges={derivedSubgraphs.complement.edges}
                maxNodes={maxComplementNodes}
                seedNode={selectedTaxa}
                controlsDisabled={interactionLocked}
              />
            </div>
          </div>

          <div className={styles.exportTables}>
            <ExportableTable
              title={`Intersection · ${firstRef?.condition || 'condition 1'} ∩ ${secondRef?.condition || 'condition 2'}`}
              description='Edges joining the same two nodes in both clusters. Edges with opposite signs keep both correlation values.'
              rows={intersectionEdgeRows(
                derivedSubgraphs.intersection.edges,
                safeFileName(firstRef?.condition),
                safeFileName(secondRef?.condition),
              )}
              fileName='intersection_edges.csv'
              testId='export-intersection'
              emptyLabel='No shared edge between the two clusters.'
            />
            <ExportableTable
              title='Complement · union without intersection'
              description='Edges present in exactly one cluster (nodes shared by both clusters are excluded).'
              rows={derivedEdgeRows(derivedSubgraphs.complement.edges, 'complement')}
              fileName='complement_edges.csv'
              testId='export-complement'
              emptyLabel='The complement is empty.'
            />
            <ExportableTable
              title='Nodes present in only one condition'
              description='For every node in exactly one cluster: whether it is missing only from the other cluster or from the other condition graph entirely.'
              rows={nodePresenceRows(
                reclusterByGraph[firstRef?.cacheKey]?.subgraph,
                reclusterByGraph[secondRef?.cacheKey]?.subgraph,
                conditionGraphs[firstRef?.cacheKey]?.graph,
                conditionGraphs[secondRef?.cacheKey]?.graph,
                firstRef?.condition || 'condition 1',
                secondRef?.condition || 'condition 2',
              )}
              fileName='nodes_only_in_one_condition.csv'
              testId='export-node-presence'
              emptyLabel='Every cluster node is shared by both conditions.'
            />
          </div>
          </>
        )}
      </div>
    </div>
  )
}
