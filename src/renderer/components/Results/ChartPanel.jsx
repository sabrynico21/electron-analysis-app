import { useMemo, useState } from 'react'
import styles from './ChartPanel.module.css'
import SubgraphPlot from './SubgraphPlot'
import { buildDerivedSubgraphs } from '../../utils/subgraphComparison'

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
      {twoConditionRefs.map((ref) => {
        const payload = conditionGraphs[ref.cacheKey]
        const graph = payload?.graph
        const stats = graph?.stats
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
        const clusterParams = clusterParamsByGraph[ref.cacheKey] || { conductance: 0.2, b: 1, teleport: 0 }

        return (
          <div key={ref.cacheKey} className={styles.regionCard}>
            <div className={styles.regionHeader}>
              <h3>{ref.condition}</h3>
              <span className={styles.regionBadge}>Condition</span>
            </div>
            {!graph && <p>Loading condition graph...</p>}
            {graph && (
              <>
                <section className={styles.metricSection}>
                  <h4 className={styles.metricTitle}>Graph Statistics</h4>
                  <div className={styles.propsGrid}>
                    <div className={styles.propBox}>
                      <span>Graph nodes</span>
                      <strong>{stats?.node_count ?? graph.nodes.length}</strong>
                    </div>
                    <div className={styles.propBox}>
                      <span>Graph edges</span>
                      <strong>{stats?.edge_count ?? graph.edges.length}</strong>
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
                    <div className={styles.propsGrid}>
                      <div className={styles.propBox}>
                        <span>Cluster nodes</span>
                        <strong>{subgraphNodes.length}</strong>
                      </div>
                      <div className={styles.propBox}>
                        <span>Cluster edges</span>
                        <strong>{subgraphEdges.length}</strong>
                      </div>
                      <div className={styles.propBox}>
                        <span>Result status</span>
                        <strong>{runStatus === 'done' ? 'Ready' : 'Waiting'}</strong>
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
                              const normalized = Number.isFinite(parsed) ? Math.max(0, Math.min(1, parsed)) : 0.2
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
                                onClusterParamsChange(ref.cacheKey, { teleport: 0 })
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
        )}
      </div>
    </div>
  )
}
