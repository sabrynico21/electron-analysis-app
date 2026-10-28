/**
 * Shared helpers for comparing two subgraphs.
 *
 * Used both by the single-analysis view (`ChartPanel`, which compares the two
 * *conditions* of one analysis) and by the history comparison view (which
 * compares two whole analyses), so the two views cannot drift apart.
 *
 * ## Intersection semantics
 *
 * Two edges are considered "the same edge" when they join the same unordered
 * pair of nodes, regardless of the correlation *sign*. The sign is then compared:
 *
 *  - same sign in both conditions  → the edge is drawn as a single solid line
 *    (red = positive, blue = negative);
 *  - opposite signs                → the edge is kept in the intersection but
 *    flagged `signConflict: true` and `conflictType` (`pos_to_neg` when the first
 *    condition is positive and the second negative, `neg_to_pos` otherwise), so
 *    the plot can render it as a dashed edge whose colour encodes the direction
 *    of the sign flip.
 *
 * Consequently the **complement** (union of the two clusters minus the
 * intersection) subtracts the *node-pair* intersection: a pair present in both
 * conditions is part of the intersection, even when the two signs disagree.
 */

/** Undirected node-pair identity of an edge, ignoring the sign. */
export function edgePairKey(edge = {}) {
  const source = String(edge.source || '')
  const target = String(edge.target || '')
  const ordered = [source, target].sort()
  return `${ordered[0]}::${ordered[1]}`
}

/** Stable identity for an edge including its sign (kept for external callers). */
export function edgeKey(edge = {}) {
  const weight = Number(edge.weight) || 0
  const sign = weight > 0 ? 'p' : weight < 0 ? 'n' : 'z'
  return `${edgePairKey(edge)}::${sign}`
}

/** Human-readable sign of a correlation weight. */
export function edgeSign(weight) {
  const numeric = Number(weight) || 0
  if (numeric > 0) return 'positive'
  if (numeric < 0) return 'negative'
  return 'neutral'
}

const hasBothEndpoints = (edge, nodeIds) => nodeIds.has(edge.source) && nodeIds.has(edge.target)

/**
 * Derive the intersection / exclusive / complement views of two subgraphs.
 * @param {{nodes?: object[], edges?: object[]}} firstSubgraph
 * @param {{nodes?: object[], edges?: object[]}} secondSubgraph
 */
export function buildDerivedSubgraphs(firstSubgraph = {}, secondSubgraph = {}) {
  const firstNodes = firstSubgraph?.nodes || []
  const secondNodes = secondSubgraph?.nodes || []
  const firstEdges = firstSubgraph?.edges || []
  const secondEdges = secondSubgraph?.edges || []

  const firstNodeIds = new Set(firstNodes.map((node) => node.id))
  const secondNodeIds = new Set(secondNodes.map((node) => node.id))
  const firstNodeMap = new Map(firstNodes.map((node) => [node.id, node]))
  const secondNodeMap = new Map(secondNodes.map((node) => [node.id, node]))

  const intersectionNodeIds = new Set([...firstNodeIds].filter((id) => secondNodeIds.has(id)))
  const intersectionNodes = [...intersectionNodeIds]
    .map((id) => firstNodeMap.get(id) || secondNodeMap.get(id))
    .filter(Boolean)

  const firstEdgeMap = new Map(firstEdges.map((edge) => [edgePairKey(edge), edge]))
  const secondEdgeMap = new Map(secondEdges.map((edge) => [edgePairKey(edge), edge]))
  const intersectionPairKeys = new Set(
    [...firstEdgeMap.keys()].filter((key) => secondEdgeMap.has(key))
  )

  const intersectionEdges = [...intersectionPairKeys].map((key) => {
    const firstEdge = firstEdgeMap.get(key)
    const secondEdge = secondEdgeMap.get(key)
    const weightA = Number(firstEdge.weight) || 0
    const weightB = Number(secondEdge.weight) || 0
    const weightTypeA = firstEdge.weight_type || 'correlation'
    const weightTypeB = secondEdge.weight_type || 'correlation'
    const metricMismatch = weightTypeA !== weightTypeB
    const weightType = metricMismatch ? 'mixed' : weightTypeA
    const signA = metricMismatch
      ? 'not comparable'
      : weightType === 'selection_stability' ? 'selected' : edgeSign(weightA)
    const signB = metricMismatch
      ? 'not comparable'
      : weightType === 'selection_stability' ? 'selected' : edgeSign(weightB)
    const signConflict = !metricMismatch && weightType !== 'selection_stability' && signA !== signB
    // Which way the sign flips, so the plot can give each direction its own
    // colour: positive→negative vs negative→positive.
    const conflictType = signConflict
      ? (signA === 'positive' ? 'pos_to_neg' : 'neg_to_pos')
      : ''
    // For a solid edge both signs agree, so `weight` only drives the stroke
    // width; for a conflicting edge the two halves use weightA/weightB.
    const magnitudeSource = Math.abs(weightA) >= Math.abs(weightB) ? weightA : weightB
    return {
      source: firstEdge.source,
      target: firstEdge.target,
      weight: signConflict ? magnitudeSource : weightA,
      weightA,
      weightB,
      pValueA: firstEdge.p_value ?? null,
      pValueB: secondEdge.p_value ?? null,
      weight_type: weightType,
      weight_type_a: weightTypeA,
      weight_type_b: weightTypeB,
      metricMismatch,
      evidence_type: firstEdge.evidence_type || secondEdge.evidence_type || 'pairwise_p_value',
      signA,
      signB,
      signConflict,
      conflictType,
    }
  })

  const firstExclusiveNodeIds = new Set([...firstNodeIds].filter((id) => !intersectionNodeIds.has(id)))
  const secondExclusiveNodeIds = new Set([...secondNodeIds].filter((id) => !intersectionNodeIds.has(id)))
  const firstExclusiveNodes = [...firstExclusiveNodeIds].map((id) => firstNodeMap.get(id)).filter(Boolean)
  const secondExclusiveNodes = [...secondExclusiveNodeIds].map((id) => secondNodeMap.get(id)).filter(Boolean)

  // "Only in <condition>": edges of that cluster that are not part of the
  // intersection and whose two endpoints are both exclusive to that condition.
  const firstExclusiveEdges = firstEdges.filter((edge) => (
    !intersectionPairKeys.has(edgePairKey(edge)) && hasBothEndpoints(edge, firstExclusiveNodeIds)
  ))
  const secondExclusiveEdges = secondEdges.filter((edge) => (
    !intersectionPairKeys.has(edgePairKey(edge)) && hasBothEndpoints(edge, secondExclusiveNodeIds)
  ))

  const unionNodeMap = new Map([...firstNodes, ...secondNodes].map((node) => [node.id, node]))
  const complementNodeIds = new Set([...unionNodeMap.keys()].filter((id) => !intersectionNodeIds.has(id)))
  const complementNodes = [...complementNodeIds].map((id) => unionNodeMap.get(id)).filter(Boolean)
  const complementEdges = [...firstEdges, ...secondEdges]
    .filter((edge) => !intersectionPairKeys.has(edgePairKey(edge)) && hasBothEndpoints(edge, complementNodeIds))

  return {
    intersection: { nodes: intersectionNodes, edges: intersectionEdges, pairKeys: intersectionPairKeys },
    exclusiveA: { nodes: firstExclusiveNodes, edges: firstExclusiveEdges, nodeIds: firstExclusiveNodeIds },
    exclusiveB: { nodes: secondExclusiveNodes, edges: secondExclusiveEdges, nodeIds: secondExclusiveNodeIds },
    complement: { nodes: complementNodes, edges: complementEdges, nodeIds: complementNodeIds },
  }
}
