/**
 * Table builders for the downloadable graph/cluster exports.
 *
 * Every builder returns a plain array of flat objects, ready to be rendered by
 * `ExportableTable` and serialized to CSV. Keeping the builders pure makes the
 * two views (single analysis and history comparison) produce byte-identical
 * tables for the same inputs.
 */

const round = (value, digits = 6) => {
  if (value === null || value === undefined) return ''
  const numeric = Number(value)
  if (!Number.isFinite(numeric)) return ''
  return Number(numeric.toFixed(digits))
}

export const signLabel = (weight) => {
  const numeric = Number(weight) || 0
  if (numeric > 0) return 'positive'
  if (numeric < 0) return 'negative'
  return 'neutral'
}

const weightColumn = (edge) => (
  edge?.weight_type === 'selection_stability'
    ? 'selection_stability'
    : edge?.weight_type === 'mixed' ? 'weight' : 'correlation'
)

/** RFC-4180 CSV serialization with a union of all row keys as header. */
export function toCsv(rows) {
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

/** Every edge of a full condition graph, with its correlation value. */
export function conditionEdgeRows(graph, condition) {
  return (graph?.edges || []).map((edge) => ({
    condition,
    source: edge.source,
    target: edge.target,
    [weightColumn(edge)]: round(edge.weight),
    weight_type: edge.weight_type || 'correlation',
    sign: edge.weight_type === 'selection_stability' ? 'selected' : signLabel(edge.weight),
    p_value: round(edge.p_value),
    evidence_type: edge.evidence_type || 'pairwise_p_value',
    method: edge.method || '',
  }))
}

/** Every edge of a cluster subgraph, with its correlation value. */
export function clusterEdgeRows(subgraph, condition) {
  return (subgraph?.edges || []).map((edge) => ({
    condition,
    source: edge.source,
    target: edge.target,
    [weightColumn(edge)]: round(edge.weight),
    weight_type: edge.weight_type || 'correlation',
    sign: edge.weight_type === 'selection_stability' ? 'selected' : signLabel(edge.weight),
    p_value: round(edge.p_value),
    evidence_type: edge.evidence_type || 'pairwise_p_value',
  }))
}

/** Edges of a derived view (exclusive / complement) of one condition. */
export function derivedEdgeRows(edges, condition) {
  return (edges || []).map((edge) => ({
    condition,
    source: edge.source,
    target: edge.target,
    [weightColumn(edge)]: round(edge.weight),
    weight_type: edge.weight_type || 'correlation',
    sign: edge.weight_type === 'selection_stability' ? 'selected' : signLabel(edge.weight),
    p_value: round(edge.p_value),
    evidence_type: edge.evidence_type || 'pairwise_p_value',
  }))
}

/**
 * Intersection edges. Each row reports the correlation of the edge in *both*
 * conditions, so a pair with opposite signs keeps both values.
 */
export function intersectionEdgeRows(intersectionEdges, labelA, labelB) {
  return (intersectionEdges || []).map((edge) => {
    const relation = edge.metricMismatch
      ? `different metrics (${edge.weight_type_a} vs ${edge.weight_type_b})`
      : !edge.signConflict
      ? 'same sign'
      : edge.conflictType === 'neg_to_pos'
        ? `negative in ${labelA} → positive in ${labelB}`
        : `positive in ${labelA} → negative in ${labelB}`
    return {
      source: edge.source,
      target: edge.target,
      [`${weightColumn(edge)}_${labelA}`]: round(edge.weightA),
      [`sign_${labelA}`]: edge.signA,
      [`p_value_${labelA}`]: round(edge.pValueA),
      [`${weightColumn(edge)}_${labelB}`]: round(edge.weightB),
      [`sign_${labelB}`]: edge.signB,
      [`p_value_${labelB}`]: round(edge.pValueB),
      relation,
    }
  })
}

/**
 * Nodes belonging to exactly one of the two clusters, with the reason they are
 * missing from the other condition:
 *  - absent from the other *cluster* only (the node is still in that
 *    condition's full graph), or
 *  - absent from the other condition's *graph* entirely.
 */
export function nodePresenceRows(firstCluster, secondCluster, firstGraph, secondGraph, labelA, labelB) {
  const clusterA = new Set((firstCluster?.nodes || []).map((node) => node.id))
  const clusterB = new Set((secondCluster?.nodes || []).map((node) => node.id))
  const graphA = new Set((firstGraph?.nodes || []).map((node) => node.id))
  const graphB = new Set((secondGraph?.nodes || []).map((node) => node.id))

  const rows = []
  const allNodes = new Set([...clusterA, ...clusterB])
  for (const node of allNodes) {
    const inClusterA = clusterA.has(node)
    const inClusterB = clusterB.has(node)
    if (inClusterA === inClusterB) continue // in both clusters (or in neither): not exclusive

    const presentIn = inClusterA ? labelA : labelB

    rows.push({
      node,
      present_in: presentIn,
      [`in_cluster_${labelA}`]: inClusterA,
      [`in_cluster_${labelB}`]: inClusterB,
      [`in_graph_${labelA}`]: graphA.has(node),
      [`in_graph_${labelB}`]: graphB.has(node)
    })
  }
  return rows
}
