/**
 * Shared helpers for comparing two subgraphs.
 *
 * Used both by the single-analysis view (`ChartPanel`, which compares the two
 * *conditions* of one analysis) and by the history comparison view (which
 * compares two whole analyses), so the two views cannot drift apart.
 */

/** Stable identity for an edge: undirected node pair + weight sign. */
export function edgeKey(edge = {}) {
  const source = String(edge.source || '')
  const target = String(edge.target || '')
  const ordered = [source, target].sort()
  const weight = Number(edge.weight) || 0
  const sign = weight > 0 ? 'p' : weight < 0 ? 'n' : 'z'
  return `${ordered[0]}::${ordered[1]}::${sign}`
}

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
  const intersectionNodes = [...intersectionNodeIds].map((id) => firstNodeMap.get(id) || secondNodeMap.get(id)).filter(Boolean)

  const firstEdgeMap = new Map(firstEdges.map((edge) => [edgeKey(edge), edge]))
  const secondEdgeMap = new Map(secondEdges.map((edge) => [edgeKey(edge), edge]))
  const intersectionEdgeKeys = new Set([...firstEdgeMap.keys()].filter((key) => secondEdgeMap.has(key)))
  const intersectionEdges = [...intersectionEdgeKeys].map((key) => firstEdgeMap.get(key) || secondEdgeMap.get(key)).filter(Boolean)

  const firstExclusiveNodeIds = new Set([...firstNodeIds].filter((id) => !intersectionNodeIds.has(id)))
  const secondExclusiveNodeIds = new Set([...secondNodeIds].filter((id) => !intersectionNodeIds.has(id)))
  const firstExclusiveNodes = [...firstExclusiveNodeIds].map((id) => firstNodeMap.get(id)).filter(Boolean)
  const secondExclusiveNodes = [...secondExclusiveNodeIds].map((id) => secondNodeMap.get(id)).filter(Boolean)

  const firstExclusiveEdges = firstEdges.filter((edge) => {
    const key = edgeKey(edge)
    return !intersectionEdgeKeys.has(key) && firstExclusiveNodeIds.has(edge.source) && firstExclusiveNodeIds.has(edge.target)
  })
  const secondExclusiveEdges = secondEdges.filter((edge) => {
    const key = edgeKey(edge)
    return !intersectionEdgeKeys.has(key) && secondExclusiveNodeIds.has(edge.source) && secondExclusiveNodeIds.has(edge.target)
  })

  const unionNodeMap = new Map([...firstNodes, ...secondNodes].map((node) => [node.id, node]))
  const complementNodeIds = new Set([...unionNodeMap.keys()].filter((id) => !intersectionNodeIds.has(id)))
  const complementNodes = [...complementNodeIds].map((id) => unionNodeMap.get(id)).filter(Boolean)
  const unionEdgeMap = new Map([...firstEdges, ...secondEdges].map((edge) => [edgeKey(edge), edge]))
  const complementEdges = [...unionEdgeMap.entries()]
    .filter(([key, edge]) => !intersectionEdgeKeys.has(key) && complementNodeIds.has(edge.source) && complementNodeIds.has(edge.target))
    .map(([, edge]) => edge)

  return {
    intersection: { nodes: intersectionNodes, edges: intersectionEdges },
    exclusiveA: { nodes: firstExclusiveNodes, edges: firstExclusiveEdges },
    exclusiveB: { nodes: secondExclusiveNodes, edges: secondExclusiveEdges },
    complement: { nodes: complementNodes, edges: complementEdges },
  }
}
