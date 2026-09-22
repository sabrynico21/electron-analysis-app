import { useEffect, useMemo, useRef, useState } from 'react'
import styles from './SubgraphPlot.module.css'
import { GROUPING_MODE, getGroupingNodeKey, getNodeDisplayName } from '../../utils/taxonomyGrouping'

const WIDTH = 700
const HEIGHT = 440
const GOLDEN_ANGLE = Math.PI * (3 - Math.sqrt(5))

function computeLayout(nodes = []) {
  if (!nodes.length) return []
  const centerX = WIDTH / 2
  const centerY = HEIGHT / 2
  const maxRadius = Math.min(WIDTH, HEIGHT) * 0.42
  const safeCount = Math.max(1, nodes.length)

  return nodes.map((node, idx) => {
    const t = Math.sqrt((idx + 0.5) / safeCount)
    const radius = t * maxRadius
    const angle = idx * GOLDEN_ANGLE
    return {
      ...node,
      x: centerX + radius * Math.cos(angle),
      y: centerY + radius * Math.sin(angle),
    }
  })
}

function buildVisibleLabelSet(laidOutNodes = []) {
  const visible = new Set()
  const minLabelDistance = 18

  laidOutNodes.forEach((node, idx) => {
    const canShow = laidOutNodes
      .slice(0, idx)
      .every((prev) => {
        const dx = node.x - prev.x
        const dy = node.y - prev.y
        return Math.sqrt((dx * dx) + (dy * dy)) > minLabelDistance
      })
    if (canShow) visible.add(node.id)
  })

  return visible
}

function getLabelAfterPrefix(value = '') {
  if (typeof value !== 'string') return ''
  const prefixIndex = value.indexOf('g_')
  if (prefixIndex < 0) return value
  const trimmed = value.slice(prefixIndex + 2).trim()
  return trimmed || value
}

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value))
}

function toSvgPoint(svgElement, clientX, clientY) {
  const rect = svgElement.getBoundingClientRect()
  if (!rect.width || !rect.height) return { x: 0, y: 0 }
  return {
    x: ((clientX - rect.left) / rect.width) * WIDTH,
    y: ((clientY - rect.top) / rect.height) * HEIGHT,
  }
}

export default function SubgraphPlot({
  nodes = [],
  edges = [],
  maxNodes = 40,
  groupingMode = GROUPING_MODE.ALL,
  seedNode = '',
  controlsDisabled = false,
}) {
  const svgRef = useRef(null)
  const [zoom, setZoom] = useState(1)
  const [pan, setPan] = useState({ x: 0, y: 0 })
  const [draggingNode, setDraggingNode] = useState(null)
  const [manualPositions, setManualPositions] = useState({})
  const [panning, setPanning] = useState(null)

  const grouped = useMemo(() => {
    if (groupingMode === GROUPING_MODE.ALL) {
      return {
        nodes,
        edges,
      }
    }

    const groupMap = new Map()
    nodes.forEach((node) => {
      const groupId = getGroupingNodeKey(node, groupingMode)
      if (!groupMap.has(groupId)) {
        groupMap.set(groupId, {
          id: groupId,
          label: getNodeDisplayName(node, groupingMode),
          kind: node.kind,
          originalCount: 0,
        })
      }
      const current = groupMap.get(groupId)
      current.originalCount += 1
    })

    const groupedNodes = Array.from(groupMap.values())
    const idToGroupId = new Map(nodes.map((node) => [node.id, getGroupingNodeKey(node, groupingMode)]))

    const groupedEdgeMap = new Map()
    edges.forEach((edge) => {
      const sourceGroup = idToGroupId.get(edge.source)
      const targetGroup = idToGroupId.get(edge.target)
      if (!sourceGroup || !targetGroup || sourceGroup === targetGroup) return

      const pairKey = [sourceGroup, targetGroup].sort().join('::')
      const numericWeight = Math.abs(Number(edge.weight) || 1)
      const existing = groupedEdgeMap.get(pairKey)
      if (existing) {
        existing.weight += numericWeight
      } else {
        groupedEdgeMap.set(pairKey, {
          source: sourceGroup,
          target: targetGroup,
          weight: numericWeight,
        })
      }
    })

    return {
      nodes: groupedNodes,
      edges: Array.from(groupedEdgeMap.values()),
    }
  }, [nodes, edges, groupingMode])

  const limitedNodes = useMemo(() => {
    const safeLimit = Math.max(1, Number(maxNodes) || 1)
    return grouped.nodes.slice(0, safeLimit)
  }, [grouped, maxNodes])

  const limitedNodeIds = useMemo(() => new Set(limitedNodes.map((node) => node.id)), [limitedNodes])

  const limitedEdges = useMemo(
    () => grouped.edges.filter((edge) => limitedNodeIds.has(edge.source) && limitedNodeIds.has(edge.target)),
    [grouped, limitedNodeIds],
  )

  const laidOutNodes = useMemo(() => computeLayout(limitedNodes), [limitedNodes])
  const finalNodes = useMemo(
    () => laidOutNodes.map((node) => ({ ...node, ...(manualPositions[node.id] || {}) })),
    [laidOutNodes, manualPositions],
  )

  const nodeMap = useMemo(() => {
    const map = new Map()
    finalNodes.forEach((node) => map.set(node.id, node))
    return map
  }, [finalNodes])

  const visibleLabels = useMemo(() => buildVisibleLabelSet(finalNodes), [finalNodes])

  const onWheelZoom = (event) => {
    event.preventDefault()
    event.stopPropagation()
    const step = event.deltaY > 0 ? -0.1 : 0.1
    setZoom((prev) => clamp(Number((prev + step).toFixed(2)), 0.5, 3.5))
  }

  useEffect(() => {
    const svgElement = svgRef.current
    if (!svgElement) return undefined

    const handleWheel = (event) => onWheelZoom(event)
    svgElement.addEventListener('wheel', handleWheel, { passive: false })

    return () => {
      svgElement.removeEventListener('wheel', handleWheel)
    }
  }, [])

  const startPan = (event) => {
    if (event.button !== 0 || !svgRef.current) return
    const point = toSvgPoint(svgRef.current, event.clientX, event.clientY)
    setPanning({ pointerX: point.x, pointerY: point.y, startPanX: pan.x, startPanY: pan.y })
  }

  const stopInteractions = () => {
    setDraggingNode(null)
    setPanning(null)
  }

  const onMouseMove = (event) => {
    if (!svgRef.current) return
    const point = toSvgPoint(svgRef.current, event.clientX, event.clientY)

    if (draggingNode) {
      const graphX = (point.x - pan.x) / zoom
      const graphY = (point.y - pan.y) / zoom
      setManualPositions((prev) => ({
        ...prev,
        [draggingNode]: { x: graphX, y: graphY },
      }))
      return
    }

    if (panning) {
      const dx = point.x - panning.pointerX
      const dy = point.y - panning.pointerY
      setPan({
        x: panning.startPanX + dx,
        y: panning.startPanY + dy,
      })
    }
  }

  const beginNodeDrag = (event, nodeId) => {
    event.stopPropagation()
    setDraggingNode(nodeId)
  }

  const resetView = () => {
    setZoom(1)
    setPan({ x: 0, y: 0 })
    setManualPositions({})
  }

  if (!limitedNodes.length) {
    return <p className={styles.empty}>No subgraph to plot yet.</p>
  }

  return (
    <div className={styles.wrapper}>
      <div className={styles.toolbar}>
        <button type='button' className={styles.toolButton} disabled={controlsDisabled} onClick={() => setZoom((prev) => clamp(Number((prev + 0.15).toFixed(2)), 0.5, 3.5))}>+</button>
        <button type='button' className={styles.toolButton} disabled={controlsDisabled} onClick={() => setZoom((prev) => clamp(Number((prev - 0.15).toFixed(2)), 0.5, 3.5))}>-</button>
        <button type='button' className={styles.toolButton} disabled={controlsDisabled} onClick={resetView}>Reset</button>
        <span className={styles.zoomInfo}>Zoom {Math.round(zoom * 100)}%</span>
      </div>
      <svg
        ref={svgRef}
        viewBox={`0 0 ${WIDTH} ${HEIGHT}`}
        className={styles.svg}
        role='img'
        aria-label='Subgraph plot'
        onWheelCapture={(event) => {
          event.preventDefault()
          event.stopPropagation()
        }}
        onMouseDown={startPan}
        onMouseMove={onMouseMove}
        onMouseUp={stopInteractions}
        onMouseLeave={stopInteractions}
      >
        <rect x='0' y='0' width={WIDTH} height={HEIGHT} className={styles.bg} />
        <g transform={`translate(${pan.x} ${pan.y}) scale(${zoom})`}>

          {limitedEdges.map((edge, idx) => {
            const source = nodeMap.get(edge.source)
            const target = nodeMap.get(edge.target)
            if (!source || !target) return null
            const numericWeight = Number(edge.weight) || 0
            const weight = Math.max(2.5, Math.min(5.5, Math.abs(numericWeight || 1) * 3))
            const edgeClassName = numericWeight > 0
              ? styles.edgePositive
              : numericWeight < 0
                ? styles.edgeNegative
                : styles.edgeNeutral
            return (
              <line
                key={`${edge.source}-${edge.target}-${idx}`}
                x1={source.x}
                y1={source.y}
                x2={target.x}
                y2={target.y}
                strokeWidth={weight}
                className={`${styles.edge} ${edgeClassName}`}
              />
            )
          })}

          {finalNodes.map((node) => {
            const isSeedNode = Boolean(seedNode) && node.id === seedNode
            const baseRadius = 8
            const nodeRadius = isSeedNode ? baseRadius + 4 : baseRadius

            return (
            <g key={node.id}>
              <circle
                cx={node.x}
                cy={node.y}
                r={nodeRadius}
                className={`${node.kind === 'bacteria' ? styles.bacteria : styles.fungi} ${isSeedNode ? styles.seedNode : ''}`}
                onMouseDown={(event) => beginNodeDrag(event, node.id)}
              />
              <title>{getLabelAfterPrefix(node.label || node.id)}</title>
              {visibleLabels.has(node.id) && (
                <text x={node.x + (isSeedNode ? 14 : 10)} y={node.y + 4} className={`${styles.label} ${isSeedNode ? styles.seedLabel : ''}`}>
                  {getLabelAfterPrefix(node.label || node.id)}
                </text>
              )}
            </g>
          )})}
        </g>
      </svg>
    </div>
  )
}
