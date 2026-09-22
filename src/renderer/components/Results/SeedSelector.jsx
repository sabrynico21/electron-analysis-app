import { useMemo } from 'react'
import styles from './SeedSelector.module.css'

function formatTaxaLabel(value = '') {
  const source = String(value || '')
  const idx = source.indexOf('g_')
  if (idx < 0) return source
  const compact = source.slice(idx + 2).trim()
  return compact || source
}

export default function SeedSelector({
  nodes = [],
  searchValue,
  selectedTaxa,
  onSearchChange,
  onTaxaSelect,
  inputId = 'taxa-search',
  disabled = false,
}) {
  const sortedNodes = useMemo(() => {
    return [...nodes].sort((left, right) => {
      const leftName = String(left.label || left.id || '').toLowerCase()
      const rightName = String(right.label || right.id || '').toLowerCase()
      return leftName.localeCompare(rightName)
    })
  }, [nodes])

  const filteredNodes = useMemo(() => {
    const term = (searchValue || '').trim().toLowerCase()
    if (!term) return sortedNodes
    return sortedNodes
      .filter((node) => {
        const displayName = String(node.label || node.id || '').toLowerCase()
        const rawName = (node.label || node.id || '').toLowerCase()
        return displayName.includes(term) || rawName.includes(term)
      })
  }, [sortedNodes, searchValue])

  const selectedNode = useMemo(
    () => sortedNodes.find((node) => (node.id || node.label) === selectedTaxa),
    [sortedNodes, selectedTaxa],
  )

  const selectedLabel = selectedNode
    ? formatTaxaLabel(selectedNode.label || selectedNode.id)
    : (selectedTaxa ? formatTaxaLabel(selectedTaxa) : 'None')

  return (
    <div className={styles.wrapper}>
      <div className={styles.topBar}>
        <div className={styles.searchBlock}>
          <label htmlFor={inputId} className={styles.label}></label>
          <input
            id={inputId}
            className={styles.input}
            type='text'
            value={searchValue}
            disabled={disabled}
            onChange={(event) => onSearchChange(event.target.value)}
            placeholder='Type to search taxa'
          />
        </div>

        <p className={styles.selectedText}>Selected taxa: <strong>{selectedLabel}</strong></p>

        <div className={styles.legend} aria-label='Taxa type legend'>
          <span className={styles.legendItem}>
            <i className={`${styles.legendSwatch} ${styles.legendBacteria}`} aria-hidden='true' />
            Bacteria
          </span>
          <span className={styles.legendItem}>
            <i className={`${styles.legendSwatch} ${styles.legendFungi}`} aria-hidden='true' />
            Fungi
          </span>
        </div>
      </div>

      <p className={styles.meta}>
        {!searchValue?.trim()
          ? `${filteredNodes.length} taxa shown`
          : `${filteredNodes.length} taxa match your search`}
      </p>

      <div className={styles.list}>
        {filteredNodes.map((node) => {
          const nodeId = node.id || node.label
          const active = selectedTaxa === nodeId
          const kindClass = node.kind === 'fungi' ? styles.itemFungi : styles.itemBacteria
          return (
            <button
              key={nodeId}
              type='button'
              className={active ? `${styles.item} ${kindClass} ${styles.active}` : `${styles.item} ${kindClass}`}
              disabled={disabled}
              onClick={() => onTaxaSelect(nodeId)}
            >
              <span>{node.label || nodeId}</span>
            </button>
          )
        })}
        {filteredNodes.length === 0 && <p className={styles.empty}>No matching nodes found.</p>}
      </div>
    </div>
  )
}
