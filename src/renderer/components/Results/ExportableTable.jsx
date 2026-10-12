import { useMemo, useState } from 'react'
import styles from './ExportableTable.module.css'
import { toCsv } from '../../utils/graphExports'

const PREVIEW_LIMIT = 200

const formatCell = (value) => {
  if (value === null || value === undefined) return ''
  if (typeof value === 'boolean') return value ? 'yes' : 'no'
  if (typeof value === 'number') return Number.isInteger(value) ? String(value) : value.toFixed(4)
  return String(value)
}

/**
 * A downloadable table.
 *
 * Renders nothing above the fold except a header: the rows are only shown when
 * the user clicks "Show table", while "Download CSV" always writes the complete
 * table through the native save dialog (main process).
 */
export default function ExportableTable({
  title,
  description = '',
  rows = [],
  fileName = 'table.csv',
  testId = '',
  emptyLabel = 'No rows to show.',
}) {
  const [open, setOpen] = useState(false)
  const [status, setStatus] = useState('')
  const [busy, setBusy] = useState(false)

  const columns = useMemo(
    () => [...new Set((rows || []).flatMap((row) => Object.keys(row || {})))],
    [rows]
  )
  const preview = (rows || []).slice(0, PREVIEW_LIMIT)
  const canDownload = (rows || []).length > 0

  const handleDownload = async () => {
    if (!canDownload) return
    setBusy(true)
    setStatus('')
    try {
      const content = toCsv(rows)
      const api = window.electronAPI
      if (!api?.saveTextFile) {
        setStatus('Saving is only available in the desktop application.')
        return
      }
      const result = await api.saveTextFile({ defaultFileName: fileName, content, format: 'csv' })
      if (result?.success) {
        setStatus(`Saved to ${result.path}`)
      } else if (result?.error && result.error !== 'cancelled') {
        setStatus(`Export failed: ${result.error}`)
      }
    } catch (error) {
      setStatus(`Export failed: ${error.message}`)
    } finally {
      setBusy(false)
    }
  }

  return (
    <div className={styles.card} data-testid={testId || undefined}>
      <div className={styles.header}>
        <div className={styles.headerText}>
          <span className={styles.title}>{title}</span>
          {description && <span className={styles.description}>{description}</span>}
        </div>
        <div className={styles.actions}>
          <span className={styles.count}>{rows.length} row{rows.length === 1 ? '' : 's'}</span>
          <button
            type='button'
            className={styles.button}
            onClick={() => setOpen((value) => !value)}
            disabled={!canDownload}
            data-testid={testId ? `${testId}-toggle` : undefined}
          >
            {open ? 'Hide table' : 'Show table'}
          </button>
          <button
            type='button'
            className={`${styles.button} ${styles.primary}`}
            onClick={handleDownload}
            disabled={!canDownload || busy}
            data-testid={testId ? `${testId}-download` : undefined}
          >
            {busy ? 'Saving…' : 'Download CSV'}
          </button>
        </div>
      </div>

      {status && <p className={styles.status}>{status}</p>}

      {open && canDownload && (
        <div className={styles.tableWrap}>
          <table className={styles.table}>
            <thead>
              <tr>{columns.map((column) => <th key={column}>{column}</th>)}</tr>
            </thead>
            <tbody>
              {preview.map((row, index) => (
                <tr key={index}>
                  {columns.map((column) => <td key={column}>{formatCell(row[column])}</td>)}
                </tr>
              ))}
            </tbody>
          </table>
          {rows.length > PREVIEW_LIMIT && (
            <p className={styles.note}>
              Showing the first {PREVIEW_LIMIT} of {rows.length} rows — the CSV contains all of them.
            </p>
          )}
        </div>
      )}

      {open && !canDownload && <p className={styles.note}>{emptyLabel}</p>}
    </div>
  )
}
