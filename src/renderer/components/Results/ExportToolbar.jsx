import { useState } from 'react'
import styles from './ExportToolbar.module.css'

const FORMATS = [
  { id: 'csv', label: 'CSV' },
  { id: 'json', label: 'JSON' },
]

export default function ExportToolbar({ jobId, disabled = false }) {
  const [status, setStatus] = useState('')
  const [busy, setBusy] = useState(false)

  const handleExport = async (format) => {
    setBusy(true)
    setStatus('')
    try {
      const result = await window.electronAPI.exportResults(jobId, format)
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
    <div className={styles.toolbar}>
      {FORMATS.map((fmt) => (
        <button
          key={fmt.id}
          type="button"
          className={styles.btn}
          disabled={disabled || busy || !jobId}
          onClick={() => handleExport(fmt.id)}
        >
          Export {fmt.label}
        </button>
      ))}
      {status && <span className={styles.status}>{status}</span>}
    </div>
  )
}
