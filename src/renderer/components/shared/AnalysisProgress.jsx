import { useRef, useEffect, useState, useMemo } from 'react'
import styles from './AnalysisProgress.module.css'

const isUsefulLog = (line = '') => {
  const text = String(line).trim()
  if (!text) return false
  if (text.startsWith('{') && text.endsWith('}')) return false
  const lowered = text.toLowerCase()
  if (lowered.startsWith('debug')) return false
  if (lowered.includes('matplotlib') || lowered.includes('font manager')) return false
  return true
}

const toReadableLog = (entry = {}) => {
  const raw = String(entry.line || '').trim()
  const cleaned = raw
    .replace(/^\[stderr\]\s*/i, '')
    .replace(/^info\s*:\s*/i, '')
    .replace(/^warning\s*:\s*/i, '')
    .replace(/^warn\s*:\s*/i, '')
    .replace(/^error\s*:\s*/i, '')
    .trim()

  const lowered = cleaned.toLowerCase()
  let level = entry.level || 'info'
  if (lowered.includes('error') || lowered.includes('failed') || lowered.includes('traceback')) level = 'error'
  else if (lowered.includes('warning') || lowered.includes('warn')) level = 'warn'

  return { level, line: cleaned }
}

export default function AnalysisProgress({ progress, logs }) {
  const [expanded, setExpanded] = useState(false)
  const logRef = useRef(null)
  const expandedLogRef = useRef(null)

  // Full log, kept from start to finish, shown in the separate window.
  const allLogs = useMemo(
    () => logs.filter((entry) => isUsefulLog(entry?.line)).map((entry) => toReadableLog(entry)),
    [logs]
  )
  // Compact preview in the inline box (last 120 lines only).
  const displayLogs = useMemo(() => allLogs.slice(-120), [allLogs])

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight
  }, [displayLogs])

  useEffect(() => {
    if (expanded && expandedLogRef.current) {
      expandedLogRef.current.scrollTop = expandedLogRef.current.scrollHeight
    }
  }, [expanded, allLogs.length])

  useEffect(() => {
    if (!expanded) return
    const onKey = (e) => {
      if (e.key === 'Escape') setExpanded(false)
    }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [expanded])

  return (
    <div className={styles.container}>
      <div className={styles.header}>
        <span className={styles.pct}>{progress}%</span>
        <button
          type="button"
          className={styles.expandBtn}
          onClick={() => setExpanded(true)}
          title="Apri il log completo in una finestra separata"
        >
          <span aria-hidden="true">⛶</span> Log completo
        </button>
      </div>
      <div className={styles.barTrack}>
        <div className={styles.barFill} style={{ width: `${progress}%` }} />
      </div>
      <div ref={logRef} className={styles.logBox}>
        {displayLogs.map((entry, i) => (
          <div key={i} className={`${styles.logLine} ${styles[entry.level || 'info']}`}>
            {entry.line}
          </div>
        ))}
        {displayLogs.length === 0 && <div className={styles.logLine}>Preparing analysis runtime...</div>}
      </div>

      {expanded && (
        <div
          className={styles.overlay}
          onClick={(e) => {
            if (e.target === e.currentTarget) setExpanded(false)
          }}
        >
          <div className={styles.modal} role="dialog" aria-modal="true" aria-label="Log di analisi completo">
            <div className={styles.modalHeader}>
              <span className={styles.modalTitle}>Log di analisi — {allLogs.length} righe</span>
              <button
                type="button"
                className={styles.closeBtn}
                onClick={() => setExpanded(false)}
                title="Chiudi (Esc)"
                aria-label="Chiudi log"
              >
                ✕
              </button>
            </div>
            <div ref={expandedLogRef} className={styles.modalLog}>
              {allLogs.map((entry, i) => (
                <div key={i} className={`${styles.logLine} ${styles[entry.level || 'info']}`}>
                  {entry.line}
                </div>
              ))}
              {allLogs.length === 0 && <div className={styles.logLine}>Nessun log disponibile.</div>}
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
