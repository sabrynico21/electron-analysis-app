import { useState, useEffect, useCallback, useRef } from 'react'
import { useAnalysisStore } from '../store/analysisStore'

/**
 * Drives a single analysis run and keeps the history archive in sync.
 *
 * `runAnalysis` accepts an optional `jobId` inside the payload: when present the
 * run overwrites that existing history entry instead of creating a new one, which
 * is what makes a stored analysis fully editable and re-runnable.
 */
export function useAnalysis() {
  const [status, setStatus] = useState('idle') // idle | running | completed | failed
  const [progress, setProgress] = useState(0)
  const [logs, setLogs] = useState([])
  const [jobId, setJobId] = useState(null)
  const [error, setError] = useState(null)
  const mountedRef = useRef(true)

  const addToHistory = useAnalysisStore((state) => state.addToHistory)

  useEffect(() => {
    mountedRef.current = true

    window.electronAPI.onAnalysisProgress((data) => {
      if (!mountedRef.current) return
      setStatus(data.status)
      setProgress(data.percent ?? 0)
    })
    window.electronAPI.onAnalysisLog((data) => {
      if (!mountedRef.current) return
      setLogs((prev) => [...prev, data])
    })

    return () => {
      mountedRef.current = false
      window.electronAPI.removeAllListeners('analysis:progress')
      window.electronAPI.removeAllListeners('analysis:log')
    }
  }, [])

  /** Clear the run state so a new or freshly opened analysis starts clean. */
  const resetRun = useCallback(() => {
    setStatus('idle')
    setProgress(0)
    setLogs([])
    setJobId(null)
    setError(null)
  }, [])

  const runAnalysis = useCallback(async (payload) => {
    setStatus('running')
    setProgress(0)
    setLogs([])
    setError(null)

    const startedAt = Date.now()
    const result = await window.electronAPI.runAnalysis(payload)

    if (result?.success) {
      setJobId(result.jobId)
      setStatus('completed')
      await addToHistory({
        id: result.jobId,
        scriptName: payload.scriptName || 'analysis.py',
        analysisName: payload.analysisName || '',
        status: 'completed',
        createdAt: payload.createdAt || startedAt,
        updatedAt: Date.now(),
        payloadSnapshot: {
          files: payload.files || {},
          params: payload.params || {},
        },
        hasResults: true,
      })
      return result
    }

    setStatus('failed')
    setError(result?.error || 'The analysis failed')
    return result || { success: false, error: 'The analysis failed' }
  }, [addToHistory])

  const cancelAnalysis = useCallback(async (cancelJobId) => {
    if (!cancelJobId) return
    await window.electronAPI.cancelAnalysis(cancelJobId)
  }, [])

  return { runAnalysis, cancelAnalysis, resetRun, status, progress, logs, jobId, error }
}
