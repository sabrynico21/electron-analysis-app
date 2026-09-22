/**
 * AnalysisPage
 *
 * Entry point reserved for **new** analyses. All the actual workflow lives in
 * `AnalysisWorkspace`, which is shared with the History views.
 *
 * This page always starts from an empty form: re-entering it means "create a new
 * analysis", so leftover files/parameters from a previous session must never
 * reappear. The `key` forces a full remount (and therefore a reset form) whenever
 * the user navigates here again.
 */
import { useLocation } from 'react-router-dom'
import AnalysisWorkspace from '../components/Analysis/AnalysisWorkspace'

export default function AnalysisPage() {
  const location = useLocation()
  const resetKey = location.key || 'new-analysis'

  return <AnalysisWorkspace key={resetKey} mode='new' />
}
