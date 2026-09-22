/**
 * HistoryPage
 *
 * Two-pane layout: a deliberately narrow list on the left and the analysis
 * content on the right.
 *
 * Opening an entry never sends the user back to the "New Analysis" workspace:
 * the stored analysis is rendered in place (fully editable) with a "Back to
 * history" button at the top. Two entries can also be compared side by side,
 * each highlighted with its own colour tint.
 */
import { useEffect, useMemo } from 'react'
import { useNavigate, useParams, useSearchParams } from 'react-router-dom'
import AnalysisWorkspace from '../components/Analysis/AnalysisWorkspace'
import ComparisonView from '../components/History/ComparisonView'
import { useAnalysisStore } from '../store/analysisStore'
import styles from '../styles/HistoryPage.module.css'

const formatDateTime = (value) => (value ? new Date(value).toLocaleString() : '—')

const formatRelative = (value) => {
  if (!value) return ''
  const diffMs = Date.now() - value
  const minutes = Math.round(diffMs / 60000)
  if (minutes < 1) return 'just now'
  if (minutes < 60) return `${minutes} min ago`
  const hours = Math.round(minutes / 60)
  if (hours < 24) return `${hours} h ago`
  const days = Math.round(hours / 24)
  if (days < 30) return `${days} d ago`
  return ''
}

export default function HistoryPage() {
  const navigate = useNavigate()
  const { jobId = '' } = useParams()
  const [searchParams, setSearchParams] = useSearchParams()

  const { history, loadHistory, removeFromHistory, rebuildHistory, isLoading, error } = useAnalysisStore()

  useEffect(() => {
    loadHistory()
  }, [loadHistory])

  const selectedEntry = useMemo(
    () => history.find((entry) => entry.id === jobId) || null,
    [history, jobId]
  )

  // `?compare=<idA>,<idB>` drives the side-by-side view.
  const compare = useMemo(() => {
    const raw = searchParams.get('compare') || ''
    const [a, b] = raw.split(',').map((value) => value.trim())
    return { a: a || '', b: b || '' }
  }, [searchParams])

  const compareEntries = useMemo(() => {
    const find = (id) => history.find((entry) => entry.id === id) || null
    return { a: find(compare.a), b: find(compare.b) }
  }, [compare.a, compare.b, history])

  const isComparing = Boolean(compareEntries.a && compareEntries.b)

  const openEntry = (id) => {
    setSearchParams({}, { replace: true })
    navigate(`/history/${id}`)
  }

  const startCompare = (id) => {
    // Fill slot A first, then B; a third click restarts the selection.
    let next
    if (!compare.a || (compare.a && compare.b)) next = { a: id, b: '' }
    else if (compare.a === id) next = { a: '', b: '' }
    else next = { a: compare.a, b: id }

    const params = {}
    const value = [next.a, next.b].filter(Boolean).join(',')
    if (value) params.compare = value
    setSearchParams(params, { replace: true })
  }

  const handleDelete = async (entry) => {
    const label = entry.analysisName || entry.scriptName || entry.id
    const ok = window.confirm(`Delete analysis "${label}"? Its results and cached graphs will be removed too.`)
    if (!ok) return

    await removeFromHistory(entry.id)
    if (entry.id === jobId) navigate('/history', { replace: true })
    if (entry.id === compare.a || entry.id === compare.b) setSearchParams({}, { replace: true })
  }

  const handleRebuild = async () => {
    const response = await rebuildHistory()
    if (response?.success) {
      window.alert(`Recovered ${response.recovered ?? 0} analyses from the results stored on disk.`)
    }
  }

  const compareHint = !compare.a
    ? ''
    : !compare.b
      ? `Select a second analysis to compare with "${compareEntries.a?.analysisName || ''}".`
      : ''

  return (
    <div className={`${styles.page} ${jobId || isComparing ? styles.detailOpen : ''}`}>
      <aside className={styles.listPane}>
        <header className={styles.listHeader}>
          <div>
            <h1 className={styles.listTitle}>History</h1>
            <p className={styles.listSubtitle}>
              {history.length} saved analys{history.length === 1 ? 'is' : 'es'} · stored in the app data folder
            </p>
          </div>
          <button type='button' className={styles.iconButton} onClick={handleRebuild} title='Recover entries from the results on disk'>
            ⟳
          </button>
        </header>

        {isLoading && <p className={styles.empty}>Loading…</p>}
        {!isLoading && error && <p className={styles.errorText}>{error}</p>}

        {!isLoading && !error && history.length === 0 && (
          <div className={styles.emptyBox}>
            <p className={styles.empty}>No analyses stored yet.</p>
            <p className={styles.emptyHint}>
              Analyses run from the “Create new analysis” page are archived here, outside the browser
              storage, so they survive updates and reinstalls.
            </p>
            <button type='button' className={styles.secondaryButton} onClick={handleRebuild}>
              Look for results on disk
            </button>
          </div>
        )}

        <ul className={styles.list} data-testid='history-list'>
          {history.map((entry) => {
            const tint = entry.id === compare.a ? 'tintA' : entry.id === compare.b ? 'tintB' : ''
            const isActive = entry.id === jobId && !isComparing
            return (
              <li
                key={entry.id}
                data-testid='history-item'
                data-job-id={entry.id}
                className={`${styles.item} ${isActive ? styles.itemActive : ''} ${tint ? styles[tint] : ''}`}
              >
                <button type='button' className={styles.itemBody} onClick={() => openEntry(entry.id)}>
                  <span className={styles.itemName}>{entry.analysisName || entry.scriptName}</span>
                  <span className={styles.itemMeta}>
                    <span>Created {formatDateTime(entry.createdAt)}</span>
                    <span>
                      Updated {formatDateTime(entry.updatedAt)}
                      {formatRelative(entry.updatedAt) ? ` · ${formatRelative(entry.updatedAt)}` : ''}
                    </span>
                  </span>
                  <span className={styles.itemFooter}>
                    <span className={`${styles.status} ${styles[entry.status] || ''}`}>{entry.status}</span>
                    <span className={styles.itemParams}>
                      {entry.payloadSnapshot?.params?.conditionA && entry.payloadSnapshot?.params?.conditionB
                        ? `${entry.payloadSnapshot.params.conditionA} vs ${entry.payloadSnapshot.params.conditionB}`
                        : 'no conditions'}
                    </span>
                  </span>
                </button>
                <div className={styles.itemActions}>
                  <button
                    type='button'
                    className={styles.actionButton}
                    onClick={() => startCompare(entry.id)}
                    title='Add this analysis to the comparison'
                  >
                    {entry.id === compare.a || entry.id === compare.b ? '− Compare' : '+ Compare'}
                  </button>
                  <button
                    type='button'
                    className={styles.deleteButton}
                    onClick={() => handleDelete(entry)}
                    title='Delete this analysis'
                  >
                    Delete
                  </button>
                </div>
              </li>
            )
          })}
        </ul>
      </aside>

      <section className={styles.detailPane} data-testid='history-detail'>
        {isComparing ? (
          <ComparisonView
            entryA={compareEntries.a}
            entryB={compareEntries.b}
            onClose={() => setSearchParams({}, { replace: true })}
            onOpen={(id) => openEntry(id)}
          />
        ) : selectedEntry ? (
          <AnalysisWorkspace
            key={selectedEntry.id}
            mode='history'
            jobId={selectedEntry.id}
            entry={selectedEntry}
            onBack={() => navigate('/history')}
            onEntryChanged={loadHistory}
          />
        ) : (
          <div className={styles.placeholder}>
            <h2 className={styles.placeholderTitle}>
              {compareHint || 'Select an analysis'}
            </h2>
            <p className={styles.placeholderText}>
              Pick an analysis from the list to reopen it here — files, parameters and results are
              fully editable. Use <strong>+ Compare</strong> on two entries to view them side by side.
            </p>
          </div>
        )}
      </section>
    </div>
  )
}
