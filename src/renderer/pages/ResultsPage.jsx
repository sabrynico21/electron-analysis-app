import { useEffect } from 'react'
import { useParams } from 'react-router-dom'
import { useResultsStore } from '../store/resultsStore'
import ChartPanel from '../components/Results/ChartPanel'
import ResultsSummary from '../components/Results/ResultsSummary'
import DataTable from '../components/Results/DataTable'
import ExportToolbar from '../components/Results/ExportToolbar'
import styles from '../styles/ResultsPage.module.css'

export default function ResultsPage() {
  const { jobId } = useParams()
  const {
    results,
    conditionGraphs,
    reclusterByGraph,
    fetchResults,
    fetchConditionGraph,
    reclusterConditionGraph,
    isLoading,
    error,
  } = useResultsStore()

  useEffect(() => {
    fetchResults(jobId)
  }, [jobId])

  useEffect(() => {
    const refs = results?.graph_pipeline?.graph_refs || []
    refs.forEach((ref) => {
      fetchConditionGraph(ref.cacheKey)
    })
  }, [results, fetchConditionGraph])

  if (isLoading) return <div className={styles.loading}>Loading results…</div>
  if (error) return <div className={styles.empty}>{error}</div>
  if (!results) return <div className={styles.empty}>No results found.</div>

  return (
    <div className={styles.page}>
      <div className={styles.header}>
        <h1>Cluster Results</h1>
        <ExportToolbar jobId={jobId} />
      </div>

      <ResultsSummary summary={results?.summary} />

      <ChartPanel
        graphRefs={results?.graph_pipeline?.graph_refs || []}
        conditionGraphs={conditionGraphs}
        reclusterByGraph={reclusterByGraph}
        onRecluster={reclusterConditionGraph}
      />

      {Array.isArray(results?.table) && results.table.length > 0 && (
        <section className={styles.section}>
          <h2 className={styles.sectionTitle}>Correlation table</h2>
          <DataTable data={results.table} />
        </section>
      )}
    </div>
  )
}
