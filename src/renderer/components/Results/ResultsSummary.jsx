import styles from './ResultsSummary.module.css'

const SUMMARY_LABELS = {
  correlation_method: 'Correlation method',
  evidence_type: 'Edge evidence',
  method_settings: 'Model settings',
  p_value_threshold: 'P-value threshold',
  p_value_threshold_applied: 'Pairwise p-value filter applied',
  weight_type: 'Edge weight type',
}

const formatValue = (value) => {
  if (value === null || value === undefined) return 'Not applicable'
  if (typeof value === 'number') return value.toFixed(4)
  if (typeof value === 'boolean') return value ? 'Yes' : 'No'
  return String(value)
}

export default function ResultsSummary({ summary }) {
  if (!summary) return null
  return (
    <div className={styles.grid}>
      {Object.entries(summary).map(([key, value]) => (
        <div key={key} className={styles.card}>
          <span className={styles.label}>{SUMMARY_LABELS[key] || key}</span>
          <span className={styles.value}>{formatValue(value)}</span>
        </div>
      ))}
    </div>
  )
}
