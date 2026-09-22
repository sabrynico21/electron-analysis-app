import styles from './ResultsSummary.module.css'

export default function ResultsSummary({ summary }) {
  if (!summary) return null
  return (
    <div className={styles.grid}>
      {Object.entries(summary).map(([key, value]) => (
        <div key={key} className={styles.card}>
          <span className={styles.label}>{key}</span>
          <span className={styles.value}>{typeof value === 'number' ? value.toFixed(4) : String(value)}</span>
        </div>
      ))}
    </div>
  )
}
