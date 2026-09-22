import styles from './DataTable.module.css'

export default function DataTable({ data }) {
  if (!data || data.length === 0) return null
  const columns = Object.keys(data[0])
  return (
    <div className={styles.wrapper}>
      <table className={styles.table}>
        <thead>
          <tr>{columns.map((col) => <th key={col}>{col}</th>)}</tr>
        </thead>
        <tbody>
          {data.map((row, i) => (
            <tr key={i}>
              {columns.map((col) => <td key={col}>{row[col]}</td>)}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}
