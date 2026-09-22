import styles from './FileList.module.css'

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1048576) return `${(bytes / 1024).toFixed(1)} KB`
  return `${(bytes / 1048576).toFixed(1)} MB`
}

export default function FileList({ files, onRemove }) {
  return (
    <ul className={styles.list}>
      {files.map((file, idx) => (
        <li key={idx} className={styles.item}>
          <span className={styles.name}>{file.name}</span>
          <span className={styles.size}>{formatBytes(file.size)}</span>
          <button className={styles.remove} onClick={() => onRemove(idx)} aria-label="Remove">✕</button>
        </li>
      ))}
    </ul>
  )
}
