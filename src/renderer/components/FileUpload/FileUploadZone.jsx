import { useCallback, useRef, useState } from 'react'
import { describePaths, extractDroppedEntries, pickFilesWithDialog } from '../../utils/fileLoading'
import styles from './FileUploadZone.module.css'

/**
 * File picker used for all four upload slots.
 *
 * Both "browse" and drag & drop end up in `utils/fileLoading`, which resolves and
 * validates absolute paths on the main process, so the two entry points cannot
 * behave differently.
 *
 * Drag & drop deliberately uses the native DOM events instead of `react-dropzone`:
 * the drop handler must read `dataTransfer.files` **synchronously** to keep the
 * real OS path (see `extractDroppedEntries`). `react-dropzone` re-derives the
 * files asynchronously, which loses the path and made the app report
 * "file not found" for every dropped file.
 *
 * @param {Object} props
 * @param {(files: import('../../utils/fileLoading').FileDescriptor[]) => Promise<void>|void} props.onFilesAdded
 * @param {(index: number) => void} [props.onRemove]
 * @param {import('../../utils/fileLoading').FileDescriptor[]} [props.files]
 * @param {boolean} [props.disabled]
 * @param {string} [props.testId] stable hook used by the end-to-end tests
 */
export default function FileUploadZone({ onFilesAdded, files = [], onRemove, disabled = false, testId }) {
  const [busy, setBusy] = useState(false)
  const [isDragActive, setIsDragActive] = useState(false)
  const [error, setError] = useState('')
  const dragDepth = useRef(0)

  const ingest = useCallback(async (accepted, rejected) => {
    setError(
      rejected?.length
        ? rejected.map((item) => `${item.name || item.path}: ${item.error}`).join(' · ')
        : ''
    )
    if (accepted.length > 0) {
      await onFilesAdded(accepted)
    }
  }, [onFilesAdded])

  const handleBrowse = useCallback(async (event) => {
    event.preventDefault()
    event.stopPropagation()

    if (disabled || busy) return

    setBusy(true)
    try {
      const { canceled, accepted, rejected } = await pickFilesWithDialog()
      if (canceled) {
        setError('')
        return
      }
      await ingest(accepted, rejected)
    } finally {
      setBusy(false)
    }
  }, [busy, disabled, ingest])

  const handleDragEnter = useCallback((event) => {
    event.preventDefault()
    event.stopPropagation()
    if (disabled || busy) return
    dragDepth.current += 1
    setIsDragActive(true)
  }, [busy, disabled])

  const handleDragOver = useCallback((event) => {
    event.preventDefault()
    event.stopPropagation()
    if (disabled || busy) return
    if (event.dataTransfer) event.dataTransfer.dropEffect = 'copy'
  }, [busy, disabled])

  const handleDragLeave = useCallback((event) => {
    event.preventDefault()
    event.stopPropagation()
    dragDepth.current = Math.max(0, dragDepth.current - 1)
    if (dragDepth.current === 0) setIsDragActive(false)
  }, [])

  const handleDrop = useCallback(async (event) => {
    event.preventDefault()
    event.stopPropagation()
    dragDepth.current = 0
    setIsDragActive(false)

    if (disabled || busy) return

    // Extract the real paths BEFORE awaiting anything: the DataTransfer is only
    // valid for the duration of this handler.
    const entries = extractDroppedEntries(event.dataTransfer)

    if (entries.length === 0) {
      setError('No usable file was found in the drop. Use “click to browse” instead.')
      return
    }

    setBusy(true)
    try {
      const { accepted, rejected } = await describePaths(entries.map((entry) => entry.path))
      await ingest(accepted, rejected)
    } finally {
      setBusy(false)
    }
  }, [busy, disabled, ingest])

  return (
    <div
      className={`${styles.zone} ${isDragActive ? styles.active : ''} ${disabled ? styles.disabled : ''}`}
      data-testid={testId}
      onClick={disabled ? undefined : handleBrowse}
      onDragEnter={handleDragEnter}
      onDragOver={handleDragOver}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
      role='button'
      tabIndex={disabled ? -1 : 0}
      onKeyDown={(event) => {
        if (event.key === 'Enter' || event.key === ' ') handleBrowse(event)
      }}
    >
      <div className={styles.inner}>
        <span className={styles.icon}>📂</span>
        {isDragActive
          ? <p>Release to add files…</p>
          : <p>{busy ? 'Reading files…' : <>Drag &amp; drop files here, or <strong>click to browse</strong></>}</p>
        }

        {files.length > 0 && (
          <ul className={styles.fileList}>
            {files.map((file, idx) => (
              <li key={`${file.path || file.name}-${idx}`} className={styles.fileItem}>
                <span className={styles.fileName} title={file.path}>{file.name}</span>
                {typeof onRemove === 'function' && (
                  <button
                    type='button'
                    className={styles.remove}
                    disabled={disabled}
                    onClick={(event) => {
                      event.stopPropagation()
                      onRemove(idx)
                    }}
                    aria-label='Remove file'
                  >
                    ×
                  </button>
                )}
              </li>
            ))}
          </ul>
        )}

        {error && <p className={styles.error} data-testid='upload-error'>{error}</p>}
      </div>
    </div>
  )
}


