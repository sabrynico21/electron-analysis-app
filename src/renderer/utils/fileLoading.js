/**
 * Unified file-acquisition pipeline.
 *
 * Both entry points — the native "browse" dialog and drag & drop — funnel through
 * this module so their behaviour can never diverge again.
 *
 * Background: after the upgrade to Electron 44 drag & drop silently produced an
 * empty path, because `File.path` (the old way of reading a dropped file's
 * location) was removed in Electron 32. Files therefore reached the analysis
 * pipeline with `path === ''` and the app reported "file not found". The
 * supported replacement is `webUtils.getPathForFile`, exposed by the preload
 * bridge as `window.electronAPI.getPathForFile`.
 *
 * Every descriptor returned here is validated on the main process (`fs:describeFiles`),
 * so consumers can trust `exists`/`isFile` without extra stat calls.
 */

/** Column of the metadata file that defines the experimental conditions. */
export const METADATA_CONDITION_COLUMN = 'Combination_treat2'

/**
 * @typedef {Object} FileDescriptor
 * @property {string}  path        Absolute path on disk.
 * @property {string}  name        Base name, shown in the UI.
 * @property {number}  size        Size in bytes (0 when unknown).
 * @property {boolean} exists      Whether the path could be stat'ed.
 * @property {boolean} isFile      Whether the path points to a regular file.
 * @property {number|null} modifiedAt
 * @property {string=} error       Reason the file was rejected, when applicable.
 */

/** Detect the column separator of a delimited text file. */
export function detectDelimiter(headerLine = '') {
  if (headerLine.includes('\t')) return '\t'
  if (headerLine.includes(';')) return ';'
  return ','
}

/** Split a text file into trimmed, non-empty lines. */
function toRows(content = '') {
  return String(content)
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter(Boolean)
}

/**
 * Read a metadata file and return the distinct conditions it defines.
 * @param {string} filePath
 * @returns {Promise<string[]>}
 */
export async function parseMetadataConditions(filePath) {
  if (!filePath) {
    throw new Error('The file path is not available. Re-add the file with drag & drop or Browse.')
  }

  const readResult = await window.electronAPI.readFile(filePath)
  if (!readResult?.success) {
    throw new Error(readResult?.error || 'Unable to read the metadata file')
  }

  const rows = toRows(readResult.content)
  if (rows.length < 2) {
    throw new Error('The metadata file is empty or contains no data rows')
  }

  const delimiter = detectDelimiter(rows[0])
  const header = rows[0].split(delimiter).map((value) => value.trim().replace(/^["']|["']$/g, ''))
  const columnIndex = header.indexOf(METADATA_CONDITION_COLUMN)

  if (columnIndex < 0) {
    throw new Error(`The metadata file must contain the "${METADATA_CONDITION_COLUMN}" column`)
  }

  const conditions = [
    ...new Set(
      rows
        .slice(1)
        .map((row) => row.split(delimiter)[columnIndex]?.trim().replace(/^["']|["']$/g, ''))
        .filter(Boolean)
    ),
  ]

  if (conditions.length < 2) {
    throw new Error(`At least 2 distinct values are required in "${METADATA_CONDITION_COLUMN}"`)
  }

  return conditions
}

/**
 * Turn a drop event's `DataTransfer` into path entries.
 *
 * **This must run synchronously inside the drop handler.** A `DataTransfer` is
 * only valid while the drop event is being dispatched: anything that reads it
 * later (which is what `react-dropzone` does, handing over re-derived `File`
 * objects) loses the backing OS path, `webUtils.getPathForFile` silently falls
 * back to a bogus relative `./<name>` path and the file is reported as missing.
 *
 * @param {DataTransfer|null} dataTransfer
 * @returns {{path: string, name: string, size: number}[]}
 */
export function extractDroppedEntries(dataTransfer) {
  if (!dataTransfer) return []

  const files = Array.from(dataTransfer.files || [])
  if (files.length === 0) return []

  const entries = []
  const seen = new Set()
  for (const file of files) {
    const resolvedPath = window.electronAPI?.getPathForFile?.(file) || ''
    if (!resolvedPath || seen.has(resolvedPath)) continue
    seen.add(resolvedPath)
    entries.push({ path: resolvedPath, name: file.name, size: file.size || 0 })
  }
  return entries
}

/**
 * Open the native file picker and return validated descriptors.
 * @returns {Promise<{canceled: boolean, accepted: FileDescriptor[], rejected: FileDescriptor[]}>}
 */
export async function pickFilesWithDialog() {
  if (!window.electronAPI?.openFileDialog) {
    return { canceled: true, accepted: [], rejected: [] }
  }

  const result = await window.electronAPI.openFileDialog({
    properties: ['openFile', 'multiSelections'],
  })

  if (result?.canceled || !Array.isArray(result?.filePaths) || result.filePaths.length === 0) {
    return { canceled: true, accepted: [], rejected: [] }
  }

  // The main process already returns validated descriptors.
  if (Array.isArray(result.files) && result.files.length > 0) {
    return splitValid(result.files)
  }

  const { accepted, rejected } = await describePaths(result.filePaths)
  return { canceled: false, accepted, rejected }
}

/**
 * Validate raw absolute paths through the main process.
 * @param {string[]} paths
 * @returns {Promise<{accepted: FileDescriptor[], rejected: FileDescriptor[]}>}
 */
export async function describePaths(paths = []) {
  const unique = [...new Set(paths.filter(Boolean))]
  if (unique.length === 0) return { accepted: [], rejected: [] }

  if (!window.electronAPI?.describeFiles) {
    // Defensive fallback: keep the raw path usable so the analysis can still run.
    return {
      accepted: unique.map((filePath) => ({
        path: filePath,
        name: filePath.split(/[/\\]/).pop() || filePath,
        size: 0,
        exists: true,
        isFile: true,
        modifiedAt: null,
      })),
      rejected: [],
    }
  }

  const described = await window.electronAPI.describeFiles(unique)
  return splitValid(described)
}

function splitValid(descriptors = []) {
  const accepted = []
  const rejected = []
  for (const descriptor of descriptors) {
    if (descriptor?.exists && descriptor?.isFile) {
      accepted.push(descriptor)
    } else {
      rejected.push({
        ...descriptor,
        error: descriptor?.error
          || (descriptor?.isDirectory
            ? 'Folders are not supported, drop the files inside it'
            : 'The file no longer exists at this path'),
      })
    }
  }
  return { accepted, rejected }
}

/** Remove duplicates (same path) while keeping the first occurrence. */
export function dedupeDescriptors(existing = [], incoming = []) {
  const seen = new Set(existing.map((file) => file.path))
  const fresh = []
  for (const file of incoming) {
    if (seen.has(file.path)) continue
    seen.add(file.path)
    fresh.push(file)
  }
  return { merged: [...existing, ...fresh], added: fresh }
}

/** Path list for the analysis payload. */
export function toPathList(files = []) {
  return files.map((file) => file?.path).filter(Boolean)
}
