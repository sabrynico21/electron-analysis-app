/**
 * Detects the Python interpreter path and resolves the FastSpar (SparCC) runtime.
 * Falls back to bundled binaries if system ones are not found.
 */
const { which } = require('shelljs')
const path = require('path')
const fs = require('fs')
const { app } = require('electron')
const Store = require('electron-store')

const store = new Store()

/**
 * Conventional Python 3 install locations, probed only when `PATH` lookup fails.
 *
 * A GUI application launched from Finder (macOS) or the Start Menu (Windows)
 * inherits a minimal environment, so `python3` is frequently not on `PATH` even
 * though it is installed. Checking these locations keeps auto-detection working
 * for users who never open a terminal.
 * @returns {string[]}
 */
function _wellKnownPythonCandidates() {
  if (process.platform === 'darwin') {
    return [
      '/opt/homebrew/bin/python3', // Homebrew on Apple Silicon
      '/usr/local/bin/python3', // Homebrew on Intel
      '/usr/bin/python3', // Xcode Command Line Tools
      '/Library/Frameworks/Python.framework/Versions/Current/bin/python3', // python.org
    ]
  }

  if (process.platform === 'win32') {
    const candidates = []
    const localAppData = process.env.LOCALAPPDATA
    if (localAppData) {
      for (const version of ['313', '312', '311', '310', '39']) {
        candidates.push(path.join(localAppData, 'Programs', 'Python', `Python${version}`, 'python.exe'))
      }
    }
    for (const version of ['313', '312', '311', '310']) {
      candidates.push(`C:\\Python${version}\\python.exe`)
      candidates.push(`C:\\Program Files\\Python${version}\\python.exe`)
    }
    return candidates
  }

  return ['/usr/bin/python3', '/usr/local/bin/python3', '/bin/python3']
}

async function getPythonPath() {
  const custom = store.get('pythonPath')
  if (custom && fs.existsSync(custom)) return custom

  // Prefer local project virtualenv when available (development checkout).
  const appPath = app.getAppPath()
  const venvCandidate = process.platform === 'win32'
    ? path.join(appPath, '.venv', 'Scripts', 'python.exe')
    : path.join(appPath, '.venv', 'bin', 'python')
  if (fs.existsSync(venvCandidate)) {
    return venvCandidate
  }

  for (const cmd of ['python3', 'python']) {
    const found = which(cmd)
    if (found) return found.toString()
  }

  for (const candidate of _wellKnownPythonCandidates()) {
    if (fs.existsSync(candidate)) return candidate
  }

  // Fail early with an actionable message instead of letting spawn() report a
  // bare ENOENT from a path that was never going to exist.
  throw new Error(
    'Python 3 was not found on this computer. Install Python 3 (python.org), ' +
    'then restart the application — or select the interpreter manually in Settings → Python.'
  )
}

function _candidateFastsparBasename() {
  return process.platform === 'win32' ? 'fastspar.exe' : 'fastspar'
}

function _bundledFastsparCandidate() {
  const appPath = app.getAppPath()
  const baseName = _candidateFastsparBasename()

  // Works in dev (project resources folder).
  const devCandidate = path.join(appPath, 'resources', 'bin', process.platform, process.arch, baseName)
  if (fs.existsSync(devCandidate)) {
    return devCandidate
  }

  // Works in packaged apps (extraResources copied under process.resourcesPath/bin).
  const packagedCandidate = path.join(process.resourcesPath || appPath, 'bin', process.platform, process.arch, baseName)
  if (fs.existsSync(packagedCandidate)) {
    return packagedCandidate
  }

  return null
}

/**
 * Guarantees that a POSIX binary carries the executable bit. Packaged builds
 * normally preserve file modes, but archives (zip), non-POSIX checkouts and some
 * installers can drop it, which would leave SparCC silently unusable.
 * @returns {{ ok: boolean, error?: string }}
 */
function _ensureExecutable(filePath) {
  if (process.platform === 'win32') return { ok: true }

  const alreadyExecutable = () => {
    try {
      fs.accessSync(filePath, fs.constants.X_OK)
      return true
    } catch {
      return false
    }
  }

  if (alreadyExecutable()) return { ok: true }

  try {
    fs.chmodSync(filePath, 0o755)
  } catch (error) {
    return { ok: false, error: error.message }
  }

  return alreadyExecutable()
    ? { ok: true }
    : { ok: false, error: 'file is not executable after chmod' }
}

function resolveFastsparPath(options = {}) {
  const notExecutable = (filePath, source) => ({
    available: false,
    path: filePath,
    source,
    error: `FastSpar binary at ${filePath} is not executable and the permission could not be set.`,
  })

  const custom = (options.customPath || '').trim()
  if (custom && fs.existsSync(custom)) {
    return _ensureExecutable(custom).ok
      ? { available: true, path: custom, source: 'settings' }
      : notExecutable(custom, 'settings')
  }

  const bundled = _bundledFastsparCandidate()
  if (bundled) {
    return _ensureExecutable(bundled).ok
      ? { available: true, path: bundled, source: 'bundled' }
      : notExecutable(bundled, 'bundled')
  }

  const fromEnv = (process.env.FASTSPAR_BIN || '').trim()
  if (fromEnv && fs.existsSync(fromEnv)) {
    return _ensureExecutable(fromEnv).ok
      ? { available: true, path: fromEnv, source: 'env' }
      : notExecutable(fromEnv, 'env')
  }

  const found = which('fastspar')
  if (found) {
    const resolved = found.toString()
    if (_ensureExecutable(resolved).ok) {
      return { available: true, path: resolved, source: 'path' }
    }
  }

  return {
    available: false,
    path: '',
    source: 'none',
    error: 'FastSpar binary not found in settings path, bundled resources, FASTSPAR_BIN, or PATH.',
  }
}

function _siblingBinaryName(name) {
  return process.platform === 'win32' ? `${name}.exe` : name
}

function resolveFastsparCompanions(fastsparPath) {
  // FastSpar's permutation p-values need fastspar_bootstrap and fastspar_pvalues
  // alongside fastspar. Resolve them relative to the resolved fastspar path.
  const result = { fastspar_bootstrap: '', fastspar_pvalues: '' }
  if (!fastsparPath) return result
  const dir = path.dirname(fastsparPath)
  for (const name of Object.keys(result)) {
    const candidate = path.join(dir, _siblingBinaryName(name))
    if (fs.existsSync(candidate)) {
      result[name] = candidate
    }
  }
  return result
}

function getSparccRuntimeStatus(options = {}) {
  const resolved = resolveFastsparPath(options)
  const companions = resolved.available
    ? resolveFastsparCompanions(resolved.path)
    : { fastspar_bootstrap: '', fastspar_pvalues: '' }
  const missingCompanions = resolved.available
    ? ['fastspar_bootstrap', 'fastspar_pvalues'].filter((name) => !companions[name])
    : []
  return {
    available: !!resolved.available && missingCompanions.length === 0,
    path: resolved.path || '',
    source: resolved.source,
    platform: process.platform,
    arch: process.arch,
    companions,
    missingCompanions,
    error: resolved.error
      || (missingCompanions.length
        ? `Missing FastSpar companion binaries: ${missingCompanions.join(', ')}`
        : null),
  }
}

module.exports = { getPythonPath, resolveFastsparPath, resolveFastsparCompanions, getSparccRuntimeStatus }
