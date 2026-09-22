#!/usr/bin/env node
/**
 * Build the frozen Python analysis backend with PyInstaller.
 *
 * The packaged application must not depend on whatever Python the user happens
 * to have installed, so the analysis pipeline (scipy + networkx + numpy and the
 * clustering package) is compiled into a single self-contained executable that
 * electron-builder ships under `resources/python-dist`.
 *
 * Usage:
 *   node scripts/build-backend.mjs                  # for the current OS
 *   node scripts/build-backend.mjs --target win32   # Windows, via wine (Linux/macOS host)
 *
 * PyInstaller cannot cross-compile, so a Windows build normally needs a Windows
 * machine. When the host is Linux or macOS and the target is `win32`, this script
 * instead installs a real Windows Python inside a wine prefix and runs PyInstaller
 * there, which produces a genuine PE64 executable.
 *
 * The output directory is always recreated, so `python-dist/analysis-backend`
 * holds exactly one platform's backend: the one you are about to package.
 *
 * Environment variables:
 *   PYTHON                 interpreter used to create the native build virtualenv
 *   BACKEND_INCREMENTAL    set to "1" to reuse the PyInstaller work directory
 *   WINE_PYTHON_VERSION    Windows Python version for the wine route (3.12.10)
 */
import { execFileSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const projectRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..')
const requirementsPath = path.join(projectRoot, 'python', 'requirements.txt')
const specPath = path.join(projectRoot, 'python', 'backend.spec')
const distPath = path.join(projectRoot, 'python-dist')
const backendDir = path.join(distPath, 'analysis-backend')
const workPath = path.join(projectRoot, '.pyinstaller-build')
const venvPath = path.join(projectRoot, '.pyinstaller-venv')
const markerPath = path.join(venvPath, '.build-fingerprint')
const incremental = process.env.BACKEND_INCREMENTAL === '1'
const winePrefix = path.join(projectRoot, '.pyinstaller-wine')
const winePythonVersion = process.env.WINE_PYTHON_VERSION || '3.12.10'

const SUPPORTED_TARGETS = ['linux', 'win32', 'darwin']

const hostPlatform = process.platform
const targetPlatform = resolveTarget()
const isWindows = targetPlatform === 'win32'
const executableName = isWindows ? 'analysis-backend.exe' : 'analysis-backend'

/**
 * Which platform we are building for: `--target <platform>`, or the host.
 * @returns {'linux'|'win32'|'darwin'}
 */
function resolveTarget() {
  const index = process.argv.indexOf('--target')
  if (index === -1) return hostPlatform
  const value = String(process.argv[index + 1] || '').trim()
  if (!SUPPORTED_TARGETS.includes(value)) {
    throw new Error(`--target must be one of ${SUPPORTED_TARGETS.join(', ')} (received "${value}")`)
  }
  return value
}

const venvPython = isWindows
  ? path.join(venvPath, 'Scripts', 'python.exe')
  : path.join(venvPath, 'bin', 'python')

function run(command, args, label) {
  console.log(`\n[backend] ${label}`)
  console.log(`[backend] $ ${command} ${args.join(' ')}`)
  execFileSync(command, args, { cwd: projectRoot, stdio: 'inherit' })
}

function candidateInterpreters() {
  const candidates = []
  if (process.env.PYTHON) candidates.push(process.env.PYTHON)
  candidates.push(
    isWindows ? path.join(projectRoot, '.venv', 'Scripts', 'python.exe') : path.join(projectRoot, '.venv', 'bin', 'python'),
  )
  candidates.push(isWindows ? 'python' : 'python3', 'python')
  return candidates
}

function resolveBasePython() {
  for (const candidate of candidateInterpreters()) {
    try {
      execFileSync(candidate, ['--version'], { stdio: 'pipe' })
      return candidate
    } catch {
      // try the next candidate
    }
  }
  throw new Error(
    'Python 3 was not found. Install it, or point the PYTHON environment variable at the interpreter to use for the build.',
  )
}

/** Requirements + PyInstaller version decide whether the venv is still valid. */
function buildFingerprint() {
  const hash = createHash('sha256')
  hash.update(fs.readFileSync(requirementsPath))
  hash.update('pyinstaller>=6.0')
  return hash.digest('hex')
}

function ensureVenv(basePython, fingerprint) {
  const needsCreation = !fs.existsSync(venvPython)
  if (needsCreation) {
    run(basePython, ['-m', 'venv', venvPath], 'Creating the isolated build environment')
  }

  const current = fs.existsSync(markerPath)
    ? fs.readFileSync(markerPath, 'utf-8').trim()
    : null
  if (!needsCreation && current === fingerprint) {
    console.log('[backend] Build environment is up to date — skipping pip install.')
    return
  }

  run(venvPython, ['-m', 'pip', 'install', '--disable-pip-version-check', '--upgrade', 'pip'], 'Upgrading pip')
  run(
    venvPython,
    ['-m', 'pip', 'install', '--disable-pip-version-check', '-r', requirementsPath, 'pyinstaller>=6.0'],
    'Installing analysis requirements and PyInstaller',
  )
  fs.writeFileSync(markerPath, fingerprint, 'utf-8')
}

function resetOutputDir() {
  // extraResources copies `python-dist` verbatim, so it must never accumulate
  // two platforms' backends at the same time.
  fs.rmSync(backendDir, { recursive: true, force: true })
  fs.mkdirSync(backendDir, { recursive: true })
}

function verifyOutput() {
  const exePath = path.join(backendDir, executableName)
  if (!fs.existsSync(exePath)) {
    throw new Error(`Build finished but ${exePath} is missing.`)
  }
  if (!isWindows) {
    fs.chmodSync(exePath, 0o755)
  }
  const bytes = fs.statSync(exePath).size
  console.log(`\n[backend] Ready: ${exePath} (${(bytes / 1024 / 1024).toFixed(1)} MB) for ${targetPlatform}`)
  console.log('[backend] The packaged app will pick this up automatically.')
}

// ---------------------------------------------------------------------------
// Native route: the target matches the host, use a local virtualenv.
// ---------------------------------------------------------------------------

function buildNative(fingerprint) {
  const basePython = resolveBasePython()
  console.log(`[backend] Using ${basePython}`)

  ensureVenv(basePython, fingerprint)

  const pyinstallerArgs = [
    '-m', 'PyInstaller',
    specPath,
    '--noconfirm',
    '--distpath', distPath,
    '--workpath', workPath,
  ]
  if (!incremental) pyinstallerArgs.push('--clean')

  run(venvPython, pyinstallerArgs, 'Freezing the analysis backend')
}

// ---------------------------------------------------------------------------
// Wine route: build a Windows executable from a Linux / macOS host.
// ---------------------------------------------------------------------------

/** Environment shared by every wine invocation. */
function wineEnv() {
  return {
    ...process.env,
    WINEPREFIX: winePrefix,
    WINEARCH: 'win64',
    // Skip the Mono/Gecko download prompts and window-driver noise.
    WINEDLLOVERRIDES: 'mscoree,mshtml=',
    WINEDEBUG: '-all',
  }
}

function wine(args, label) {
  console.log(`\n[backend] ${label}`)
  console.log(`[backend] $ wine ${args.join(' ')}`)
  execFileSync('wine', args, { cwd: projectRoot, stdio: 'inherit', env: wineEnv() })
}

/** Windows path (C:\..., Z:\...) for a host path, as wine sees it. */
function toWindowsPath(hostPath) {
  return execFileSync('winepath', ['-w', hostPath], { env: wineEnv() }).toString().trim()
}

// python.org installs into "Python3<minor>" (e.g. 3.12 -> Python312).
const WINDOWS_PYTHON_DIR = `C:\\Program Files\\Python${winePythonVersion.split('.')[0]}${winePythonVersion.split('.')[1]}`
const WINDOWS_PYTHON = `${WINDOWS_PYTHON_DIR}\\python.exe`

function ensureWineAvailable() {
  for (const tool of ['wine', 'winepath']) {
    try {
      execFileSync('which', [tool], { stdio: 'pipe' })
    } catch {
      throw new Error(
        `${tool} is not installed. Install it first (on Arch: "sudo pacman -S wine"), or build the Windows package on a Windows machine.`,
      )
    }
  }
}

function ensureWinePrefix() {
  if (fs.existsSync(path.join(winePrefix, 'drive_c'))) return
  console.log(`[backend] Creating the wine prefix at ${winePrefix}`)
  execFileSync('wineboot', ['-u'], { cwd: projectRoot, stdio: 'inherit', env: wineEnv() })
}

function ensureWindowsPython(fingerprint) {
  const installed = (() => {
    try {
      return execFileSync('wine', [WINDOWS_PYTHON, '--version'], { env: wineEnv(), stdio: 'pipe' })
        .toString()
        .trim()
    } catch {
      return null
    }
  })()

  if (!installed) {
    const installerName = `python-${winePythonVersion}-amd64.exe`
    const cacheDir = path.join(os.homedir(), '.cache', 'electron-analysis-app')
    fs.mkdirSync(cacheDir, { recursive: true })
    const installerPath = path.join(cacheDir, installerName)

    if (!fs.existsSync(installerPath)) {
      const url = `https://www.python.org/ftp/python/${winePythonVersion}/${installerName}`
      console.log(`[backend] Downloading ${url}`)
      execFileSync('curl', ['-fL', '--retry', '3', '-o', installerPath, url], { stdio: 'inherit' })
    }

    wine(
      [
        installerPath,
        '/quiet',
        'InstallAllUsers=1',
        'PrependPath=1',
        'Include_test=0',
        'Include_launcher=0',
        `TargetDir=${WINDOWS_PYTHON_DIR}`,
      ],
      `Installing Windows Python ${winePythonVersion} into the wine prefix`,
    )
  } else {
    console.log(`[backend] ${installed} is already available in the wine prefix`)
  }

  const markerFile = path.join(winePrefix, 'backend-build-fingerprint')
  const current = fs.existsSync(markerFile) ? fs.readFileSync(markerFile, 'utf-8').trim() : null
  if (current === fingerprint) {
    console.log('[backend] Windows build environment is up to date — skipping pip install.')
    return
  }

  const requirementsWindows = toWindowsPath(requirementsPath)
  wine([WINDOWS_PYTHON, '-m', 'pip', 'install', '--disable-pip-version-check', '--upgrade', 'pip'], 'Upgrading pip')
  wine(
    [WINDOWS_PYTHON, '-m', 'pip', 'install', '--disable-pip-version-check', '-r', requirementsWindows, 'pyinstaller>=6.0'],
    'Installing analysis requirements and PyInstaller (Windows)',
  )
  fs.writeFileSync(markerFile, fingerprint, 'utf-8')
}

function buildViaWine(fingerprint) {
  ensureWineAvailable()
  ensureWinePrefix()
  ensureWindowsPython(fingerprint)

  const args = [
    WINDOWS_PYTHON, '-m', 'PyInstaller',
    toWindowsPath(specPath),
    '--noconfirm',
    '--distpath', toWindowsPath(distPath),
    '--workpath', toWindowsPath(path.join(workPath, `${targetPlatform}-${process.arch}`)),
  ]
  if (!incremental) args.push('--clean')

  wine(args, 'Freezing the analysis backend for Windows (inside wine)')
}

function main() {
  if (targetPlatform !== hostPlatform && !(targetPlatform === 'win32' && (hostPlatform === 'linux' || hostPlatform === 'darwin'))) {
    throw new Error(
      `Cannot build the ${targetPlatform} backend on a ${hostPlatform} host. ` +
      'PyInstaller has no cross-compiler: run this on the matching OS, or let the wine route handle win32.',
    )
  }

  const fingerprint = buildFingerprint()
  resetOutputDir()

  if (targetPlatform === 'win32' && hostPlatform !== 'win32') {
    console.log(`[backend] Building for win32 on ${hostPlatform} through wine`)
    buildViaWine(fingerprint)
  } else {
    buildNative(fingerprint)
  }

  verifyOutput()
}

try {
  main()
} catch (error) {
  console.error(`\n[backend] Build failed: ${error.message}`)
  process.exitCode = 1
}
