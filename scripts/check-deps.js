#!/usr/bin/env node
/**
 * Pre-build check: reports how analyses will run.
 *
 * A packaged install ships the frozen analysis backend, so Python is not needed.
 * In a development checkout the backend is usually absent and the app falls back
 * to a system interpreter, which is what this script verifies.
 */
const { execSync } = require('child_process')
const path = require('path')
const fs = require('fs')

const projectRoot = path.resolve(__dirname, '..')
const backendExe = process.platform === 'win32' ? 'analysis-backend.exe' : 'analysis-backend'
const bundledBackend = path.join(projectRoot, 'python-dist', 'analysis-backend', backendExe)

const preferredPython = process.platform === 'win32'
  ? path.join(projectRoot, '.venv', 'Scripts', 'python.exe')
  : path.join(projectRoot, '.venv', 'bin', 'python')

const pythonCmd = fs.existsSync(preferredPython) ? preferredPython : 'python3'

function check(cmd, label) {
  try {
    const version = execSync(`"${cmd}" --version`, { stdio: 'pipe' }).toString().trim()
    console.log(`✅ ${label}: ${version}`)
    return true
  } catch {
    console.error(`❌ ${label} not found — install it or set the path in Settings.`)
    return false
  }
}

if (fs.existsSync(bundledBackend)) {
  const sizeMb = (fs.statSync(bundledBackend).size / 1024 / 1024).toFixed(1)
  console.log(`✅ Bundled analysis engine: ${bundledBackend} (${sizeMb} MB)`)
  console.log('   Python is NOT required at runtime — the packaged app is self-contained.')
  console.log('\n✅ All dependencies found.')
} else {
  console.log(`ℹ️  Bundled analysis engine not built yet (${bundledBackend} is missing).`)
  console.log('   Run "npm run build:backend" to create it — required before packaging.')
  console.log('   Development runs fall back to a system Python interpreter, checked below.\n')

  const ok = [check(pythonCmd, 'Python')]
  let pythonPackagesOk = true

  if (ok[0]) {
    try {
      execSync(`"${pythonCmd}" -c "import scipy, networkx, numpy"`, { stdio: 'pipe' })
      console.log('✅ Python packages: scipy, numpy, networkx')
    } catch {
      pythonPackagesOk = false
      console.error('❌ Missing Python packages: scipy, numpy and/or networkx')
      console.error('   Run: python3 -m pip install -r python/requirements.txt')
    }
  }

  ok.push(pythonPackagesOk)

  if (ok.every(Boolean)) {
    console.log('\n✅ All dependencies found.')
  } else {
    console.log('\n⚠️  Some runtimes are missing. The app will still run but those engines will fail.')
  }
}
