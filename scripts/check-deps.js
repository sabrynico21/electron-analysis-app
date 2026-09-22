#!/usr/bin/env node
/**
 * Pre-build check: verifies Python and the required analysis packages are available.
 */
const { execSync } = require('child_process')
const path = require('path')
const fs = require('fs')

const projectRoot = path.resolve(__dirname, '..')
const preferredPython = process.platform === 'win32'
  ? path.join(projectRoot, '.venv', 'Scripts', 'python.exe')
  : path.join(projectRoot, '.venv', 'bin', 'python')

const pythonCmd = fs.existsSync(preferredPython) ? preferredPython : 'python3'

function check(cmd, label) {
  try {
    const version = execSync(`${cmd} --version`, { stdio: 'pipe' }).toString().trim()
    console.log(`✅ ${label}: ${version}`)
    return true
  } catch {
    console.error(`❌ ${label} not found — install it or set the path in settings.`)
    return false
  }
}

const ok = [check(pythonCmd, 'Python')]
let pythonPackagesOk = true

if (ok[0]) {
  try {
    execSync(`${pythonCmd} -c \"import scipy, networkx, numpy\"`, { stdio: 'pipe' })
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
