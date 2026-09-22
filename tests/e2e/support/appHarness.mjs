/**
 * Launch a real Electron instance for the end-to-end suites.
 *
 * The application is started exactly like a user would start it (`electron .`,
 * loading the built renderer from `dist/renderer`) but with an isolated
 * `--user-data-dir`, so the tests are repeatable and never touch the developer's
 * own archive.
 *
 * The page is driven over the Chrome DevTools Protocol, which is what makes it
 * possible to test things Playwright cannot express here — most importantly real
 * drag & drop, which needs `Input.dispatchDragEvent` with actual files on disk.
 */
import { spawn } from 'node:child_process'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'

export const REPO_ROOT = path.resolve(import.meta.dirname, '..', '..', '..')
const ELECTRON_BIN = path.join(REPO_ROOT, 'node_modules', 'electron', 'dist', 'electron')
const RENDERER_ENTRY = path.join(REPO_ROOT, 'dist', 'renderer', 'index.html')

let nextPort = 9400

export class CdpClient {
  constructor(ws) {
    this.ws = ws
    this.sequence = 0
    this.pending = new Map()
    ws.addEventListener('message', (event) => {
      const message = JSON.parse(event.data)
      const resolve = this.pending.get(message.id)
      if (resolve) {
        this.pending.delete(message.id)
        resolve(message)
      }
    })
  }

  send(method, params = {}, timeoutMs = 60000) {
    return Promise.race([
      new Promise((resolve) => {
        const id = ++this.sequence
        this.pending.set(id, resolve)
        this.ws.send(JSON.stringify({ id, method, params }))
      }),
      new Promise((_, reject) => setTimeout(() => reject(new Error(`CDP timeout: ${method}`)), timeoutMs)),
    ])
  }

  /** Evaluate an async expression in the page and return its value. */
  async evaluate(expression, timeoutMs) {
    const message = await this.send('Runtime.evaluate', {
      expression: `(async () => { ${expression} })()`,
      awaitPromise: true,
      returnByValue: true,
    }, timeoutMs)
    const exception = message.result?.exceptionDetails
    if (exception) {
      throw new Error(exception.exception?.description || exception.text || 'evaluate failed')
    }
    return message.result?.result?.value
  }

  /** Evaluate a statement list that must `return` a boolean, polling until it is true. */
  async waitFor(predicate, { timeoutMs = 20000, intervalMs = 120 } = {}) {
    const result = await this.evaluate(`
      const deadline = Date.now() + ${timeoutMs}
      while (Date.now() < deadline) {
        if (${predicate}) return true
        await new Promise((resolve) => setTimeout(resolve, ${intervalMs}))
      }
      return false
    `, timeoutMs + 5000)
    return result === true
  }
}

async function findPageTarget(port, timeoutMs) {
  const deadline = Date.now() + timeoutMs
  let lastError = null
  while (Date.now() < deadline) {
    try {
      const response = await fetch(`http://127.0.0.1:${port}/json/list`)
      const targets = await response.json()
      const page = targets.find((target) => target.type === 'page')
      if (page) return page
    } catch (error) {
      lastError = error
    }
    await new Promise((resolve) => setTimeout(resolve, 200))
  }
  throw new Error(`The application did not expose a page on port ${port}${lastError ? `: ${lastError.message}` : ''}`)
}

/**
 * Start the application. Call `close()` when done.
 *
 * @param {{ userDataDir: string, cwd?: string }} options
 */
export async function launchApp({ userDataDir, cwd = REPO_ROOT }) {
  if (!fs.existsSync(ELECTRON_BIN)) {
    throw new Error(`Electron binary not found at ${ELECTRON_BIN}. Run "npm install" first.`)
  }
  if (!fs.existsSync(RENDERER_ENTRY)) {
    throw new Error(
      `The renderer bundle is missing (${RENDERER_ENTRY}). Run "npm run build:renderer" before the e2e suite.`
    )
  }

  const port = nextPort++
  // The child output goes to a file rather than a pipe: Electron's helper
  // processes inherit the stdio descriptors, and a pipe kept open by an orphan
  // helper would stop the test runner from ever exiting.
  const logPath = path.join(userDataDir, 'app-output.log')
  const logFd = fs.openSync(logPath, 'a')
  const child = spawn(
    ELECTRON_BIN,
    ['.', `--remote-debugging-port=${port}`, `--user-data-dir=${userDataDir}`, '--no-sandbox'],
    { cwd, stdio: ['ignore', logFd, logFd] }
  )
  fs.closeSync(logFd)
  child.unref()

  const readOutput = () => {
    try {
      return fs.readFileSync(logPath, 'utf-8')
    } catch {
      return ''
    }
  }

  let client
  let ws
  try {
    const target = await findPageTarget(port, 30000)
    ws = new WebSocket(target.webSocketDebuggerUrl)
    await Promise.race([
      new Promise((resolve, reject) => {
        ws.addEventListener('open', resolve)
        ws.addEventListener('error', reject)
      }),
      new Promise((_, reject) => setTimeout(() => reject(new Error('CDP socket did not open')), 10000)),
    ])
    client = new CdpClient(ws)
    await client.send('Runtime.enable')
    await client.send('Page.enable')
    // Required before synthetic drag events can be dispatched.
    await client.send('Input.setInterceptDrags', { enabled: true })
  } catch (error) {
    child.kill('SIGKILL')
    throw new Error(`${error.message}\n--- application output ---\n${readOutput()}`)
  }

  return {
    client,
    port,
    readOutput,
    /** Relaunch the page so a fresh renderer bundle / clean state is used. */
    async reload() {
      await client.send('Page.reload', { ignoreCache: true })
      await new Promise((resolve) => setTimeout(resolve, 2500))
    },
    /** Drop real files from disk onto the element matching the CSS selector. */
    async dropFiles(selector, filePaths) {
      const box = await client.evaluate(`
        const element = document.querySelector(${JSON.stringify(selector)})
        if (!element) return null
        const rect = element.getBoundingClientRect()
        return { x: Math.round(rect.left + rect.width / 2), y: Math.round(rect.top + rect.height / 2) }
      `)
      if (!box) throw new Error(`Cannot drop files: "${selector}" was not found`)

      const data = { items: [], files: filePaths, dragOperationsMask: 1 }
      await client.send('Input.dispatchDragEvent', { type: 'dragEnter', x: box.x, y: box.y, data })
      await client.send('Input.dispatchDragEvent', { type: 'dragOver', x: box.x, y: box.y, data })
      await client.send('Input.dispatchDragEvent', { type: 'drop', x: box.x, y: box.y, data })
    },
    async close() {
      try {
        ws.close()
      } catch {
        // already closed
      }
      if (child.exitCode !== null || child.signalCode !== null) return
      child.kill('SIGTERM')
      await new Promise((resolve) => {
        const timer = setTimeout(() => {
          child.kill('SIGKILL')
          resolve()
        }, 5000)
        child.on('exit', () => {
          clearTimeout(timer)
          resolve()
        })
      })
      // Give the helper processes a moment to release the user-data directory.
      await new Promise((resolve) => setTimeout(resolve, 250))
    },
  }
}

/** Create a temporary workspace directory for one suite. */
export function makeTempDir(prefix = 'analysis-app-e2e-') {
  return fs.mkdtempSync(path.join(os.tmpdir(), prefix))
}

/** Set a <select> value the way React expects (native setter + change event). */
export const SET_SELECT_HELPER = `
  const setSelectValue = (element, value) => {
    Object.getOwnPropertyDescriptor(window.HTMLSelectElement.prototype, 'value').set.call(element, value)
    element.dispatchEvent(new Event('change', { bubbles: true }))
  }
`
