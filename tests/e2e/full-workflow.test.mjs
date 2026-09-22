/**
 * End-to-end: a complete analysis lifecycle.
 *
 *   drop real files -> configure -> run (Python) -> the run is archived
 *   -> reopen it from the history -> re-run it in place -> compare it.
 *
 * This is the only suite that executes the Python backend, so it is skipped
 * unless explicitly requested:
 *
 *   npm run test:e2e:full
 *
 * It needs a working Python 3 interpreter. If the project virtualenv
 * (`<repo>/.venv`) already provides scipy/numpy/networkx nothing is installed;
 * otherwise the app creates a private environment inside the temporary
 * user-data directory, which requires internet access the first time.
 */
import assert from 'node:assert/strict'
import fs from 'node:fs'
import path from 'node:path'
import { after, before, describe, it } from 'node:test'

import { createFixtures, CONDITION_A, CONDITION_B } from './support/fixtures.mjs'
import { SET_SELECT_HELPER, launchApp, makeTempDir } from './support/appHarness.mjs'

const ANALYSIS_TIMEOUT_MS = 240_000

/** The Python backend is only exercised when the suite is explicitly enabled. */
const enabled = process.env.E2E_FULL === '1'

describe('full analysis lifecycle', { skip: enabled ? false : 'set E2E_FULL=1 (npm run test:e2e:full) to run the Python-backed suite' }, () => {
  let tempDir
  let app
  let fixtures
  let input

  before(async () => {
    tempDir = makeTempDir()
    fixtures = createFixtures(tempDir)
    input = fixtures.inputDir
    app = await launchApp({ userDataDir: fixtures.userDataDir })
    await app.reload()
    assert.equal(await app.client.waitFor(`document.querySelectorAll("[class*='zone']").length === 4`), true)
  })

  after(async () => {
    if (app) await app.close()
    if (tempDir) fs.rmSync(tempDir, { recursive: true, force: true })
  })

  /** Poll until no analysis is in progress and the result section is on screen. */
  const waitForCompletion = async () => {
    const finished = await app.client.waitFor(`
      [...document.querySelectorAll('h2')].some((heading) => /Cluster Results/i.test(heading.textContent)) &&
      ![...document.querySelectorAll('button')].some((button) => /Running/.test(button.textContent))
    `, { timeoutMs: ANALYSIS_TIMEOUT_MS })
    if (!finished) {
      const errors = await app.client.evaluate(`
        return [...document.querySelectorAll("[class*='error']")].map((node) => node.textContent).join(' | ')
      `)
      const output = app.readOutput().split('\n').slice(-12).join('\n')
      assert.fail(`the analysis did not complete. Errors: ${errors}\nLast application output:\n${output}`)
    }
  }

  it('runs an analysis and archives it', { timeout: ANALYSIS_TIMEOUT_MS + 60_000 }, async () => {
    const before = (await app.client.evaluate(`return (await window.electronAPI.listHistory()).entries.length`))

    await app.dropFiles("[data-testid='upload-zone-bacteria']", [path.join(input, 'abundance.csv')])
    await app.dropFiles("[data-testid='upload-zone-bacteria-metadata']", [path.join(input, 'metadata.tsv')])
    assert.equal(await app.client.waitFor(`document.querySelector('#conditionA').options.length > 2`), true)

    const configured = await app.client.evaluate(`
      ${SET_SELECT_HELPER}
      const setValue = (id, value) => {
        const element = document.getElementById(id)
        const prototype = element.tagName === 'SELECT'
          ? window.HTMLSelectElement.prototype
          : window.HTMLInputElement.prototype
        Object.getOwnPropertyDescriptor(prototype, 'value').set.call(element, value)
        element.dispatchEvent(new Event('input', { bubbles: true }))
        element.dispatchEvent(new Event('change', { bubbles: true }))
      }
      setValue('analysis-name', 'E2E lifecycle')
      setValue('datasetType', 'bacteria')
      setValue('correlationMethod', 'spearman')
      setValue('groupingMode', 'all')
      setValue('conditionA', ${JSON.stringify(CONDITION_A)})
      setValue('conditionB', ${JSON.stringify(CONDITION_B)})
      await new Promise((resolve) => setTimeout(resolve, 400))
      const button = [...document.querySelectorAll('button')].find((item) => /Start Analysis/.test(item.textContent))
      return {
        name: document.querySelector('#analysis-name').value,
        datasetType: document.querySelector('#datasetType').value,
        conditionA: document.querySelector('#conditionA').value,
        conditionB: document.querySelector('#conditionB').value,
        canStart: Boolean(button) && !button.disabled,
      }
    `)
    assert.equal(configured.canStart, true, `the run button is disabled: ${JSON.stringify(configured)}`)

    await app.client.evaluate(`
      [...document.querySelectorAll('button')].find((item) => /Start Analysis/.test(item.textContent)).click()
      return true
    `)
    await waitForCompletion()

    const archived = await app.client.evaluate(`
      const entries = (await window.electronAPI.listHistory()).entries
      return entries[0]
    `)
    assert.equal(archived.analysisName, 'E2E lifecycle')
    assert.equal(archived.id !== undefined, true)
    assert.ok(Object.keys(archived.summary || {}).length > 3, 'the result summary must be stored with the entry')
    assert.equal((archived.payloadSnapshot.files.bacteria || []).length, 1)
    assert.equal((archived.payloadSnapshot.files.bacteria_metadata || []).length, 1)
    assert.ok(archived.createdAt > 0 && archived.updatedAt >= archived.createdAt)

    const after = (await app.client.evaluate(`return (await window.electronAPI.listHistory()).entries.length`))
    assert.equal(after, before + 1, 'the run must add exactly one archive entry')

    // The stored analysis must reopen with its results on screen.
    await app.client.evaluate(`location.hash = '#/history/${archived.id}'; return true`)
    assert.equal(await app.client.waitFor(`document.querySelector("[data-testid='back-to-history']") !== null`, { timeoutMs: 30000 }), true)
    assert.equal(
      await app.client.waitFor(`[...document.querySelectorAll('h2')].some((h) => /Cluster Results/i.test(h.textContent))`, { timeoutMs: 60000 }),
      true
    )

    // Re-running must update that entry in place, not duplicate it.
    await app.client.evaluate(`
      [...document.querySelectorAll('button')].find((item) => /Re-run analysis/.test(item.textContent)).click()
      return true
    `)
    await waitForCompletion()

    const rerun = await app.client.evaluate(`
      const entries = (await window.electronAPI.listHistory()).entries
      return { count: entries.length, entry: entries.find((item) => item.id === '${archived.id}') }
    `)
    assert.equal(rerun.count, after, 're-running must not create a new entry')
    assert.equal(rerun.entry.runCount, 2, 'the run counter must be incremented')
    assert.equal(rerun.entry.createdAt, archived.createdAt, 'the creation time must be preserved')

    // Cleaning up must remove the entry.
    const deleted = await app.client.evaluate(`
      const response = await window.electronAPI.deleteHistory('${archived.id}')
      return response.success === true && !response.entries.some((item) => item.id === '${archived.id}')
    `)
    assert.equal(deleted, true)
  })
})
