/**
 * End-to-end: the file-loading pipeline.
 *
 * This is the regression suite for the reported bugs:
 *   - a file chosen with the native picker must have its contents read so that
 *     section 2 can offer the conditions;
 *   - a file dropped from the file manager must keep its real path (Electron 32
 *     removed `File.path`; `webUtils.getPathForFile` must be used *synchronously*
 *     inside the drop handler, otherwise the path degrades to `./<name>` and the
 *     file is reported as missing).
 *
 * Files are dropped here from real paths through the DevTools protocol, which is
 * the only way to reproduce a genuine OS drop.
 *
 *   npm run test:e2e        (after `npm run build:renderer`)
 */
import assert from 'node:assert/strict'
import fs from 'node:fs'
import path from 'node:path'
import { after, before, describe, it } from 'node:test'

import { createFixtures, CONDITION_A, CONDITION_B } from './support/fixtures.mjs'
import { launchApp, makeTempDir } from './support/appHarness.mjs'

const ZONE_BACTERIA = "[data-testid='upload-zone-bacteria']"
const ZONE_METADATA = "[data-testid='upload-zone-bacteria-metadata']"

/** Text shown in the bacteria-metadata upload column (zone + error message). */
const metadataColumnText = `document.querySelector(${JSON.stringify(ZONE_METADATA)}).parentElement.innerText`

describe('file loading (drag & drop and native picker)', () => {
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
  })

  after(async () => {
    if (app) await app.close()
    if (tempDir) fs.rmSync(tempDir, { recursive: true, force: true })
  })

  /** Reset to an empty new-analysis form. */
  const resetForm = async () => {
    await app.client.evaluate(`location.hash = '#/history'; return true`)
    await app.client.waitFor(`document.querySelector("[data-testid='history-list']") !== null`)
    await app.client.evaluate(`location.hash = '#/'; return true`)
    assert.equal(await app.client.waitFor(`document.querySelector(${JSON.stringify(ZONE_BACTERIA)}) !== null`), true)
  }

  it('resolves the real path of a dropped file', async () => {
    await resetForm()
    await app.dropFiles(ZONE_BACTERIA, [path.join(input, 'abundance.csv')])

    assert.equal(
      await app.client.waitFor(`document.querySelectorAll("[class*='fileItem']").length === 1`),
      true,
      'the dropped file was not registered'
    )
    const state = await app.client.evaluate(`
      const zone = document.querySelector(${JSON.stringify(ZONE_BACTERIA)})
      return {
        name: zone.querySelector("[class*='fileName']")?.textContent || '',
        error: zone.querySelector("[data-testid='upload-error']")?.textContent || '',
      }
    `)
    assert.equal(state.name, 'abundance.csv')
    assert.equal(state.error, '', 'the drop must not report a path problem')
  })

  it('populates the section-2 selects from the dropped metadata contents', async () => {
    await app.dropFiles(ZONE_METADATA, [path.join(input, 'metadata.tsv')])

    assert.equal(
      await app.client.waitFor(`document.querySelector('#conditionA').options.length > 2`, { timeoutMs: 20000 }),
      true,
      'the condition selects were not filled from the metadata file'
    )
    const conditions = await app.client.evaluate(`
      const options = [...document.querySelector('#conditionA').options].map((option) => option.value).filter(Boolean)
      return { options, datasetTypeEnabled: !document.querySelector('#datasetType').disabled }
    `)
    assert.deepEqual(conditions.options.sort(), [CONDITION_A, CONDITION_B].sort())
  })

  it('de-duplicates a file added twice', async () => {
    await app.dropFiles(ZONE_BACTERIA, [path.join(input, 'abundance.csv')])
    await new Promise((resolve) => setTimeout(resolve, 500))
    const count = await app.client.evaluate(`return document.querySelectorAll("[class*='fileItem']").length`)
    assert.equal(count, 2, 'the same path must not be added twice (1 abundance + 1 metadata expected)')
  })

  it('reads a semicolon-delimited metadata file too', async () => {
    await resetForm()
    await app.dropFiles(ZONE_METADATA, [path.join(input, 'metadata_semicolon.csv')])
    assert.equal(
      await app.client.waitFor(`document.querySelector('#conditionA').options.length === 3`),
      true,
      `expected the ";" delimiter to be detected, got ${await app.client.evaluate(`return document.querySelector('#conditionA').options.length`)} options`
    )
  })

  it('reports a metadata file without the required column', async () => {
    await resetForm()
    await app.dropFiles(ZONE_METADATA, [path.join(input, 'metadata_no_column.tsv')])
    assert.equal(
      await app.client.waitFor(`/Combination_treat2/.test(${metadataColumnText} || '')`),
      true,
      'the missing column was not reported'
    )
    const options = await app.client.evaluate(`return document.querySelector('#conditionA').options.length`)
    assert.equal(options, 1, 'the selects must stay empty when the metadata cannot be read')
  })

  it('reports a metadata file with a single condition', async () => {
    await resetForm()
    await app.dropFiles(ZONE_METADATA, [path.join(input, 'metadata_no_conditions.tsv')])
    assert.equal(
      await app.client.waitFor(`/At least 2 distinct values/i.test(${metadataColumnText} || '')`),
      true,
      `the single-condition metadata was not reported (column text: ${await app.client.evaluate(`return ${metadataColumnText}`)})`
    )
  })

  it('removes a file and clears the conditions', async () => {
    await resetForm()
    await app.dropFiles(ZONE_METADATA, [path.join(input, 'metadata.tsv')])
    assert.equal(await app.client.waitFor(`document.querySelector('#conditionA').options.length > 2`), true)

    await app.client.evaluate(`
      const zone = document.querySelector(${JSON.stringify(ZONE_METADATA)})
      zone.querySelector("[class*='remove']").click()
      return true
    `)
    assert.equal(
      await app.client.waitFor(`
        document.querySelectorAll("[class*='fileItem']").length === 0 &&
        document.querySelector('#conditionA').options.length === 1
      `),
      true,
      'removing the metadata file must clear the selects'
    )
  })

  it('starts a new analysis from a blank form again', async () => {
    await app.dropFiles(ZONE_BACTERIA, [path.join(input, 'abundance.csv')])
    await app.dropFiles(ZONE_METADATA, [path.join(input, 'metadata.tsv')])
    assert.equal(await app.client.waitFor(`document.querySelector('#conditionA').options.length > 2`), true)
    await app.client.evaluate(`
      const field = document.querySelector('#analysis-name')
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value').set
      setter.call(field, 'Half filled analysis')
      field.dispatchEvent(new Event('input', { bubbles: true }))
      return true
    `)

    await resetForm()
    const fresh = await app.client.evaluate(`
      return {
        name: document.querySelector('#analysis-name').value,
        files: document.querySelectorAll("[class*='fileItem']").length,
        conditions: document.querySelector('#conditionA').options.length,
        results: [...document.querySelectorAll('h2')].some((heading) => /Cluster Results/i.test(heading.textContent)),
      }
    `)
    assert.equal(fresh.name, '', 'the analysis name must be reset')
    assert.equal(fresh.files, 0, 'the uploaded files must be cleared')
    assert.equal(fresh.conditions, 1, 'the conditions must be cleared')
    assert.equal(fresh.results, false, 'no results should be shown')
  })
})
