/**
 * End-to-end: history archive, reopening a stored analysis, the comparison view
 * and the scroll architecture.
 *
 * Runs against an isolated user-data directory (see `support/fixtures.mjs`), so
 * it is repeatable and needs no prior analysis.
 *
 *   npm run test:e2e        (after `npm run build:renderer`)
 */
import assert from 'node:assert/strict'
import fs from 'node:fs'
import { after, before, describe, it } from 'node:test'

import { createFixtures, CONDITION_A, CONDITION_B, IDS, TIMESTAMPS } from './support/fixtures.mjs'
import { SET_SELECT_HELPER, launchApp, makeTempDir } from './support/appHarness.mjs'

const LIST_ITEM = "[data-testid='history-item']"

describe('history archive and comparison view', () => {
  let tempDir
  let app
  let fixtures

  before(async () => {
    tempDir = makeTempDir()
    fixtures = createFixtures(tempDir)
    app = await launchApp({ userDataDir: fixtures.userDataDir })
    await app.reload()
    assert.equal(
      await app.client.waitFor(`document.querySelectorAll("[class*='zone']").length === 4`, { timeoutMs: 25000 }),
      true,
      'the application did not finish rendering'
    )
  })

  after(async () => {
    if (app) await app.close()
    if (tempDir) fs.rmSync(tempDir, { recursive: true, force: true })
  })

  it('opens on a reset "Create new analysis" form', async () => {
    const state = await app.client.evaluate(`
      return {
        hash: location.hash,
        title: document.querySelector("[data-testid='workspace-title']")?.textContent?.trim() || '',
        name: document.querySelector('#analysis-name')?.value ?? null,
        uploaded: document.querySelectorAll("[class*='fileItem']").length,
        conditions: document.querySelector('#conditionA')?.options.length ?? null,
      }
    `)
    // HashRouter leaves the hash empty while on the default route.
    assert.ok(state.hash === '' || state.hash === '#/', `unexpected route: ${state.hash}`)
    assert.equal(state.name, '')
    assert.equal(state.uploaded, 0)
    assert.equal(state.conditions, 1, 'only the placeholder option should be present')
    assert.match(state.title, /create new analysis/i)
  })

  it('never scrolls the window and keeps the navigation sidebar fixed', async () => {
    const result = await app.client.evaluate(`
      const sidebar = document.querySelector('nav')
      const before = sidebar.getBoundingClientRect().top
      const scroller = document.querySelector('main')
      scroller.scrollTop = scroller.scrollHeight
      await new Promise((resolve) => setTimeout(resolve, 120))
      return {
        windowScroll: window.scrollY,
        documentOverflows: document.documentElement.scrollHeight > window.innerHeight + 1,
        sidebarMoved: Math.abs(sidebar.getBoundingClientRect().top - before) > 1,
        contentOverflowY: getComputedStyle(scroller).overflowY,
      }
    `)
    assert.equal(result.windowScroll, 0, 'the window itself must not scroll')
    assert.equal(result.documentOverflows, false)
    assert.equal(result.sidebarMoved, false, 'the sidebar must stay in place')
    assert.equal(result.contentOverflowY, 'auto')
  })

  it('lists every stored analysis with its name, creation and update time', async () => {
    await app.client.evaluate(`location.hash = '#/history'; return true`)
    assert.equal(await app.client.waitFor(`document.querySelectorAll("${LIST_ITEM}").length === 3`), true)

    const items = await app.client.evaluate(`
      return [...document.querySelectorAll("${LIST_ITEM}")].map((item) => ({
        id: item.dataset.jobId,
        text: item.innerText,
      }))
    `)
    // Newest update first: the fixture timestamps make the order deterministic.
    assert.deepEqual(items.map((item) => item.id), [IDS.rich, IDS.plain, IDS.missingFiles])
    assert.match(items[0].text, /Rich analysis/)
    assert.match(items[0].text, /Created/)
    assert.match(items[0].text, /Updated/)
    const expectedDate = new Date(TIMESTAMPS.richAnalysis).toLocaleString()
    assert.ok(items[0].text.includes(expectedDate), `expected the creation date "${expectedDate}" in "${items[0].text}"`)
  })

  it('keeps the history list narrow so the analysis has room', async () => {
    const width = await app.client.evaluate(`
      const list = document.querySelector("[data-testid='history-list']").closest('aside')
      return Math.round(list.getBoundingClientRect().width)
    `)
    assert.ok(width <= 340, `the list pane is ${width}px wide`)
  })

  it('survives wiping the browser profile (the archive lives on disk)', async () => {
    await app.client.evaluate(`window.localStorage.clear(); return true`)
    await app.reload()
    await app.client.evaluate(`location.hash = '#/history'; return true`)
    assert.equal(
      await app.client.waitFor(`document.querySelectorAll("${LIST_ITEM}").length === 3`),
      true,
      'the history was lost with the browser storage'
    )
  })

  it('reopens a stored analysis in the history view, not in the new-analysis page', async () => {
    await app.client.evaluate(`location.hash = '#/history/${IDS.rich}'; return true`)
    assert.equal(
      await app.client.waitFor(`document.querySelector("[data-testid='back-to-history']") !== null`, { timeoutMs: 30000 }),
      true
    )
    // The stored analysis is hydrated asynchronously.
    assert.equal(
      await app.client.waitFor(`document.querySelector('#analysis-name')?.value === 'Rich analysis'`, { timeoutMs: 20000 }),
      true,
      'the stored analysis was not restored'
    )

    const state = await app.client.evaluate(`
      return {
        hash: location.hash,
        back: document.querySelector("[data-testid='back-to-history']")?.textContent?.trim() || '',
        title: document.querySelector("[data-testid='workspace-title']")?.textContent?.trim() || '',
        name: document.querySelector('#analysis-name')?.value ?? null,
      }
    `)
    assert.match(state.hash, /^#\/history\//)
    assert.match(state.back, /back to history/i)
    assert.equal(state.title, 'Rich analysis')
    assert.equal(state.name, 'Rich analysis')
  })

  it('reads the metadata file contents to populate the condition selects', async () => {
    assert.equal(
      await app.client.waitFor(`document.querySelector('#conditionA').options.length > 2`, { timeoutMs: 25000 }),
      true,
      'the condition selects were not populated from the metadata file'
    )
    const conditions = await app.client.evaluate(`
      const options = (id) => [...document.querySelector(id).options].map((option) => option.value).filter(Boolean)
      return {
        a: options('#conditionA'),
        b: options('#conditionB'),
        valueA: document.querySelector('#conditionA').value,
        valueB: document.querySelector('#conditionB').value,
        files: [...document.querySelectorAll("[class*='fileName']")].map((node) => node.textContent),
      }
    `)
    assert.deepEqual(conditions.a.sort(), [CONDITION_A, CONDITION_B].sort())
    assert.equal(conditions.valueA, CONDITION_A, 'the stored condition A must be restored, not pruned')
    assert.equal(conditions.valueB, CONDITION_B, 'the stored condition B must be restored, not pruned')
    assert.ok(conditions.files.length >= 2, 'the input file names should be listed')
  })

  it('leaves a stored analysis fully editable', async () => {
    const editable = await app.client.evaluate(`
      const name = document.querySelector('#analysis-name')
      const method = document.querySelector('#correlationMethod')
      const rerun = [...document.querySelectorAll('button')].find((button) => /Re-run analysis/i.test(button.textContent))
      return {
        nameDisabled: name.disabled,
        methodDisabled: method.disabled,
        hasRerun: Boolean(rerun) && !rerun.disabled,
        hasSave: [...document.querySelectorAll('button')].some((button) => /save changes/i.test(button.textContent)),
      }
    `)
    assert.equal(editable.nameDisabled, false)
    assert.equal(editable.methodDisabled, false)
    assert.equal(editable.hasRerun, true)
    assert.equal(editable.hasSave, true)
  })

  it('warns when the stored input files are no longer on disk', async () => {
    await app.client.evaluate(`location.hash = '#/history/${IDS.missingFiles}'; return true`)
    assert.equal(
      await app.client.waitFor(`/[0-9]+ files? from this analysis could not be found/i.test(document.body.innerText)`, { timeoutMs: 25000 }),
      true,
      'no warning about the missing files'
    )
  })

  it('offers no way to create a new analysis from the history page', async () => {
    const offenders = await app.client.evaluate(`
      const scope = document.querySelector('main') || document.body
      return [...scope.querySelectorAll('a, button')]
        .map((element) => element.textContent.trim().toLowerCase())
        .filter((text) => text.includes('new analysis'))
    `)
    assert.deepEqual(offenders, [])
  })

  it('shows no export buttons in the analysis interface', async () => {
    const exports = await app.client.evaluate(`
      return [...document.querySelectorAll('button')]
        .map((button) => button.textContent.trim().toLowerCase())
        .filter((text) => text.includes('export'))
    `)
    assert.deepEqual(exports, [])
  })

  describe('comparison', () => {
    before(async () => {
      await app.client.evaluate(`location.hash = '#/history?compare=${IDS.rich},${IDS.plain}'; return true`)
      assert.equal(
        await app.client.waitFor(`document.querySelectorAll("[data-testid='graph-select-a-conditions']").length === 1`, { timeoutMs: 30000 }),
        true
      )
      await new Promise((resolve) => setTimeout(resolve, 500))
    })

    it('renders two columns with distinct tints and the same height', async () => {
      const geometry = await app.client.evaluate(`
        const a = document.querySelector("[data-testid='compare-column-a']")
        const b = document.querySelector("[data-testid='compare-column-b']")
        const styleA = getComputedStyle(a)
        const styleB = getComputedStyle(b)
        const status = document.querySelector("[data-testid='compare-status-a']")
        const open = document.querySelector("[data-testid='compare-open-a']")
        const columnBox = a.getBoundingClientRect()
        const statusBox = status.getBoundingClientRect()
        const openBox = open.getBoundingClientRect()
        return {
          heightA: Math.round(columnBox.height),
          heightB: Math.round(b.getBoundingClientRect().height),
          differentTint: styleA.backgroundColor !== styleB.backgroundColor && styleA.borderLeftColor !== styleB.borderLeftColor,
          statusBeforeOpen: statusBox.left < openBox.left,
          statusOnSameRow: Math.abs(statusBox.top - openBox.top) < 4,
          statusAtTopRight: Math.round(columnBox.right - statusBox.right) < 160 && statusBox.top - columnBox.top < 40,
          statusText: status.textContent.trim(),
        }
      `)
      assert.equal(geometry.differentTint, true, 'the two analyses must keep distinct tints')
      assert.equal(geometry.heightA, geometry.heightB, 'the two columns must have the same height')
      assert.equal(geometry.statusBeforeOpen, true)
      assert.equal(geometry.statusOnSameRow, true)
      assert.equal(geometry.statusAtTopRight, true)
      assert.equal(geometry.statusText, 'completed')
    })

    it('gives condition graphs and clusters their own fixed-size area', async () => {
      const areas = await app.client.evaluate(`
        const height = (id) => Math.round(document.querySelector("[data-testid='" + id + "']").getBoundingClientRect().height)
        const plots = (id) => document.querySelectorAll("[data-testid='" + id + "'] svg[role='img']").length
        return {
          panels: ['conditions', 'clusters'].map((kind) => Boolean(document.querySelector("[data-testid='graph-panel-a-" + kind + "']"))),
          selectCount: document.querySelectorAll("[data-testid='compare-column-a'] select").length,
          heights: [
            height('graph-viewport-a-conditions'), height('graph-viewport-a-clusters'),
            height('graph-viewport-b-conditions'), height('graph-viewport-b-clusters'),
          ],
          plotsA: plots('graph-viewport-a-conditions') + plots('graph-viewport-a-clusters'),
        }
      `)
      assert.deepEqual(areas.panels, [true, true], 'condition graphs and clusters need separate areas')
      assert.equal(areas.selectCount, 2, 'each area uses exactly one select switcher')
      assert.equal(new Set(areas.heights).size, 1, `the graph areas must have one fixed height: ${areas.heights}`)
      assert.ok(areas.heights[0] >= 250)
      assert.equal(areas.plotsA, 2, 'exactly one graph per area must be rendered')
    })

    it('shows the selected option in full and switches the rendered graph', async () => {
      const info = await app.client.evaluate(`
        ${SET_SELECT_HELPER}
        const select = document.querySelector("[data-testid='graph-select-a-conditions']")
        const label = (element) => element.options[element.selectedIndex].textContent
        const nodeSet = () => [...document.querySelector("[data-testid='graph-viewport-a-conditions'] svg[role='img']")
          .querySelectorAll('title')].map((title) => title.textContent).sort().join('|')
        const before = nodeSet()
        const beforeLabel = label(select)
        setSelectValue(select, select.options[select.options.length - 1].value)
        await new Promise((resolve) => setTimeout(resolve, 400))
        return {
          optionCount: select.options.length,
          truncated: select.scrollWidth > select.clientWidth + 1,
          beforeLabel,
          afterLabel: label(select),
          changed: before !== nodeSet(),
          singlePlot: document.querySelectorAll("[data-testid='graph-viewport-a-conditions'] svg[role='img']").length === 1,
        }
      `)
      assert.ok(info.optionCount >= 2, 'both conditions should be offered')
      assert.equal(info.truncated, false, 'the selected label must be fully visible')
      assert.equal(info.singlePlot, true)
      assert.equal(info.changed, true, `the graph did not change (${info.beforeLabel} -> ${info.afterLabel})`)
    })

    it('keeps the graph details collapsed until asked for', async () => {
      const details = await app.client.evaluate(`
        const node = document.querySelector("[data-testid='graph-details-a-clusters']")
        const openBefore = node.open
        node.open = true
        await new Promise((resolve) => setTimeout(resolve, 120))
        return { openBefore, openAfter: node.open, text: node.innerText }
      `)
      assert.equal(details.openBefore, false, 'the details must start collapsed')
      assert.equal(details.openAfter, true)
      assert.ok(details.text.length > 0)
    })

    it('enlarges a graph in an overlay and closes it with Escape', async () => {
      const overlay = await app.client.evaluate(`
        const inline = document.querySelector("[data-testid='graph-viewport-a-conditions'] svg[role='img']")
        const inlineWidth = Math.round(inline.getBoundingClientRect().width)
        document.querySelector("[data-testid='graph-expand-a-conditions']").click()
        await new Promise((resolve) => setTimeout(resolve, 300))
        const modal = document.querySelector("[data-testid='graph-overlay-a-conditions']")
        if (!modal) return { opened: false }
        const enlarged = modal.querySelector("svg[role='img']")
        return { opened: true, inlineWidth, enlargedWidth: Math.round(enlarged.getBoundingClientRect().width) }
      `)
      assert.equal(overlay.opened, true, 'the expand button did not open the overlay')
      assert.ok(overlay.enlargedWidth > overlay.inlineWidth, `enlarged graph is ${overlay.enlargedWidth}px vs ${overlay.inlineWidth}px`)

      assert.equal(
        await app.client.evaluate(`
          window.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }))
          await new Promise((resolve) => setTimeout(resolve, 300))
          return document.querySelector("[data-testid='graph-overlay-a-conditions']") === null
        `),
        true
      )
    })

    it('shows no difference counter and no correlation-value table', async () => {
      const found = await app.client.evaluate(`
        return {
          counter: document.querySelector("[data-testid='comparison-differences']") !== null,
          counterText: /All differences/i.test(document.body.innerText),
          correlationTable: /Correlation table/i.test(document.body.innerText),
        }
      `)
      assert.equal(found.counter, false)
      assert.equal(found.counterText, false)
      assert.equal(found.correlationTable, false)
    })

    it('scrolls the comparison without moving either sidebar', async () => {
      const result = await app.client.evaluate(`
        const nav = document.querySelector('nav')
        const list = document.querySelector("[data-testid='history-list']").closest('aside')
        const detail = document.querySelector("[data-testid='history-detail']")
        const navTop = nav.getBoundingClientRect().top
        const listTop = list.getBoundingClientRect().top
        detail.scrollTop = detail.scrollHeight
        await new Promise((resolve) => setTimeout(resolve, 200))
        return {
          scrolled: detail.scrollTop > 0 || detail.scrollHeight <= detail.clientHeight,
          navMoved: Math.abs(nav.getBoundingClientRect().top - navTop) > 1,
          listMoved: Math.abs(list.getBoundingClientRect().top - listTop) > 1,
        }
      `)
      assert.equal(result.navMoved, false, 'the main sidebar moved')
      assert.equal(result.listMoved, false, 'the history list moved')
    })
  })
})
