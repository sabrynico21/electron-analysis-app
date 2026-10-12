import { useState, useEffect } from 'react'
import styles from '../styles/SettingsPage.module.css'

/**
 * Plain-language labels for how FastSpar was located. The raw values coming from
 * the main process are technical, so a non-technical user gets a friendly phrase.
 */
const FASTSPAR_SOURCE_LABELS = {
  settings: 'selected by you',
  bundled: 'included with this app',
  env: 'set by the FASTSPAR_BIN variable',
  path: 'found automatically on this computer',
  none: 'not found',
}

function StatusBadge({ tone, children }) {
  const toneClass = tone === 'ok' ? styles.badgeOk : tone === 'warn' ? styles.badgeWarn : styles.badgeBad
  const icon = tone === 'ok' ? '✓' : tone === 'warn' ? '!' : '✕'
  return (
    <span className={`${styles.badge} ${toneClass}`}>
      <span className={styles.badgeIcon} aria-hidden="true">{icon}</span>
      {children}
    </span>
  )
}

export default function SettingsPage() {
  const [settings, setSettings] = useState({ pythonPath: '', fastsparPath: '', theme: 'dark' })
  const [saved, setSaved] = useState(false)
  const [sparccStatus, setSparccStatus] = useState(null)
  const [checkingSparcc, setCheckingSparcc] = useState(false)
  const [runtimeStatus, setRuntimeStatus] = useState(null)
  const [detectMessage, setDetectMessage] = useState(null)

  useEffect(() => {
    loadEverything()
  }, [])

  const loadEverything = async () => {
    const loaded = await window.electronAPI.getSettings()
    setSettings((prev) => ({ ...prev, ...(loaded || {}) }))
    setRuntimeStatus(await window.electronAPI.getAnalysisRuntimeStatus())
    await refreshSparccStatus()
  }

  const refreshSparccStatus = async () => {
    setCheckingSparcc(true)
    try {
      const status = await window.electronAPI.getSparccStatus()
      setSparccStatus(status)
      return status
    } finally {
      setCheckingSparcc(false)
    }
  }

  const handleSave = async () => {
    await window.electronAPI.saveSettings(settings)
    setSaved(true)
    setTimeout(() => setSaved(false), 2500)
    setRuntimeStatus(await window.electronAPI.getAnalysisRuntimeStatus())
    await refreshSparccStatus()
  }

  const handleBrowse = async (key) => {
    const result = await window.electronAPI.openFileDialog({ properties: ['openFile'] })
    if (!result.canceled) setSettings((s) => ({ ...s, [key]: result.filePaths[0] }))
  }

  /**
   * Look for FastSpar on this computer and, when found, remember its location in
   * the override field. When nothing is found the analysis still works — it just
   * uses Spearman instead of SparCC.
   */
  const handleAutoDetect = async () => {
    setDetectMessage(null)
    const status = await refreshSparccStatus()
    if (status?.available && status.path) {
      setSettings((s) => ({ ...s, fastsparPath: status.path }))
      setDetectMessage({ tone: 'ok', text: `FastSpar found and selected: ${status.path}` })
    } else {
      setDetectMessage({
        tone: 'warn',
        text: 'FastSpar was not found on this computer. SparCC analyses will automatically run with Spearman instead — everything else keeps working.',
      })
    }
  }

  const engineReady = runtimeStatus?.kind === 'bundled'
  const sparccReady = !!sparccStatus?.available
  const everythingReady = engineReady && sparccReady

  return (
    <div className={styles.page}>
      <header className={styles.header}>
        <h1>Settings</h1>
        <p className={styles.helpText}>
          Most people never need to change anything here. Use this page only if one of the checks
          below shows a warning.
        </p>
      </header>

      {/* At-a-glance summary */}
      <section className={`${styles.summary} ${everythingReady ? styles.summaryOk : styles.summaryWarn}`}>
        <div className={styles.summaryIcon} aria-hidden="true">{everythingReady ? '✓' : '!'}</div>
        <div className={styles.summaryBody}>
          <h2 className={styles.summaryTitle}>
            {everythingReady ? 'Everything is ready' : 'The app is ready, with one optional note'}
          </h2>
          <p className={styles.summaryText}>
            {engineReady
              ? 'This app includes its own analysis engine, so nothing extra has to be installed.'
              : 'This copy runs analyses with the Python already installed on this computer. It works fine; the official installer for your system would not need Python at all.'}
            {!sparccReady && ' SparCC is optional and currently unavailable — the app will use Spearman instead, which needs no extra software.'}
          </p>
        </div>
      </section>

      {/* Analysis engine */}
      <section className={styles.card}>
        <div className={styles.cardHead}>
          <div>
            <h2 className={styles.cardTitle}>Analysis engine</h2>
            <p className={styles.cardHint}>The part of the app that performs the statistical calculations.</p>
          </div>
          <StatusBadge tone={engineReady ? 'ok' : 'warn'}>
            {engineReady ? 'Built-in engine' : 'Uses Python on this PC'}
          </StatusBadge>
        </div>

        <div className={styles.detailList}>
          <div className={styles.detailRow}>
            <span className={styles.detailLabel}>How analyses run</span>
            <span className={styles.detailValue}>
              {engineReady
                ? 'Together with the app (nothing to install)'
                : 'Using the Python interpreter found on this computer'}
            </span>
          </div>
          <div className={styles.detailRow}>
            <span className={styles.detailLabel}>Computer type</span>
            <span className={styles.detailValue}>
              {runtimeStatus ? `${runtimeStatus.platform} · ${runtimeStatus.arch}` : 'Checking…'}
            </span>
          </div>
          <div className={styles.detailRow}>
            <span className={styles.detailLabel}>Engine location</span>
            <span className={styles.detailValue}>
              <code className={styles.path}>
                {runtimeStatus?.bundledPath || 'not present in this copy of the app'}
              </code>
            </span>
          </div>
        </div>

        {!engineReady && (
          <p className={styles.note}>
            This is normal when the app is started from its source code, or when the engine was
            built for a different operating system. Installed versions of the app always include
            the right engine for your computer.
          </p>
        )}
      </section>

      {/* SparCC / FastSpar */}
      <section className={styles.card}>
        <div className={styles.cardHead}>
          <div>
            <h2 className={styles.cardTitle}>SparCC correlation</h2>
            <p className={styles.cardHint}>
              An extra correlation method that needs a small helper program called FastSpar.
            </p>
          </div>
          <StatusBadge tone={sparccReady ? 'ok' : 'warn'}>
            {sparccReady ? 'Available' : 'Not available'}
          </StatusBadge>
        </div>

        <p className={styles.note}>
          FastSpar is included on Windows and Linux (64-bit). If it is missing, the app simply uses
          the standard Spearman method instead — no analysis is ever blocked.
        </p>

        {sparccStatus?.path && (
          <div className={styles.detected}>
            <span className={styles.detectedLabel}>Detected automatically</span>
            <code className={styles.path}>{sparccStatus.path}</code>
            <span className={styles.detectedSource}>
              {FASTSPAR_SOURCE_LABELS[sparccStatus.source] || sparccStatus.source}
            </span>
          </div>
        )}

        <div className={styles.actions}>
          <button type="button" className={styles.button} onClick={handleAutoDetect} disabled={checkingSparcc}>
            {checkingSparcc ? 'Looking…' : 'Find FastSpar automatically'}
          </button>
          <button type="button" className={styles.buttonGhost} onClick={refreshSparccStatus} disabled={checkingSparcc}>
            Check again
          </button>
        </div>

        {detectMessage && (
          <p className={`${styles.inlineMessage} ${detectMessage.tone === 'ok' ? styles.inlineOk : styles.inlineWarn}`}>
            {detectMessage.text}
          </p>
        )}

        {sparccStatus && !sparccStatus.available && sparccStatus.error && (
          <p className={styles.note}>{sparccStatus.error}</p>
        )}
      </section>

      {/* Advanced (collapsed by default) */}
      <details className={styles.advanced}>
        <summary className={styles.advancedSummary}>Advanced options</summary>
        <p className={styles.cardHint}>
          Leave these empty unless support asks you to change them. Empty means “find it automatically”.
        </p>

        <div className={styles.field}>
          <label className={styles.fieldLabel} htmlFor="pythonPath">Python interpreter</label>
          <p className={styles.fieldHint}>
            Only used when the built-in engine is missing. If empty, the app finds Python on its own.
          </p>
          <div className={styles.fieldRow}>
            <input
              id="pythonPath"
              value={settings.pythonPath || ''}
              onChange={(e) => setSettings((s) => ({ ...s, pythonPath: e.target.value }))}
              placeholder="Automatic"
            />
            <button type="button" className={styles.buttonGhost} onClick={() => handleBrowse('pythonPath')}>
              Browse…
            </button>
          </div>
        </div>

        <div className={styles.field}>
          <label className={styles.fieldLabel} htmlFor="fastsparPath">FastSpar program</label>
          <p className={styles.fieldHint}>
            Only needed if you installed FastSpar yourself and the automatic detection does not find it.
          </p>
          <div className={styles.fieldRow}>
            <input
              id="fastsparPath"
              value={settings.fastsparPath || ''}
              onChange={(e) => setSettings((s) => ({ ...s, fastsparPath: e.target.value }))}
              placeholder="Automatic"
            />
            <button type="button" className={styles.buttonGhost} onClick={() => handleBrowse('fastsparPath')}>
              Browse…
            </button>
          </div>
        </div>

        {runtimeStatus && !engineReady && runtimeStatus.expectedBundledPath && (
          <div className={styles.field}>
            <span className={styles.fieldLabel}>Expected engine location</span>
            <p className={styles.fieldHint}>
              The built-in engine was looked for here (useful when reporting a problem).
            </p>
            <code className={styles.path}>{runtimeStatus.expectedBundledPath}</code>
          </div>
        )}
      </details>

      <div className={styles.saveRow}>
        <button className={styles.saveButton} onClick={handleSave}>
          {saved ? 'Saved ✓' : 'Save settings'}
        </button>
        {saved && <span className={styles.savedHint}>Your settings have been saved.</span>}
      </div>
    </div>
  )
}
