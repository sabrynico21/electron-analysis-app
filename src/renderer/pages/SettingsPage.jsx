import { useState, useEffect } from 'react'
import styles from '../styles/SettingsPage.module.css'

export default function SettingsPage() {
  const [settings, setSettings] = useState({ pythonPath: '', fastsparPath: '', theme: 'dark' })
  const [saved, setSaved] = useState(false)
  const [sparccStatus, setSparccStatus] = useState(null)
  const [checkingSparcc, setCheckingSparcc] = useState(false)

  useEffect(() => {
    window.electronAPI.getSettings().then((loaded) => {
      setSettings((prev) => ({ ...prev, ...(loaded || {}) }))
    })
    refreshSparccStatus()
  }, [])

  const refreshSparccStatus = async () => {
    setCheckingSparcc(true)
    try {
      const status = await window.electronAPI.getSparccStatus()
      setSparccStatus(status)
    } finally {
      setCheckingSparcc(false)
    }
  }

  const handleSave = async () => {
    await window.electronAPI.saveSettings(settings)
    setSaved(true)
    setTimeout(() => setSaved(false), 2000)
    refreshSparccStatus()
  }

  const handleBrowse = async (key) => {
    const result = await window.electronAPI.openFileDialog({ properties: ['openFile'] })
    if (!result.canceled) setSettings((s) => ({ ...s, [key]: result.filePaths[0] }))
  }

  return (
    <div className={styles.page}>
      <h1>Settings</h1>
      <p className={styles.helpText}>
        This section lets you configure the runtimes used during analyses. If the fields are empty,
        the app tries to auto-detect the binaries from the PATH.
      </p>

      <fieldset className={styles.group}>
        <legend>Python Interpreter</legend>
        <p className={styles.groupHint}>
          Set Python only if you want to force a specific interpreter. If the field is empty, the app
          uses auto-detect (the project virtualenv or Python from the PATH).
        </p>
        <div className={styles.row}>
          <label>Python</label>
          <input
            value={settings.pythonPath || ''}
            onChange={(e) => setSettings((s) => ({ ...s, pythonPath: e.target.value }))}
            placeholder="Auto-detect"
          />
          <button onClick={() => handleBrowse('pythonPath')}>Browse</button>
        </div>
      </fieldset>

      <fieldset className={styles.group}>
        <legend>SparCC Runtime</legend>
        <p className={styles.groupHint}>
          FastSpar is only required when you select SparCC correlation. With Spearman you can leave this field empty.
        </p>
        <div className={styles.row}>
          <label>FastSpar path</label>
          <input
            value={settings.fastsparPath || ''}
            onChange={(e) => setSettings((s) => ({ ...s, fastsparPath: e.target.value }))}
            placeholder="Optional override (auto-detect if empty)"
          />
          <button onClick={() => handleBrowse('fastsparPath')}>Browse</button>
        </div>

        <div className={styles.statusRow}>
          <button type="button" onClick={refreshSparccStatus} disabled={checkingSparcc}>
            {checkingSparcc ? 'Checking...' : 'Check SparCC Runtime'}
          </button>
          {sparccStatus?.available ? (
            <span className={styles.statusOk}>Available</span>
          ) : (
            <span className={styles.statusError}>Not available</span>
          )}
        </div>

        {sparccStatus && (
          <div className={styles.statusDetails}>
            <div><strong>Platform:</strong> {sparccStatus.platform} / {sparccStatus.arch}</div>
            <div><strong>Resolved binary:</strong> {sparccStatus.path || 'none'}</div>
            {!sparccStatus.available && sparccStatus.error && (
              <div><strong>Error:</strong> {sparccStatus.error}</div>
            )}
          </div>
        )}
      </fieldset>

      <button className={styles.saveButton} onClick={handleSave}>
        {saved ? 'Saved ✓' : 'Save Settings'}
      </button>
    </div>
  )
}
