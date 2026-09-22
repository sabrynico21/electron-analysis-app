import { useEffect } from 'react'
import { Routes, Route } from 'react-router-dom'
import Sidebar from './components/shared/Sidebar'
import AnalysisPage from './pages/AnalysisPage'
import ResultsPage from './pages/ResultsPage'
import HistoryPage from './pages/HistoryPage'
import SettingsPage from './pages/SettingsPage'
import styles from './styles/App.module.css'

export default function App() {
  // A file dropped outside an upload zone would otherwise make Electron navigate
  // the window to that file and lose the whole application. Blocks the default
  // for the entire window; the upload zones handle their own drops.
  useEffect(() => {
    const block = (event) => event.preventDefault()
    window.addEventListener('dragover', block)
    window.addEventListener('drop', block)
    return () => {
      window.removeEventListener('dragover', block)
      window.removeEventListener('drop', block)
    }
  }, [])

  return (
    <div className={styles.appLayout}>
      <Sidebar />
      <main className={styles.content}>
        <Routes>
          {/* Reserved for creating brand new analyses: the form always starts empty. */}
          <Route path="/" element={<AnalysisPage />} />
          <Route path="/history" element={<HistoryPage />} />
          {/* Reopens a stored analysis next to the list, with a back button on top. */}
          <Route path="/history/:jobId" element={<HistoryPage />} />
          <Route path="/results/:jobId" element={<ResultsPage />} />
          <Route path="/settings" element={<SettingsPage />} />
        </Routes>
      </main>
    </div>
  )
}
