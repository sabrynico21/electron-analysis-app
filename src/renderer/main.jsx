import React from 'react'
import ReactDOM from 'react-dom/client'
import { HashRouter } from 'react-router-dom'
import App from './App'
import './styles/global.css'

// The application always opens on the "Create new analysis" page with an empty
// form. A hash left over from a previous session (or a dev-server reload) must
// not bring the user back to a half-filled workspace.
if (window.location.hash && window.location.hash !== '#/') {
  try {
    window.history.replaceState(null, '', '#/')
  } catch {
    // `file://` URLs can reject replaceState; a plain assignment always works.
    window.location.hash = '#/'
  }
}

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <HashRouter>
      <App />
    </HashRouter>
  </React.StrictMode>
)
