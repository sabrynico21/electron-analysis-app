const { app, BrowserWindow, dialog } = require('electron')
const path = require('path')
const fs = require('fs')
const log = require('electron-log')
const { setupIpcHandlers } = require('./ipc/handlers')
const { createMenu } = require('./menu')

log.initialize({ preload: true })

let mainWindow

const DEV_SERVER_URL = process.env.VITE_DEV_SERVER_URL || 'http://localhost:5173'

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1280,
    height: 800,
    minWidth: 900,
    minHeight: 600,
    icon: path.join(__dirname, '../../resources/icons/icon.png'),
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
    },
  })

  // The Vite dev server is used only when development mode is explicitly
  // requested (`npm run dev`). `npm start` and packaged builds load the built
  // renderer bundle from dist/renderer.
  if (!app.isPackaged && process.env.NODE_ENV === 'development') {
    mainWindow.loadURL(DEV_SERVER_URL)
    mainWindow.webContents.openDevTools()
  } else {
    const indexPath = path.join(__dirname, '../../dist/renderer/index.html')
    if (fs.existsSync(indexPath)) {
      mainWindow.loadFile(indexPath)
    } else {
      // Fail loudly: a missing renderer bundle must never show a blank window.
      log.error(`Renderer bundle not found at ${indexPath}`)
      dialog.showErrorBox(
        'Installation incomplete',
        `The application interface could not be found:\n${indexPath}\n\n` +
          'This build is missing the compiled renderer. Rebuild the app with "npm run build".'
      )
    }
  }

  createMenu(mainWindow)
}

app.whenReady().then(() => {
  createWindow()
  setupIpcHandlers()
  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow()
  })
})

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit()
})
