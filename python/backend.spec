# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller build description for the frozen analysis backend.

Produces a self-contained `analysis-backend` executable (one directory, not a
single file, so start-up stays fast and antivirus scanners do not unpack a
huge binary on every launch). The result is shipped through electron-builder's
`extraResources`, which is what lets the packaged app run analyses without any
Python installation on the user's computer.

Build it with `npm run build:backend` rather than calling PyInstaller directly:
the wrapper creates an isolated build virtualenv and installs the requirements.
"""

from pathlib import Path

PYTHON_DIR = Path(SPECPATH).resolve()
CLUSTERING_DIR = PYTHON_DIR / "clustering"

analysis = Analysis(
    [str(PYTHON_DIR / "backend.py")],
    pathex=[str(PYTHON_DIR), str(CLUSTERING_DIR)],
    binaries=[],
    datas=[],
    # The clustering package is imported at function level from backend.py and
    # uses flat imports (`from bacteria_graph import ...`), so make sure the
    # modules are collected even though they are not top-level imports.
    hiddenimports=["bacteria_graph", "bacteria_pagerank", "pagerank_nibble"],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # Nothing in the analysis pipeline needs a GUI toolkit, a notebook stack or
    # a data-frame library; excluding them keeps the payload reasonable.
    excludes=[
        "tkinter",
        "matplotlib",
        "PyQt5",
        "PyQt6",
        "PySide2",
        "PySide6",
        "IPython",
        "jupyter",
        "notebook",
        "pytest",
        "sphinx",
        "pandas",
    ],
    noarchive=False,
)

pyz = PYZ(analysis.pure)

exe = EXE(
    pyz,
    analysis.scripts,
    [],
    exclude_binaries=True,
    name="analysis-backend",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    analysis.binaries,
    analysis.datas,
    strip=False,
    upx=False,
    name="analysis-backend",
)
