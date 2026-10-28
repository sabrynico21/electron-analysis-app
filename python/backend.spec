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
REPO_ROOT = PYTHON_DIR.parent
CLUSTERING_DIR = PYTHON_DIR / "clustering"
# SPIEC-EASI lives in the standalone library so the app and the library share a
# single source instead of a duplicated copy.
SPIEC_EASI_SRC = REPO_ROOT / "spiec-easi-python" / "src"
CCLASSO_SRC = REPO_ROOT / "cclasso-python" / "src"
SPRING_SRC = REPO_ROOT / "spring-python" / "src"
SPRING_IPOL = SPRING_SRC / "spring_microbiome" / "latentcor" / "data" / "ipol_22.npz"
THIRD_PARTY_LICENSES = [
    (str(REPO_ROOT / "cclasso-python" / "LICENSE"), "licenses/cclasso"),
    (str(REPO_ROOT / "spiec-easi-python" / "LICENSE"), "licenses/spieceasi"),
    (str(REPO_ROOT / "spring-python" / "LICENSE"), "licenses/spring"),
]

analysis = Analysis(
    [str(PYTHON_DIR / "backend.py")],
    pathex=[
        str(PYTHON_DIR),
        str(CLUSTERING_DIR),
        str(SPIEC_EASI_SRC),
        str(CCLASSO_SRC),
        str(SPRING_SRC),
    ],
    binaries=[],
    # SPRING's approximate latentcor path loads its pinned interpolation table
    # through importlib.resources at runtime.
    datas=[
        (str(SPRING_IPOL), "spring_microbiome/latentcor/data"),
        *THIRD_PARTY_LICENSES,
    ],
    # The clustering package is imported at function level from backend.py and
    # uses flat imports (`from bacteria_graph import ...`), so make sure the
    # modules are collected even though they are not top-level imports.
    hiddenimports=[
        "bacteria_graph",
        "bacteria_pagerank",
        "pagerank_nibble",
        "selection_ladder",
        "spieceasi",
        "spieceasi._validation",
        "spieceasi.cli",
        "spieceasi.lambda_path",
        "spieceasi.mb",
        "spieceasi.pipeline",
        "spieceasi.preprocess",
        "spieceasi.stars",
        "spieceasi.transforms",
        "cclasso",
        "cclasso._core",
        "cclasso.validation",
        "spring_microbiome",
        "spring_microbiome.spring",
        "spring_microbiome.mclr",
        "spring_microbiome.huge_kmb",
        "spring_microbiome.pulsar",
        "spring_microbiome.latentcor",
        "spring_microbiome.latentcor.latentcor",
        "spring_microbiome.latentcor.bridge",
        "spring_microbiome.latentcor.interp",
        "spring_microbiome.latentcor.kendall",
        "spring_microbiome.latentcor.nearpd",
        "spring_microbiome.huge.mb",
    ],
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
