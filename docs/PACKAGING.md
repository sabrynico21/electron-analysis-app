# Packaging — compiled MB core, wheels and the frozen backend

Two artifacts ship this project:

1. the **Python package** `spiec-easi-python` (a wheel containing the optional
   compiled MB core `spieceasi._mb_core`), and
2. the **frozen backend** `analysis-backend`, a PyInstaller bundle that the
   Electron app ships under `resources/python-dist` so end users need no Python.

## 1. The compiled MB core

`src/spieceasi/_mb_core.cpp` is a near-literal translation of `huge_core.cpp`
`mb` / `mb_partial_residual` / `mb_refine_active`. Only the solver is compiled;
CLR, correlation, lambda path, StARS, selection and refit stay in Python and are
shared by both implementations.

```bash
cd spiec-easi-python
python -m pip install --no-build-isolation -e .   # needs a C++17 compiler + pybind11
python -c "import spieceasi.mb_solver as m; print(m.active_solver())"   # -> cpp
```

Compile flags (`setup.py`): `-O2 -ffp-contract=off -fno-lto` on POSIX
(`/O2 /fp:precise` on MSVC); **no** `-ffast-math`, **no** `-march=native`. They
deliberately stay close to R's own build; `docs/PERFORMANCE.md` records the
flags R uses and the measured effect of `-O3` (0–1.2 %, not adopted) and of the
row-major S layout (adopted, bit-identical).

### Fallback — the extension is optional

If `spieceasi._mb_core` cannot be imported (no compiler, unsupported platform,
wheel without the binary), `spieceasi.mb_solver` falls back to the pure-Python
reference solver in `mb_solver_reference.py`, emits a `RuntimeWarning` naming the
import error, and `active_solver()` returns `"python"`. The desktop app logs

```
SPIEC-EASI solver: python
WARNING: the compiled MB core (spieceasi._mb_core) is not available in this
build; falling back to the pure-Python reference solver, which is roughly two
orders of magnitude slower.
```

so a silent performance regression is impossible. A CI job removes the compiled
module after installing and re-runs the determinism tests on that path.

## 2. Wheels

`[tool.cibuildwheel]` in `pyproject.toml` builds CPython 3.10–3.14 wheels for
Linux, macOS and Windows. Its `test-command` asserts that
`spieceasi.mb_solver.HAVE_CPP_SOLVER` is true **inside the built wheel**, so a
wheel silently missing the binary fails the build. macOS wheels are repaired
with `delocate-wheel --require-archs {delocate_archs}`.

`.github/workflows/ci.yml` runs cibuildwheel on `ubuntu-latest`, `macos-13`
(Intel), `macos-14` (Apple silicon) and `windows-latest`. PyInstaller cannot
cross-compile and neither can cibuildwheel for macOS/Windows, so those wheels are
produced by CI only — they cannot be built from this Linux development machine.

## 3. Frozen backend (PyInstaller)

```bash
npm run build:backend          # or: node scripts/build-backend.mjs [--target win32]
```

* **One directory, not one file** (`EXE(exclude_binaries=True)` + `COLLECT`): the
  one-file variant unpacks itself into a temporary directory at every launch,
  which costs seconds and trips antivirus scanners.
* `python/backend.spec` lists `spieceasi.mb_solver`, `spieceasi.mb_solver_reference`,
  `spieceasi._mb_core` and `spieceasi.runtime` in `hiddenimports`, so the compiled
  core is packaged explicitly instead of being discovered by luck.
* Regression tests (slow, skipped when the bundle is absent):
  * `test_frozen_backend_runs_spieceasi_with_app_defaults` — the packaged binary
    runs a SPIEC-EASI analysis with the app defaults, exit code 0, the thread
    backend (never a loky process pool, which would re-execute the binary with
    `-B`), and a payload byte-identical to `n_jobs=1`;
  * `test_frozen_backend_uses_the_compiled_mb_core` — the same run must log
    `SPIEC-EASI solver: cpp`, i.e. the bundle ships *and uses* the compiled core
    rather than the fallback.

## 4. OpenMP (not implemented)

An OpenMP-parallel loop over the `d` nodewise regressions is the natural next
step: R's `huge` already does it and gains ~2.1× at p=100 (`docs/PERFORMANCE.md`
has the measurement). It is **not** implemented, and no code depends on it. The
proposal — bit-identity gates, per-platform runtime handling and the
anti-oversubscription rule — is written up under "Future work: OpenMP" in
`docs/PERFORMANCE.md`.

## 5. Startup time of the frozen backend

See `docs/PERFORMANCE.md` for the measured breakdown (import profile, one-folder
launch, effect of the import changes). Summary: the frozen executable needs about
+5.3 s over the source run before any computation starts; the decomposition and
what of it was removed is documented there.
