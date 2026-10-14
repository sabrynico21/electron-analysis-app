# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

- **Default joblib backend for StARS: `threads` (was `processes`).** The compiled
  MB core releases the GIL and the subsamples are aggregated in a fixed order, so
  a thread pool is both safe and faster (2.2×–3.2× at 4–8 workers on the
  benchmark cases). Results are bit-identical for any worker count; `prefer`
  still accepts `"processes"` for callers that need it. This also removes the
  failure mode that broke frozen executables (loky re-executing the binary).
- **Default worker count now comes from the physical cores**, not from
  `min(4, cpu_count)`: `spieceasi.runtime.default_n_jobs()` reads the physical
  core count where the platform exposes it, falls back to `logical // 2`
  (assume SMT) and then to 1. The desktop app uses it; `spieceasiJobs` (or
  `n_jobs`) overrides it. On the development machine (4 physical / 8 logical) the
  p=100 pipeline goes from 21.6 s serial to 8.7 s at 4 workers.
- **Build backend.** `hatchling` → `setuptools` + `pybind11`, so the compiled MB
  core can be built. `[tool.cibuildwheel]` builds wheels for CPython 3.10–3.14 on
  Linux, macOS and Windows; its `test-command` fails if the compiled core is
  missing from the wheel.
- **SPIEC-EASI / StARS subsample fraction (application default).**
  `python/analysis.py` no longer hard-codes `subsample_fraction=0.8`. When the
  user does not specify a fraction, the library now follows `pulsar::.ratcheck`
  exactly, as the R SpiecEasi reference does:
  - `n_samples <= 144` → ratio `0.8`, size `floor(0.8 * n)`;
  - `n_samples > 144` → ratio `10*sqrt(n)/n`, size `floor(10*sqrt(n))`.

  **Behaviour change:** for datasets with more than 144 samples the StARS
  subsamples are now *smaller* than before (e.g. `n = 200` → 141 instead of 160;
  `n = 300` → 173 instead of 240). Results for such datasets may therefore differ
  from previous releases. An explicit value can still be forced through the
  `spieceasiSubsampleFraction` parameter, which is passed through unchanged.

### Added

- **Compiled MB core (`spieceasi._mb_core`).** A pybind11 extension translating
  `huge_core.cpp` `mb` / `mb_partial_residual` / `mb_refine_active`. It is an
  accelerator, not a requirement: `spieceasi.mb_solver_reference` keeps the
  pure-Python solver as the reference implementation and as the automatic
  fallback (a `RuntimeWarning` is emitted when the compiled module cannot be
  imported). See `docs/AUDIT.md` §7 for the declared deviations.

  **Speed-up on the measured cases** (serial, median of 3, full pipeline):
  97.8 s → 0.57 s (`ds2`, p=30), 453.4 s → 3.99 s (p=50), 1565.8 s → 22.06 s
  (p=100); versus R that is 0.27× / 0.90× / 1.57×. The non-slow test suite goes
  from 662 s to 18.66 s. **The measured outputs are unchanged**: selected lambda
  and edge count are identical to the pure-Python solver and to R on all
  benchmark cases.

  **Accuracy declaration.** On the tested cases the extension is *not*
  bit-identical to the pure-Python reference nor to R: the maximum coefficient
  difference is 9.470e-14 on 244 random problems (1.652e-13 versus the released
  goldens) and the discovered support is identical in every non-degenerate case.
  In a deliberately degenerate family (relative collinearity 1e-9) the lasso
  solution is not unique, the two implementations and R land in different points
  of the same solution set (max difference 1.881e-01, 25/244 problems, 40 lambda
  blocks), and R does not single out either of them. This is reported as-is; no
  tolerance was widened to hide it.
- `spieceasi.stars.infer_subsample_size`: public helper that resolves the StARS
  subsample size with the `pulsar::.ratcheck` rule, now used internally and
  covered by tests.
- `spieceasi.runtime`: `physical_core_count()` and `default_n_jobs()`, used by the
  desktop app for the SPIEC-EASI worker count.
- Golden fixtures are content-addressed: `golden_extended/CHECKSUMS.sha256`
  records the SHA-256 of the *decompressed* content of every fixture, and
  `solver_adj_path.tsv` is stored gzipped (44 MB → 2.0 MB, 22.2×). The data is
  byte-identical — the tests load `.gz` transparently and `tests/test_golden_checksums.py`
  fails if any content changes.
- Extended R/Python equivalence golden set (`tests/golden_extended/`, 23 cases)
  with `tests/test_equivalence_extended.py`.
- `docs/PERFORMANCE.md` with measured MB-solver timings and known performance
  limitations.

### Fixed

- **The packaged (PyInstaller) backend crashed on every SPIEC-EASI run.**
  `stars_select` called `joblib.Parallel(n_jobs=…)` with the default *loky*
  (process) backend. Inside a frozen executable loky re-executes the binary with
  `-B`, which the CLI rejected (`unknown command: -B`, exit code 2). The app never
  sends `spieceasiJobs`, so it always used `min(4, cpu_count)` and therefore
  always crashed. Fixed by adding an explicit `prefer` argument to the StARS
  subset loop (`"processes"` by default, so normal installs are unchanged) and by
  selecting `prefer="threads"` with `n_jobs=1` for frozen builds. `backend.py`
  additionally calls `multiprocessing.freeze_support()` and reports
  interpreter-style arguments explicitly instead of the generic
  "unknown command".

### Notes

- No golden file was modified and no comparison tolerance was widened.
