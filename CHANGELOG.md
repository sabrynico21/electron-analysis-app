# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

- **SPIEC-EASI backend replaced.** The standalone Python library in
  `spiec-easi-python/` was rebuilt from scratch as a faithful port of the R
  reference (`SpiecEasi` 1.99.0 + `huge` 2.0.1 + `pulsar` 0.3.13), `method="mb"`
  only. The previous implementation — with its pybind11 `_mb_core` accelerator,
  joblib `prefer="threads"` backend and `spieceasi.mb_solver` /
  `spieceasi.runtime` helpers — was removed. `python/analysis.py` now calls the
  new public API `spieceasi.spiec_easi(...)`.
- **Worker count.** The StARS subsample loop now runs in a process pool
  (`spieceasiJobs`); the subsamples are independent and aggregated in a fixed
  order, so results are identical for any worker count. The joblib
  `threads`-vs-`processes` backend choice and the compiled core are gone.
- **`spieceasiMaxIter` is no longer used.** The port follows R's fixed MB solver
  iteration limit (`MAX_ITER = 10000`); the previous per-case cap was a deviation
  from R and has been removed. The `spieceasiNLambdas` / `spieceasiNSubsamples`
  reductions in the tiny/high-dimensional regime (`n_samples <= 10` and
  `n_taxa >= 200`) are kept.

### Removed

- References to the deleted library internals (`spieceasi._mb_core`,
  `spiec-easi-python/docs/AUDIT.md`, `spiec-easi-python/docs/PERFORMANCE.md`,
  the extension-vs-reference sweep) from `README.md`, `CHANGELOG.md` and
  `python/backend.spec`.

### Notes

- Equivalence with R is documented in `spiec-easi-python/docs/EQUIVALENCE.md`
  (what is identical, what is identical within a numerical tolerance, what cannot
  be); deviations are in `spiec-easi-python/docs/DIFFERENCES.md`.
- No comparison tolerance was widened to make a test pass.
