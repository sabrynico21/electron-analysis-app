# Changelog

All notable changes to this project are documented here.
This project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Step 0 porting analysis (`docs/PORTING_NOTES.md`).
- Repository skeleton: `src/cclasso`, `tests/`, `docs/`, CI workflow.
- R golden-fixture generator (`tests/r_reference/generate_golden.R`) + 8 committed
  cases under `tests/golden/` (with `MANIFEST.sha256`).
- Core port: `cclasso`, `cclasso_sub`, `cv_loss_cclasso`, `boot_cclasso`
  (line-mapped to `R/cclasso.R`).
- Input validation (`cclasso.validation`), CLI (`python -m cclasso` / `cclasso`).
- Tests: equivalence vs R, mathematical properties, CLI contract, validation.
- Docs: `DIFFERENCES.md`, `EQUIVALENCE.md`, `PACKAGING.md`.
- Full `LICENSE` (LGPL-3.0-or-later): official LGPLv3 + GPLv3 texts from gnu.org.
- Four synthetic R golden cases closing the coverage gap: `counts=False`
  (`synth_frac_*`) and non-default `k_cv`/`n_boot`/`lam_int` (`synth_*`).
- QMP pinned at `tests/r_reference/data/QMP.rda` so fixture regeneration works
  from a clean checkout.
- Dataset provenance fields (`source_package`, version, license, citation,
  retrieval date, `redistribution_cleared`) stamped into every public
  `params.json`; audit in `../docs/DATASET_PROVENANCE.md`.
- `VENDORED.md` and `VENDOR_HISTORY.bundle`: the package is now vendored into
  the root repository (no longer an untracked nested Git repo).

### Verified
- R↔Python equivalence: `lams`/`lambda` **bit-identical**, ADMM iteration counts
  **exact**, outputs within ~1e-13 (contract 1e-6). See `docs/EQUIVALENCE.md`.
- 49 tests pass; the equivalence suite now also runs the synthetic
  `counts=False` / varied-parameter cases.

### Notes
- Ported from the R implementation **v2.0** (`CCLasso-master/R/cclasso.R`).
- `v1.0/` and SparCC are out of scope.

