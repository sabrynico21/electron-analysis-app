# SPRING stress sweep — 20 Oct 2026, seed 20261020

Predeclared seeded stress suite required by
`docs/EQUIVALENCE_ACCEPTANCE_CONTRACT.md` ("predeclared seeded stress suite of
at least 1,000 small cases per method family, weighted toward sparse, singular,
and threshold-near inputs").

## Command

```bash
cd spring-python
.venv/bin/python tests/r_reference/run_status_stress.py \
    /tmp/spring-stress-1000 --cases 1000 --seed 20261020
```

Runner: `tests/r_reference/run_status_stress.py`; R side:
`tests/r_reference/stress_status_R.R` (sources `SPRING-master/R/helpers.R`).

## Environment

| item | value |
|---|---|
| date | 2026-10-20 |
| R | 4.6.1 |
| Python | 3.14 (venv `spring-python/.venv`) |
| NumPy | see `environment.json` |
| platform | Linux x86-64, serial run |
| settings | `quantitative=FALSE` (mclr), `nlambda=5`, `lambda.min.ratio=1e-2` (explicit default path), `thresh=0.1`, `subsample.ratio=0.8`, `rep.num=4`, `Rmethod="approx"`, `use.nearPD=TRUE`, `nu=0.001`, `ratio=0.9` |

## Predeclared design

Six equally represented strata (`case_index % 6`): moderate sparsity (~15%
zeros), sparse (~65% zeros), rank-deficient (constant positive feature) + ~35%
zeros, low library sizes + ~45% zeros, `p > n`, and threshold-near weak signal
+ ~30% zeros. `n ∈ [12, 30]`, `p ∈ [2, 8]` (p > n stratum: `p ∈ [n+1, n+4]`,
capped at 12). Sample rows with zero total are repaired with one positive count
(both languages reject all-zero rows). Per-case seeds and StARS subsample
indices come from one NumPy master stream and are injected into both languages.

**Bridge-inversion exclusion (predeclared).** A single pair outside latentcor's
interpolation region (`|K| > 0.9 * bound22(z1, z2)`) costs ~20-25 s in Python,
because the deterministic 4-D normal CDF is evaluated inside Brent's method;
one case with such pairs in its replicas took 102 s. That region is measured
separately by `compare_approx_original.R` and the QMP goldens, so candidates
whose full data *or any replica* contains a cutoff pair are rejected and
re-drawn from the same master stream (max 20 attempts). This run rejected 78
candidates and executed 0 cases with cutoff pairs, as recorded in
`summary.json`.

## Result

| metric | value |
|---|---:|
| cases | 1,000 |
| status (`finite`/`nonfinite`/`error`), Python | 967 / 0 / 33 |
| status (`finite`/`nonfinite`/`error`), R | 967 / 0 / 33 |
| status mismatches | **0** |
| `opt_index` mismatches | **0** |
| refit edge-set mismatches | **1** (boundary-classified) |
| unexplained mismatches | **0** |
| runner exit code | **0** |

The single edge-set difference (`case-00531`) removes/adds edge `0:2` with
`min_abs_coef = 0` in one implementation, i.e. the predeclared soft-threshold
tie (`|beta| <= 1e-6`); it is recorded in `boundary_preview` and is not an
unexplained failure.

Error cases (all status-matched, mclr/sparse strata only: 25 × mode 1,
8 × mode 3): Python raises `ValueError: Matrix seems negative semi-definite`
(30) or `numpy.linalg.LinAlgError: Eigenvalues did not converge` (3), while R
raises `simpleError: missing value where TRUE/FALSE needed` on the same cases
(sampled `case-00007`, `case-00517`, `case-00691`, `case-00763`). The failure
*mechanism text* differs because Python's `near_pd` raises on a non-finite /
degenerate latent correlation while R's `nearPD` propagates `NA` into an `if`;
the failure *category* (error, no result) is identical, which is what this
sweep compares. These are shared degenerate-input failures, not port
divergences, and remain a limitation for such inputs.

## Timing (this machine, serial)

| phase | seconds | per case |
|---|---:|---:|
| Python (classification + bridge screening) | 342.5 | mean 0.34 s, median 0.07 s, p95 1.57 s |
| R classification | 110.0 | mean 0.11 s |
| **total wall clock** | **475.3 (≈ 7 min 56 s)** | |

Screening (cheap Kendall-τ pass over full data and four replicas) is included
in the Python time. Without it, a single cutoff pair would add ~20-25 s per
occurrence; the `case-00001`-style datasets cost ~100 s each.

## Files

- `summary.json` — machine-readable run result, including the bridge screening
  counters.
- `python_status.csv`, `r_status.csv` — per-case status, error class, selected
  index, refit edges, coefficients at the selected lambda, StARS margin and
  merge frequencies.
- `environment.json` — Python/NumPy/platform/R versions, CPU count.
- `cases.sha256` — SHA-256 of every archived case file (3,000 entries).
