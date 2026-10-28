# CCLasso stress sweep — 20 Oct 2026, seed 20261019

Predeclared seeded stress suite required by
`docs/EQUIVALENCE_ACCEPTANCE_CONTRACT.md` ("predeclared seeded stress suite of
at least 1,000 small cases per method family, weighted toward sparse, singular,
and threshold-near inputs").

## Command

```bash
cd cclasso-python
.venv/bin/python tests/r_reference/run_status_stress.py \
    /tmp/cclasso-stress-1000 --cases 1000 --seed 20261019
```

Runner: `tests/r_reference/run_status_stress.py` (already present before this
run); R side: `tests/r_reference/stress_status_R.R`, which sources the pinned
`tests/r_reference/reference/R/cclasso.R` and replays the injected bootstrap
indices and null draws.

## Environment

| item | value |
|---|---|
| date | 2026-10-20 |
| R | 4.6.1 |
| Python | 3.14.7 (`.venv`) |
| NumPy / SciPy | 2.5.4 / 1.18.1 |
| cclasso | 0.1.0 |
| platform | Linux x86-64, serial run |
| settings | `pseudo=0.5`, `k_cv=3`, `lam_int=[1e-4, 1]`, `k_max=3`, `n_boot=3` |

## Predeclared design

Cases are drawn from one NumPy master stream (`--seed`); per case,
`case_index % 4` selects: moderate sparsity (~15% zeros), high sparsity (~65%
zeros), constant-feature/rank-deficient design with ~35% zeros, and low library
sizes with ~45% zeros. `n ∈ [8, 18]`, `p ∈ [2, min(9, n+3)]`. Bootstrap indices
and the null normal matrix are drawn once and written 1-based for R, so both
languages consume identical randomness. Case tree: 1,000 directories, 4,000
files, checksummed in `cases.sha256` (cases themselves stay in the run output
directory and are regenerable from the seed).

## Result

| metric | value |
|---|---:|
| cases | 1,000 |
| status (`finite`/`nonfinite`/`error`), Python | 525 / 394 / 81 |
| status (`finite`/`nonfinite`/`error`), R | 525 / 394 / 81 |
| status mismatches | **0** |
| non-finite output fields | `cor_w` and `p_vals` in all 394 cases; `var_w` never |
| Python error class | `ValueError` (81) |
| R error class | `simpleError` (81) |

By generator stratum (`finite`/`nonfinite`/`error`):

| mode | design | Python | R |
|---|---|---:|---:|
| 0 | moderate sparsity | 97 / 118 / 35 | 97 / 118 / 35 |
| 1 | ~65% zeros | 141 / 99 / 10 | 141 / 99 / 10 |
| 2 | constant feature / rank-deficient | 171 / 71 / 8 | 171 / 71 / 8 |
| 3 | low library sizes | 116 / 106 / 28 | 116 / 106 / 28 |

Runner exit code: **0** (no status mismatch).

## Interpretation

- The port reproduces the reference **failure category on every case**, which
  is what the contract requires for the failure/control-flow gate: no
  Python-only failure and no Python-only non-finite output in this domain.
- The 394 non-finite cases are the shared negative-bootstrap-covariance
  behaviour already documented (`docs/STATISTICAL_RELEASE_READINESS.md`,
  `docs/NUMERICAL_EQUIVALENCE.md`): Python emits NaN where R emits NA, at the
  same positions. They are counted as matched categories, **not** as valid
  estimates.
- The 81 errors are the shared missing-logical condition (sampled case
  `case-00006`): Python raises `ValueError: missing value where TRUE/FALSE
  needed`, R raises the corresponding `simpleError`. Class names differ by
  language; the status classification is identical.
- This sweep is supporting evidence over one seeded design, not proof over all
  inputs, and it does not change the readiness verdict: CCLasso convergence and
  the shared non-finite failure remain documented limitations.

## Timing (this machine, serial)

| phase | seconds |
|---|---:|
| Python classification (1,000 cases) | ~144 |
| R classification (1,000 cases) | ~88 |
| **total wall clock** | **232 (≈ 3 min 52 s)** |
| 20-case smoke test (incl. startup) | 6 |

Cost is dominated by the per-case ADMM/bootstrap work in both languages; a
10,000-case sweep would be roughly 35–40 min at this scale.

## Files

- `summary.json` — machine-readable run result.
- `python_status.csv`, `r_status.csv` — per-case status, error class and
  non-finite counts per output.
- `environment.json` — Python/NumPy/SciPy/R versions and platform.
- `cases.sha256` — SHA-256 of every archived case file (4,000 entries).
