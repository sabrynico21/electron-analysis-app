# R reference fixtures

These fixtures pin the Python port to the R implementation **CCLasso v2.0**.

## Layout

- `reference/R/cclasso.R` — verbatim copy of `CCLasso-master/R/cclasso.R`
  (the source of truth; do not edit).
- `reference/DESCRIPTION` — verbatim copy of the R package `DESCRIPTION`.
- `generate_golden.R` — generator script.
- `../golden/` — real-data fixtures and explicit R-failure cases, **committed** so
  Python tests never need R.

## Regenerating

From the `cclasso-python/` directory:

```bash
Rscript tests/r_reference/generate_golden.R tests/golden
```

Environment used for the committed fixtures (see `../golden/sessionInfo.txt`):

- R 4.6.1 (2026-06-24), x86-64 Linux.
- The pinned `cclasso.R` algorithm requires base R. Public-data extraction also
  uses `SpiecEasi` 1.99.0 and `phyloseq` 1.56.0; these are not algorithm
  dependencies. (`SparCC.R` needs `gtools`, but SparCC is out of scope.)
- Data sources: QMP is **vendored** at `tests/r_reference/data/QMP.rda`
  (SHA-256 `603ddba2...ff1e8b`) so regeneration works from a clean checkout;
  set `QMP_RDA` to override. American Gut Round 1 comes from
  `SpiecEasi::amgut1.filt`, and Caporaso et al. GlobalPatterns from
  `phyloseq::GlobalPatterns`. Public-data preprocessing is shared from
  `../scripts/public_microbiome_cases.R`. Provenance and redistribution status
  are audited in `../../../docs/DATASET_PROVENANCE.md`; public fixtures are
  marked `redistribution_cleared: false` until terms are confirmed.

After regenerating, update the manifest:

```bash
cd tests/golden && find . -type f ! -name MANIFEST.sha256 -print0 \
  | sort -z | xargs -0 sha256sum > MANIFEST.sha256
```

## Cases

| case | dimensions | selection |
|---|---:|---|
| `qmp_n106_p20` | 106 × 20 | all QMP samples; 20 least-prevalent taxa |
| `qmp_n20_p30` | 20 × 30 | real p > n subset; 30 taxa across prevalence ranks |
| `amgut_r1_n80_p20_seed3101..3103` | 80 × 20 | deterministic sample/taxon subset; 3 independent R bootstrap seeds |
| `amgut_r1_n20_p30_seed3101..3103` | 20 × 30 | real p > n subset; 3 independent R bootstrap seeds |
| `globalpatterns_n26_p15_seed3101..3103` | 26 × 15 | most-prevalent varying taxa; 3 independent R bootstrap seeds |
| `globalpatterns_n26_p40_seed3101..3103` | 26 × 40 | real p > n subset; 3 independent R bootstrap seeds |
| `synth_frac_n37_p8_k3_b4_seed4101` | 37 × 8 | `counts=False`; deterministic logistic-normal fractions |
| `synth_frac_n45_p8_k5_b3_seed4106` | 45 × 8 | `counts=False`; non-default `k_cv=5` |
| `synth_frac_n34_p8_k3_b5_lam_seed4103` | 34 × 8 | `counts=False`; non-default `lam_int=(1e-3,1e-1)` |
| `synth_counts_n36_p8_k4_b7_lam_seed4104` | 36 × 8 | `counts=True`; non-default `k_cv=4`, `n_boot=7`, `lam_int=(1e-5,1)` |

The separate `tests/golden/failures/` tree records reference-side failure
outcomes. `globalpatterns_n26_p40_spread_seed3102` is a deliberately sparse
prevalence-spread stress input where R errors during CCLasso's ADMM loop; Python
is tested to raise the corresponding missing-logical error rather than return a
non-finite estimate. It is not counted as a successful equivalence case.

Each case contains `taxa.txt` and `sample_ids.txt` to preserve the selected
features and observations, in addition to the numeric matrices below.

## What is saved per case (`tests/golden/<case>/`)

| file | content |
|---|---|
| `input.csv` | raw input matrix `n×p` (before pseudo-count / log) |
| `params.json` | dataset and `counts,pseudo,k_cv,lam_int,k_max,n_boot,seed,n,p`; public cases also carry `source`, `source_package`, `source_package_version`, `source_license`, `source_citation`, `retrieval_date`, `redistribution_cleared` |
| `boot_indices.csv` | bootstrap row indices `n×n_boot`, **1-based** (R order) |
| `null_norm.csv` | raw `rnorm(n*p)` draws `n×p` for the artificial null sample |
| `x_log.csv` | matrix after pseudo-count + log |
| `vx2.csv` | covariance of `x_log` (p×p, ddof=1) |
| `wd.csv` | diagonal loss weights (length p) |
| `lams.csv` | λ values **in log10** visited by the golden section, in order |
| `fvals.csv` | CV loss for each `lams` entry, same order |
| `info_cv.json` | `{k}` — number of CV evaluations (`k+2`) |
| `lambda.json` | `{lambda, lam_lo, lam_hi}` |
| `var_w.csv` | variance estimates (length p) |
| `cor_w.csv` | final correlation matrix (p×p) |
| `p_vals.csv` | p-values (p×p) |
| `admm_iters.csv` | ADMM iteration count for every `cclasso_sub` call, in call order |
| `admm_relative_errors.csv` | final relative residual for every ADMM call, including capped calls |

All numeric values are written with `sprintf("%.17g", x)` (near full double precision).

## Random-index injection

`boot_indices.csv` and `null_norm.csv` are produced by replaying the R RNG with the
same `seed` used for the real run. `boot_cclasso` consumes RNG only via
`sample()` then `rnorm()`, so the replay is exact; the generator **asserts** this for
every case. The Python port injects these two arrays, which makes the random parts
bit-comparable across languages.

## Replaying the known NumPy-seeded non-finite case

The determinism property test uses a seeded synthetic input and had been observed
to produce non-finite correlations after one bootstrap covariance acquired a
negative diagonal. The following commands regenerate its exact input and NumPy
draws, then replay those same values through the pinned R source (base R only):

```bash
.venv/bin/python tests/r_reference/make_numpy_nan_replay.py /tmp/cclasso-nan-replay
Rscript tests/r_reference/replay_numpy_nan.R \
  /tmp/cclasso-nan-replay /tmp/cclasso-nan-replay
```

The generated data are diagnostic artifacts, not synthetic golden fixtures. The
R report writes the bootstrap variance matrix and per-output finite/non-finite
counts so the cross-language failure status can be checked directly.

## Seeded failure-status stress runner

Run the registered 1,000-case CCLasso sweep (small count dimensions, varying
zero prevalence/library sizes, constants, and singular designs) with:

```bash
.venv/bin/python tests/r_reference/run_status_stress.py \
    /tmp/cclasso-stress-1000 --cases 1000 --seed 20261019
```

The runner creates the Python cases and injected draws once, runs Python and the
pinned R source on those same inputs, and writes per-language status CSVs plus a
summary. Any disagreement in `finite`/`nonfinite`/`error` classification or the
count of non-finite output entries exits nonzero. The complete case archive is
retained at the output path for investigation.

Recorded run: `../stress/2026-10-20-seed20261019/` (1,000 cases, seed
20261019, 232 s wall clock; 525 finite / 394 nonfinite / 81 error in both
languages, 0 mismatches, exit code 0).
