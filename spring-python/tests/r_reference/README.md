# R reference / golden fixtures

This folder contains the SPRING-specific R reference and fixture generators.
They generate the goldens in `../golden/`, which the Python test suite consumes
so **running the tests never requires R**.

## Self-contained regeneration

`reference/` vendors everything the generator needs from the R side:

| Path | Provenance | License |
|---|---|---|
| `reference/R/` | `SPRING-master/R/` (SPRING 0.99.0, `IrinaStatsLab/SPRING`) | GPL-3 |
| `reference/data/QMP.rda` | `SPRING-master/data/QMP.rda` | GPL-3 |
| `reference/public_microbiome_cases.R` | workspace `scripts/public_microbiome_cases.R` | GPL-3 |
| `reference/DESCRIPTION`, `reference/NAMESPACE` | `SPRING-master/` | GPL-3 |

`generate_golden.R` prefers these vendored paths, so regenerating the fixtures
does **not** depend on an external `SPRING-master` checkout. `SPRING_R_HOME` can
still override `reference/R/`. The only external requirements are the R packages
listed in `../golden/sessionInfo.txt` (and `SpiecEasi`/`phyloseq` for the public
datasets).

## Files

| File | Purpose |
|---|---|
| `generate_golden.R` | Generates every golden case (inputs, ALL intermediate outputs, exact subsample indices, `sessionInfo()`). |
| `run_status_stress.py` | Predeclared seeded stress sweep (status + discrete outputs, R/Python). |
| `stress_status_R.R` | R classifier used by the stress runner on the archived cases. |
| `compare_approx_original.R` | Quantifies the `approx`/`original` latentcor sensitivity per source/dimension family. |

## Golden configuration coverage

The `generate_golden.R` case list deliberately spans settings that the original
fixture set did not:

| Family | Settings | What it proves |
|---|---|---|
| `*_seed310x` public cases | `nlambda=3`, `rep.num=2`, `thresh=0.8`, R-derived numeric λ | broad, cheap parity across three seeds |
| `qmp_defaults_n106_p8`, `amgut_r1_n80_p8_defaults` | SPRING defaults: `nlambda=20`, `rep.num=20`, `thresh=0.1`, default λ path | parity at the real defaults, not just a reduced StARS run |
| `original_latentcor_n60_p3` (dedicated `kind="latentcor"` fixture) | `Rmethod="original"` | the numeric bridge-inversion branch against R |
| `*_dataseq` (American Gut, GlobalPatterns) | `lambdaseq="data-specific"` on `mclr` data | full Python `latentcor → λ-max → grid` chain on compositional input |
| `qmp_n106_p20`, `qmp_n20_p30` | `lambdaseq="data-specific"`, `p>n` included | data-specific path on quantitative QMP incl. `p>n` |


## Seeded failure-status stress sweep

Run the registered 1,000-case SPRING sweep (small count dimensions in six
predeclared strata: moderate sparsity, ~65% zeros, rank-deficient/constant
feature, low library sizes, `p > n`, and threshold-near weak signal) with:

```bash
.venv/bin/python tests/r_reference/run_status_stress.py \
    /tmp/spring-stress-1000 --cases 1000 --seed 20261020
```

The runner draws one seed per case and one StARS-index draw seed per case from
a NumPy master stream, archives every case (`counts.csv`,
`subsample_indices.csv`, `params.json`), classifies each case in Python, then
runs `stress_status_R.R` (which sources `SPRING-master/R/helpers.R`) on the same
counts and injected subsample indices. Both sides use `quantitative=FALSE`,
`nlambda=5`, the explicit default lambda path (`lambda.min.ratio=1e-2`),
`thresh=0.1`, `subsample.ratio=0.8`, `rep.num=4`, `Rmethod="approx"`,
`use.nearPD=TRUE`. The data-specific lambda path is **not** used here; its known
bridge-CDF lambda-max difference stays covered by the QMP goldens.

It writes `python_status.csv`, `r_status.csv`, `summary.json` and
`environment.json`. A case is an unexplained failure when the
`finite`/`nonfinite`/`error` status differs, or when the selected index or
refit edge set differs without meeting the predeclared boundary rules:
`|stability - thresh| <= 1e-3` for index ties, and for edge differences either
a selection frequency within 0.05 of 0.5 or a soft-threshold-boundary
coefficient (`|beta| <= 1e-6`). Any unexplained difference exits nonzero. The
complete case archive is retained at the output path.

**Bridge-inversion exclusion (predeclared).** A pair outside latentcor's
interpolation region (`|K| > 0.9 * bound22(z1, z2)`) costs ~20-25 s in Python
(deterministic 4-D normal CDF inside Brent's method), so candidates with any
such pair on the full data or any injected replica are re-drawn from the same
master stream (max 20 attempts). That region is measured separately by
`compare_approx_original.R` and the QMP goldens; `summary.json` records the
number of rejected candidates and verifies that no executed case contains a
cutoff pair.

Recorded run: `../stress/2026-10-20-seed20261020/` (1,000 cases, seed
20261020, 475 s wall clock, 78 candidates redrawn; 967 finite / 0 nonfinite /
33 error in both languages, 0 status mismatches, 0 selected-index mismatches,
1 boundary-classified single-edge difference, exit code 0).

## Regenerating the goldens

You need a working R with the SPRING dependencies plus `SpiecEasi` 1.99.0 and
`phyloseq` 1.56.0 for public-data extraction. On the machine where the current
goldens were produced:

```
R 4.6.1 (2026-06-24), x86_64-pc-linux-gnu
latentcor 2.0.2, huge 2.0.1, pulsar 0.3.13,
pcaPP 2.0.7, fMultivar 4031.84, mnormt 2.1.2, Matrix 1.7.5, Rcpp 1.1.1.1.1
```

`SPRING-master` must be present at `<repo>/../SPRING-master` (or set the
`SPRING_R_HOME` environment variable to its `R/` directory). The QMP dataset is
loaded from `../SPRING-master/data/QMP.rda` by default; set `QMP_RDA` to override
that path. Then, from the package root (`spring-python/`):

```bash
Rscript tests/r_reference/generate_golden.R tests/golden
```

The full `sessionInfo()` of the generating session is written to
`tests/golden/sessionInfo.txt` and in each case directory. **Always** keep it
with the goldens.
QMP cases retain the data-specific lambda path. For the two independent public
sources, the R-derived numeric lambda sequence is injected into Python so those
fixtures isolate MB/StARS parity from the documented bridge-CDF difference in
lambda-max calculation. The Python/R data-specific path remains covered by QMP.

To quantify the `approx`/`original` latentcor sensitivity on a representative
case from each source/dimension family (without running the MB/StARS pipeline):

```bash
Rscript tests/r_reference/compare_approx_original.R \
    tests/golden /tmp/spring-bridge-sensitivity
```

This reports pairwise deviations separately for table-interpolation and
bridge-inversion regions. `approx` and `original` are distinct computational
methods, so their difference is not a numerical-noise floor; see
`docs/NUMERICAL_EQUIVALENCE.md` for the measured interpretation.

After regenerating, refresh the checksum manifest:

```bash
find tests/golden -type f ! -name MANIFEST.sha256 -print0 | sort -z \
  | xargs -0 sha256sum > tests/golden/MANIFEST.sha256
```

## The interpolation tables (`ipol_22.npz`)

`latentcor`'s default (`method = "approx"`) latent correlation uses precomputed
multilinear interpolation tables shipped inside `latentcor/R/sysdata.rda`. To
guarantee bit-identical interpolation *inputs*, we extract the table SPRING
actually uses (`ipol_22`, truncated/truncated) and ship it with the Python
package:

```bash
Rscript scripts/extract_ipol_tables.R /tmp/ipol_raw
.venv/bin/python scripts/pack_ipol22.py /tmp/ipol_raw \
    src/spring_microbiome/latentcor/data/ipol_22.npz
```

`scripts/extract_ipol_tables.R` also writes `ipol_meta.json` recording
`latentcor_version` and `r_version`, so the provenance of the table is explicit.

## What is written per case

At 17 significant digits (round-trip safe):

- `input.csv` — raw counts (n×p).
- `params.json` — every SPRING parameter + `opt.index` + the numeric `lambdaseq`.
- `qdat.csv` — the matrix actually fed to pulsar (counts or `mclr` output).
- `lambdaseq.csv`, `summary.csv`, `opt_index.txt`.
- `full_K.csv`, `full_Rpointwise.csv`, `full_R.csv` — full-data latentcor
  intermediates (Kendall τ, pre-`nearPD`, post-`nearPD`+`nu`).
- `merge_XX.csv` — StARS selection-frequency matrix per λ (X = 1..nlambda).
- `beta_full_XX.csv`, `adjacency_full_XX.csv` — full-data MB coefficients /
  adjacency per λ.
- `adjacency_refit.csv` — `fit$refit$stars`, the final adjacency.
- `subsample_indices.csv` — 1-based R subsample indices (rep.num rows).
- `subsample_R_XX.csv` — latentcor `R` for each subsample.
- `subsample_paths_long.csv` — nonzeros of every per-subsample per-λ adjacency
  (`k,lambda,row,col`).
- `subsample_beta_long.csv` — nonzeros of every per-subsample per-λ MB
  coefficient matrix (`k,lambda,row,col,value`).

## Cases

| Case | n×p | zeros | λ path | notes |
|---|---|---|---|---|
| `qmp_n106_p20` | 106×20 | 61.51% | data-specific | all samples, quantitative QMP counts |
| `qmp_n20_p30` | 20×30 | 32.0% | data-specific | real `p > n` subset, quantitative counts |
| `amgut_r1_n80_p20_seed3101..3103` | 80×20 | 29.50% | R-injected | American Gut Round 1, 3 StARS seeds |
| `amgut_r1_n20_p30_seed3101..3103` | 20×30 | 31.33% | R-injected | American Gut Round 1, `p > n`, 3 StARS seeds |
| `globalpatterns_n26_p15_seed3101..3103` | 26×15 | 0% | R-injected | GlobalPatterns prevalent taxa, 3 StARS seeds |
| `globalpatterns_n26_p40_seed3101..3103` | 26×40 | 0% | R-injected | GlobalPatterns, `p > n`, 3 StARS seeds |

The source QMP matrix has 106 samples × 91 taxa and 27.10% zeros. The smaller
fixtures retain sample/taxon names in `sample_ids.txt` and `taxa.txt`; the
selection procedure is in `generate_golden.R`. The unit/CLI tests continue to
exercise the compositional MCLR path and `Rmethod="original"` without relying
on synthetic golden files.

> **Note on `use.nearPD = FALSE`:** no golden case uses it. With
> `use.nearPD = FALSE` the pointwise latent correlation can be non-PSD, and
> `huge.mb` (2.0.1) then raises *"Covariance input must be positive
> semidefinite."*. This is faithful R behaviour (SPRING would crash too), so it
> is documented in `docs/DIFFERENCES.md` rather than tested as a golden.
