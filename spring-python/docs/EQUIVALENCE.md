# Equivalence with the R implementation

This document states exactly **what is guaranteed identical**, what is identical
only **within numerical tolerance**, and what **cannot** be identical, between
this Python port and the R implementation. It is intentionally blunt: nothing is
hidden.

Reference R stack (the machine that produced the goldens):

```
R 4.6.1 (2026-06-24), x86_64-pc-linux-gnu
latentcor 2.0.2, huge 2.0.1, pulsar 0.3.13 (GitHub zdk123 @ 2bc391f),
pcaPP 2.0.7, fMultivar 4031.84, mnormt 2.1.2, Matrix 1.7.5
```

Python stack: NumPy 2.5.4, SciPy 1.18.1 on CPython 3.14 (Linux x86-64).

The current equivalence suite (`tests/test_equivalence.py`, `-m slow`) replays
golden cases from QMP, American Gut Round 1 and GlobalPatterns, plus a synthetic
latent-correlation fixture for `Rmethod="original"`. It covers
`p < n`, `p > n`, and distinct sparsity/library-size profiles. The case matrix
is chosen so that the guarantees below are not limited to one configuration:

| Coverage | Cases | Settings |
|---|---|---|
| SPRING defaults | `qmp_defaults_n106_p8`, `amgut_r1_n80_p8_defaults` | `nlambda=20`, `rep.num=20`, `thresh=0.1`, default λ path; one quantitative (`QMP`) and one compositional (`mclr`) |
| `Rmethod="original"` | `original_latentcor_n60_p3` (dedicated latentcor fixture) | numeric bridge inversion for every pair, directly compared with R; a full StARS run is impractical at ~20 s/inverted pair |
| data-specific λ on compositional public data | `amgut_r1_n80_p20_dataseq`, `globalpatterns_n26_p15_dataseq` | full Python `latentcor → λ-max → grid` chain, not an injected R grid |
| data-specific λ on quantitative data | `qmp_n106_p20`, `qmp_n20_p30` | incl. `p > n` |
| broad reduced-StARS sweep | six `*_seed310x` public cases | `nlambda=3`, `rep.num=2`, `thresh=0.8`, R-derived numeric λ |

Public `*_seed310x` cases use the R-derived numeric λ sequence to isolate
MB/StARS parity; the `*_dataseq` cases and QMP retain the data-specific path
check. The former synthetic golden fixtures have been replaced; historical
measurements below that refer to synthetic cases are retained for context, not
claimed as current fixtures.

Current QMP measurements: the selected StARS index and binary final adjacency
match exactly in both cases. In `qmp_n20_p30`, the data-specific lambda grid
differs by `4.02e-8` absolute (`4.87e-8` relative), and the maximum pointwise
full-data coefficient difference is `3.28e-5`; support is identical. This
real p>n case is therefore checked by support and each implementation's KKT
residuals rather than the strict pointwise coefficient tolerance used in the
well-conditioned case. Its nearPD-regularized matrix is full-rank but has
condition number about `8.1e3`. The largest observed
active KKT residual is `8.86e-5`, within huge's `thol=1e-4`.

## 1. Guaranteed *identical* (up to the RNG-injection described in §4)

If you inject the same StARS subsample indices, the following match **exactly**
(0 differing entries) or to machine precision:

| Stage | R | Guarantee |
|---|---|---|
| mclr transformation | `mclr` | ≤ 2e-15 abs diff |
| Kendall τ matrix `K` | `latentcor$K` | ≤ 2e-16 abs diff |
| pointwise latent correlation `Rpointwise` | `latentcor$Rpointwise` | ≤ 1e-15 abs diff |
| `nearPD` + `nu` shrinkage → `R` | `latentcor$R` | ≤ 3e-15 abs diff |
| λ grid | `lambdaseq` | tight for fixed paths; data-specific p>n case uses the measured tolerance in §2 |
| per-subsample MB **edge sets** | `hugeKmb$path` | **identical sets** (0 symmetric difference) |
| MB coefficients | `hugeKmb$beta` | ≤ 3e-13 abs diff in well-conditioned cases; p>n uses support/KKT checks |
| StARS selection frequencies `merge` | `pulsar$stars$merge` | ≤ 3e-16 abs diff |
| StARS instability curve `summary` | `pulsar$stars$summary` | ≤ 3e-16 abs diff |
| optimal index `opt.index` | `pulsar$stars$opt.index` | **exact integer** |
| final adjacency (`refit`) | `refit$stars` | **exact binary matrix** |

The identical edge sets and exact `opt.index`/adjacency on the committed cases
are the strongest result. Bridge-inversion cases are additionally assessed by
support, KKT, and objective invariants where pointwise coefficient equality is
not achieved; see the measured limits below.

## 2. Identical only within a documented tolerance

| Node | Why | Observed |
|---|---|---|
| `Rmethod = "original"` λ grid | the λ grid is derived from the latent correlation, whose inversion uses `optimize` (R) vs the ported R Brent algorithm | ≤ 7e-6 relative; `opt.index` and adjacency still match exactly |
| `Rmethod = "original"` MB coefficients | they inherit the λ-grid difference (β is Lipschitz in λ) | ≤ 1.3e-6 absolute / 3.3e-6 relative |
| QMP p>n data-specific λ grid | SciPy and R bridge inversions differ slightly in the zero-rich rank-deficient subset | ≤ 4.02e-8 absolute / 4.87e-8 relative; test tolerance `rtol=1e-7, atol=1e-9` |
| data-specific λ grid when the full-data latent correlation contains a bridge-inversion pair | λ-max is an entry of the matrix computed through the bridge inversion, so it inherits the documented bridge-CDF difference | GlobalPatterns p15: 1.27e-6 relative; test tolerance `rtol=5e-6, atol=1e-8` |
| QMP p>n pointwise MB coefficients | the input correlation is rank-deficient and the optimizer may have multiple minimizers | max absolute difference `3.28e-5`; support exact and both KKT residuals ≤ `8.86e-5` |
| latent correlation for pairs with `\|K\| > ratio·bound` | these do **not** use the interpolation table: they invert the bridge function, which needs the 4-variate normal CDF. R uses `mnormt::pmnorm` (Genz `sadmvn`, `abseps=1e-6, releps=0`); Python uses `scipy.stats.multivariate_normal.cdf`. On the mildly ill-conditioned 4×4 bridge matrices that occur here the two CDFs differ by up to ~8e-5. **We arbitrated this with independent 4-D Gauss–Legendre quadrature** (`n=45` per dimension): SciPy agrees with the quadrature to ~3e-6, whereas `mnormt` is off by ~7.6e-5 on the worst matrix (min eigenvalue ≈ 3e-3); `mnormt`'s value does not move when its tolerance is tightened from 1e-6 to 1e-12, so the error is structural. **We therefore deliberately do not reproduce R's value here**: our result is the more accurate one. The bridge objective is very flat near its root, so the inverted correlation shifts by ≤ 5e-5 absolute (~2e-4 relative) relative to R | `Rpointwise`: strict (≤ 3e-15) on interpolated pairs, ≤ 3.2e-5 on bridge-inversion pairs; `R` (post-`nearPD`, which mixes all entries): ≤ 9.4e-6 |
| determinism of the latent correlation | SciPy's `multivariate_normal.cdf` uses **randomized QMC** and is non-deterministic by default (its value changes call to call). The port passes a fresh `default_rng(0)` **on every call** with `abseps=releps=1e-10`, `maxpts=1e6`, making each evaluation reproducible to ~1e-7 and independent of call order | fully deterministic |
| `mnormt::pmnorm` vs `scipy.stats.multivariate_normal.cdf` | different quadrature constants | ≤ 1e-12; only reaches the result for pairs with `|K| > ratio·bound` |
| `pcaPP::cor.fk` vs our τ-b | different sort/accumulation order | ≤ 1e-16 |
| `Matrix::nearPD` (pure R) vs our port | BLAS vs LAPACK eigen routines | ≤ 3e-15 |
| MB coordinate descent | summation order (BLAS `dsyrk` vs NumPy) | coefficients ≤ 3e-13 |
| BLAS/LAPACK generally | library differences | 1e-12 … 1e-15 |

## 3. Cannot be identical

1. **Random number generation.** R's `set.seed`/`sample` and NumPy's
   `Generator` are different generators; the same integer seed yields different
   subsamples. Consequences:
   - *With* `subsample_indices` injected → full equivalence (§1).
   - *Without* injection → **statistically** equivalent only. The selected
     `opt.index` can differ because the StARS instability curve is evaluated on
     different subsamples.
2. **Boundary ties in the Lasso path.** A coefficient exactly at the
   soft-threshold boundary can be selected on one side and dropped on the other.
   Our coordinate descent fixes the iteration order and tolerances to match
   `huge` 2.0.1, so this was **not observed** in any golden case; if it ever
   occurs it must be reported, not hidden.
3. **`use.nearPD = FALSE` on non-PSD data.** When the pointwise latent
   correlation is not positive semidefinite, `huge.mb` refuses the input
   (`"Covariance input must be positive semidefinite."`). This is faithful R
   behaviour (SPRING crashes too) and is therefore not tested as a golden case.
4. **`synthData_from_ecdf` is not ported** (out of scope by decision): it depends
   on R's RNG and on `rootSolve::uniroot.all`, so only statistical equivalence
   would be possible. Use `docs/PACKAGING.md`-documented simple generators for
   simulation studies instead.

## 4. How to obtain identical results

```python
import numpy as np
from spring_microbiome import spring

# indices must be 0-based (R is 1-based): convert with `- 1`
indices = [np.asarray(ix) - 1 for ix in r_indices]  # one array per replicate
result = spring(counts, subsample_indices=indices, seed=10010)
```

The golden generator writes R's indices to `subsample_indices.csv`, which is
exactly how the equivalence tests are able to reproduce R bit-for-bit.

## 5. Scope of the evidence

The measurements above were produced on **Linux x86-64** with the versions
listed at the top. Different BLAS/LAPACK builds or different R/CRAN versions can
shift the tolerance-level quantities; the *discrete* results (edge sets,
`opt.index`) are robust to those shifts because they are determined by
inequalities with wide margins in the tested cases.

## 6. Reproducibility of the fixtures

The fixture directory is **self-contained**:

- `tests/r_reference/reference/R/` vendors the SPRING R source (0.99.0);
- `tests/r_reference/reference/data/QMP.rda` vendors the QMP data;
- `tests/r_reference/reference/public_microbiome_cases.R` vendors the public-data
  selector (American Gut via `SpiecEasi`, GlobalPatterns via `phyloseq`).

`generate_golden.R` prefers these paths, so regenerating the goldens does not
require an external `SPRING-master` checkout. Regenerating all fixtures from the
vendored sources reproduced every pre-existing numeric artefact **bit-for-bit**;
only `params.json` gained provenance metadata (source/zero-fraction/library-size
fields). `tests/golden/MANIFEST.sha256` records a SHA-256 for every fixture file.
The whole `spring-python/` tree is tracked by the workspace Git repository, and
the goldens are committed so the test-suite never needs R.
