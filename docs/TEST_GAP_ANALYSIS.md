# Test coverage and remaining gaps

## Scope

The three Python ports are checked against their R references with committed
QMP real-data fixtures. The two fixtures per port are derived from the public
QMP stool microbiome data in `SPRING-master/data/QMP.rda` (106 samples × 91
taxa). They cover all samples with a zero-rich 20-taxon subset and a 20-sample ×
30-taxon `p > n` subset. Taxon/sample names and source parameters are retained
with the fixtures. Small synthetic cases remain generated in memory by unit and
property tests where they exercise known mathematical invariants; synthetic
golden datasets have been replaced.

## Coverage matrix

| Area | Existing coverage | Added/expanded coverage | Remaining limitation |
|---|---|---|---|
| CCLasso R equivalence | R-generated intermediates, bootstrap indices, null draws, CV path, final variances/correlations/p-values, exact ADMM call counts | QMP plus American Gut and GlobalPatterns; three R seeds per added case family; per-call residual parity; matching R/Python error fixture | All 70 QMP ADMM solves hit the 200-iteration cap; GlobalPatterns also has shared non-finite/error outcomes. Parity is not convergence |
| CCLasso edge/stability | Input validation, scale invariance, PSD/symmetry, p-value bounds, determinism | Minimum `p=2`, common count scaling by `1e100`, Student-t CDF tails; exact seeded non-finite and missing-logical R failure replays | Shared R failures/non-finite outputs remain unusable estimates and must be excluded or explicitly reported |
| SPIEC-EASI R equivalence | CLR, correlation, lambda grid, MB path, StARS/refit, RNG-injected subsamples | QMP, American Gut and GlobalPatterns; three R seeds per added case family; objective/KKT checks remain | The added fixture suite passes; `glasso`, `slr`, `coat`, `ising`, and `poisson` remain intentionally unported |
| SPIEC-EASI edge/stability | Zero handling, known-graph recovery, validation, degenerate Lasso KKT/objective checks | One-lambda/equal-endpoint lambda paths and duplicate-feature singular-design test | R/NumPy BLAS boundary ties can differ by ulps at `lambda.max`; this remains explicitly diagnosed rather than hidden |
| SPRING R equivalence | MCLR, Kendall tau, latent correlation, MB paths, StARS/refit, injected R indices | QMP, American Gut and GlobalPatterns; three R seeds per added case family; lambda-path injection isolates public-data MB/StARS comparisons | One GlobalPatterns inversion pair shifts computed lambda-max by about `1.3e-6` relative; public tests inject R's path. Coefficients are invariant-checked on inversion cases |
| SPRING edge/stability | `nearPD`, interpolation, ties, bridge determinism, original inversion, validation, CLI, synthetic recovery | Explicit `use_nearPD=False` pointwise behavior and unsupported latentcor type rejection | Only `types="tru"` is implemented because that is the type SPRING itself invokes; `con`/`bin`/`ter`/`ord` latentcor modes are not ported |

## Test sufficiency assessment

Before the data expansion, the largest gap was lack of independent real
microbiome inputs in Python/R equivalence tests. QMP, American Gut Round 1, and
GlobalPatterns now cover distinct study sources, varied zero/library-size
profiles, and `p < n`/`p > n`. Additional unit cases target the minimal CCLasso
dimension, large finite count scales, Student-t tails, singular SPIEC-EASI
designs, lambda-path endpoints, and SPRING branches for `nearPD` and latentcor
type validation.

The tests are extensive for the implemented `mb` and CCLasso paths but are not a
proof of universal correctness. They do not exhaust every possible count scale,
BLAS/LAPACK implementation, operating system, or R package version. In
particular, outputs close to an optimizer threshold and CCLasso's observed ADMM
non-convergence must be interpreted using the numerical-equivalence contract in
`docs/NUMERICAL_EQUIVALENCE.md`.

## CI wiring note

The root `.github/workflows/ci.yml` still references the removed compiled
`spieceasi._mb_core`/`mb_solver` path and test modules that are absent from the
current SPIEC-EASI package. The standalone `spiec-easi-python/` is also an
ignored nested Git repository, not a submodule in the Electron app checkout.
Repairing the root workflow requires a repository-integration decision (vendor
the package, add a submodule, or delegate to the nested repository's CI), so the
stale workflow was documented rather than rewritten to point at files that a
root checkout cannot obtain. The standalone SPIEC-EASI repository has its own
workflow under `spiec-easi-python/.github/workflows/ci.yml`.
