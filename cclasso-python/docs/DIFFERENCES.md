# Differences from R (deliberate, documented)

This file records every place where the Python port intentionally deviates from the
R source. Anything not listed here should be considered a bug in the port.

## D1. Parameter name `lam_int` (not `lam_min_ratio`)
The R source (`R/cclasso.R:36`) uses `lam_int`, a length-2 **log10** interval
(default `c(1e-4, 1)`). The task brief suggested `lam_min_ratio`; the source wins,
so the public parameter is `lam_int`.

## D2. CV folds are deterministic
`cv_loss_cclasso` (`R/cclasso.R:211-233`) uses contiguous training/test blocks
`(n_b*(k-1)+1):(n_b*k)` with `n_b = floor(n/k_cv)`. There is **no randomness** in the
cross-validation. A `cv_folds=None` parameter is accepted for API symmetry with the
brief, but when `None` it reproduces R's deterministic blocks.

## D3. Warnings go to stderr, never stdout
R's `cat()` writes "Optimal lambda is near boundary" and the `cclasso_sub` maximum-
iteration warning to **stdout**. The CLI contract forbids extra stdout, so the Python
port writes these to **stderr**. The textual content is otherwise analogous.

## D4. `v1.0/` is not ported
`CCLasso-master/v1.0/R/cclasso.R` is a different algorithm (fixed λ grid
`10^seq(0,-8,-0.01)`, no golden section, `nearPD` post-processing, `kmax=5000`,
`x.tol=1e-6`). Only **v2.0** is ported.

## D5. SparCC not ported
`SparCC.count` / `SparCC.frac` (`R/SparCC.R`) are out of scope for this library.

## D6. RNG-dependent sums require injected indices
Bootstrap indices (`sample`) and the artificial-null draws (`rnorm`) use R's RNG,
which cannot be reproduced by NumPy. They are exposed via `boot_indices` / `null_norm`.
See `docs/EQUIVALENCE.md`.

## D7. `cors_boot`/`lower.tri` ordering
R extracts lower-triangular elements in **column-major** order; NumPy's
`tril_indices` is row-major. The port uses an explicit column-major index vector so
intermediate arrays align with the R golden fixtures.

## D8. Non-convergence / NaN reproduced (not "fixed")
When the ADMM reaches `k_max=200` without meeting `x_tol=1e-4`, R returns the last
iterate even if it yields a non-PSD covariance. The port does the same (it raises a
`RuntimeWarning`). An exact replay of the known NumPy-seeded case confirms one shared
negative bootstrap variance and matching non-finite positions in `cor_w`/`p_vals`
(R represents these as `NA`, Python as `NaN`). These estimates remain unusable;
preserving them is not a recommendation to analyze them. See the replay harness in
`tests/r_reference/README.md`.

An additional prevalence-spread GlobalPatterns seed made R stop inside the ADMM
loop because its stopping expression became a missing logical value. Python had
previously treated `NaN > x_tol` as false and returned non-finite estimates.
The loop now explicitly propagates a NaN stopping condition as the corresponding
`ValueError`, matching R's failure category; finite iteration arithmetic and
operation order are unchanged. This regression is covered by the R failure
fixture.

## D9. Extra private diagnostic fields
`cclasso()` returns private keys `_admm_iters`, `_admm_relative_errors`, and
`_bootstrap_variances` for equivalence/stability diagnostics. They are not part of the
R output and are **stripped** by the CLI. All public output keys mirror the R list.
