# Numerical equivalence: R and Python

## Short answer

Byte-for-byte equality is not a general promise for two independently implemented
numerical programs. It is achievable for particular outputs on a pinned machine
and environment, and should be required for discrete outputs when the test inputs
are well away from decision boundaries. It is not a sound universal requirement
for floating-point intermediates, optimizer coefficients, or results across
platforms.

The current regression fixtures compare all three Python ports with their R
references on real QMP microbiome data. R's random draws are saved and injected
into Python, so RNG differences do not obscure deterministic algorithm
comparisons. The fixtures and regeneration instructions are under each port's
`tests/r_reference/` directory.

## Why floating-point results can differ

1. **RNG streams differ.** R's `set.seed()`/`sample()`/`rnorm()` and NumPy's
   generators are different algorithms. The same integer seed does not produce
   the same resamples or bootstrap draws. The equivalence tests therefore inject
   R's saved draws. Without injection, only statistical—not numeric—equivalence
   is expected.
2. **Floating-point addition is order-dependent.** IEEE-754 arithmetic rounds
   each operation. Matrix reductions performed in different orders can differ by
   a few ulps even when both implementations use binary64. BLAS vendor, threading,
   vectorization, compiler, and fused multiply-add choices can change that order.
3. **BLAS/LAPACK implementations are not identical.** R and NumPy/SciPy may call
   different routines for matrix products, eigenproblems, and decompositions.
   Eigenvectors are also sign-indeterminate, and bases in repeated-eigenvalue
   spaces are not unique.
4. **Special functions and numerical integration differ.** Student-t and
   multivariate-normal CDFs, logarithms, exponentials, and bridge inversions can
   use different approximations, stopping criteria, or quadrature. SPRING's
   latent-correlation bridge inversion is especially sensitive in the
   non-interpolated region; its R reference uses `mnormt`, while Python uses
   SciPy's multivariate-normal CDF.
5. **Iterative algorithms amplify tiny differences.** ADMM, coordinate descent,
   Brent optimization, and StARS compare values to thresholds. A small numeric
   perturbation may change an iteration count, selected coefficient, or lambda
   index when the result lies close to a stopping or soft-threshold boundary.
6. **Some optimization solutions are not unique.** With rank-deficient designs
   (notably `p > n`), different coefficient vectors can minimize the same Lasso
   objective. In those cases, support, objective, and KKT conditions are more
   meaningful comparisons than pointwise coefficient equality.

## Comparison contract

- Compare deterministic transforms and integer/control-flow results exactly or
  at tight, explicitly stated tolerances when the implementations follow the same
  arithmetic path.
- Require exact edge sets, selected indices, and binary adjacency only for
  fixtures whose decisions are not numerically ambiguous. Boundary cases must be
  reported and analyzed, not silently hidden behind a larger tolerance.
- Compare continuous solver outputs using measured absolute/relative errors and
  algorithmic invariants (objective and KKT checks for non-unique solutions).
- Record the OS, R/Python versions, numerical-library versions, and RNG injection
  conditions whenever publishing a maximum deviation. A tolerance measured on
  one machine is not a cross-platform guarantee.

## Current real-data observations

The source dataset is the public QMP dataset distributed with SPRING:
`SPRING-master/data/QMP.rda` (106 samples × 91 taxa; 27.10% zero entries; no
all-zero rows or taxa); `cclasso-python` vendors a pinned copy at
`cclasso-python/tests/r_reference/data/QMP.rda`. The expanded fixtures add
American Gut Project Round 1
and GlobalPatterns / Caporaso et al. 2011. Their deterministic selection rules
and profile metadata are in `scripts/public_microbiome_cases.R`; provenance and
redistribution status are recorded per fixture and audited in
`DATASET_PROVENANCE.md` (both are marked `redistribution_cleared: false`, so the
matrices are CI artifacts only). The fixtures include:

- `qmp_n106_p20`: all 106 samples and 20 low-prevalence taxa (61.51% zeros).
- `qmp_n20_p30`: 20 observed samples and 30 taxa (`p > n`; 32.0% zeros).
- American Gut Round 1: `80×20` (29.50% zeros) and `20×30` (`p > n`, 31.33% zeros), each with three R RNG seeds.
- GlobalPatterns: `26×15` and `26×40` (`p > n`), each with three R RNG seeds. The core cases select the most prevalent varying taxa (0% zeros); a separate prevalence-spread case is retained as a CCLasso failure stress fixture.
- Synthetic coverage cases: `counts=False` (three cases) and non-default `k_cv`/`n_boot`/`lam_int` (four cases), generated deterministically without third-party data.

On the pinned Linux x86-64 environment (R 4.6.1, NumPy 2.5.4, SciPy 1.18.1):

- **CCLasso:** both QMP cases match the R golden outputs within the existing
  contract, with finite reported variances, correlations, and p-values. The R
  algorithm and Python port both hit CCLasso's 200-iteration ADMM limit on all
  70 recorded QMP solves; the largest observed relative residual warning was
  about `1.97e-2`. This is a solver-convergence limitation shared with the R
  reference, not evidence that the estimate is fully converged.
- **SPIEC-EASI:** transform, correlation, lambda path, MB paths, StARS selection,
  and refit tests pass on both QMP subsets. In the p>n case the nodewise
  design has rank 19/30; coefficients match to `1.12e-14` absolute and support
  is exact. The maximum measured active KKT residual is `4.09e-5` (inactive
  violation is zero), and objective agreement is within `2.1e-14` relative.
  These invariants are checked in addition to pointwise coefficients.
- **SPRING:** the selected index and final adjacency agree exactly on both QMP
  subsets. For `qmp_n20_p30`, the data-derived lambda path differs by at most
  `4.02e-8` absolute (`4.87e-8` relative); MB support agrees exactly while the
  largest pointwise coefficient difference is `3.28e-5`. Both p>n solutions
  satisfy the solver's KKT tolerance (largest measured active residual
   `8.86e-5`), so that case is tested by support and KKT checks instead of
   pointwise coefficient equality.

Across the expanded committed fixture set, the equivalence suites pass: CCLasso
**49**, SPIEC-EASI **145**, and SPRING **140** tests (measured 20 Oct 2026).
CCLasso's GlobalPatterns prevalence-spread stress case is an explicit R/Python
matching failure status; standard public cases also compare identical
non-finite masks when the R bootstrap emits NA.

SPIEC-EASI passes discrete-output and coefficient comparisons on the selected
public cases. SPRING's selected index/adjacency and StARS outputs also pass. For
GlobalPatterns, the single pair taking bridge inversion in each selected case
differs between Python and R `approx` by `1.42e-5` (`p=15`) and `1.10e-5`
(`p=40`), shifting the data-specific lambda maximum by `1.27e-6` and `1.33e-6`
relative. Those SPRING public cases therefore inject the R-derived numeric
lambda path to isolate MB/StARS parity; QMP still exercises Python's end-to-end
data-specific path. On the bridge-sensitive cases, support is exact, active KKT
residuals stay below `9.1e-5`, and scaled objective gaps stay below `8.3e-7`;
coefficient differences reach `2.63e-4` (`p=15`) and `8.81e-4` (`p=40`). These
are invariant-based comparisons under slightly different covariance inputs,
not pointwise coefficient equality.

### SPRING bridge sensitivity follow-up

The R `approx` versus `original` comparison was measured on the committed QMP
`qdat` matrices with latentcor 2.0.2. It is **not** a numerical-noise floor:
`approx` uses a precomputed interpolation table for most pairs, while `original`
numerically inverts the bridge for every pair, so their difference includes the
intentional approximation error.

| Case | Pairs using interpolation / inversion in `approx` | Python vs R `approx`, interpolation max abs | Python vs R `approx`, inversion max abs | R `approx` vs `original`, interpolation max abs | R `approx` vs `original`, inversion max abs |
|---|---:|---:|---:|---:|---:|
| `qmp_n106_p20` | 190 / 0 | `2.22e-16` | n/a | `4.53e-3` | n/a |
| `qmp_n20_p30` | 434 / 1 | `3.33e-16` | `1.004e-5` | `7.39e-3` | `0` |
| `amgut_r1_n80_p20` | 190 / 0 | `2.78e-16` | n/a | `1.03e-2` | n/a |
| `amgut_r1_n20_p30` | 435 / 0 | `4.44e-16` | n/a | `5.35e-3` | n/a |
| `globalpatterns_n26_p15` | 104 / 1 | `1.22e-15` | `1.416e-5` | `9.09e-3` | `0` |
| `globalpatterns_n26_p40` | 779 / 1 | `1.22e-15` | `1.099e-5` | `9.16e-3` | `0` |

Thus, on these fixtures Python reproduces R's interpolation table path to
roundoff. In the single QMP p>n pair that invokes inversion, Python/R `approx`
differs by `1.004e-5`; R's `approx` delegates that pair to the same bridge
inversion used by `original`, so their outputs agree there. The large
`approx`/`original` differences in the table arise on interpolated pairs and
are intentional method differences, not a numerical-noise estimate. The
comparison can be regenerated with
`spring-python/tests/r_reference/compare_approx_original.R`.

In both GlobalPatterns cases, `approx` delegates that pair to the same R
inversion routine used by `original`, so the R-method difference on that pair is
zero; Python/R still differs by about `1.1–1.4e-5` because the numerical CDF
implementations differ. The much larger `9.1e-3` `approx`/`original` difference
is on table-interpolated pairs. These are method sensitivities, not acceptable
cross-language error allowances. The public SPRING fixtures inject R's numeric
lambda path and use solver invariants on the bridge-sensitive full-data fit.

These measurements support the stated contracts for the tested data and
environment only; they do not prove equality for every dataset or platform.

The known NumPy-seeded CCLasso non-finite case has now been replayed through the
pinned R v2.0 source with the exact same input, bootstrap indices, and artificial-
null draws. Both implementations produce one negative bootstrap variance at the
same taxon/replicate, finite mean `var_w`, and the exact same eight non-finite
positions in `cor_w` and `p_vals` (Python represents them as NaN; R represents
them as NA). Finite output entries agree within `1.02e-13` for `var_w`,
`5.23e-14` for `cor_w`, and `8.19e-14` for `p_vals`; the per-entry bootstrap
variance matrices differ by at most `2.31e-13` absolute (`1.99e-10` relative).
This closes the question of whether that particular failure is Python-specific:
it is shared R behavior following non-converged ADMM solves. It remains a
scientific-use limitation, not an acceptable basis for silently changing the
Python algorithm.

Diagnostic-only recording of the terminal ADMM relative residual was also added
to both the R golden generator and Python result diagnostics. Across the two QMP
fixtures, all 36 and 34 `cclasso_sub` calls respectively reach the 200-iteration
cap, with iteration counts agreeing exactly. The maximum terminal-residual
difference is `3.96e-14` absolute (`6.33e-12` relative). This provides stronger
evidence of shared CCLasso solver behavior, but also confirms that these QMP
solutions are cap-limited rather than converged to the requested tolerance. The
regression check uses `rtol=1e-10, atol=1e-13`, roughly an order of magnitude
tighter than the largest observed relative discrepancy while allowing for
rounding on the pinned environment. See
`EQUIVALENCE_ACCEPTANCE_CONTRACT.md` for the broader pass/fail policy.

## Reproduction

Run each port's tests from its package directory. The committed QMP goldens make
Python-only runs independent of R. To regenerate, use the corresponding
`tests/r_reference/generate_golden.R`; it requires R and the reference packages
listed in that port's `tests/r_reference/README.md`. By default, the R scripts
load QMP from `../SPRING-master/data/QMP.rda`; set `QMP_RDA` to override that
path.
