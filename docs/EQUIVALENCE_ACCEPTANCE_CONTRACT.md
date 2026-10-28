# R/Python equivalence acceptance contract

**Status: preregistered before the follow-up investigations in this change.**
This contract fixes the scope and decision rules before adding expanded public-
dataset equivalence fixtures. The initial CCLasso and SPRING diagnostic results
are recorded below and in `NUMERICAL_EQUIVALENCE.md`; they do not satisfy the
broader dataset/stress-suite gates. If a criterion needs to change, record the
reason and date before inspecting the affected result, and retain the original
criterion in the report.

## Claim scope

The intended claim is method-specific and limited to the code paths implemented
in these ports, compared with the named R reference versions, under the pinned
Linux environment recorded in each fixture manifest. It is not a claim of
universal mathematical equality, cross-platform identity, or equivalence of
unported options.

- **CCLasso:** CCLasso R 2.0 API and behavior represented by the in-repository
  `reference/R/cclasso.R`.
- **SPIEC-EASI:** the implemented CLR + Meinshausen-Bühlmann + StARS/refit path.
  `glasso`, `slr`, `coat`, `ising`, and `poisson` are out of scope.
- **SPRING:** the implemented SPRING MB + StARS path and `latentcor` type
  `tru`; `con`, `bin`, `ter`, and `ord` are out of scope.

The Electron app, PyInstaller packaging, and app exports are not part of the
scientific equivalence claim.

## Data and test design

1. Retain the two existing QMP fixtures as the initial reference set.
2. Add at least **two independent public microbiome datasets** after verifying
   source, license, access method, sample/taxon annotations, and preprocessing.
   Prefer datasets distributed with the R reference packages when their
   provenance and usage terms permit. Record any filtering, aggregation, and
   sample selection as executable scripts; do not hand-edit matrices.
   Initial local candidates are inventoried in `PUBLIC_DATASET_INVENTORY.md`:
   American Gut Round 1 and GlobalPatterns. American Gut Round 2 is a useful
   sensitivity case but does not count as a separate study from Round 1. HMP
   16S remains a reserve candidate pending access to its `hmp2` package.
3. Across datasets, cover sparse and less-sparse tables, variable library sizes,
   `p < n`, and `p > n`. Use a compact, predeclared design of representative
   cases rather than an unreported post-hoc subset of a large grid.
4. For each stochastic R path, inject the same saved R draws/indices into Python.
   Use at least three independently generated R RNG seeds for each stochastic
   case family. Preserve inputs, parameters, injected randomness, outputs, and
   environment manifests with checksums.
   For SPRING's added public-data cases, inject the R-derived numeric lambda
   path to isolate MB/StARS parity; retain QMP data-specific cases as the test of
   Python's end-to-end lambda-max/path calculation. Report these scopes
   separately.
5. Include deterministic edge cases (ties, constant/near-constant taxa,
   singular designs, threshold-near cases) as explicitly labeled stress cases;
   do not count them as independent real datasets.

## Comparison and pass/fail rules

There is no universal tolerance that makes a numerical port "academic quality".
Use the following rules; justify every numerical tolerance against the algorithm,
reference precision, and observed error distribution. Do not widen tolerances
after a failing result simply to obtain a pass.

### Discrete and control-flow outputs

- Require exact equality for shapes, selected lambda/index, supports, and binary
  adjacency on cases whose decision margins are demonstrably larger than the
  measured floating-point discrepancy.
- For threshold-near cases, report the decision margin and both outputs. A
  difference is classified as boundary-sensitive only if an independent
  perturbation/resolution analysis supports that explanation; otherwise it is an
  unexplained failure.
- Require the same convergence/failure category: converged, reached iteration
  cap, returned non-finite output, or raised an error. No Python-only failure or
  Python-only non-finite output is acceptable in the claimed domain.

### Continuous outputs

- Compare deterministic transforms and intermediate matrices with tight
  absolute/relative tolerances established from the pinned reference run.
- Compare solver outputs using the existing method-specific tolerances, plus
  support, objective, and KKT checks for non-unique/rank-deficient solutions.
  Report max absolute and relative errors and their locations for every case.
- For methods involving numerical special functions, separate regions using
  different computational paths (e.g., SPRING interpolation vs bridge inversion)
  and report the error by region. Compare the Python/R discrepancy with R's own
  `approx`/`original` sensitivity where applicable; do not present that
  sensitivity as proof if the implementations use different definitions.

### Failure and stochastic behavior

- Replay every known failure with identical inputs and injected random arrays in
  both languages. Classify shared failures separately from port divergence.
- Run a predeclared seeded stress suite of at least 1,000 small cases per method
  family, weighted toward sparse, singular, and threshold-near inputs. Any
  method-specific numerical failure must be triaged and retained as a regression
  fixture; the sweep is supporting evidence, not proof over all inputs.
- Native-RNG distribution comparisons are supplementary. They cannot replace
  injected-RNG deterministic comparisons where injection is possible.

## Environment and reproducibility

Record OS/architecture, R version and reference package versions, Python version,
NumPy/SciPy/scikit-learn versions, BLAS/LAPACK vendor and thread settings, and
whether R-generated randomness was injected. Pin the tested environment and
make fixture regeneration and Python-only test commands executable from a clean
checkout. No cross-platform equivalence claim is made because only one
environment is available.

The current workspace reproduction entry point is
`bash scripts/reproduce_equivalence.sh <new-output-directory>`. It regenerates
all three R fixture trees in a fresh destination, runs their Python equivalence
suites, writes R/Python environment manifests, and checksums each result tree.
Set `PY_CCLASSO`, `PY_SPIEC`, and `PY_SPRING` when the three packages use
different Python environments. `cclasso-python/` is vendored into the root
repository and pins its own `QMP.rda`, so its component reproduces from a clean
checkout. The combined command still depends on `spiec-easi-python/` (ignored
nested Git repository) and `spring-python/` (untracked), which need the same
vendor/submodule integration before a clean checkout can run all three.

## Readiness decision

A port is **citable as equivalent within the stated scope** only when:

1. every declared supported path has passing deterministic comparisons over the
   predeclared public-data cases and stress cases;
2. there are no unexplained differences in outputs, convergence category, or
   failure status;
3. boundary-sensitive exceptions and shared algorithm limitations are
   enumerated with measurements rather than hidden in aggregate scores;
4. a clean checkout can reproduce the fixtures and reports; and
5. paper wording names the R version, Python port version, tested cases, and
   environment, without implying untested inputs/platforms are covered.

If a gate is not met, the result is **not yet established**—not automatically a
bad implementation. The paper may make narrower claims for the methods/cases
that pass, provided the scope and exclusions are explicit. Statistical
equivalence under independent RNG streams is a separate claim requiring a
predeclared distributional analysis and should not be conflated with numerical
equivalence.

## First investigation results and remaining work

- **CCLasso replay completed:** the exact NumPy-seeded case is non-finite in both
  Python and pinned R v2.0 at the same bootstrap replicate. On QMP fixtures,
  iteration counts match exactly and maximum terminal-residual difference is
  `3.96e-14` absolute. The shared non-finite behavior remains a limitation; do
  not silently change the default algorithm.
- **SPRING bridge sensitivity measured:** interpolation and inversion regions
  have been separated on the two QMP fixtures. Python reproduces the R
  interpolation path to roundoff, while one p>n pair uses inversion. The R
  `approx`/`original` difference is an algorithmic sensitivity comparison, not a
  general error floor. See `NUMERICAL_EQUIVALENCE.md` for the measurements.
- **Still open:** add at least two eligible independent public datasets, execute
  the predeclared representative test design and seeded stress suite, and
  produce the clean-checkout environment/fixture reproduction. SPIEC-EASI keeps
  exact discrete-output plus objective/KKT requirements for the implemented MB
  path.

No solver or numerical algorithm is to be optimized or changed as part of these
investigations. Instrumentation may record diagnostics only if it does not alter
the operation order or returned scientific values.
