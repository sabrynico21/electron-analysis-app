# Statistical-method release readiness

**Verdict: not yet ready for an unqualified scientific-paper release.** The
ports now have reproducible R-equivalence evidence on QMP and two additional
public study sources, including multiple seeded stochastic cases. That still
does not establish universal correctness, cross-platform behavior, or
publication-grade performance over the declared input domain. See
`NUMERICAL_EQUIVALENCE.md` and `TEST_GAP_ANALYSIS.md` for the measured details.

## Method-specific blockers

- **CCLasso:** QMP and the added public cases match their R outputs/statuses,
  but all 70 recorded ADMM solves across the two QMP fixtures reach the
  200-iteration cap (largest observed relative residual warning about
  `1.97e-2`). The public GlobalPatterns cases also expose shared R/Python
  non-finite bootstrap results and a shared missing-logical error in the
  prevalence-spread stress case. A separate seeded
  bootstrap can produce non-finite correlations/p-values after a negative
  bootstrap covariance diagonal. Exact replay now confirms the pinned R source
  produces the same failure at the same bootstrap replicate; this is shared
  algorithm behavior, not evidence of a Python-only mismatch. It remains a
  limitation: affected outputs must not be treated as valid estimates. The
  predeclared 1,000-case seeded stress sweep has been executed (seed 20261019,
  evidence in `cclasso-python/tests/stress/2026-10-20-seed20261019/`): both
  languages classify 525 finite / 394 nonfinite / 81 error with **zero status
  mismatches**, so the failure-category gate is met on this design while the
  shared non-finite and error behaviour remains a documented limitation, not a
  port divergence. The earlier fixture-coverage gap is closed: four synthetic
  goldens now exercise `counts=False` and non-default `k_cv`, `n_boot` and
  `lam_int`, and a test asserts that coverage persists. Public-data provenance
  is recorded per fixture and audited in `DATASET_PROVENANCE.md`; the American
  Gut Round 1 and GlobalPatterns matrices remain `redistribution_cleared:
  false`, so they are CI artifacts only until terms are confirmed.
- **SPIEC-EASI:** R equivalence is strong for the tested MB path, but StARS
  selects network edges; it does not provide pairwise p-values. App output now
  preserves StARS evidence instead of inventing `p=0/1`. Report the network
  selection rule and tuning settings, not a p-value interpretation. The
  predeclared 1,000-case seeded stress sweep has now been executed for the MB
  path (seed 20261020, evidence in
  `spiec-easi-python/tests/stress/2026-10-20-seed20261020/`): failure categories
  agree in 1,000/1,000 cases, the selected index differs in 21/1,000 (all p=2
  extreme-index ties), and the refit adjacency differs by exactly one edge in
  256/1,000, always the maximum-`|cor|` pair at `lambda.max` - the boundary tie
  already documented for `huge` MB. The sweep therefore supports the existing
  scoped claim and adds no new failure class; it does not by itself justify a
  universal-equivalence claim.
- **SPRING:** QMP data-specific-path tests and expanded public-data MB/StARS
  comparisons pass. On GlobalPatterns, one pair per selected case takes numeric
  bridge inversion; Python/R pointwise differences are about `1.1–1.4e-5`, which
  shift the computed lambda maximum by about `1.3e-6` relative. Public-data
  tests inject the R-derived path to isolate MB/StARS; only QMP tests the
  Python-derived data-specific path. Bridge-sensitive coefficients are checked
  by exact support, KKT, and objective invariants, not pointwise equality.
  SPRING's StARS network selection is not a pairwise p-value test. The
  predeclared 1,000-case seeded stress sweep has been executed (seed 20261020,
  evidence in `spring-python/tests/stress/2026-10-20-seed20261020/`): failure
  categories agree in 1,000/1,000 cases, the selected index agrees in
  1,000/1,000, and the refit adjacency differs in exactly one case where the
  differing edge has a soft-threshold-boundary coefficient. Cases whose full
  data or replicas enter latentcor's bridge-inversion region were excluded by
  a predeclared deterministic redraw (78 candidates); that region is covered by
  the QMP goldens and `compare_approx_original.R`, not by this sweep.

## Release options

1. **Recommended near-term: experimental/qualified release.** Keep all three
   choices available, retain the method-specific evidence labels, expose the
   exact settings and software versions in exported results, and clearly state
   that CCLasso may fail to converge on difficult data. Require users to inspect
   warnings and non-finite values before interpreting results.
2. **Before making paper-wide correctness claims:** QMP plus two independent
   public study sources now have fixtures, including `p < n`, `p > n`, and
   varied sparsity/library-size profiles. Public-data provenance and
   redistribution status are recorded per fixture and audited in
   `DATASET_PROVENANCE.md`; American Gut Round 1 and GlobalPatterns are marked
   `redistribution_cleared: false`, so do not publish or ship their matrices
   until permission is obtained or they are replaced with synthetic data.
   Broaden dimensions/sparsity before claiming
   more than these measured cases. All three predeclared 1,000-case stress
   sweeps are now recorded (`cclasso-python/tests/stress/2026-10-20-seed20261019/`,
   `spiec-easi-python/tests/stress/2026-10-20-seed20261020/`,
   `spring-python/tests/stress/2026-10-20-seed20261020/`); they are supporting
   evidence for the scoped claims, not proof over untested inputs. The SPRING
   sweep excludes the bridge-inversion region by design; that region remains a
   measured sensitivity, not a swept case family.
3. **For CCLasso:** investigate a higher iteration limit or a convergence
   safeguard, but treat this as an algorithm change. Compare the altered solver
   to the R reference and report convergence criteria/results separately rather
   than silently weakening parity tolerances.
4. **For numerical reproducibility:** pin the R/Python/BLAS environments for
   reference runs, preserve input and parameter hashes, and report exact versus
   tolerance-level outputs separately. Avoid claiming byte identity across
   platforms.
5. **Licensing review is required before distributing the combined app.** The
   app is MIT-licensed, CCLasso's Python port is LGPL-3.0-or-later, and the
   SPRING Python port is GPL-3.0-only. Bundling these into the Electron/PyInstaller
   product may create source-distribution and licensing obligations; do not
   assume the existing app license alone is sufficient. SPIEC-EASI's standalone
   port is GPL-3.0-only as well.
   Options to resolve this are: (a) release the combined application under a
   GPL-compatible license and meet the corresponding source/notice obligations;
   (b) obtain alternative licensing from the relevant rights holders; or
   (c) keep the MIT app and methods as separately distributed components only
   after counsel confirms that architecture and distribution model are compliant.
   Merely moving a GPL method to another process is not, by itself, a legal
   conclusion.
6. **Pin the method sources for reproducible builds.** `cclasso-python/` is now
   **vendored** as ordinary tracked files in the root repository (see
   `cclasso-python/VENDORED.md`), including a pinned `QMP.rda` copy at
   `cclasso-python/tests/r_reference/data/QMP.rda`, so its fixtures and tests
   are reproducible from a clean checkout. `spring-python/` remains untracked
   and `spiec-easi-python/` is still an ignored nested Git repository; a clean
   root-repo checkout still needs those two vendored or pinned as submodules
   before the combined Electron integration is reproducible.

7. **Current port-level follow-up:** the CCLasso NumPy-seeded non-finite case
   now reproduces in R with identical injected inputs/draws, and QMP terminal
   ADMM residuals match across languages to at most `3.96e-14` absolute. This
   resolves the Python-only divergence question but does not make affected
   CCLasso estimates valid. SPRING bridge-path sensitivities are measured in
   `NUMERICAL_EQUIVALENCE.md`; their small fixture count is not a universal
   error bound.

## Integration semantics

The method selector now distinguishes pairwise correlation methods from network
selection methods. Spearman, SparCC, and CCLasso have pairwise p-value evidence;
SPIEC-EASI and SPRING use StARS-selected edges and do not. For SPIEC-EASI the
edge score is StARS selection frequency; for SPRING it is the signed latent
correlation on the StARS-selected network. The p-value threshold is disabled for
the two network methods, and exported network edges have an empty/null p-value.

This prevents the interface from equating network selection with statistical
significance, but it does not substitute for independent scientific validation.

## Method-equivalence follow-up

- The exact seeded CCLasso bootstrap failure now replays in R with identical
  input and draws. Both implementations have the same negative bootstrap
  variance and non-finite result positions; Python emits NaN while R emits NA.
  The failure is shared, but those estimates remain unusable.
- Across the two QMP cases, all 70 CCLasso ADMM solves hit the 200-iteration
  cap. R/Python iteration counts are exact and terminal residuals differ by at
  most `3.96e-14` absolute. This confirms parity of cap behavior, not convergence
  of the estimates.
- SPRING bridge sensitivity has been split by interpolation/inversion path. The
  `approx`/`original` difference is not a numerical-noise floor; see
  `NUMERICAL_EQUIVALENCE.md` for measured results.
- QMP, American Gut Round 1, and GlobalPatterns fixture families now exist, with
  three R RNG seeds per added stochastic case family. HMP 16S remains an
  optional reserve dataset. Provenance and redistribution status for every
  public source are audited in `DATASET_PROVENANCE.md`; American Gut Round 1 and
  GlobalPatterns are recorded as **not cleared** for redistribution and are
  retained only as CI fixtures. Four synthetic goldens now cover `counts=False`
  and non-default `k_cv`/`n_boot`/`lam_int`.
- SPIEC-EASI's expanded public-data suite passes 145 tests; SPRING passes 140;
  CCLasso passes 49 including matched R failure-status cases and the new
  `counts=False`/varied-parameter synthetic cases (measured 20 Oct 2026). These
  results are recorded in each package's golden manifest/checksum tree.
- A fixture-reproduction script and R/Python environment manifests have been
  added. `cclasso-python/` is vendored into the root repository and regenerates
  its fixtures from a clean checkout using the pinned `QMP.rda`. The combined
  three-method reproduction still depends on the ignored nested SPIEC-EASI and
  untracked SPRING source trees, which must be vendored or pinned before the
  whole app is reproducible.
- No solver/algorithm optimization or behavior change was made.

## Integration checks in this change

- Python method/backend tests: **6 passed, 1 skipped** (the clustering-edge
  null-p-value test is skipped in the CCLasso test venv because that venv lacks
  NetworkX; it was separately smoke-tested in the app venv). The passing tests
  include real calls to all three methods and a run-level check that StARS
  selections bypass p-value filtering.
- Electron Jest unit tests: **18 passed**; ESLint and the Vite renderer build
  passed.
- The Electron E2E suite is **not green** in this workspace: existing file
  deduplication/history expectations fail and several CDP calls time out. The
  method-selector option assertion was added, but the history suite did not
  reach it. An isolated Linux PyInstaller build and frozen-executable smoke test
  for all three methods passed from `/tmp`; the app's `python-dist` artifact was
  not replaced, and a complete electron-builder installer was not built.
