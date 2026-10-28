# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

- **Statistical ports are now git submodules.** `cclasso-python/`,
  `spiec-easi-python/` and `spring-python/` are standalone git repositories
  (released independently) integrated as submodules instead of vendored copies
  / a nested ignored repo. The app pins one commit per port; clone with
  `git clone --recurse-submodules` (or run `git submodule update --init`).
  Packaging is unaffected: PyInstaller still freezes the ports into
  `python-dist`, so releases remain self-contained.
- **SPIEC-EASI backend replaced.** The standalone Python library in
  `spiec-easi-python/` was rebuilt from scratch as a faithful port of the R
  reference (`SpiecEasi` 1.99.0 + `huge` 2.0.1 + `pulsar` 0.3.13), `method="mb"`
  only. The previous implementation — with its pybind11 `_mb_core` accelerator,
  joblib `prefer="threads"` backend and `spieceasi.mb_solver` /
  `spieceasi.runtime` helpers — was removed. `python/analysis.py` now calls the
  new public API `spieceasi.spiec_easi(...)`.
- **Worker count.** The StARS subsample loop now runs in a process pool
  (`spieceasiJobs`); the subsamples are independent and aggregated in a fixed
  order, so results are identical for any worker count. The joblib
  `threads`-vs-`processes` backend choice and the compiled core are gone.
- **`spieceasiMaxIter` is no longer used.** The port follows R's fixed MB solver
  iteration limit (`MAX_ITER = 10000`); the previous per-case cap was a deviation
  from R and has been removed. The `spieceasiNLambdas` / `spieceasiNSubsamples`
  reductions in the tiny/high-dimensional regime (`n_samples <= 10` and
  `n_taxa >= 200`) are kept.
- **CCLasso "input model" option removed.** CCLasso now always receives the same
  abundance table as every other correlation method and applies its own
  0.5 pseudo-count closure internally (`counts=True`). The old
  "strictly positive fractions" mode rejected zero-containing tables and had no
  counterpart in the other methods, so it is gone (`cclassoCounts` is no longer
  sent by the renderer or read by `python/analysis.py`).
- **SPRING "data model" option removed.** SPRING now always receives the same
  abundance table as every other method and applies its compositional transform
  (MCLR) internally (`quantitative=False`); `springQuantitative` is no longer
  sent by the renderer or read by `python/analysis.py`.
- **SPRING threshold label unified.** The shared `correlationThreshold` field is
  now labelled "Correlation filter · min |r|" for SPRING too (its edge weight is
  a latent correlation). SPIEC-EASI keeps "StARS edge stability · minimum"
  because its edge weight is a selection frequency, not a correlation. The graph
  provenance label for SPRING is now `minimum_absolute_correlation`; the
  distinction remains recorded in `weight_type` ("latent_correlation").
- **Cluster selection is now a five-step ladder.**
  `python/clustering/selection_ladder.py` replaces the "first prefix below φ"
  rule and the vote-count fallback rank: strict → quality (mass certificate
  dropped) → multi-scale (volume floor halved, reusing the same PageRank
  vector) → relaxed φ (×1.25 … ×2.0 rungs, teleport frozen) → normalised
  penalty (`score/φ`, octaves below the volume floor, small mass bonus) with a
  hard size floor → honest seed singleton when nothing is admissible. The step
  that produced the cluster is reported as `ladder_step` / `ladder_rule`, with
  `phi_effective`, `volume_min_used`, `penalty_breakdown`, `push_iterations`
  and `push_converged` added to the metadata. The methods specification is
  `python/clustering/tests/docs/SELECTION_LADDER.md`.
- **Sweep-order consistency.** The sweep orders nodes by `p(v)/d(v)` using the
  same degree measure as the volumes and the conductance (previously the order
  always used unweighted degrees, even in weighted mode); ties are broken by
  node id, so the order is a deterministic function of the vector.
- **Volume window capped at half.** `volume_max` is `⌊vol(V)/2⌋` instead of
  `⌊2·vol(V)/3⌋`: conductance is complement-invariant, so prefixes beyond half
  of the volume only duplicate smaller sets, and every admissible prefix has
  conductance in `[0, 1]`.
- **Size bounds enforced.** `min_cut_size` (default 5 taxa — previously
  accepted and discarded with `del`) is a floor of *every* ladder step,
  `max_cut_size` an optional ceiling, and every candidate must contain the
  seed.
- **Sweep evaluation is incremental.** Each prefix's cut is updated in
  `O(deg)` instead of recomputing the whole boundary: on network 101 the sweep
  plus the ladder cost 2.4 ms against 667 ms of push. The redundant inner
  `pagerank_nibble` best-cut computation (discarded by its caller) is gone.
- **Clustering defaults are now the reference experiment values.** The app
  sends `phi = 0.4`, the dimension-based scale `b = compute_default_b(n, m)`
  (10 on both bundled networks, derived from a target cluster of ~10% of the
  taxa) and a constant `alpha = 0.05` — the `Default` configuration of the
  published tuning/comparison reports (strict acceptance 100%, conductance
  0.17–0.29). The legacy values (`phi = 0.2`, `b = 4`,
  `alpha = phi²/(225·ln(100√m)) ≈ 2e-5`, visible in the UI as `0.000021`)
  are retired: on network 101 they never pass strictly (the ladder answers
  with step 4) and the tiny teleport consumed the whole iteration budget.
  The same recipes live in `python/clustering/clustering_profiles.py` and in
  `AnalysisWorkspace.jsx` (`computeDefaultB` / `computeDefaultTeleportAlpha`);
  `bacteria_pagerank` function defaults follow, and the experiment harness
  keeps passing its base values explicitly, so the published reports stay
  reproducible.
- **Experiment harness decoupled from profiles.** `CONFIGS` in
  `tests/stability_experiment.py` pass their base values (`alpha=0.05`,
  `b=10`, `phi=0.4`) as explicit overrides (`EXPERIMENT_DEFAULTS`), so the
  published experiments remain reproducible with a single parameter set; the
  report generators' prose follows the new code and `rebuild_reports.py`
  re-renders the Markdown from the unchanged historical JSON.
- **Taxa-filtering thresholds now use the same scale as the model in combined
  mode.** With `datasetType="both"` the *mean relative abundance* and
  *variability* filters closed each sample on its **own domain** total, while
  the correlation methods close the fused `bac::`/`fun::` matrix on bacteria +
  fungi **together** — so the minority domain looked ~100x more abundant to the
  filter than to the model (and its variance ~10^4x, i.e. those thresholds were
  effectively that much looser for it). `_apply_taxa_filters` now takes an
  optional `sample_totals` and `run()` passes the combined per-sample total
  whenever both domains are loaded (a sample absent from one domain contributes
  0, exactly as in the model); a sample missing from the caller's totals falls
  back to its domain total. *Prevalence* is scale-free and is unaffected. The
  report records the choice in `relative_abundance_denominator`
  (`combined_bacteria_fungi` | `domain`) together with `scope`
  (`all_samples_in_abundance_table` — the thresholds stay deliberately global
  so both compared conditions start from the same candidate taxon set),
  `removed_share_percent` and `aggressive_removal`, and a warning is logged
  when a threshold removes more than 60% of a domain.

### Fixed

- **Clustering no longer crashes the app on degenerate inputs.** The backend
  could emit `Infinity` in its JSON stream — the honest Python value of an
  undefined conductance on the seed-singleton fallback (a graph with fewer
  than `min_cut_size` taxa, an edge-less condition graph, or a non-positive
  `phi`) — and `JSON.parse` in the Electron side rejects the whole response
  with "Unexpected token 'I'". `run_clustering.py` now sanitises the payload
  (`Infinity`/`NaN` → `null`, `allow_nan=False` as a hard guard), clamps a
  non-positive `targetConductance` back to the app default and records
  `phi_clamped_from`; and the renderer treats an *empty* conductance/`b`/
  teleport field as "unset" instead of sending `Number('') === 0` (which
  produced `phi=0`/`alpha=0` runs), falling back to `null` so the backend
  applies its own defaults. Verified end-to-end through both backend modes
  (frozen and source) plus a 12-seed fuzz run on network 101.
- **Benign solver output no longer shown as an error.** The run log painted any
  line containing the substring `error` in red, so CCLasso's ADMM diagnostics
  (`Maximum Iteration: 200 && Relative error: …`) and the `Optimal lambda is near
  boundary` notice looked like failures. Warning levels are now decided by
  explicit failure markers, and `python/analysis.py` captures CCLasso's
  `RuntimeWarning`s and re-emits them as readable `CCLasso solver note: …` lines
  instead of raw tracebacks on stderr.
- **The taxa-filtering panel no longer claims "every taxon is used".** With the
  master toggle (or all individual filters) off, the readout promised "every
  taxon is used" / "the full taxon set", yet every backend still discards
  constant and all-zero taxa (zero variance cannot produce a correlation) and
  SPRING also drops all-zero samples. `describeTaxaFiltering` and the panel now
  say that explicitly, state the global scope of the thresholds (all samples of
  the abundance table, every condition — not only the selected ones) and, in
  combined mode, which denominator the mean/variability thresholds use; the
  prevalence read-out is labelled as the metadata-count estimate it is.
- **Method-side removals are now in the results payload.** Taxa dropped as
  constant inside a condition (and SPRING's empty samples) existed only as log
  lines, so `taxa_filtering.<domain>.before/after` disagreed with the number of
  taxa the model actually saw, with no way to reconcile them from the results.
  Each backend drop is now recorded with its method, condition/scope, item kind
  and item list and exposed as `taxa_filtering.meta.method_removals` (the exact
  model input size stays in `graph.stats.method_diagnostics.n_taxa`).

### Removed

- References to the deleted library internals (`spieceasi._mb_core`,
  `spiec-easi-python/docs/AUDIT.md`, `spiec-easi-python/docs/PERFORMANCE.md`,
  the extension-vs-reference sweep) from `README.md`, `CHANGELOG.md` and
  `python/backend.spec`.
- **The `recommended` parameter profile.** Gone from code
  (`clustering_profiles.recommended_parameters`, the `RECOMMENDED_*` constants,
  the `run_clustering.py` branch) and from documentation (`ALGORITHM.md` §3 and
  the payload examples). A payload still asking for
  `parameterProfile: "recommended"` is rejected with an explicit error; the
  values live on as the experiments' reference parameters
  (`EXPERIMENT_DEFAULTS`), so the published reports stay reproducible.
- **The `selection` option and the `first` / `min_conductance` rules**, both in
  `pagerank_cluster` and in the payload (renderer and `AnalysisManager` no
  longer send `parameterProfile` / `selection`): the five-step ladder is the
  only selection rule.

### Notes

- Equivalence with R is documented in `spiec-easi-python/docs/EQUIVALENCE.md`
  (what is identical, what is identical within a numerical tolerance, what cannot
  be); deviations are in `spiec-easi-python/docs/DIFFERENCES.md`.
- No comparison tolerance was widened to make a test pass.
