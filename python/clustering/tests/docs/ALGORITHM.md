# How clustering parameters, the selection ladder and signed edges work

This document explains, in plain language, the parameters of the
PageRank-Nibble clustering, the five-step ladder that selects the
returned cluster, and the signed-edge handling. The shipped defaults are
the reference values of the published experiments; the legacy app values
remain reproducible through explicit overrides.

The relevant code (in `python/clustering/`, one directory up):

| File | Purpose |
|---|---|
| `selection_ladder.py` | sweep scan and the five-step selection ladder |
| `clustering_profiles.py` | the `app` parameter set, overrides |
| `clustering_advanced.py` | signed score, advanced clustering entry |
| `bacteria_pagerank.py` | adaptive entry (parameters + ladder) |
| `bacteria_graph.py` | stores both `weight = |correlation|` and the sign |
| `run_clustering.py` | reads the options from the JSON payload |

For the experiments that produced the numbers (design, columns, how to
interpret them) see [METHODOLOGY.md](METHODOLOGY.md); the full
specification of the selection ladder - criteria, steps, parameters,
complexity - is [SELECTION_LADDER.md](SELECTION_LADDER.md).

## 1. The knobs, in one sentence each

* **`phi`** – the quality bar. A cluster is accepted in step 1 of the
  ladder only if its score (conductance, or signed score) is below
  `phi`. Lower = tighter clusters; if nothing clears the bar the ladder
  relaxes it and reports the bar it actually used.
* **`alpha`** – the teleport probability of the random walk. Small
  values produce larger, smoother clusters and *many more* push
  iterations (the iteration budget is consumed at rate `1/alpha`); large
  values produce smaller, more local clusters.
* **`b`** – the size scale. It sets the minimum cluster volume
  (`volume_min = 2^(b-1)`) and the accuracy of the PageRank
  approximation (`epsilon = 1/(2^b · 48 · ⌈log2 m⌉)`). Larger `b` =
  larger minimum size and more accurate vector. The ladder can halve
  the floor at step 3 without recomputing anything.
* **`epsilon`** – derived from `b` and the number of edges `m`; never
  set by hand.
* **`volume_min` / `volume_max`** – the allowed interval for the sum of
  degrees of the candidate cluster: `(2^(b-1), ⌊vol(V)/2⌋]`. The
  ceiling is half the total volume because conductance is invariant
  under complementation.
* **`tau`** – the minimum PageRank mass a candidate must carry to pass
  the *strict* test (`1/(48·⌈log2 m⌉)`); it is a certificate of the
  paper's proof and is dropped at ladder step 2 (see
  [SELECTION_LADDER.md](SELECTION_LADDER.md) §4 for why).
* **`min_cut_size` / `max_cut_size`** – hard size bounds of *every*
  ladder step (default: at least 5 taxa).

Cluster *selection* is not a parameter any more: the five-step ladder
(strict → quality → multi-scale → relaxed phi → penalty → honest
singleton) is the only rule, and the step that produced the answer is
reported as `ladder_step` / `ladder_rule` in the metadata.

## 2. Parameter set `app` (the only set, reference values)

```
phi    = 0.4
b      = compute_default_b(n, m)      # 10 on both bundled networks
alpha  = 0.05                         # constant, independent of phi
weights   = |correlation|
```

These are the values the Electron UI sends by default and the reference
values of the published experiments (the `Default` configuration of
RESULTS_TUNING / RESULTS_COMPARISON: strict acceptance 100 %, cluster
conductance 0.17-0.29 on the bundled networks):

* `b` is derived from the graph dimensions alone: target cluster ~10 % of
  the taxa, converted into a target volume with the average degree and
  rounded to the paper's scale;
* `alpha` is a constant because the paper's
  `phi^2 / (225·ln(100·sqrt(m)))` is a worst-case bound for the theorem,
  not an operating point: it yields ~2e-5, which barely localises the
  walk and burns the whole iteration budget.

The retired legacy defaults (`phi = 0.2`, `b = 4`, `alpha ≈ 2e-5`) remain
reproducible through explicit overrides; note that with `phi = 0.2`
network 101 (best conductance ≈ 0.28) never passes strictly and is served
by ladder step 4 (`phi_effective = 0.30`).

Every value can still be overridden from the UI. In the Python payload:

```json
{
  "targetConductance": 0.35,
  "b": 10,
  "teleportAlpha": 0.05
}
```

## 3. The five-step selection ladder

For every sweep prefix the ladder evaluates three admission criteria –
volume window, PageRank mass, quality below `phi` – and then, in
order:

1. **strict**: all three hold → return the lowest-conductance prefix
   (ties go to the smallest prefix);
2. **quality**: volume and quality hold, mass dropped;
3. **multi-scale**: the volume floor is halved (once or twice),
   reusing the same PageRank vector;
4. **relaxed phi**: the bar is raised `×1.25 … ×2.0` inside the
   original window, reporting the rung used;
5. **penalty**: the admissible prefix with the lowest normalised
   penalty (`score/phi`, octaves below the floor, small mass bonus) is
   returned together with the breakdown.

If nothing is admissible, the seed alone is returned
(`seed_singleton_fallback`) - an explicit "no cluster found" instead of
a degenerate micro-cluster. Full specification, formulas, parameter
defaults and measured costs: [SELECTION_LADDER.md](SELECTION_LADDER.md).

## 4. Signed edge weights

### Why

The graph stores correlations in `[-1, 1]`. About **43% (101) / 46%
(rs4-101)** of the edges are negative. The default method uses
`weight = |correlation|`, so two strongly anti-correlated taxa attract
each other exactly like two positively correlated ones.

With signs, the rule is weak balance: **positive edges should stay
inside a cluster, negative edges should cross it**.

### What the signed method does, step by step

1. **Store the sign.** `BacteriaGraph` already keeps the signed
   correlation in the edge data (`correlation`); `weight` remains
   `|correlation|` for backward compatibility. No migration is needed.
2. **Positive edges drive the walk.** The random walk runs on a copy of
   the graph containing only positive edges (`A+ = max(A, 0)`). All
   probabilities stay non-negative, so the original PageRank code can
   be reused.
3. **Negative edges never add rank.** They are not part of the
   transition matrix, so a negative link cannot pull a taxon into the
   cluster.
4. **Negative edges are counted as violations in the score.** For a
   candidate set `S`:

   ```
   violations(S) = positive edges leaving S
                 + lambda * negative edges inside S

   signed_score(S) = violations(S) / min(positive volume of S,
                                         positive volume of V \ S)
   ```

   Positive volume is the sum of positive edge weights incident to the
   nodes. `lambda` (default **1.0**) controls how strongly negative
   edges inside the cluster are penalised; `lambda = 0` ignores signs
   in scoring.
5. **The score replaces conductance in the ladder.** With
   `signed = true`, every ladder step compares `signed_score` instead
   of conductance (acceptance, relaxation rungs, and the step-5
   penalty); the metadata reports both `score` (signed) and
   `conductance` (plain `|r|`).
6. **Report the sign composition.** The metadata of a signed run
   contains `cut_pos`, `neg_inside`, `neg_cut`, `positive_inside`,
   `violations` and `vol_pos`. `results/RESULTS_COMPARISON.md` shows
   the signed score for the Signed and Both configurations; the raw
   JSON also contains it for the unsigned runs.

### Options in the payload

```json
{
  "signed": true,
  "signLambda": 1.0
}
```

### Practical caveat

SparCC/SpiecEasi negative correlations can be noisy or compositional.
If a negative edge is unreliable, lower `lambda` (e.g. 0.5) or filter
weak negatives before building the graph. A false negative edge should
not be allowed to split a genuine cluster.

## 5. How to read the results

Every cluster carries the step that produced it:

| `ladder_step` | meaning | key metadata |
|---|---|---|
| 1 | all three criteria hold | `strategy = target_conductance`, `conditions_met_count = 3` |
| 2 | mass certificate dropped | `failed_conditions = [prob_ok]` |
| 3 | volume floor halved | `volume_min_used < volume_min` |
| 4 | quality bar relaxed | `phi_effective > phi` |
| 5 | graded penalty fallback | `penalty_breakdown` |
| 0 | no admissible cluster | singleton, `seed_singleton_fallback` |

`conditions_met_count` and `failed_conditions` always refer to the
*strict* criteria, so "strict success rate" keeps its meaning across
steps. The stability metrics (size, volume, chance-corrected support,
consensus, the random conductance baseline) and the rules for telling a
non-trivial cluster from a lucky one are documented once in
[METHODOLOGY.md](METHODOLOGY.md). The reports live in `results/`; the
repository [README.md](../README.md) lists them together with the exact
commands that regenerate them.

## 6. Which parameter values should I start with?

* **`b = 10`, `phi = 0.4`** are the shipped defaults and the values of
  the reference experiments (strict acceptance 100 %, conductance
  0.17-0.29 on the bundled networks). `phi` only controls acceptance; `b`
  is a size dial above the argmin's volume: on 101, `b = 12` triples the
  cluster (50 → 171 taxa) and weakens cross-seed consistency, so stay
  with the dimension-based default unless you *want* larger clusters.
* **`alpha` is the only parameter worth tuning per condition.** From
  `results/RESULTS_TUNING.md`: `0.01` wins on rs4-101 (both weight
  modes) and on 101/unweighted, `0.10` on 101/weighted (nested support
  above chance 0.247 → 0.729); `0.05` is a good neutral default. The
  direction flips between networks, so pick it per network, not once
  for all.
* **`phi = 0.2` (the legacy default) is tight for network 101**: the
  ladder will regularly answer with step 4 (`phi_effective = 0.30`).
  That is correct, transparent behaviour - but if you want strict
  acceptance, stay at `phi = 0.4`.
* **Signed edges** help on rs4-101 but reduce cross-seed consistency on
  101. Do not enable them globally without checking the
  `support above chance` values in the Findings sections of the
  reports.

## 7. Naming in the UI / payload

`run_clustering.py` accepts:

* `targetConductance`, `b`, `teleportAlpha` – explicit parameter
  overrides (also accepted as `parameterOverrides`).
* `signed`: `true` / `false` (default `false`).
* `signLambda`: default `1.0`.
* `weighted`: `true` / `false`.

Removed: `parameterProfile` (only `app` exists; a payload asking for
`recommended` is rejected with an explicit error) and `selection`
(the ladder is the only rule).

The Electron side (`AnalysisManager._runPythonClustering`) forwards
these options, so a UI switch only needs to send the same keys.

## 8. Quick reference

| | value (bundled networks) |
|---|---|
| `phi` | 0.4 |
| `b` | 10 (dimension-based; both bundled networks) |
| `alpha` | 0.05 (constant) |
| `epsilon` | 1.57e-6 (at b=10) |
| `volume_min` / `volume_max` | 512 / ⌊vol(V)/2⌋ (on 101 weighted: 512 / 3320) |
| `tau` | 1.6e-3 |
| selection | five-step ladder (`min_cut_size = 5`) |
| edge weights | `|r|`, or signed if `signed=true` |
