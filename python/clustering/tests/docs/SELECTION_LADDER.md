# The five-step selection ladder for PageRank-Nibble

*This document specifies the cluster-selection stage of the clustering
pipeline in `python/clustering/`. It replaces the previous
first-prefix / minimum-conductance rules and is intended to be quotable
as-is in a methods section: every criterion, relaxation and default is
defined below with its symbol, its value and its justification.*

Related documents: [ALGORITHM.md](ALGORITHM.md) (parameters, signed
edges, payload), [METHODOLOGY.md](METHODOLOGY.md) (the stability
experiments), [../../results/README.md](../../results/README.md)
(reports and raw data).

---

## 1. Problem statement

Given an undirected weighted graph *G = (V, E)* (bacteria taxa as
nodes, `weight = |correlation|` as edge weights), a seed taxon
*v₀ ∈ V* and a quality bar *φ*, the algorithm must return a cluster
*S ∋ v₀* and an honest account of how good *S* is. Local spectral
partitioning answers this by computing an approximate PageRank vector
from *v₀*, ordering the nodes by their *density* *p(v)/d(v)* and
examining the nested prefixes (the *sweep*) of that order. The
remaining question — **which prefix should be returned?** — is answered
by the five-step ladder implemented in `selection_ladder.py`.

The contract is that the pipeline *always* returns a cut, and that the
metadata states precisely which of the three admission criteria
(Sec. 4) the returned cut satisfies. A degenerate answer (a
three-taxa clique, or "half the microbiome") is never returned
silently: it is either excluded by construction or reported through
`ladder_step` and `failed_conditions`.

## 2. Notation

| Symbol | Meaning | Code |
|---|---|---|
| *d(v)* | degree of *v*; weighted degree `Σ|correlation|` in weighted mode, edge count otherwise | `graph.get_deg(weight=...)` |
| *vol(S)* | `Σ_{v∈S} d(v)` | prefix `volume` |
| *cut(S)* | total weight of edges with exactly one endpoint in *S* | prefix `cut` |
| *φ(S)* | conductance = `cut(S) / min(vol(S), vol(V∖S))` | `compute_conductance` |
| *score(S)* | `φ(S)`, or the signed score when signed clustering is enabled (Sec. 3.3) | prefix `score` |
| *p(v)* | approximate PageRank mass on *v* (a sub-probability vector: `Σp ≤ 1`) | `PageRankResult.scores` |
| *m*, *n* | number of edges, number of nodes | — |
| *b* | scale parameter (integer, powers of two) | `LadderParameters.b` |

All quantities below are computed in the same degree measure: in
weighted mode the sweep order, the volumes and the conductance all use
weighted degrees; in unweighted mode all of them use edge counts.

## 3. Stage 1 — the approximate PageRank push

One linear-approximation push (Andersen, Chung & Lang, 2006) computes
*p* from the seed with teleport probability *α* and residual threshold
*ε*, where

```
ε(b, m) = 1 / (2^b · 48 · ⌈log₂ m⌉)          (paper's accuracy bound)
α       = teleport probability                (parameter; see Sec. 7)
```

On the two bundled networks (`m = 5121` / `5223`, `⌈log₂m⌉ = 13`):
*ε(4) = 1.00·10⁻⁴*, *ε(10) = 1.57·10⁻⁶*, *ε(13) = 1.96·10⁻⁷*.

The push is capped at `max_iterations = 100 000`. The cap is reported,
not hidden: `metadata.push_iterations` and `metadata.push_converged`
record whether the residual actually converged. The iteration count is
governed by *α*, because each pop absorbs only `α·r` of the residual
mass (Sec. 9).

### 3.1 Signed mode

With `signed = true` the walk runs on the positive sub-graph
(`A⁺ = max(A, 0)`), so all transition probabilities stay non-negative,
and the *score* replaced by the conductance everywhere in the ladder is

```
score(S) = (cut⁺(S) + λ · negInside(S)) / min(vol⁺(S), vol⁺(V∖S))
```

with `λ = signLambda` (default 1.0). Positive edges outside the cluster
and negative edges inside it are both violations. The *reported*
`conductance` key stays the plain `|correlation|` conductance of the
returned set; the signed value is reported as `score`.

## 4. The three admission criteria

For every sweep prefix *S* the ladder evaluates:

1. **Volume** — `2^(b−1) < vol(S) ≤ ⌊vol(V)/2⌋`.
   The lower bound is the paper's scale guarantee: it forbids
   degenerate tiny prefixes. The upper bound is **half** of the total
   volume (the previous implementation used `2/3`): conductance is
   invariant under complementation — `φ(S) = φ(V∖S)` — so every prefix
   beyond half of the volume merely duplicates the value of a smaller
   set, and the half-cap additionally bounds the conductance of any
   admissible prefix to `[0, 1]`.
2. **Mass** — `p(S) > τ` with `τ = 1 / (48·⌈log₂m⌉)`. On both bundled
   networks `τ = 0.001603` against a total mass `Σp ≤ 1`.
3. **Quality** — `score(S) < φ`.

Two structural guards apply to **all** steps:

* **seed membership** — a prefix is admissible only if it contains the
  seed (the sweep order is not guaranteed to start with it);
* **size floor** — `|S| ≥ min_cut_size` (default **5** taxa); an upper
  bound `|S| ≤ max_cut_size` applies when the caller sets it. The
  admissible range may legitimately be empty, in which case the ladder
  falls through to the singleton (step 0).

The mass certificate deserves a remark. It is a device of the paper's
proof, and it is *near-vacuous* when the push converges normally: `τ`
is `0.0016` while the vector sums to about one, so essentially any
prefix passes. It only *binds* in the opposite regime — when *α* is so
small (`α = 2·10⁻⁵`, the retired legacy default) that the push hits the
iteration cap before absorbing the mass, and `Σp` itself falls below
`τ`. In other words the certificate never separates good from bad
clusters; it only ever reports "the walk did not finish". This is
exactly why step 2 drops it from the acceptance test and keeps it as
metadata.

## 5. Stage 2 — the sweep scan

```
order   ← nodes of p sorted by p(v)/d(v) descending,
          ties broken lexicographically by node id (deterministic)
vol ← 0; cut ← 0; mass ← 0
for v in order:
    mass  ← mass + p(v)
    cut   ← cut + d(v) − 2·w(v, S)      # S = prefix before adding v
    vol   ← vol + d(v)
    if vol > ⌊vol(V)/2⌋: stop            # nested prefixes only grow
    record (|S|, vol, cut, score = cut/min(vol, vol(V)−vol), mass)
```

Properties:

* **Incremental evaluation.** The cut update `d(v) − 2·w(v, S)` makes
  the whole scan `O(Σ_v d(v))`; the previous implementation recomputed
  every cut from scratch at `O(K·Σ_v d(v))`. Measured on network 101:
  push 667 ms vs. **scan + selection 2.4 ms** (197 recorded prefixes,
  193 admissible) — the selection stage is ≤ 0.4 % of the runtime.
* **Single pass, five answers.** The scan stores, per prefix, the four
  numbers (`vol, cut, score, mass`) that *all* five steps read. Steps
  1–5 are linear filters over this array; none of them re-runs the
  push.
* **Early stop.** Scanning halts at the volume ceiling; the remaining
  prefixes are out of the window for every step.

## 6. The five steps

The steps are tried in order; the first that yields a candidate wins.

### Step 1 — strict (rule `strict`)

`volume ∧ mass ∧ quality` all hold. The returned cut is **the
lowest-score admissible candidate, with ties resolved towards the
smallest prefix**: the first candidate whose score is within
`epsilon_tie` (default **10 %**) of the best score of the set.

*Why not "the first prefix below φ" (previous behaviour)?* It stops on
the descending limb of the conductance curve and leaves quality on the
table, and it is knife-edge sensitive: a candidate moving from 0.399 to
0.401 between two queue orders changes the answer, which directly
hurts repeatability.

*Why not a raw `argmin`?* With a flat conductance curve a raw minimum
drifts arbitrarily far to the right — the classic "argmin selects a
huge set" failure. The relative tie rule returns the smallest of
equivalent prefixes and only moves outward when the score genuinely
improves by more than `epsilon_tie`.

Outcome metadata: `strategy = "target_conductance"`,
`conditions_met_count = 3`, `satisfied = true`.

### Step 2 — quality without mass (rule `quality_without_mass`)

`volume ∧ quality`, mass dropped (Sec. 4). A prefix with the right size
and a conductance below *φ* is a valid cluster even when the walk
truncated its mass.

### Step 3 — multi-scale (rule `multi_scale`)

If no prefix *inside* the window satisfies the quality bar, the volume
floor is halved up to `scale_relaxations = 2` times — floors
`2^(b−2)` then `2^(b−3)` — and the search is repeated at each floor
(strictest first), with the mass certificate still dropped.

**No second push is executed.** The scale `b` enters the algorithm
only through *ε* and the floor; a finer *ε* (that of the original `b`)
is always admissible for a coarser scale, so the *same* PageRank
vector and the *same* scan serve every floor. The floor actually used
is reported as `volume_min_used`, the number of halvings as
`scale_relaxation`.

### Step 4 — relaxed quality bar (rule `relaxed_phi`)

If even the relaxed floors admit nothing, the quality bar is raised
geometrically — `φ·{1.25, 1.5, 1.75, 2.0}`, capped at `2φ` — inside
the **original** volume window, returning at the first rung that
admits a prefix (again the smallest within tie-distance of the best).

The teleport probability is **frozen** during relaxation: only the
acceptance threshold moves. (With the retired legacy parameter set `α`
was computed from `φ`; and in general recomputing `α` per rung would
force a new push per rung and change the vector itself — the rung
would no longer be comparable.)

The rung is fully reported: `phi_effective` (the bar actually met),
`phi_relax_factor` (which rung). `strategy` becomes
`"target_conductance_fallback"` and `failed_conditions` still names
`phi_ok` — the strict bar was *not* met.

### Step 5 — best effort (rule `penalty_best_effort`)

If no step above fires, every admissible prefix is graded by a
dimensionless penalty and the minimum is returned:

```
penalty(S) = score(S)/φ                                  # quality, dominant
           + 0.5 · max(0, log₂(volume_min / vol(S)))      # octaves below the floor
           − 0.1 · min(p(S), 1)                           # mass, mild bonus
```

* `score/φ` is *cardinal*: 1.0 sits exactly on the bar, 1.4 is a 40 %
  overshoot. The previous rank counted satisfied criteria, which
  treated a 2.5 % miss and a 3× miss identically.
* The volume term is measured in **octaves** — the same unit `b` is
  defined in — so one halving outside the floor is priced at 0.5
  quality units; inside the window the term is exactly zero.
* The mass term can never exceed 0.1: the weakest certificate is only
  a tie-breaker.
* The `min_cut_size` floor still applies, so a dense three-taxon pair
  with near-zero conductance cannot win (it is excluded *before*
  scoring, not merely penalised).

The full breakdown is stored in `penalty_breakdown` (`quality`,
`volume_octaves`, `mass_bonus`, `penalty`), so the reason a cluster
was accepted is inspectable, not implicit.

*Worked example (network 101, φ = 0.05, α = 0.05, b = 10).* Nothing
satisfies φ ≈ 0.05 — the best achievable conductance from this seed is
≈ 0.26. Step 5 returns a 61-taxon prefix with conductance 0.2645,
inside the volume window: `penalty = 5.290 + 0 − 0.048 = 5.242`, i.e.
the candidate is dominated by the quality overshoot and is not penalised
for size at all.

### Step 0 — no admissible cluster (rule `no_admissible_cluster`)

If the admissible set is empty (e.g. the graph has fewer than
`min_cut_size` taxa, or every prefix's penalty is non-finite), the seed
is returned alone with `strategy = "seed_singleton_fallback"` and all
three criteria listed in `failed_conditions`. This is the only honest
"nothing found" answer; the previous implementation could instead
surface a sub-floor micro-cluster through the fallback rank.

## 7. Parameters

### 7.1 Algorithm parameters (parameter set `app`, the only set left)

| Parameter | Symbol | Default | Source / justification |
|---|---|---|---|
| quality bar | *φ* | `0.4` (`targetConductance`) | value of the reference experiments (strict acceptance 100 % on the bundled networks); the legacy `0.2` is reproducible via overrides |
| scale | *b* | `compute_default_b(n, m)` → **10** on both bundled networks | dimension-based: target cluster ~10 % of the taxa converted to a volume; sets `volume_min = 2^(b−1)` and *ε* |
| teleport | *α* | **0.05**, constant | selected by the reference experiments; the paper's bound `φ²/(225·ln(100√m))` gives ~2·10⁻⁵ and burns the iteration budget (Sec. 9) |
| accuracy | *ε* | `1/(2^b·48·⌈log₂m⌉)` → 1.57·10⁻⁶ at b=10 | never set by hand |
| weights | — | `weight = |correlation|` | app semantics; `weighted=false` uses edge counts |

The former `recommended` profile (φ = 0.4, b = 10, α = 0.05,
dimension-based *b*) and the `first` / `min_conductance` selection
rules were removed as *variants*; their values became the shipped
defaults of the single parameter set (the legacy app values —
φ = 0.2, b = 4, α ≈ 2·10⁻⁵ — are reproducible through explicit
overrides). The experiment harness still passes its base values
explicitly (`EXPERIMENT_DEFAULTS` in `tests/stability_experiment.py`),
so the published experiments remain reproducible even if the product
defaults ever move again.

### 7.2 Ladder parameters (`selection_ladder.py`)

| Parameter | Default | Meaning |
|---|---|---|
| `min_cut_size` | **5** | universal size floor (admissibility of every step) |
| `max_cut_size` | `None` | optional upper size bound (payload `parameterOverrides`) |
| `scale_relaxations` | **2** | step 3: floors `2^(b−2)`, `2^(b−3)` |
| `phi_relax_factors` | **(1.25, 1.5, 1.75, 2.0)** | step 4 rungs, capped at `2φ` |
| `epsilon_tie` | **0.10** | step 1/2/3/4: smallest prefix within 10 % of the best score |
| `penalty_volume_weight` | **0.5** | step 5: quality units per octave below the floor |
| `penalty_mass_weight` | **0.10** | step 5: maximum mass bonus |
| `max_iterations` | 100 000 | push budget; reported via `push_converged` |

## 8. Metadata and strategy taxonomy

Every run emits the legacy keys the reports and the UI consume
(`strategy`, `conditions_met_count`, `failed_conditions`, `constraints`,
`volume`, `conductance`, `p_at_vol`, `ranked_nodes`, ...) plus the
ladder keys:

| Step | `ladder_step` | `ladder_rule` | `strategy` | distinctive keys |
|---|---|---|---|---|
| strict | 1 | `strict` | `target_conductance` | `phi_effective = φ` |
| drop mass | 2 | `quality_without_mass` | `target_conductance_fallback` | `failed_conditions ∋ prob_ok` |
| multi-scale | 3 | `multi_scale` | `target_conductance_fallback` | `scale_relaxation`, `volume_min_used < 2^(b−1)` |
| relaxed φ | 4 | `relaxed_phi` | `target_conductance_fallback` | `phi_effective > φ`, `phi_relax_factor` |
| best effort | 5 | `penalty_best_effort` | `target_conductance_fallback` | `penalty_breakdown` |
| none | 0 | `no_admissible_cluster` | `seed_singleton_fallback` | `size = 1`, all criteria failed |

`conditions_met_count` and `failed_conditions` are always computed
**against the strict criteria** (original window, `τ`, `φ`), regardless
of which step produced the answer — so `satisfied = (count == 3)` keeps
its meaning and the stability reports' `strict_success_rate` counts
exactly step-1 outcomes. Signed runs additionally carry `score`,
`sign_lambda`, `cut_pos`, `neg_inside`, `neg_cut`, `positive_inside`,
`violations`, `vol_pos`.

## 9. Complexity and runtime budget

| Component | Cost | Measured (network 101, weighted, seed *Chloroplast*) |
|---|---|---|
| push | `O(iterations · d̄)`, capped at 10⁵ | α=0.1: 12 620 it / 287 ms · α=0.05: 26 708 it / 657 ms · α=0.01: 10⁵ it (capped) / 2 308 ms · α=2·10⁻⁵: 10⁵ it (capped) / 2 550 ms |
| sweep scan | `O(Σ d(v))` | ≈ 2 ms (197 prefixes) |
| ladder (5 steps) | `O(K)` per step | included in the 2.4 ms below |
| **scan + selection total** | | **2.4 ms** |

The selection stage cannot blow the budget: it is a handful of linear
filters over one array. The wall-clock ceiling is set by the push alone
and is *hard*, through `max_iterations`. Two regimes follow:

* **Converged vector** (`α ≳ 0.03` on these graphs): total run ≈
  0.3–0.7 s.
* **Iteration-capped vector** (`α ≲ 0.01`): ≈ 2.5 s, with
  `push_converged = false` recorded in the metadata. The sweep and the
  ladder are unaffected (they only need relative densities), but the
  mass criterion is likely to fail — which step 2 is designed for.

Reference behaviours on the bundled networks (weighted, seed
*Chloroplast*):

| Parameter set | Network | Outcome |
|---|---|---|
| reference (φ=0.4, b=10, α=0.05) | 101 | step 1, 50 taxa, φ = 0.2833, vol = 830 ∈ (512, 3320] |
| reference | rs4-101 | step 1, 127 taxa, φ = 0.1764, vol = 3077 ∈ (512, 3291] |
| legacy (φ=0.2, b=4, α=2·10⁻⁵) | 101 | step 4, 52 taxa, φ = 0.2752, `phi_effective = 0.30`, `push_converged = false` |
| legacy | rs4-101 | step 1, 128 taxa, φ = 0.165 |
| φ = 0.4, b = 13 (floor 4096 > any real prefix) | 101 | step 3, 171 taxa, `volume_min_used = 2048`, φ = 0.2973 |
| φ = 0.05 | 101 | step 5, 61 taxa, φ = 0.2645, penalty 5.242 |
| signed, φ=0.4, b=10, α=0.05 | 101 | step 1, 133 taxa, signed score 0.1972 (plain conductance 0.4890) |
| any | 3-node synthetic graph | step 0, singleton, `seed_singleton_fallback` |

Determinism: with `random_seed = None` repeated runs are bit-identical;
five runs with randomized queue orders (`random_seed = 42…46`,
reference parameters, 101) produced five identical 50-taxon clusters
(pairwise Jaccard 1.000).

## 10. Changes relative to the previous implementation

1. **Removed**: the `recommended` parameter profile
   (`clustering_profiles.recommended_parameters`, payload key
   `parameterProfile: "recommended"`), the `selection = "first" /
   "min_conductance"` rules, and the payload keys `parameterProfile`
   and `selection`. Selection is exclusively the ladder; a payload that
   still asks for `recommended` is rejected with an explicit error.
2. **Corrected**: the sweep order now divides by the *same* degree
   measure used for the volumes and the conductance (the previous order
   always used unweighted degrees, even in weighted mode).
3. **Tightened**: volume ceiling `2/3 → 1/2` of `vol(V)`
   (complement-invariance argument, Sec. 4).
4. **Enforced**: `min_cut_size` (previously accepted and discarded via
   `del`) is now a floor of *every* step, and seed membership is
   required of every candidate.
5. **Replaced**: the vote-count fallback rank
   `(conditions_met, −conductance, mass, −size, −index)` by the graded
   penalty of step 5, and the "first prefix below φ" rule by the
   min-score rule with the `epsilon_tie` tie-break.
6. **Accelerated**: the `O(K·E)` boundary recomputation per run was
   replaced by the incremental scan; the double evaluation (the inner
   `pagerank_nibble` computed a best cut that the caller discarded) is
   gone.
7. **Reported**: `ladder_step`, `ladder_rule`, `phi_effective`,
   `volume_min_used`, `penalty_breakdown`, `push_iterations`,
   `push_converged`.
8. **Defaults**: the shipped app defaults are now the reference
   experiment values (φ = 0.4, dimension-based `b` → 10, constant
   α = 0.05); the legacy values (φ = 0.2, b = 4, α ≈ 2·10⁻⁵ from the
   paper's bound) are retired and reproducible through explicit
   overrides. The UI treats an empty parameter field as "unset" and the
   backend clamps a non-positive `φ` back to the default, so the ladder
   always receives a usable bar.

The raw JSON files in `results/` (dated 2026-10-26) predate this
refactor and still record the removed profile as run provenance; they
are historical data and were not edited. The Markdown reports can be
re-rendered from them with `rebuild_reports.py` — the numbers do not
change, only the static prose.

## 11. Reproducibility

```bash
# unit + smoke tests (includes the ladder tests)
.venv/bin/python -m unittest discover -s python/clustering/tests -p "test_*.py" -v

# stability experiments (base values: EXPERIMENT_DEFAULTS)
.venv/bin/python python/clustering/tests/tuning_experiment.py --jobs 8
.venv/bin/python python/clustering/tests/comparison_experiment.py --reuse-unsigned --jobs 8

# re-render the reports from existing JSON without re-running anything
.venv/bin/python python/clustering/tests/rebuild_reports.py --findings
```

`LadderSelectionTests` in `tests/test_stability.py` pins one test per
step (strict, multi-scale, relaxed-φ, penalty, singleton, size bound,
metadata schema).

## 12. Reference

Andersen, R., Chung, F., & Lang, K. (2006). *Local Graph Partitioning
using PageRank Vectors.* Proceedings of the 47th Annual IEEE Symposium
on Foundations of Computer Science (FOCS).
<https://www.cs.cmu.edu/~15859n/RelatedWork/local_partitioning_full.pdf>
