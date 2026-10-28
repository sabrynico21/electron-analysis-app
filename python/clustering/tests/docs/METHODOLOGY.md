# Methodology: the stability experiments

This document is the single definition of the two experiments, of every
column shown in the reports, and of how to interpret the values. The
generated reports ([RESULTS_TUNING.md](../results/RESULTS_TUNING.md) and
[RESULTS_COMPARISON.md](../results/RESULTS_COMPARISON.md)) link here instead
of repeating the legend. The algorithm itself (parameters, signed
edges, payload format) is documented in [ALGORITHM.md](ALGORITHM.md);
the cluster-selection stage is specified in
[SELECTION_LADDER.md](SELECTION_LADDER.md).

## 1. The two experiments

Every configuration is evaluated with the same two experiments, on both the
weighted graph (`|correlation|` edges) and the unweighted graph (every edge
has weight 1).

### Experiment 1 - repeatability

Run the clustering **N times from the same seed taxon** (`Chloroplast`),
changing only the residual-queue order. The algorithm is deterministic
without that change, so ordering is the only source of variation. `N` is the
`--runs` argument (100 in the published reports).

*Answers:* does one seed with one parameter set always return the same
cluster?

### Experiment 2 - nested consistency

1. Compute the deterministic cluster of the seed - the **base cluster**.
2. Re-run the algorithm **once from every member of that cluster**.
3. Compare the resulting clusters with each other and with the base cluster.

The run count of Experiment 1 does **not** apply here: the number of runs
equals the size of the base cluster, because every member must get the
chance to recover the community. (`--max-nested-nodes` caps this only for
smoke tests and is never used for the published reports.)

*Answers:* is the community recovered from its own members, whichever member
is used as the seed? A real community is; a deterministic artefact of a
lucky seed is not.

## 2. Columns of the report tables

| Column | Meaning |
|---|---|
| Runs | Runs used in that row, **one value per tested parameter value (or per configuration), separated by `/`**: fixed `N` in Experiment 1 (identical across values by construction); in Experiment 2 the size of each value's base cluster, so the counts differ whenever the base clusters differ (e.g. a larger `b` yields a larger base cluster and therefore more re-seeded runs). |
| Mean cluster size (taxa) | Mean number of taxa in the returned cluster; there is one cluster per run. |
| Jaccard (mean / min / max) | Overlap `\|A ∩ B\| / \|A ∪ B\|` between two runs; 1 means identical clusters. |
| Intersection (n / %) | Taxa present in **every** run, and that number as a share of the mean cluster size (Experiment 1) or of the base cluster (Experiment 2). |
| Union | Taxa present in at least one run. |
| Conductance (mean / min / max) | Cut quality of the cluster, lower is better. Plain `\|correlation\|` conductance is shown for every configuration; the signed configurations are accepted on the signed score, which is stored in the JSON. |
| Identical clusters | Number of runs sharing the most frequent cluster; equals the run count when every run produced the same cluster. |
| Consensus >=80% | Taxa found by at least 80% of the runs. |
| Asymmetry (Experiment 2 only) | Share of ordered pairs `(u, v)` where `u`'s cluster contains `v` but `v`'s cluster does not contain `u`. |

Additional fields (support, signed scores, per-run rows, parameters,
thresholds, condition pass rates) are stored in the raw JSON under
`results/`; the "Findings" sections of the reports quote the relevant ones.

The ranking columns of the tuning report (`Best parameters`, `Best score`,
`Default score`, `Gap`) are not part of these tables: they are defined next
to that table in
[RESULTS_TUNING.md](../results/RESULTS_TUNING.md#best-value-per-condition).

## 3. Support above chance

Raw support favours large clusters: a cluster covering most of the graph
trivially recovers most taxa from any seed, so raw support cannot be
compared across conditions of different sizes. The reports therefore quote
**support above chance** in the Findings sections:

* `support_null = k / n` is the per-taxon support expected from random
  re-runs of a cluster of size `k` in a graph of `n` taxa;
* `support above chance = (support - support_null) / (1 - support_null)`.

It is 0 for "no better than random re-seeds" and 1 for "every member
recovers the cluster".

## 4. Reproducibility

* Experiment 1 uses the run seeds `42, 43, ...` (recorded as
  `ordering_seeds` in the JSON); the seed taxon and all parameters are fixed.
* Experiment 2 is fully deterministic (no random seeds): the base cluster is
  recomputed, then every member is used as a seed.
* The optional smoke-test subsample of Experiment 2 uses the fixed seed
  `20241023`.
* Every report is generated from the raw JSON, which records the graph, the
  parameters, the derived thresholds and every per-run result, so any table
  in the reports can be recomputed from `results/*.json` alone.
