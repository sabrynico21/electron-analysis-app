# Parameter tuning: default values and one-at-a-time sweep

Generated: `2026-10-27T13:29:14.388349+00:00`  
Command: `python tuning_experiment.py --network both --runs 100 --jobs 8 --output-dir python/clustering/tests/results --summary-md python/clustering/tests/results/RESULTS_TUNING.md`

## Design

- **Experiment 1 - repeatability**: 100 runs from the same seed taxon, changing only the residual-queue order (the only source of variation: the algorithm is otherwise deterministic).
- **Experiment 2 - nested consistency**: the deterministic cluster of the seed (the *base cluster*) is re-computed, and the algorithm is then re-run from **every member** of that cluster.
- Every condition is run on the weighted graph (`|correlation|` edges) and on the unweighted graph (every edge has weight 1).

Column definitions, the metrics and the interpretation rules live in [docs/METHODOLOGY.md](../docs/METHODOLOGY.md); this report only lists the values.

## Defaults and sweep ranges

| Parameter | Default | Values tested |
|---|---|---|
| alpha | 0.05 | 0.01 / 0.05 / 0.1 |
| b | 10 | 8 / 10 / 12 |
| phi | 0.4 | 0.2 / 0.3 / 0.4 |

Each parameter row changes one parameter at a time; the other two stay at their default. `alpha` controls how local the random walk is, `b` sets the minimum cluster volume, and `phi` is the acceptance threshold.

## Reading note

Each parameter row lists its values in the order of the `Values` column; `/` separates the values, and inside a mean/min/max cell each `(mean/min/max)` group belongs to one value. The `default` row is the reference; the `alpha`, `b` and `phi` rows repeat the default value inside their sequence.

## Experiment 1 - Repeatability (100 runs per value)

| Condition | Parameter | Values | Runs | Mean cluster size (taxa) | Jaccard (mean / min / max) | Intersection (n / %) | Union | Conductance (mean / min / max) | Identical clusters | Consensus >=80% |
|---|---|---|---|---|---|---|---|---|---|---|
| rs4-101 / weighted | default | - | 100 | 126.9 | (0.998/0.984/1.000) | 125 (99%) | 127 | (0.177/0.176/0.183) | 86 | 127 |
| rs4-101 / weighted | alpha | 0.01 / 0.05 / 0.1 | 100 / 100 / 100 | 127.4 / 126.9 / 130.0 | (0.995/0.984/1.000) / (0.998/0.984/1.000) / (0.998/0.985/1.000) | 125 (98%) / 125 (99%) / 129 (99%) | 128 / 127 / 131 | (0.170/0.169/0.171) / (0.177/0.176/0.183) / (0.189/0.189/0.190) | 53 / 86 / 94 | 127 / 127 / 130 |
| rs4-101 / weighted | b | 8 / 10 / 12 | 100 / 100 / 100 | 126.3 / 126.9 / 127.0 | (0.992/0.976/1.000) / (0.998/0.984/1.000) / (1.000/1.000/1.000) | 124 (98%) / 125 (99%) / 127 (100%) | 127 / 127 / 127 | (0.179/0.175/0.183) / (0.177/0.176/0.183) / (0.176/0.176/0.176) | 45 / 86 / 100 | 125 / 127 / 127 |
| rs4-101 / weighted | phi | 0.2 / 0.3 / 0.4 | 100 / 100 / 100 | 126.9 / 126.9 / 126.9 | (0.998/0.984/1.000) / (0.998/0.984/1.000) / (0.998/0.984/1.000) | 125 (99%) / 125 (99%) / 125 (99%) | 127 / 127 / 127 | (0.177/0.176/0.183) / (0.177/0.176/0.183) / (0.177/0.176/0.183) | 86 / 86 / 86 | 127 / 127 / 127 |
| rs4-101 / unweighted | default | - | 100 | 126.0 | (1.000/1.000/1.000) | 126 (100%) | 126 | (0.183/0.183/0.183) | 100 | 126 |
| rs4-101 / unweighted | alpha | 0.01 / 0.05 / 0.1 | 100 / 100 / 100 | 140.8 / 126.0 / 125.0 | (0.993/0.972/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 137 (97%) / 126 (100%) / 125 (100%) | 143 / 126 / 125 | (0.164/0.163/0.166) / (0.183/0.183/0.183) / (0.193/0.193/0.193) | 64 / 100 / 100 | 140 / 126 / 125 |
| rs4-101 / unweighted | b | 8 / 10 / 12 | 100 / 100 / 100 | 126.0 / 126.0 / 126.0 | (1.000/0.992/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 126 (100%) / 126 (100%) / 126 (100%) | 127 / 126 / 126 | (0.183/0.181/0.183) / (0.183/0.183/0.183) / (0.183/0.183/0.183) | 97 / 100 / 100 | 126 / 126 / 126 |
| rs4-101 / unweighted | phi | 0.2 / 0.3 / 0.4 | 100 / 100 / 100 | 126.0 / 126.0 / 126.0 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 126 (100%) / 126 (100%) / 126 (100%) | 126 / 126 / 126 | (0.183/0.183/0.183) / (0.183/0.183/0.183) / (0.183/0.183/0.183) | 100 / 100 / 100 | 126 / 126 / 126 |
| 101 / weighted | default | - | 100 | 50.0 | (1.000/1.000/1.000) | 50 (100%) | 50 | (0.283/0.283/0.283) | 100 | 50 |
| 101 / weighted | alpha | 0.01 / 0.05 / 0.1 | 100 / 100 / 100 | 51.4 / 50.0 / 48.0 | (0.989/0.962/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 50 (97%) / 50 (100%) / 48 (100%) | 52 / 50 / 48 | (0.278/0.275/0.283) / (0.283/0.283/0.283) / (0.306/0.306/0.306) | 52 / 100 / 100 | 51 / 50 / 48 |
| 101 / weighted | b | 8 / 10 / 12 | 100 / 100 / 100 | 50.0 / 50.0 / 171.0 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 50 (100%) / 50 (100%) / 171 (100%) | 50 / 50 / 171 | (0.283/0.283/0.283) / (0.283/0.283/0.283) / (0.297/0.297/0.297) | 100 / 100 / 100 | 50 / 50 / 171 |
| 101 / weighted | phi | 0.2 / 0.3 / 0.4 | 100 / 100 / 100 | 50.0 / 50.0 / 50.0 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 50 (100%) / 50 (100%) / 50 (100%) | 50 / 50 / 50 | (0.283/0.283/0.283) / (0.283/0.283/0.283) / (0.283/0.283/0.283) | 100 / 100 / 100 | 50 / 50 / 50 |
| 101 / unweighted | default | - | 100 | 51.0 | (1.000/1.000/1.000) | 51 (100%) | 51 | (0.290/0.290/0.290) | 100 | 51 |
| 101 / unweighted | alpha | 0.01 / 0.05 / 0.1 | 100 / 100 / 100 | 53.0 / 51.0 / 53.0 | (0.982/0.963/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 52 (98%) / 51 (100%) / 53 (100%) | 55 / 51 / 53 | (0.277/0.274/0.278) / (0.290/0.290/0.290) / (0.297/0.297/0.297) | 67 / 100 / 100 | 52 / 51 / 53 |
| 101 / unweighted | b | 8 / 10 / 12 | 100 / 100 / 100 | 51.0 / 51.0 / 77.0 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 51 (100%) / 51 (100%) / 77 (100%) | 51 / 51 / 77 | (0.290/0.290/0.290) / (0.290/0.290/0.290) / (0.301/0.301/0.301) | 100 / 100 / 100 | 51 / 51 / 77 |
| 101 / unweighted | phi | 0.2 / 0.3 / 0.4 | 100 / 100 / 100 | 51.0 / 51.0 / 51.0 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 51 (100%) / 51 (100%) / 51 (100%) | 51 / 51 / 51 | (0.290/0.290/0.290) / (0.290/0.290/0.290) / (0.290/0.290/0.290) | 100 / 100 / 100 | 51 / 51 / 51 |

## Experiment 2 - Nested consistency (one run per base-cluster member)

| Condition | Parameter | Values | Runs | Mean cluster size (taxa) | Jaccard (mean / min / max) | Intersection (n / %) | Union | Conductance (mean / min / max) | Identical clusters | Consensus >=80% | Asymmetry |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rs4-101 / weighted | default | - | 127 | 130.7 | (0.387/0.000/0.963) | 0 (0%) | 368 | (0.232/0.143/0.361) | 1 | 26 | 15% |
| rs4-101 / weighted | alpha | 0.01 / 0.05 / 0.1 | 127 / 127 / 130 | 131.8 / 130.7 / 130.3 | (0.413/0.007/0.975) / (0.387/0.000/0.963) / (0.361/0.000/0.939) | 0 (0%) / 0 (0%) / 0 (0%) | 368 / 368 / 368 | (0.224/0.142/0.356) / (0.232/0.143/0.361) / (0.240/0.152/0.362) | 1 / 1 / 1 | 32 / 26 / 17 | 15% / 15% / 15% |
| rs4-101 / weighted | b | 8 / 10 / 12 | 125 / 127 / 127 | 130.2 / 130.7 / 134.0 | (0.391/0.000/0.954) / (0.387/0.000/0.963) / (0.396/0.007/0.954) | 0 (0%) / 0 (0%) / 0 (0%) | 368 / 368 / 368 | (0.232/0.143/0.361) / (0.232/0.143/0.361) / (0.232/0.143/0.361) | 1 / 1 / 1 | 26 / 26 / 28 | 15% / 15% / 15% |
| rs4-101 / weighted | phi | 0.2 / 0.3 / 0.4 | 127 / 127 / 127 | 134.7 / 131.9 / 130.7 | (0.397/0.000/0.963) / (0.391/0.000/0.963) / (0.387/0.000/0.963) | 0 (0%) / 0 (0%) / 0 (0%) | 368 / 368 / 368 | (0.228/0.143/0.350) / (0.231/0.143/0.361) / (0.232/0.143/0.361) | 1 / 1 / 1 | 30 / 29 / 26 | 16% / 15% / 15% |
| rs4-101 / unweighted | default | - | 126 | 133.8 | (0.431/0.000/0.957) | 0 (0%) | 368 | (0.233/0.151/0.371) | 1 | 36 | 12% |
| rs4-101 / unweighted | alpha | 0.01 / 0.05 / 0.1 | 141 / 126 / 125 | 137.8 / 133.8 / 130.9 | (0.444/0.000/0.979) / (0.431/0.000/0.957) / (0.405/0.000/0.947) | 0 (0%) / 0 (0%) / 0 (0%) | 368 / 368 / 368 | (0.222/0.147/0.370) / (0.233/0.151/0.371) / (0.243/0.151/0.374) | 1 / 1 / 1 | 42 / 36 / 30 | 11% / 12% / 12% |
| rs4-101 / unweighted | b | 8 / 10 / 12 | 126 / 126 / 126 | 133.5 / 133.8 / 135.8 | (0.431/0.000/0.957) / (0.431/0.000/0.957) / (0.436/0.009/0.957) | 0 (0%) / 0 (0%) / 0 (0%) | 368 / 368 / 368 | (0.233/0.151/0.371) / (0.233/0.151/0.371) / (0.233/0.151/0.371) | 1 / 1 / 1 | 36 / 36 / 37 | 12% / 12% / 12% |
| rs4-101 / unweighted | phi | 0.2 / 0.3 / 0.4 | 126 / 126 / 126 | 137.8 / 133.9 / 133.8 | (0.441/0.000/0.964) / (0.432/0.000/0.957) / (0.431/0.000/0.957) | 0 (0%) / 0 (0%) / 0 (0%) | 368 / 368 / 368 | (0.229/0.151/0.349) / (0.233/0.151/0.371) / (0.233/0.151/0.371) | 1 / 1 / 1 | 37 / 36 / 36 | 13% / 12% / 12% |
| 101 / weighted | default | - | 50 | 99.1 | (0.488/0.170/0.943) | 26 (52%) | 383 | (0.270/0.231/0.295) | 1 | 39 | 20% |
| 101 / weighted | alpha | 0.01 / 0.05 / 0.1 | 52 / 50 / 48 | 106.0 / 99.1 / 97.5 | (0.477/0.175/0.957) / (0.488/0.170/0.943) / (0.481/0.159/0.907) | 29 (56%) / 26 (52%) / 23 (48%) | 382 / 383 / 384 | (0.263/0.224/0.288) / (0.270/0.231/0.295) / (0.276/0.235/0.308) | 1 / 1 / 1 | 40 / 39 / 37 | 18% / 20% / 21% |
| 101 / weighted | b | 8 / 10 / 12 | 50 / 50 / 171 | 99.0 / 99.1 / 161.6 | (0.489/0.170/0.943) / (0.488/0.170/0.943) / (0.327/0.020/0.963) | 26 (52%) / 26 (52%) / 0 (0%) | 382 / 383 / 387 | (0.270/0.231/0.295) / (0.270/0.231/0.295) / (0.283/0.208/0.358) | 1 / 1 / 1 | 39 / 39 / 1 | 20% / 20% / 9% |
| 101 / weighted | phi | 0.2 / 0.3 / 0.4 | 50 / 50 / 50 | 115.1 / 99.1 / 99.1 | (0.458/0.170/0.945) / (0.488/0.170/0.943) / (0.488/0.170/0.943) | 26 (52%) / 26 (52%) / 26 (52%) | 385 / 383 / 383 | (0.264/0.231/0.295) / (0.270/0.231/0.295) / (0.270/0.231/0.295) | 1 / 1 / 1 | 39 / 39 / 39 | 17% / 20% / 20% |
| 101 / unweighted | default | - | 51 | 90.5 | (0.502/0.162/0.933) | 26 (51%) | 385 | (0.275/0.222/0.298) | 1 | 37 | 25% |
| 101 / unweighted | alpha | 0.01 / 0.05 / 0.1 | 53 / 51 / 53 | 92.1 / 90.5 / 92.1 | (0.503/0.161/0.956) / (0.502/0.162/0.933) / (0.479/0.139/0.918) | 26 (49%) / 26 (51%) / 18 (34%) | 383 / 385 / 387 | (0.269/0.223/0.289) / (0.275/0.222/0.298) / (0.282/0.229/0.342) | 1 / 1 / 1 | 39 / 37 / 37 | 25% / 25% / 28% |
| 101 / unweighted | b | 8 / 10 / 12 | 51 / 51 / 77 | 90.3 / 90.5 / 121.7 | (0.503/0.162/0.933) / (0.502/0.162/0.933) / (0.419/0.025/0.962) | 26 (51%) / 26 (51%) / 4 (5%) | 385 / 385 / 387 | (0.275/0.223/0.299) / (0.275/0.222/0.298) / (0.276/0.212/0.346) | 1 / 1 / 1 | 37 / 37 / 50 | 25% / 25% / 11% |
| 101 / unweighted | phi | 0.2 / 0.3 / 0.4 | 51 / 51 / 51 | 102.8 / 90.5 / 90.5 | (0.468/0.162/0.933) / (0.502/0.162/0.933) / (0.502/0.162/0.933) | 26 (51%) / 26 (51%) / 26 (51%) | 385 / 385 / 385 | (0.270/0.222/0.298) / (0.275/0.222/0.298) / (0.275/0.222/0.298) | 1 / 1 / 1 | 39 / 37 / 37 | 23% / 25% / 25% |

## Best value per condition

This table picks, for every network and weight mode, the `alpha` and `b` that `comparison_experiment.py` will use (the `Best` and `Both` configurations of the comparison report) and compares them with the default parameters. The comparison is a single number per experiment run, the **score**, defined below.

### How the score is calculated

Each swept value produces one Experiment 1 report and one Experiment 2 report. Each report is reduced to a size-corrected value in `[0, 1]` and the two are averaged:

**score = 0.5 x (repeatability above chance + support above chance)**

- **repeatability above chance** (Experiment 1): the mean pairwise Jaccard `J` of the runs, corrected for cluster size. With `k` the mean cluster size and `n` the taxa of the graph, two random clusters of that size overlap with `J_null = k / (2n - k)`, so the corrected value is `(J - J_null) / (1 - J_null)`. 0 = "no more overlap than two random clusters", 1 = "every run returned the same cluster".
- **support above chance** (Experiment 2): the fraction of base-cluster members that recover a taxon, corrected for cluster size. Random re-seeds recover it with `support_null = k / n`, so the corrected value is `(support - support_null) / (1 - support_null)`. 0 = "no better than random re-seeds", 1 = "every member recovered the cluster".

Both inputs are corrected for size on purpose: raw overlap and raw support both grow with the cluster, so the largest cluster would win any ranking. Averaging the two views balances "the same seed repeats its cluster" against "every member of that cluster finds it back": a value scores high only when it is stable under both.

### Eligibility

- Only values whose Experiment 1 **strict acceptance rate** (the share of runs that met the volume, probability and conductance conditions) is at least 90% compete; anything else is ineligible and cannot win the ranking.
- `phi` is not ranked: it changes acceptance, not the cluster shape. Only `alpha` and `b` vary in the ranking, so the winning `phi` is always its default value.

| Condition | Best parameters | Best score | Default score | Gap |
|---|---|---|---|---|
| rs4-101 / weighted | alpha=0.01, b=10, phi=0.4 | 0.727 | 0.709 | 0.018 |
| rs4-101 / unweighted | alpha=0.05, b=12, phi=0.4 | 0.755 | 0.754 | 0.001 |
| 101 / weighted | alpha=0.05, b=8, phi=0.4 | 0.932 | 0.932 | 0.000 |
| 101 / unweighted | alpha=0.05, b=8, phi=0.4 | 0.914 | 0.914 | 0.000 |

### Reading the table

- **Best parameters**: the winning `alpha` and `b` for this condition, with `phi` at its default.
- **Best score**: the score of that winner; these are the values consumed by the comparison report.
- **Default score**: the score of the defaults (`alpha=0.05`, `b=10`, `phi=0.4`) - the reference for the same condition.
- **Gap**: `Best score - Default score`. A positive gap means the sweep beat the defaults, a gap near 0 means the defaults sit on a plateau, and `-` means no value reached the 90% eligibility bar.

## Related files

- `results/tuning_*.json` - raw data of this sweep.
- [RESULTS_COMPARISON.md](RESULTS_COMPARISON.md) - the best parameters of this sweep vs the defaults, plus the signed-edge variants.
- [README.md](../README.md) - repository layout and the exact commands used to produce this report.
