# Comparison: defaults, tuned parameters, signed and both

Generated: `2026-10-27T13:52:08.458137+00:00`  
Command: `python comparison_experiment.py --runs 100 --jobs 8 --reuse-unsigned --output-dir python/clustering/tests/results --summary-md python/clustering/tests/results/RESULTS_COMPARISON.md`

## Design

- **Experiment 1 - repeatability**: 100 runs from the same seed taxon, changing only the residual-queue order (the only source of variation: the algorithm is otherwise deterministic).
- **Experiment 2 - nested consistency**: the deterministic cluster of the seed (the *base cluster*) is re-computed, and the algorithm is then re-run from **every member** of that cluster.
- Every condition is run on the weighted graph (`|correlation|` edges) and on the unweighted graph (every edge has weight 1).

Column definitions, the metrics and the interpretation rules live in [docs/METHODOLOGY.md](../docs/METHODOLOGY.md); this report only lists the values.

The tuned parameters come from `tuning_20261027T135033Z.json`.

## Configurations

Every table row follows the order **Default / Best / Signed / Both**. Each cell repeats that order, separated by `/`.

| Condition | Default | Best | Signed | Both |
|---|---|---|---|---|
| rs4-101 / weighted | alpha=0.05, b=10, phi=0.4 | alpha=0.01, b=10, phi=0.4 | alpha=0.05, b=10, phi=0.4 | alpha=0.01, b=10, phi=0.4 |
| rs4-101 / unweighted | alpha=0.05, b=10, phi=0.4 | alpha=0.05, b=12, phi=0.4 | alpha=0.05, b=10, phi=0.4 | alpha=0.05, b=12, phi=0.4 |
| 101 / weighted | alpha=0.05, b=10, phi=0.4 | alpha=0.05, b=8, phi=0.4 | alpha=0.05, b=10, phi=0.4 | alpha=0.05, b=8, phi=0.4 |
| 101 / unweighted | alpha=0.05, b=10, phi=0.4 | alpha=0.05, b=8, phi=0.4 | alpha=0.05, b=10, phi=0.4 | alpha=0.05, b=8, phi=0.4 |

`Best` uses the best `alpha` and `b` from the tuning report; `Both` is the signed version with those same parameters.

## Reading note

Every value column is ordered **Default / Best / Signed / Both**: each item inside a cell is separated by `/`, and inside a mean/min/max cell each `(mean/min/max)` group belongs to one configuration. The conductance column is the plain `|correlation|` conductance for every configuration; Signed and Both are accepted on the signed score, which is stored in the JSON.

## Experiment 1 - Repeatability (100 runs per configuration)

| Condition | Runs | Mean cluster size (taxa) | Jaccard (mean / min / max) | Intersection (n / %) | Union | Conductance (mean / min / max) | Identical clusters | Consensus >=80% |
|---|---|---|---|---|---|---|---|---|
| rs4-101 / weighted | 100 / 100 / 100 / 100 | 126.9 / 127.4 / 80.0 / 114.0 | (0.998/0.984/1.000) / (0.995/0.984/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 125 (99%) / 125 (98%) / 80 (100%) / 114 (100%) | 127 / 128 / 80 / 114 | (0.177/0.176/0.183) / (0.170/0.169/0.171) / (0.534/0.534/0.534) / (0.523/0.523/0.523) | 86 / 53 / 100 / 100 | 127 / 127 / 80 / 114 |
| rs4-101 / unweighted | 100 / 100 / 100 / 100 | 126.0 / 126.0 / 75.0 / 173.0 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 126 (100%) / 126 (100%) / 75 (100%) / 173 (100%) | 126 / 126 / 75 / 173 | (0.183/0.183/0.183) / (0.183/0.183/0.183) / (0.592/0.592/0.592) / (0.511/0.511/0.511) | 100 / 100 / 100 / 100 | 126 / 126 / 75 / 173 |
| 101 / weighted | 100 / 100 / 100 / 100 | 50.0 / 50.0 / 133.0 / 133.1 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) / (0.998/0.971/1.000) | 50 (100%) / 50 (100%) / 133 (100%) / 133 (100%) | 50 / 50 / 133 / 137 | (0.283/0.283/0.283) / (0.283/0.283/0.283) / (0.489/0.489/0.489) / (0.489/0.481/0.489) | 100 / 100 / 100 / 97 | 50 / 50 / 133 / 133 |
| 101 / unweighted | 100 / 100 / 100 / 100 | 51.0 / 51.0 / 127.0 / 127.0 | (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) / (1.000/1.000/1.000) | 51 (100%) / 51 (100%) / 127 (100%) / 127 (100%) | 51 / 51 / 127 / 127 | (0.290/0.290/0.290) / (0.290/0.290/0.290) / (0.471/0.471/0.471) / (0.471/0.471/0.471) | 100 / 100 / 100 / 100 | 51 / 51 / 127 / 127 |

## Experiment 2 - Nested consistency (one run per base-cluster member)

| Condition | Runs | Mean cluster size (taxa) | Jaccard (mean / min / max) | Intersection (n / %) | Union | Conductance (mean / min / max) | Identical clusters | Consensus >=80% | Asymmetry |
|---|---|---|---|---|---|---|---|---|---|
| rs4-101 / weighted | 127 / 127 / 80 / 114 | 130.7 / 131.8 / 94.8 / 101.2 | (0.387/0.000/0.963) / (0.413/0.007/0.975) / (0.578/0.025/1.000) / (0.581/0.007/1.000) | 0 (0%) / 0 (0%) / 1 (1%) / 1 (1%) | 368 / 368 / 291 / 307 | (0.232/0.143/0.361) / (0.224/0.142/0.356) / (0.515/0.457/0.582) / (0.519/0.462/0.661) | 1 / 1 / 2 / 3 | 26 / 32 / 60 / 57 | 15% / 15% / 19% / 28% |
| rs4-101 / unweighted | 126 / 126 / 75 / 173 | 133.8 / 135.8 / 90.1 / 163.7 | (0.431/0.000/0.957) / (0.436/0.009/0.957) / (0.631/0.195/1.000) / (0.480/0.000/1.000) | 0 (0%) / 0 (0%) / 22 (29%) / 0 (0%) | 368 / 368 / 235 / 367 | (0.233/0.151/0.371) / (0.233/0.151/0.371) / (0.570/0.486/0.621) / (0.505/0.411/0.564) | 1 / 1 / 2 / 3 | 36 / 37 / 58 / 45 | 12% / 12% / 20% / 17% |
| 101 / weighted | 50 / 50 / 133 / 133 | 99.1 / 99.0 / 112.4 / 112.4 | (0.488/0.170/0.943) / (0.489/0.170/0.943) / (0.326/0.000/0.966) / (0.326/0.000/0.966) | 26 (52%) / 26 (52%) / 0 (0%) / 0 (0%) | 383 / 382 / 380 / 380 | (0.270/0.231/0.295) / (0.270/0.231/0.295) / (0.508/0.429/0.631) / (0.508/0.429/0.631) | 1 / 1 / 1 / 1 | 39 / 39 / 2 / 2 | 20% / 20% / 12% / 12% |
| 101 / unweighted | 51 / 51 / 127 / 127 | 90.5 / 90.3 / 109.7 / 109.7 | (0.502/0.162/0.933) / (0.503/0.162/0.933) / (0.334/0.000/0.966) / (0.334/0.000/0.966) | 26 (51%) / 26 (51%) / 0 (0%) / 0 (0%) | 385 / 385 / 377 / 377 | (0.275/0.222/0.298) / (0.275/0.223/0.299) / (0.491/0.413/0.609) / (0.491/0.413/0.609) | 1 / 1 / 1 / 1 | 37 / 37 / 2 / 2 | 25% / 25% / 14% / 14% |

## Related files

- `results/comparison_*.json` - raw data of this comparison.
- [RESULTS_TUNING.md](RESULTS_TUNING.md) - the sweep that selected the tuned parameters.
- [README.md](../README.md) - repository layout and the exact commands used to produce this report.
