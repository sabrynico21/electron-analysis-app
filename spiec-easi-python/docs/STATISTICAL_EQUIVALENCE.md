# Statistical equivalence with R (Phase E)

R's StARS subsample indices are **not** injected here: each implementation
draws its own subsamples, so the runs cannot be identical. What is compared
is the *distribution* over independent seeds of:

- the selected lambda index,
- the StARS instability curve (raw and cummax), averaged over seeds,
- the edge count and the per-edge selection frequency.

Reference runs: `tests/r_reference/generate_stats_runs.R` ->
`tests/golden_stats/<case>/r_runs.json`. Python side:
`tests/test_statistical_equivalence.py` (marked `slow`).

## Summary

| case | shape | seeds | R selected index | mean edges R / Python | raw curve max|diff| | cummax curve max|diff| | distinct edges R / Python |
|---|---|---|---|---|---|---|---|
| stats_p10_n60_s401 | 10 x 60 | 30 | {'1': 30} | 0.00 / 0.00 | 0.0106 | 0.0106 | 0 / 0 |
| stats_p10_n100_s402 | 10 x 100 | 30 | {'1': 30} | 0.00 / 0.00 | 0.0096 | 0.0021 | 0 / 0 |

## Detail
### stats_p10_n60_s401

- seeds: 30 (pulsar seed 1..30 on the R side,
  `random_state` 1..30 on the Python side)
- R selected-index counts: `{'1': 30}`
- mean edge count: R 0.000 (sd 0.000), Python 0.000 (sd 0.000)
- distinct edges seen: R 0, Python 0
- per-edge selection-frequency Pearson r: n/a (degenerate)
- mean raw instability curve, max |R - Python|: 0.0106 (mean 0.0014)
- mean cummax instability curve, max |R - Python|: 0.0106 (mean 0.0011)

| lambda index | R raw | Python raw | R cummax | Python cummax |
|---|---|---|---|---|
| 0 | 0.1309 | 0.1273 | 0.1309 | 0.1273 |
| 1 | 0.3337 | 0.3231 | 0.3337 | 0.3231 |
| 2 | 0.4699 | 0.4695 | 0.4699 | 0.4695 |
| 3 | 0.3956 | 0.3862 | 0.4699 | 0.4695 |
| 4 | 0.1143 | 0.1132 | 0.4699 | 0.4695 |
| 5 | 0.0244 | 0.0236 | 0.4699 | 0.4695 |
| 6 | 0.0030 | 0.0023 | 0.4699 | 0.4695 |
| 7 | 0.0001 | 0.0004 | 0.4699 | 0.4695 |
| 8 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 9 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 10 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 11 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 12 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 13 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 14 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 15 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 16 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 17 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 18 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
| 19 | 0.0000 | 0.0000 | 0.4699 | 0.4695 |
### stats_p10_n100_s402

- seeds: 30 (pulsar seed 1..30 on the R side,
  `random_state` 1..30 on the Python side)
- R selected-index counts: `{'1': 30}`
- mean edge count: R 0.000 (sd 0.000), Python 0.000 (sd 0.000)
- distinct edges seen: R 0, Python 0
- per-edge selection-frequency Pearson r: n/a (degenerate)
- mean raw instability curve, max |R - Python|: 0.0096 (mean 0.0010)
- mean cummax instability curve, max |R - Python|: 0.0021 (mean 0.0008)

| lambda index | R raw | Python raw | R cummax | Python cummax |
|---|---|---|---|---|
| 0 | 0.1025 | 0.1023 | 0.1025 | 0.1023 |
| 1 | 0.3356 | 0.3335 | 0.3356 | 0.3335 |
| 2 | 0.3336 | 0.3271 | 0.3429 | 0.3422 |
| 3 | 0.2514 | 0.2418 | 0.3429 | 0.3422 |
| 4 | 0.0572 | 0.0585 | 0.3429 | 0.3422 |
| 5 | 0.0021 | 0.0019 | 0.3429 | 0.3422 |
| 6 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 7 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 8 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 9 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 10 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 11 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 12 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 13 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 14 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 15 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 16 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 17 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 18 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |
| 19 | 0.0000 | 0.0000 | 0.3429 | 0.3422 |

## Interpretation and limits

- For this synthetic data family StARS selects `lambda = lambda.max` for every
  seed in both implementations, so the networks are empty (or contain the single
  `argmax|cor|` edge) and the selected-index statistic is degenerate. This is a
  property of the data, not of the implementation.
- The discriminating statistics are therefore the instability curves, which agree
  up to a few hundredths in raw/cummax instability at every lambda.
- This is **not** a claim of general statistical equivalence: it holds for these
  datasets, these 30 seeds and this data family only.

## Expanded stochastic-baseline scenarios (`tests/test_r_variability_expanded_scenarios.py`)

The six original scenes all select `lambda = lambda.max` (`py index 0`) at every
seed in both stacks, so the cross-implementation distance was dominated by the
block that is numerically indeterminate even in R (`docs/AUDIT.md` §2.E): Python
returns the empty model, R keeps the single `argmax|cor|` edge. The test now:

* applies the **same narrow, explicit exception** used by step 9 of the extended
  equivalence test — only when both implementations selected index 0 and only
  when the entire difference is that single `argmax|cor|` pair. The declared
  margins (`STOCHASTIC_BASELINE_TOLERANCE`) are unchanged.
* adds three higher-signal scenes (two chain-like, one hub-like) that select an
  **interior** lambda at every seed; on those the comparison is strict and an
  assertion proves the exception is not needed.
* reports the diagnostics and fails if fewer than
  `MIN_INTERIOR_LAMBDA_SELECTIONS` (scene, seed) curves are informative, so the
  test cannot silently degrade into comparing empty models.

Measured run (Python with the compiled MB core, R 4.6.1, 9 scenes × 4 seeds):

| | value |
|---|---|
| (scene, seed) curves | 36 (3 of them strict scenes) |
| curves selecting `lambda.max` in Python | 24/36 |
| curves selecting an interior lambda in Python | 12/36 (minimum required 3) |
| `lambda.max` exceptions used | 23 (each verified to differ only by the `argmax\|cor\|` pair) |
| cross-implementation medians | Jaccard 0.000000, degree MAE 0.000000, lambda-rank gap 0.000000 (margins 0.08 / 0.20 / 0.10) |
| on the three strict scenes | Python and R select the *same* interior index at all four seeds and the edge sets agree exactly (35/35, 41/41, 15/15 edges) |

The remaining non-excepted curve is a `chain_medium` seed where R selects an
index other than 0 while Python selects 0; it does not affect any median.
Again: this is a statement about these scenes, seeds and implementations only.

