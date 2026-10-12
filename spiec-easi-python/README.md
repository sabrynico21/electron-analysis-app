# spiec-easi-python

Python implementation of SPIEC-EASI (Kurtz et al., 2015), focused on the
MB variant (Meinshausen-Buhlmann neighborhood selection) with StARS model
selection.

## What this library does

Given a count matrix (samples x taxa), it:

1. Applies CLR preprocessing for compositional data.
2. Fits sparse neighborhood models (Lasso) over a lambda path.
3. Uses StARS to select the regularization level.
4. Builds an undirected sparse association network.

The output contains adjacency, weighted edges, selected lambda, instability
curves, and a ready-to-use NetworkX graph.

## Installation

```bash
pip install -e .
```

For development:

```bash
pip install -e .[dev]
```

## Quickstart

```python
import pandas as pd
from spieceasi import spiec_easi_mb

# rows: samples, columns: taxa
counts = pd.read_csv("counts.csv", index_col=0)

result = spiec_easi_mb(
		counts,
		n_lambdas=25,
		n_subsamples=50,
		stars_threshold=0.05,
)

print("selected lambda:", result.selected_lambda)
print("adjacency shape:", result.adjacency.shape)
print("edges:", result.graph.number_of_edges())
```

## Main function

```python
spiec_easi_mb(data, **kwargs)
```

Returns a `SpiecEasiMBResult` dataclass with:

- `adjacency`: binary symmetric matrix (selected edges)
- `weights`: symmetric edge weights
- `graph`: NetworkX graph
- `selected_lambda`: StARS-selected regularization value
- `alpha_path`, `instability`, `monotone_instability`
- `directed_coefficients`: per-neighborhood regression coefficients
- `taxa_names`

## Parameters (brief)

### Preprocessing

| Parameter | Default | Meaning |
| --- | --- | --- |
| `pseudocount` | `1.0` | Value used for zero handling in CLR preprocessing. |
| `zero_method` | `"additive"` | Zero strategy before CLR (`"additive"` or `"multiplicative"`). |

### MB regularization path

| Parameter | Default | Meaning |
| --- | --- | --- |
| `n_lambdas` | `25` | Number of regularization values to evaluate. |
| `lambda_min_ratio` | `1e-2` | Ratio between smallest and largest lambda in auto path. |
| `alpha_path` | `None` | Optional explicit lambda path. If set, overrides auto path creation. |

### StARS model selection

| Parameter | Default | Meaning |
| --- | --- | --- |
| `n_subsamples` | `50` | Number of subsamples used by StARS. |
| `subsample_size` | `None` | Absolute subsample size (mutually exclusive with `subsample_fraction`). |
| `subsample_fraction` | `None` | Relative subsample size in `(0,1)`. |
| `stars_threshold` | `0.05` | Instability threshold used to choose lambda. |

### MB fitting controls

| Parameter | Default | Meaning |
| --- | --- | --- |
| `sym_rule` | `"or"` | Symmetrization rule for neighborhood edges (`"or"` or `"and"`). |
| `standardize` | `True` | Whether predictors are standardized before each Lasso fit. |
| `max_iter` | `10000` | Maximum Lasso iterations. |
| `tol` | `1e-4` | Lasso optimization tolerance. |
| `random_state` | `0` | Random seed for reproducibility. |
| `n_jobs` | `-1` | Parallel jobs for StARS subsample fitting. |
| `require_samples_gt_taxa` | `True` | Raise when `n_samples <= n_taxa`; set `False` to allow the (potentially unstable) high-dimensional regime. |
| `stars_progress_callback` | `None` | Optional callback `(step, total, alpha, instability)` called after each lambda in StARS. |

## Practical tuning guidance

- Start from defaults, then tune `n_lambdas`, `n_subsamples`, and `stars_threshold`.
- Increase `n_subsamples` for more stable StARS estimates (higher cost).
- Reduce `n_lambdas` and/or `n_subsamples` for faster exploratory runs.
- In very high-dimensional settings (taxa >> samples), either pre-filter low-prevalence
	taxa before fitting or pass `require_samples_gt_taxa=False` to allow the run.

## Development

```bash
ruff check .
black --check .
mypy src
pytest -m "not r_validation"
```

## Test suite

All tests use synthetic compositional data
(`tests/synthetic_data.py`); no external dataset is required.

```bash
# Fast suite (no R required)
pytest -m "not r_validation"

# Full suite, including R parity checks
pytest
```

| Test area | File(s) | Simple description (what it checks) | Current outcome |
| --- | --- | --- | --- |
| Transforms | `tests/test_transforms.py` | Checks data preprocessing is correct: CLR is mathematically consistent, zero handling does not break values, and invalid parameters are rejected. | Pass |
| Edge cases | `tests/test_edge_cases.py` | Checks the code stops early on problematic input (for example only one taxon, all-zero taxon, or too few samples). | Pass |
| Graph | `tests/test_graph.py` | Checks conversion from matrices to a graph keeps the expected nodes, edges, and labels. | Pass |
| Neighborhood | `tests/test_neighborhood.py` | Checks neighborhood selection creates a valid symmetric network and that strict symmetrization (`and`) is not denser than permissive (`or`). | Pass |
| Pipeline | `tests/test_pipeline.py` | Checks the full Python pipeline: recovers meaningful edges on strong-signal synthetic data and stays reproducible with fixed seed. | Pass |
| StARS | `tests/test_stars.py` | Checks model selection logic: deterministic behavior with fixed seed, valid instability curve, and correct selection rule. | Pass |
| R parity components (opt-in) | `tests/test_r_components.py` | Checks single building blocks against R (CLR, lambda grid, per-neighborhood support behavior, and weight symmetrization). | Pass |
| R parity StARS isolation (opt-in) | `tests/test_r_stars_isolation.py` | Checks the StARS instability curve against R using the same lambda path and the same subsamples. | Pass |
| R parity end-to-end (opt-in) | `tests/test_r_validation.py` | Checks full graph parity against R (nodes, edge overlap, degree profile, and edge-weight agreement with strict thresholds). | Pass |

Status (October 2026): **27 passed** (19 fast + 8 R parity).

Notes:

- Edge-recovery assertions use the strong-signal generator because on weak,
  high-variance compositional data even the R reference implementation recovers
  no ground-truth edges; the Python and R graphs agree almost exactly on such
  data (shared-edge F1 ≈ 0.95), so parity is the meaningful check there.
- Result equivalences use tolerances that absorb solver differences
  (`sklearn`/`huge`/`glmnet`) while enforcing agreement of support, StARS curves
  and selected lambda.

## Validation against R

An optional validation path compares this implementation against the R SpiecEasi
package through rpy2.

- Script: `examples/validate_against_r.py`
- Test modules: `tests/test_r_validation.py`, `tests/test_r_components.py`,
  `tests/test_r_stars_isolation.py`

Those checks are marked with `@pytest.mark.r_validation` and are skipped in
standard CI environments where R is unavailable. Running them requires R with
the `SpiecEasi`, `pulsar`, `glmnet` and `rpy2` dependencies.

### What to show in a Python vs R graph comparison

To communicate parity clearly (also for non-experts), report these metrics:

| Graph aspect | Metric to report | Meaning | Typical strict target in this repo |
| --- | --- | --- | --- |
| Nodes | Node count match | Checks both graphs have the same number of taxa/nodes. | exact match |
| Nodes | Degree Spearman, Degree MAE | Compares node connectivity profile (how many links each node has). Spearman is diagnostic; MAE is asserted. | MAE <= 0.50 |
| Edges | Precision, Recall, F1, Jaccard | Compares which edges are present/absent in Python vs R. | Precision >= 0.80, Recall >= 0.95, F1 >= 0.85, Jaccard >= 0.75 |
| Edges | Shared edge count | Number of common edges used for weight comparison. | >= 15 |
| Edge weights | Pearson correlation | Linear agreement of weights on shared edges. | >= 0.95 |
| Edge weights | Spearman correlation | Rank agreement of edge strength ordering. | >= 0.90 |
| Edge weights | MAE | Average absolute weight difference on shared edges. | <= 0.10 |
| Edge weights | Sign agreement | Fraction of shared edges with matching sign (+/-). | >= 0.95 |

These checks are implemented in `tests/test_r_validation.py` and validated on
multiple strong-signal synthetic datasets.

### Run the strict R-parity checks

```bash
cd spiec-easi-python
.venv/bin/python -m pytest tests/test_r_validation.py -q
```
