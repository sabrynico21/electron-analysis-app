"""Phase E — statistical equivalence between R and Python (no index injection).

Unlike ``test_equivalence*.py``, R's StARS subsample indices are **not**
injected: each implementation draws its own subsamples, so the two StARS runs
cannot be identical. This module therefore compares *distributions* over
``n_seeds`` independent seeds:

* the selected-lambda index,
* the StARS instability curve (raw and cummax), averaged over seeds,
* the edge count and the per-edge selection frequency.

Reference runs come from ``tests/r_reference/generate_stats_runs.R``
(``tests/golden_stats/``).

Note on the data family: for these synthetic counts StARS selects
``lambda = lambda.max`` in every seed, so the selected index is degenerate. The
discriminating statistics are the instability curves.

The thresholds asserted here are deliberate *sanity* bounds for a distributional
check, not claims of bit-level equivalence; the measured values are also written
to ``docs/STATISTICAL_EQUIVALENCE.md``.

Heavy: run with ``-m slow``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

GOLDEN_STATS = Path(__file__).parent / "golden_stats"
CASE_IDS = ["stats_p10_n60_s401", "stats_p10_n100_s402"]

# --- sanity thresholds (distributional, intentionally loose) ----------------
MAX_MEAN_RAW_CURVE_ABS_DIFF = 0.05    # mean raw-instability curve, per lambda
MAX_MEAN_MONO_CURVE_ABS_DIFF = 0.05   # mean cummax-instability curve, per lambda
MAX_MEAN_EDGE_REL_DIFF = 0.50         # |mean_py - mean_r| / max(mean_r, 1)
MIN_EDGE_FREQ_PEARSON = 0.50          # per-edge selection-frequency correlation

pytestmark = pytest.mark.skipif(
    not (GOLDEN_STATS / CASE_IDS[0] / "r_runs.json").is_file(),
    reason=(
        "Statistical reference not found. Generate with "
        "`Rscript tests/r_reference/generate_stats_runs.R`."
    ),
)


def _read_counts(path: Path) -> np.ndarray:
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    return np.asarray(
        [[float(v) for v in ln.split("\t")[1:]] for ln in lines[1:]], dtype=float
    )


def _edge_key(pair) -> tuple[str, str]:
    return (str(pair[0]), str(pair[1]))


def _frequencies(runs: list[dict]) -> dict[tuple[str, str], float]:
    keys = {_edge_key(e) for run in runs for e in run["edges"]}
    return {
        key: sum(1 for run in runs if any(_edge_key(e) == key for e in run["edges"]))
        / len(runs)
        for key in sorted(keys)
    }


def _python_runs(counts: np.ndarray, n_seeds: int) -> list[dict]:
    """Run the Python pipeline once per seed, drawing its own subsamples."""
    from spieceasi import spiec_easi_mb

    taxa = [f"f{i + 1:02d}" for i in range(counts.shape[1])]
    runs: list[dict] = []
    for seed in range(1, n_seeds + 1):
        result = spiec_easi_mb(
            counts,
            n_lambdas=20,
            lambda_min_ratio=1e-3,
            n_subsamples=20,
            stars_threshold=0.05,
            sym_rule="or",
            random_state=seed,
            require_samples_gt_taxa=False,
            n_jobs=-1,
        )
        adj = result.adjacency.astype(int)
        edges = [
            [taxa[i], taxa[j]]
            for i in range(adj.shape[0])
            for j in range(i + 1, adj.shape[1])
            if adj[i, j] != 0
        ]
        runs.append(
            {
                "seed": seed,
                "n_edges": len(edges),
                "edges": edges,
                "raw_instability": np.asarray(result.instability, dtype=float).tolist(),
                "monotone_instability": np.asarray(
                    result.monotone_instability, dtype=float
                ).tolist(),
            }
        )
    return runs


def _curve(runs: list[dict], key: str) -> np.ndarray:
    return np.asarray([run[key] for run in runs], dtype=float)


def _summary(case: str, r_runs: list[dict], py_runs: list[dict]) -> dict:
    r_raw = _curve(r_runs, "raw_instability")
    py_raw = _curve(py_runs, "raw_instability")
    r_mono = _curve(r_runs, "monotone_instability")
    py_mono = _curve(py_runs, "monotone_instability")

    r_raw_mean, py_raw_mean = r_raw.mean(axis=0), py_raw.mean(axis=0)
    r_mono_mean, py_mono_mean = r_mono.mean(axis=0), py_mono.mean(axis=0)

    r_edges = np.array([run["n_edges"] for run in r_runs], dtype=float)
    py_edges = np.array([run["n_edges"] for run in py_runs], dtype=float)

    r_freq, py_freq = _frequencies(r_runs), _frequencies(py_runs)
    keys = sorted(set(r_freq) | set(py_freq))
    rv = np.array([r_freq.get(k, 0.0) for k in keys])
    pv = np.array([py_freq.get(k, 0.0) for k in keys])
    pearson = (
        float(np.corrcoef(rv, pv)[0, 1])
        if rv.std() > 0 and pv.std() > 0
        else float("nan")
    )

    return {
        "case": case,
        "n_seeds": len(r_runs),
        "raw_curve_max_abs_diff": float(np.abs(py_raw_mean - r_raw_mean).max()),
        "raw_curve_mean_abs_diff": float(np.abs(py_raw_mean - r_raw_mean).mean()),
        "mono_curve_max_abs_diff": float(np.abs(py_mono_mean - r_mono_mean).max()),
        "mono_curve_mean_abs_diff": float(np.abs(py_mono_mean - r_mono_mean).mean()),
        "r_opt_index_counts": {
            str(i): sum(1 for run in r_runs if run["opt_index"] == i)
            for i in sorted({run["opt_index"] for run in r_runs})
        },
        "py_opt_index_counts": "see docs (Python index derived from selected_lambda)",
        "r_mean_edges": float(r_edges.mean()),
        "py_mean_edges": float(py_edges.mean()),
        "r_sd_edges": float(r_edges.std(ddof=1)) if len(r_edges) > 1 else 0.0,
        "py_sd_edges": float(py_edges.std(ddof=1)) if len(py_edges) > 1 else 0.0,
        "n_distinct_edges_r": len(r_freq),
        "n_distinct_edges_py": len(py_freq),
        "edge_freq_pearson": pearson,
        "r_raw_mean_curve": r_raw_mean.tolist(),
        "py_raw_mean_curve": py_raw_mean.tolist(),
        "r_mono_mean_curve": r_mono_mean.tolist(),
        "py_mono_mean_curve": py_mono_mean.tolist(),
    }


@pytest.mark.slow
@pytest.mark.parametrize("case", CASE_IDS)
def test_statistical_equivalence(case: str, request) -> None:
    ref = json.loads((GOLDEN_STATS / case / "r_runs.json").read_text())
    counts = _read_counts(GOLDEN_STATS / case / "counts.tsv")
    py_runs = _python_runs(counts, int(ref["n_seeds"]))
    report = _summary(case, ref["runs"], py_runs)

    cache = getattr(request.config, "_stat_eq_reports", None)
    if cache is None:
        cache = request.config._stat_eq_reports = {}
    cache[case] = report

    assert report["raw_curve_max_abs_diff"] <= MAX_MEAN_RAW_CURVE_ABS_DIFF, (
        f"{case}: mean raw instability curve differs by "
        f"{report['raw_curve_max_abs_diff']:.4f}"
    )
    assert report["mono_curve_max_abs_diff"] <= MAX_MEAN_MONO_CURVE_ABS_DIFF, (
        f"{case}: mean monotone instability curve differs by "
        f"{report['mono_curve_max_abs_diff']:.4f}"
    )
    rel = abs(report["py_mean_edges"] - report["r_mean_edges"]) / max(
        report["r_mean_edges"], 1.0
    )
    assert rel <= MAX_MEAN_EDGE_REL_DIFF, (
        f"{case}: mean edge count R={report['r_mean_edges']:.2f} "
        f"Python={report['py_mean_edges']:.2f}"
    )
    pearson = report["edge_freq_pearson"]
    assert pearson >= MIN_EDGE_FREQ_PEARSON or np.isnan(pearson), (
        f"{case}: per-edge frequency Pearson r={pearson:.3f}"
    )
