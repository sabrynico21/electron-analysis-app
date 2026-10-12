"""Compare Python SPIEC-EASI MB output against the R SpiecEasi reference.

Usage
-----
python examples/validate_against_r.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from spieceasi.pipeline import SpiecEasiMBResult, spiec_easi_mb


def generate_shared_dataset(path: Path, *, seed: int = 31) -> pd.DataFrame:
    """Generate a synthetic dataset once and reuse it for both Python and R."""
    if path.exists():
        return pd.read_csv(path)

    rng = np.random.default_rng(seed)
    n_samples = 90
    n_taxa = 8
    depth = 2_000

    precision = np.eye(n_taxa)
    for i in range(n_taxa - 1):
        precision[i, i + 1] = -0.25
        precision[i + 1, i] = -0.25
    precision += 0.2 * np.eye(n_taxa)

    covariance = np.linalg.inv(precision)
    latent = rng.multivariate_normal(np.zeros(n_taxa), covariance, size=n_samples)
    probs = np.exp(latent)
    probs /= np.sum(probs, axis=1, keepdims=True)
    counts = np.vstack([rng.multinomial(depth, p) for p in probs])

    df = pd.DataFrame(counts, columns=[f"taxon_{i}" for i in range(n_taxa)])
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


def edge_set(A: np.ndarray) -> set[tuple[int, int]]:
    idx = np.triu_indices(A.shape[0], k=1)
    return {(int(i), int(j)) for i, j in zip(*idx, strict=True) if A[i, j] != 0}


@dataclass(frozen=True)
class EdgeMetrics:
    precision: float
    recall: float
    f1: float
    jaccard: float


def edge_metrics(
    r_edges: set[tuple[int, int]],
    py_edges: set[tuple[int, int]],
) -> EdgeMetrics:
    tp = len(r_edges.intersection(py_edges))
    fp = len(py_edges.difference(r_edges))
    fn = len(r_edges.difference(py_edges))

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    jaccard = tp / len(r_edges.union(py_edges)) if r_edges.union(py_edges) else 1.0
    return EdgeMetrics(precision=precision, recall=recall, f1=f1, jaccard=jaccard)


def run_r_reference(
    counts_df: pd.DataFrame,
    *,
    nlambda: int,
    lambda_min_ratio: float,
    n_subsamples: int,
    threshold: float,
) -> tuple[np.ndarray, float, np.ndarray | None]:
    """Run SpiecEasi::spiec.easi in R through rpy2."""
    from rpy2 import robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    robjects.r("""
        suppressMessages(library(SpiecEasi))

        run_spieceasi_mb <- function(
            count_df, nlambda, lambda_min_ratio, n_subsamples, threshold
        ) {
          X <- as.matrix(count_df)
          fit <- spiec.easi(
            X,
            method = "mb",
            sel.criterion = "stars",
            nlambda = nlambda,
            lambda.min.ratio = lambda_min_ratio,
            pulsar.params = list(rep.num = n_subsamples, thresh = threshold)
          )

          adj <- as.matrix(getRefit(fit))
          lambda_path <- as.numeric(getLambda(fit))
          opt_lambda <- as.numeric(getOptLambda(fit))

          stars_curve <- NULL
          if (!is.null(fit$select) && !is.null(fit$select$stars)) {
            stars_curve <- as.numeric(fit$select$stars)
          }

          list(
            adjacency = adj,
            selected_lambda = opt_lambda,
            lambda_path = lambda_path,
            stars_curve = stars_curve
          )
        }
        """)
    run_fn = robjects.globalenv["run_spieceasi_mb"]

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_result = run_fn(
            counts_df,
            nlambda,
            lambda_min_ratio,
            n_subsamples,
            threshold,
        )

    r_adj = np.asarray(r_result.rx2("adjacency"), dtype=int)
    r_lambda = float(np.asarray(r_result.rx2("selected_lambda"))[0])

    try:
        stars_curve = np.asarray(r_result.rx2("stars_curve"), dtype=float)
        if stars_curve.size == 0:
            stars_curve = None
    except Exception:
        stars_curve = None

    return r_adj, r_lambda, stars_curve


def weight_correlation(r_adj: np.ndarray, py: SpiecEasiMBResult) -> float:
    both = np.logical_and(r_adj != 0, py.adjacency != 0)
    upper = np.triu_indices_from(both, k=1)
    mask = both[upper]
    if not np.any(mask):
        return float("nan")

    r_w = r_adj[upper][mask].astype(float)
    py_w = py.weights[upper][mask]
    return float(np.corrcoef(r_w, py_w)[0, 1])


def print_diagnostics(
    py_result: SpiecEasiMBResult,
    r_adj: np.ndarray,
    r_stars_curve: np.ndarray | None,
) -> None:
    py_edge_counts = np.sum(py_result.adjacency, axis=1)
    r_edge_counts = np.sum(r_adj, axis=1)
    print("Per-taxon neighborhood sizes (Python vs R):")
    for i, (p_count, r_count) in enumerate(
        zip(py_edge_counts, r_edge_counts, strict=True)
    ):
        print(f"  taxon_{i}: {int(p_count)} vs {int(r_count)}")

    print("\nStARS monotone instability (Python):")
    for a, s in zip(py_result.alpha_path, py_result.monotone_instability, strict=True):
        print(f"  alpha={a:.6g} instability={s:.6f}")

    if r_stars_curve is not None:
        print("\nStARS curve (R, if exposed by object):")
        for i, s in enumerate(r_stars_curve):
            print(f"  idx={i} instability={float(s):.6f}")


def main() -> None:
    root = Path(__file__).resolve().parent
    data_path = root / "data" / "synthetic_counts.csv"
    counts_df = generate_shared_dataset(data_path)
    print(f"Using shared dataset: {data_path}")

    nlambda = 12
    n_subsamples = 12
    threshold = 0.05
    lambda_min_ratio = 1e-2

    py_result = spiec_easi_mb(
        counts_df,
        n_lambdas=nlambda,
        lambda_min_ratio=lambda_min_ratio,
        n_subsamples=n_subsamples,
        subsample_fraction=0.8,
        stars_threshold=threshold,
        sym_rule="or",
        random_state=17,
        n_jobs=1,
    )

    r_adj, r_lambda, r_stars_curve = run_r_reference(
        counts_df,
        nlambda=nlambda,
        lambda_min_ratio=lambda_min_ratio,
        n_subsamples=n_subsamples,
        threshold=threshold,
    )

    metrics = edge_metrics(edge_set(r_adj), edge_set(py_result.adjacency))
    corr = weight_correlation(r_adj, py_result)

    print("\nComparison summary")
    print(f"  Precision: {metrics.precision:.3f}")
    print(f"  Recall:    {metrics.recall:.3f}")
    print(f"  F1:        {metrics.f1:.3f}")
    print(f"  Jaccard:   {metrics.jaccard:.3f}")
    print(f"  R lambda:  {r_lambda:.6g}")
    print(f"  PY lambda: {py_result.selected_lambda:.6g}")
    print(f"  Weight correlation (shared edges): {corr:.3f}")

    passed = metrics.f1 > 0.85
    print(f"\nPASS: {passed} (threshold: F1 > 0.85)")

    if not passed:
        print("\nDetailed diagnostics due to divergence:")
        print_diagnostics(py_result, r_adj, r_stars_curve)


if __name__ == "__main__":
    main()
