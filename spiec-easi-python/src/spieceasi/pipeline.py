"""End-to-end SPIEC-EASI MB pipeline."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import networkx as nx
import numpy as np

from ._types import ArrayLike2D
from .graph import to_networkx_graph
from .neighborhood import (
    MBNeighborhoodResult,
    SymRule,
    make_alpha_path,
    mb_path_neighborhood_selection,
)
from .stars import StARSResult, stars_select
from .transforms import ZeroMethod, clr_transform, coerce_taxa_names


@dataclass(frozen=True)
class SpiecEasiMBResult:
    """Output bundle from ``spiec_easi_mb``.

    Attributes
    ----------
    adjacency
        Symmetric binary adjacency matrix.
    weights
        Symmetric weighted matrix derived from MB coefficients.
    graph
        networkx graph with node labels and edge weight attributes.
    selected_lambda
        Selected regularization parameter from StARS.
    alpha_path
        Full alpha path evaluated by StARS.
    instability
        Raw instability per alpha.
    monotone_instability
        Monotonized instability curve used for model selection.
    directed_coefficients
        Directed MB coefficient matrix prior to symmetrization.
    taxa_names
        Taxa labels aligned to matrix columns.
    """

    adjacency: np.ndarray
    weights: np.ndarray
    graph: nx.Graph
    selected_lambda: float
    alpha_path: np.ndarray
    instability: np.ndarray
    monotone_instability: np.ndarray
    directed_coefficients: np.ndarray
    taxa_names: list[str]


def spiec_easi_mb(
    data: ArrayLike2D,
    *,
    pseudocount: float = 1.0,
    zero_method: ZeroMethod = "additive",
    sym_rule: SymRule = "or",
    n_lambdas: int = 20,
    lambda_min_ratio: float = 1e-3,
    alpha_path: np.ndarray | None = None,
    n_subsamples: int = 20,
    subsample_size: int | None = None,
    subsample_fraction: float | None = None,
    subsamples: Sequence[np.ndarray] | None = None,
    stars_threshold: float = 0.05,
    standardize: bool = True,
    max_iter: int = 10_000,
    tol: float = 1e-4,
    random_state: int = 0,
    n_jobs: int = -1,
    prefer: str = "processes",
    require_samples_gt_taxa: bool = True,
    stars_progress_callback: Callable[[int, int, float, float], None] | None = None,
) -> SpiecEasiMBResult:
    """Run SPIEC-EASI MB with StARS model selection.

    Parameters
    ----------
    data
        Count matrix with samples as rows and taxa as columns.
    pseudocount
        Zero-handling parameter passed to CLR transformation.
    zero_method
        Zero-handling strategy: ``'additive'`` or ``'multiplicative'``.
    sym_rule
        MB neighborhood symmetrization rule: ``'or'`` or ``'and'``.
    n_lambdas
        Number of lambda values if ``alpha_path`` is not provided.
    lambda_min_ratio
        Ratio defining minimum lambda relative to lambda max.
    alpha_path
        Optional custom decreasing lambda/alpha path.
    n_subsamples
        Number of StARS subsamples.
    subsample_size
        Optional absolute subsample size for StARS.
    subsample_fraction
        Optional relative subsample size for StARS.
    subsamples
        Optional fixed 0-based subsample index arrays for StARS. When
        provided, internal random subsampling is bypassed.
    stars_threshold
        Instability threshold for StARS selection.
    standardize
        Whether to standardize predictors in each neighborhood regression.
    max_iter
        Max Lasso iterations.
    tol
        Lasso tolerance.
    random_state
        Random seed used in StARS subsampling and deterministic regressions.
    n_jobs
        Parallel worker count for StARS.
    require_samples_gt_taxa
        If ``True`` (default), raise when ``n_samples <= n_taxa``. Set to
        ``False`` to allow the (potentially unstable) high-dimensional regime.
    stars_progress_callback
        Optional callback invoked after each alpha is evaluated by StARS, with
        signature ``(step, total, alpha, instability)``.

    Returns
    -------
    SpiecEasiMBResult
        Adjacency, weighted graph, selected lambda, and diagnostics.
    """
    clr_result = clr_transform(data, pseudocount=pseudocount, zero_method=zero_method)
    X_clr = clr_result.clr
    taxa_names = coerce_taxa_names(clr_result.taxa_names, X_clr.shape[1])

    n_samples, n_taxa = X_clr.shape
    if require_samples_gt_taxa and n_samples <= n_taxa:
        msg = (
            "The sample count must exceed the number of taxa for stable MB "
            f"regressions (n_samples={n_samples}, n_taxa={n_taxa})."
        )
        raise ValueError(msg)

    if alpha_path is None:
        alpha_path = make_alpha_path(
            X_clr,
            n_lambdas=n_lambdas,
            lambda_min_ratio=lambda_min_ratio,
        )
    else:
        alpha_path = np.asarray(alpha_path, dtype=float)
        if alpha_path.ndim != 1 or alpha_path.size < 2:
            msg = "alpha_path must be a 1D array with at least 2 values."
            raise ValueError(msg)
        if np.any(alpha_path <= 0):
            msg = "alpha_path must contain strictly positive values."
            raise ValueError(msg)
        if np.any(np.diff(alpha_path) > 0):
            msg = (
                "alpha_path must be non-increasing (huge requires a "
                "non-increasing lambda path)."
            )
            raise ValueError(msg)

    stars_result: StARSResult = stars_select(
        X_clr,
        alpha_path=alpha_path,
        n_subsamples=n_subsamples,
        subsample_size=subsample_size,
        subsample_fraction=subsample_fraction,
        subsamples=subsamples,
        threshold=stars_threshold,
        sym_rule=sym_rule,
        standardize=standardize,
        max_iter=max_iter,
        tol=tol,
        random_state=random_state,
        n_jobs=n_jobs,
        prefer=prefer,
        on_alpha_done=stars_progress_callback,
    )

    final_mb: MBNeighborhoodResult = mb_path_neighborhood_selection(
        X_clr,
        alpha_path=alpha_path,
        selected_index=stars_result.selected_index,
        sym_rule=sym_rule,
        max_iter=max_iter,
        tol=tol,
    )

    graph = to_networkx_graph(final_mb.adjacency, final_mb.weight_matrix, taxa_names)

    return SpiecEasiMBResult(
        adjacency=final_mb.adjacency,
        weights=final_mb.weight_matrix,
        graph=graph,
        selected_lambda=stars_result.selected_alpha,
        alpha_path=stars_result.alpha_path,
        instability=stars_result.instability,
        monotone_instability=stars_result.monotone_instability,
        directed_coefficients=final_mb.directed_coefficients,
        taxa_names=taxa_names,
    )
