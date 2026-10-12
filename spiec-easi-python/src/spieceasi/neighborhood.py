"""Meinshausen-Buhlmann neighborhood selection for SPIEC-EASI.

The nodewise Lasso is solved by :mod:`spieceasi.mb_solver`, a faithful
re-implementation of ``huge::mb()`` (warm-started lambda path, sequential
strong rule, huge's coordinate order and stopping criteria). Only the
huge/SpiecEasi-congruent solver is supported: R's ``method='mb'`` uses huge's
coordinate descent, **not** glmnet.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .mb_solver import (
    huge_correlation,
    solve_mb_path_huge_scale,
    symmetrize_path,
)

SymRule = Literal["or", "and"]
CoefficientScale = Literal["standardized", "original"]


# The nodewise Lasso solver (cold single-alpha and warm-started full path)
# lives in :mod:`spieceasi.mb_solver`.


@dataclass(frozen=True)
class MBNeighborhoodResult:
    """Result of MB neighborhood selection.

    Attributes
    ----------
    adjacency
        Binary symmetric adjacency matrix.
    weight_matrix
        Symmetric weighted matrix derived from directed Lasso coefficients.
    directed_coefficients
        Directed coefficient matrix where entry ``(j, k)`` is the coefficient of
        taxon ``k`` in the regression predicting taxon ``j``.
    """

    adjacency: np.ndarray
    weight_matrix: np.ndarray
    directed_coefficients: np.ndarray


def _require_huge_scale(
    standardize: bool, coefficient_scale: CoefficientScale
) -> None:
    if not standardize or coefficient_scale != "standardized":
        msg = (
            "Only the huge/SpiecEasi-congruent MB solver is supported "
            "(standardize=True, coefficient_scale='standardized'). The "
            "sklearn/glmnet-style branch was removed: R's method='mb' uses "
            "huge's coordinate descent, not glmnet."
        )
        raise NotImplementedError(msg)


def _fit_single_neighborhood(
    X: np.ndarray,
    target_idx: int,
    *,
    alpha: float,
    standardize: bool,
    max_iter: int,
    tol: float,
    random_state: int | None,
    coefficient_scale: CoefficientScale = "standardized",
) -> np.ndarray:
    """Coefficients of predictors for a single node (single-alpha solve).

    A one-point path is solved, so the result is bit-identical to what
    :func:`mb_neighborhood_selection` returns (no warm start is possible with
    a single lambda). ``random_state`` is accepted for compatibility only.
    """
    _require_huge_scale(standardize, coefficient_scale)
    directed_path = solve_mb_path_huge_scale(
        X, np.array([alpha], dtype=float), max_iter=max_iter, tol=tol
    )
    return directed_path[0, target_idx, :].copy()


def _result_from_directed(
    directed: np.ndarray, sym_rule: SymRule
) -> MBNeighborhoodResult:
    coef_abs = np.abs(directed)
    if sym_rule == "or":
        adjacency = np.logical_or(coef_abs > 0, coef_abs.T > 0)
    elif sym_rule == "and":
        adjacency = np.logical_and(coef_abs > 0, coef_abs.T > 0)
    else:
        msg = f"Unsupported sym_rule '{sym_rule}'."
        raise ValueError(msg)
    adjacency = adjacency.astype(int)
    np.fill_diagonal(adjacency, 0)
    weight_matrix = 0.5 * (directed + directed.T)
    weight_matrix *= adjacency
    np.fill_diagonal(weight_matrix, 0.0)
    return MBNeighborhoodResult(
        adjacency=adjacency,
        weight_matrix=weight_matrix,
        directed_coefficients=directed,
    )


def mb_neighborhood_selection(
    X: np.ndarray,
    *,
    alpha: float,
    sym_rule: SymRule = "or",
    standardize: bool = True,
    max_iter: int = 10_000,
    tol: float = 1e-4,
    random_state: int | None = None,
    coefficient_scale: CoefficientScale = "standardized",
) -> MBNeighborhoodResult:
    """Infer an undirected network via MB neighborhood selection (one alpha).

    Solves a one-point lambda path, hence identical to
    :func:`mb_path_neighborhood_selection` at that single point.
    ``random_state`` is accepted for backward compatibility and unused (the
    solver is deterministic).
    """
    _require_huge_scale(standardize, coefficient_scale)
    if alpha <= 0:
        msg = "alpha must be strictly positive."
        raise ValueError(msg)
    directed_path = solve_mb_path_huge_scale(
        X, np.array([alpha], dtype=float), max_iter=max_iter, tol=tol
    )
    return _result_from_directed(directed_path[0], sym_rule)


def mb_path_neighborhood_selection(
    X: np.ndarray,
    *,
    alpha_path: np.ndarray,
    selected_index: int,
    sym_rule: SymRule = "or",
    max_iter: int = 10_000,
    tol: float = 1e-4,
) -> MBNeighborhoodResult:
    """Warm-started solve of the whole MB lambda path; network at ``selected_index``.

    Matches R's ``spiec.easi``/``pulsar`` refit, which takes ``path[[opt.index]]``
    from a single warm-started solve over the full path (``huge::mb``).
    """
    alpha_path = np.asarray(alpha_path, dtype=float)
    directed_path = solve_mb_path_huge_scale(
        X, alpha_path, max_iter=max_iter, tol=tol
    )
    return _result_from_directed(directed_path[int(selected_index)], sym_rule)


def alpha_max_for_mb(X: np.ndarray) -> float:
    """Compute a global alpha upper-bound for MB Lasso regressions.

    Matches ``huge`` ``lambda.max`` semantics on the standardized (correlation)
    scale: the largest magnitude off-diagonal correlation
    ``max_{j != k} |corr(X_j, X_k)|``.

    The value is taken from the *same* correlation matrix the solver builds
    (:func:`spieceasi.mb_solver.huge_correlation`). Computing it via a different
    BLAS path (e.g. a separate nodewise ``X_other.T @ y`` product) can differ by
    1 ULP, which is enough for the max-correlation coefficient to survive the
    soft threshold at ``lambda == alpha_max`` and leave a spurious ~1e-16 edge in
    the nominal "empty" top model. Using the solver's own ``S`` makes
    ``alpha_max == max|S_off|`` exactly, so the top model is empty as in R.
    """
    S = huge_correlation(X)
    off_diagonal = S.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    alpha_max = float(np.max(np.abs(off_diagonal)))
    return max(alpha_max, 1e-6)


def make_alpha_path(
    X: np.ndarray,
    *,
    n_lambdas: int = 25,
    lambda_min_ratio: float = 1e-2,
) -> np.ndarray:
    """Create a decreasing alpha path suitable for MB neighborhood selection."""
    if n_lambdas < 2:
        msg = "n_lambdas must be at least 2."
        raise ValueError(msg)
    if not 0 < lambda_min_ratio < 1:
        msg = "lambda_min_ratio must be in (0, 1)."
        raise ValueError(msg)

    alpha_max = alpha_max_for_mb(X)
    alpha_min = alpha_max * lambda_min_ratio
    path = np.geomspace(alpha_max, alpha_min, num=n_lambdas)
    # Pin the endpoints exactly. `geomspace` computes exp(log(alpha_max)), which
    # can round 1 ULP *below* alpha_max; at such a lambda the max-correlation
    # coefficient survives the soft threshold and the sparsest model is not
    # empty. R's ``lambda.max`` guarantees an empty top model, so reproduce that.
    path[0] = alpha_max
    path[-1] = alpha_min
    return path
