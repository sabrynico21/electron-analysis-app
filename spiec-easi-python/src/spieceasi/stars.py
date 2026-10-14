"""StARS regularization selection for SPIEC-EASI MB."""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from joblib import Parallel, delayed

from .mb_solver import solve_mb_path_huge_scale, symmetrize_path
from .neighborhood import SymRule, _require_huge_scale


@dataclass(frozen=True)
class StARSResult:
    """Outputs from StARS selection.

    Attributes
    ----------
    selected_alpha
        Chosen regularization parameter.
    alpha_path
        Candidate alpha values in decreasing order.
    instability
        Raw instability per alpha.
    monotone_instability
        Monotonized instability curve used for thresholding.
    selected_index
        0-based index of ``selected_alpha`` in ``alpha_path``.
    """

    selected_alpha: float
    alpha_path: np.ndarray
    instability: np.ndarray
    monotone_instability: np.ndarray
    selected_index: int


def _infer_subsample_adjacency_path(
    X_sub: np.ndarray,
    *,
    alpha_path: np.ndarray,
    sym_rule: SymRule,
    max_iter: int,
    tol: float,
) -> np.ndarray:
    """Symmetrized adjacency path ``(n_lambda, d, d)`` for one subsample.

    One warm-started solve of the whole lambda path, mirroring ``huge::mb`` /
    the pulsar ``fun`` applied to ``data[ind, ]``. A per-lambda cold solve can
    reach a different near-optimal support in ill-conditioned regimes.
    """
    directed_path = solve_mb_path_huge_scale(
        X_sub, alpha_path, max_iter=max_iter, tol=tol
    )
    return symmetrize_path(directed_path, sym_rule)


def _instability_from_adj_sum(
    adj_sum: np.ndarray, n_subsamples: int
) -> np.ndarray:
    """StARS instability per lambda from summed adjacency paths.

    Matches ``pulsar::stars.stability``: ``4*sum(P*(1-P))/(p*(p-1))`` with P
    symmetric, computed here as ``4*mean`` over the strict upper triangle of
    ``P*(1-P)`` (the full-matrix sum is twice the upper-triangle sum and the
    factor 4 absorbs it).
    """
    P = adj_sum / float(n_subsamples)
    edge = 4.0 * P * (1.0 - P)
    d = edge.shape[-1]
    if d < 2:
        return np.zeros(edge.shape[0], dtype=float)
    iu = np.triu_indices(d, k=1)
    return edge[:, iu[0], iu[1]].mean(axis=1)


def _stars_index(monotone_instability: np.ndarray, threshold: float) -> int:
    """Replicate ``pulsar:::.starsind`` exactly.

    ``opt.index = max(first index with summary >= thresh - 1, 0)``; if no lambda
    reaches the threshold, the last (least regularized) index is returned with
    the same warning pulsar emits.
    """
    crossed = np.where(monotone_instability >= threshold)[0]
    if crossed.size:
        return max(int(crossed[0]) - 1, 0)
    warnings.warn(
        "Optimal lambda may be outside the supplied values: no lambda reached "
        "the StARS instability threshold; selecting the least-regularized lambda.",
        RuntimeWarning,
        stacklevel=2,
    )
    return int(monotone_instability.size - 1)


def _validate_decreasing(alpha_path: np.ndarray) -> None:
    if np.any(np.diff(alpha_path) > 0):
        msg = (
            "alpha_path must be non-increasing (huge requires a non-increasing "
            "lambda path)."
        )
        raise ValueError(msg)


def _instability_from_adjacencies(adjacencies: list[np.ndarray]) -> float:
    if len(adjacencies) == 0:
        msg = "At least one adjacency matrix is required."
        raise ValueError(msg)

    A = np.stack(adjacencies, axis=0)
    p_select = np.mean(A, axis=0)
    # Match pulsar::stars.stability: 4 * sum(P * (1-P)) / (p * (p-1)),
    # where sum is over the full symmetric matrix. Averaging on the strict
    # upper triangle therefore uses 4 * P * (1-P) per undirected edge.
    edge_instability = 4.0 * p_select * (1.0 - p_select)

    upper = np.triu_indices(edge_instability.shape[0], k=1)
    if upper[0].size == 0:
        return 0.0
    return float(np.mean(edge_instability[upper]))


def infer_subsample_size(
    n_samples: int,
    *,
    subsample_size: int | None = None,
    subsample_fraction: float | None = None,
) -> int:
    """Resolve the StARS subsample size following ``pulsar::.ratcheck``.

    An explicit ``subsample_size`` (or ``subsample_fraction``) always wins. When
    neither is given the pulsar default is used: ratio ``0.8`` for ``n <= 144``,
    otherwise ratio ``10*sqrt(n)/n``; in both cases the size is
    ``floor(ratio * n)`` as in ``pulsar`` (``floor(n * subsample.ratio)``).
    """
    if subsample_size is not None:
        return int(subsample_size)
    if subsample_fraction is None:
        # Match pulsar::.ratcheck defaults.
        if n_samples > 144:
            return int(np.floor(10 * np.sqrt(n_samples)))
        return int(np.floor(0.8 * n_samples))
    if not 0 < subsample_fraction < 1:
        msg = "subsample_fraction must be in (0, 1)."
        raise ValueError(msg)
    return int(np.floor(subsample_fraction * n_samples))


def stars_select(
    X: np.ndarray,
    *,
    alpha_path: np.ndarray,
    n_subsamples: int = 20,
    subsample_size: int | None = None,
    subsample_fraction: float | None = None,
    threshold: float = 0.05,
    sym_rule: SymRule = "or",
    standardize: bool = True,
    max_iter: int = 10_000,
    tol: float = 1e-4,
    random_state: int = 0,
    n_jobs: int = -1,
    prefer: str = "threads",
    on_alpha_done: Callable[[int, int, float, float], None] | None = None,
    subsamples: Sequence[np.ndarray] | None = None,
) -> StARSResult:
    """Select alpha via StARS for MB neighborhood selection.

    Each subsample's full lambda path is solved in one warm-started pass (as
    ``huge::mb`` does), matching R's pulsar pipeline.

    Notes
    -----
    Candidate alphas must be non-increasing (sparser to denser). When
    ``subsamples`` is provided, those 0-based index arrays are used directly and
    random subsampling is bypassed.

    ``prefer`` selects the joblib backend used to run the independent subsamples.
    The default is ``"threads"``: the compiled MB core releases the GIL, the
    subsamples are aggregated in a fixed order, and a thread pool cannot
    re-execute the interpreter — a *process* pool (``"processes"``, the legacy
    behaviour) breaks frozen executables, where loky relaunches the binary with
    its own bootstrap arguments. Results are bit-identical either way.
    """
    _require_huge_scale(standardize, "standardized")
    n_samples = X.shape[0]
    alpha_path = np.asarray(alpha_path, dtype=float)
    _validate_decreasing(alpha_path)
    if n_subsamples < 2 and subsamples is None:
        msg = "n_subsamples must be at least 2 for StARS."
        raise ValueError(msg)
    if threshold <= 0 or threshold >= 1:
        msg = "threshold must be in (0, 1)."
        raise ValueError(msg)

    if subsample_size is not None and subsample_fraction is not None:
        msg = "Specify only one of subsample_size or subsample_fraction."
        raise ValueError(msg)

    if subsamples is not None:
        if len(subsamples) < 2:
            msg = "subsamples must contain at least 2 index arrays for StARS."
            raise ValueError(msg)

        subsample_indices: list[np.ndarray] = []
        for i, idx in enumerate(subsamples):
            arr = np.asarray(idx, dtype=int)
            if arr.ndim != 1:
                msg = f"subsamples[{i}] must be a 1D array of indices."
                raise ValueError(msg)
            if arr.size <= 1:
                msg = f"subsamples[{i}] must contain at least 2 indices."
                raise ValueError(msg)
            if np.any(arr < 0) or np.any(arr >= n_samples):
                msg = (
                    f"subsamples[{i}] contains out-of-range values for "
                    f"n_samples={n_samples}."
                )
                raise ValueError(msg)
            if np.unique(arr).size != arr.size:
                msg = f"subsamples[{i}] contains duplicate indices."
                raise ValueError(msg)
            subsample_indices.append(arr)
    else:
        subsample_size = infer_subsample_size(
            n_samples,
            subsample_size=subsample_size,
            subsample_fraction=subsample_fraction,
        )

        if subsample_size <= 1 or subsample_size > n_samples:
            msg = (
                "Invalid subsample_size; must satisfy 1 < subsample_size <= n_samples "
                f"(got {subsample_size}, n_samples={n_samples})."
            )
            raise ValueError(msg)

        rng = np.random.default_rng(random_state)
        subsample_indices = [
            rng.choice(n_samples, size=subsample_size, replace=False)
            for _ in range(n_subsamples)
        ]

    # One warm-started path solve per subsample (as huge::mb / pulsar fun).
    # The subsamples are independent, so the result is aggregated in a fixed
    # order and is identical for every backend/n_jobs.
    adj_paths = Parallel(n_jobs=n_jobs, prefer=prefer)(
        delayed(_infer_subsample_adjacency_path)(
            X[idx, :],
            alpha_path=alpha_path,
            sym_rule=sym_rule,
            max_iter=max_iter,
            tol=tol,
        )
        for idx in subsample_indices
    )
    adj_sum = np.sum(np.stack(adj_paths, axis=0), axis=0)  # (n_lambda, d, d)

    instability = _instability_from_adj_sum(adj_sum, len(subsample_indices))
    total_alphas = int(alpha_path.size)
    for i in range(total_alphas):
        if on_alpha_done is not None:
            on_alpha_done(
                i + 1, total_alphas, float(alpha_path[i]), float(instability[i])
            )

    monotone_instability = np.maximum.accumulate(instability)
    selected_idx = _stars_index(monotone_instability, threshold)

    return StARSResult(
        selected_alpha=float(alpha_path[selected_idx]),
        alpha_path=alpha_path,
        instability=instability,
        monotone_instability=monotone_instability,
        selected_index=selected_idx,
    )
