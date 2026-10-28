"""Meinshausen-Bühlmann neighborhood selection (port of ``huge`` 2.0.1).

Faithful re-implementation of the coordinate-descent solver in
``huge/src/huge_core.cpp`` (``huge::mb``, L679-812) including the sequential
strong rule, the exact update order, ``thol = 1e-4``, ``MAX_ITER = 10000`` and
the strict soft-thresholding rule. ``scikit-learn`` is deliberately NOT used:
its Lasso has different standardisation and stopping rules.

R sources:
    huge/R/*.R         (huge.mb wrapper, lambda handling, symmetrisation)
    huge/src/huge_core.cpp  (huge::mb coordinate descent)
"""

from __future__ import annotations

import numpy as np

__all__ = ["huge_mb", "symmetrize_path"]

_THOL = 1e-4
_MAX_ITER = 10000


def _soft(x: float, thr: float) -> float:
    """``threshold_l1`` from huge_core.cpp (strict inequalities)."""
    if x > thr:
        return x - thr
    if x < -thr:
        return x + thr
    return 0.0


def _partial_residual(S, m, j, idx_a, size_a, w0) -> float:
    """``mb_partial_residual``."""
    r = S[m, j]
    for k in range(size_a):
        idx = idx_a[k]
        r -= S[j, idx] * w0[idx]
    return r


def _refine_active(S, m, ilambda, idx_a, size_a, w0, w1) -> None:
    """``mb_refine_active``: cycle coordinate descent until convergence."""
    gap_int = 1.0
    it = 0
    while gap_int > _THOL and it < _MAX_ITER:
        change_sum = 0.0
        coef_sum = 0.0
        for j in range(size_a):
            wi = idx_a[j]
            r = S[m, wi] + w0[wi]
            for k in range(size_a):
                r -= S[wi, idx_a[k]] * w0[idx_a[k]]
            w1[wi] = _soft(r, ilambda)
            coef_sum += abs(w1[wi])
            change_sum += abs(w1[wi] - w0[wi])
            w0[wi] = w1[wi]
        gap_int = change_sum / coef_sum if coef_sum > 0 else 0.0
        it += 1


def _validate_correlation_psd(correlation: np.ndarray) -> None:
    """``.huge_validate_correlation_psd`` (chol-based, with tolerance shift)."""
    d = correlation.shape[0]
    spectral_bound = max(1.0, float(np.abs(correlation).sum(axis=1).max()))
    tolerance = 100.0 * np.finfo(float).eps * max(1, d) * spectral_bound
    shifted = correlation.copy()
    shifted[np.diag_indices_from(shifted)] += tolerance
    try:
        np.linalg.cholesky(shifted)
    except np.linalg.LinAlgError as exc:  # pragma: no cover - error path
        raise ValueError("Covariance input must be positive semidefinite.") from exc


def _as_correlation(S: np.ndarray) -> np.ndarray:
    """Covariance -> correlation conversion used for symmetric input."""
    x = np.asarray(S, dtype=float)
    covariance = x.copy()
    diff = x != x.T
    covariance[diff] = x[diff] / 2.0 + x.T[diff] / 2.0
    diagonal = np.diag(covariance)
    if np.any(~np.isfinite(diagonal)) or np.any(diagonal <= 0):
        raise ValueError("Covariance input must have a positive finite diagonal.")
    inv_sd = 1.0 / np.sqrt(diagonal)
    correlation = covariance * inv_sd[:, None] * inv_sd[None, :]
    if np.any(~np.isfinite(correlation)):
        raise ValueError("Covariance input cannot produce a finite correlation matrix.")
    if np.any(np.abs(correlation) > 1 + 1e-8):
        raise ValueError("Covariance input is not a valid covariance matrix.")
    correlation = np.clip(correlation, -1.0, 1.0)
    np.fill_diagonal(correlation, 1.0)
    _validate_correlation_psd(correlation)
    return correlation


def huge_mb(S, lambda_) -> np.ndarray:
    """Solve the MB nodewise lasso path.

    Parameters
    ----------
    S : array-like of shape (d, d)
        Symmetric covariance/correlation matrix.
    lambda_ : array-like
        Strictly positive, non-increasing regularization path.

    Returns
    -------
    numpy.ndarray
        ``beta`` with shape ``(nlambda, d, d)``; ``beta[l, k, m]`` is the
        coefficient of predictor ``k`` in the regression of node ``m`` at
        ``lambda_[l]``.
    """
    Sc = _as_correlation(S)
    lam = np.asarray(lambda_, dtype=float)
    if lam.ndim != 1 or lam.size < 1:
        raise ValueError("lambda must be a one-dimensional numeric sequence.")
    if np.any(~np.isfinite(lam)) or np.any(lam <= 0):
        raise ValueError("lambda must contain strictly positive values.")
    if lam.size > 1 and np.any(lam[1:] > lam[:-1]):
        raise ValueError('lambda must be non-increasing for method = "mb".')

    d = Sc.shape[0]
    nlambda = lam.size
    beta = np.zeros((nlambda, d, d), dtype=float)

    for m in range(d):
        w0 = np.zeros(d)
        w1 = np.zeros(d)
        idx_i = np.ones(d, dtype=bool)
        idx_i[m] = False
        idx_a: list[int] = []
        size_a = 0
        grad_abs = np.abs(Sc[m, :]).copy()
        strong = np.zeros(d, dtype=bool)

        for i in range(nlambda):
            ilambda = lam[i]
            strong_thr = 2.0 * lam[i] - lam[i - 1] if i > 0 else 2.0 * lam[i]
            strong = idx_i & (grad_abs > strong_thr)

            iter_ext = 0
            certified = False
            while not certified and iter_ext < _MAX_ITER:
                gap_ext = 1
                while gap_ext != 0 and iter_ext < _MAX_ITER:
                    size_a_prev = size_a
                    for j in range(d):
                        if idx_i[j] and strong[j]:
                            r = _partial_residual(Sc, m, j, idx_a, size_a, w0)
                            grad_abs[j] = abs(r)
                            w1[j] = _soft(r, ilambda)
                            if w1[j] != 0.0:
                                idx_a.append(j)
                                size_a += 1
                                idx_i[j] = False
                            w0[j] = w1[j]
                    gap_ext = size_a - size_a_prev

                    _refine_active(Sc, m, ilambda, idx_a, size_a, w0, w1)

                    # drop zeros from the active set; keep sweeping them
                    kept = []
                    for j in idx_a:
                        if w1[j] == 0.0:
                            idx_i[j] = True
                            strong[j] = True
                        else:
                            kept.append(j)
                    size_a = len(kept)
                    idx_a = kept
                    iter_ext += 1

                # certification pass over every inactive coordinate
                violations = 0
                for j in range(d):
                    if idx_i[j]:
                        r = _partial_residual(Sc, m, j, idx_a, size_a, w0)
                        grad_abs[j] = abs(r)
                        w1[j] = _soft(r, ilambda)
                        if w1[j] != 0.0:
                            idx_a.append(j)
                            size_a += 1
                            idx_i[j] = False
                            strong[j] = True
                            violations += 1
                        w0[j] = w1[j]
                certified = violations == 0

            for j in idx_a:
                beta[i, j, m] = w1[j]

    return beta


def symmetrize_path(beta: np.ndarray, sym: str = "or") -> np.ndarray:
    """Reproduce ``huge.mb`` path symmetrisation.

    ``path = abs(beta)``; ``sym="or"`` -> ``sign(path + path.T)`` (unsigned
    0/1 adjacency); ``sym="and"`` -> ``sign(path * path.T)``.
    """
    if sym not in ("or", "and"):
        raise ValueError('sym must be exactly one of "or" or "and".')
    path = np.abs(beta)
    if sym == "or":
        return np.sign(path + np.transpose(path, (0, 2, 1)))
    return np.sign(path * np.transpose(path, (0, 2, 1)))
