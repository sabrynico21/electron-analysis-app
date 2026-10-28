"""Faithful port of ``Matrix::nearPD`` (Higham 2002).

R source: ``Matrix/R/nearPD.R`` (Matrix 1.7-5, pure R implementation).
The port keeps the exact iteration order, projections, tolerances and defaults
used by ``latentcor`` (``corr = TRUE`` and all other arguments default).
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = ["near_pd"]


def _norm_inf(M: np.ndarray) -> float:
    """R ``norm(M, "I")``: maximum absolute row sum."""
    return float(np.abs(M).sum(axis=1).max())


def _eigh_desc(M: np.ndarray):
    """``eigen(M, symmetric = TRUE)``: decreasing eigenvalues."""
    d, Q = np.linalg.eigh(M)
    return d[::-1], Q[:, ::-1]


def _symmetrize(M: np.ndarray) -> np.ndarray:
    """Exactly symmetric matrix via one triangle (mirrors BLAS ``dsyrk``)."""
    upper = np.triu(M)
    return upper + np.triu(M, 1).T


def near_pd(
    x,
    corr: bool = False,
    keepDiag: bool = False,
    base_matrix: bool = False,
    do2eigen: bool = True,
    doSym: bool = False,
    doDykstra: bool = True,
    eig_tol: float = 1e-6,
    conv_tol: float = 1e-7,
    posd_tol: float = 1e-8,
    maxit: int = 100,
    conv_norm_type: str = "I",
) -> np.ndarray:
    """Nearest positive-definite (correlation) matrix. Returns ``mat``."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or x.shape[0] != x.shape[1]:
        raise ValueError("near_pd expects a square matrix.")
    if not np.array_equal(x, x.T):
        x = (x + x.T) / 2.0

    n = x.shape[0]
    diagX0 = np.diag(x).copy() if keepDiag else None

    D_S = np.zeros_like(x) if doDykstra else None
    X = x.copy()
    it = 0
    converged = False
    conv = np.inf

    while it < maxit and not converged:
        Y = X.copy()
        R = (Y - D_S) if doDykstra else None

        d, Q = _eigh_desc(R if doDykstra else Y)
        p = d > eig_tol * d[0]
        if not np.any(p):
            raise ValueError("Matrix seems negative semi-definite")
        Qp = Q[:, p]
        dp = d[p]
        # R: tcrossprod(Q * rep(d[p], each = nrow(Q)), Q)
        X = _symmetrize((Qp * dp) @ Qp.T)

        if doDykstra:
            D_S = X - R

        if doSym:
            X = (X + X.T) / 2.0
        if corr:
            np.fill_diagonal(X, 1.0)
        elif keepDiag:
            np.fill_diagonal(X, diagX0)

        conv = _norm_inf(Y - X) / _norm_inf(Y)
        it += 1
        converged = conv <= conv_tol

    if not converged:
        warnings.warn(
            f"'nearPD()' did not converge in {it} iterations", stacklevel=2
        )

    if do2eigen:
        d, Q = _eigh_desc(X)
        eps = posd_tol * abs(d[0])
        if d[n - 1] < eps:
            d = np.where(d < eps, eps, d)
            o_diag = np.diag(X).copy()
            X = _symmetrize((Q * d) @ Q.T)
            D = np.sqrt(np.maximum(eps, o_diag) / np.diag(X))
            X = D[:, None] * X * D[None, :]
        if corr:
            np.fill_diagonal(X, 1.0)
        elif keepDiag:
            np.fill_diagonal(X, diagX0)

    return X
