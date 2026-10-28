"""Latent-correlation estimation (port of ``latentcor::latentcor``).

SPRING only ever calls ``latentcor(..., types = "tru")`` (via ``hugeKmb``), so
every variable is truncated/zero-inflated and every off-diagonal pair is
``comb = "22"``. Only that path is implemented; any other type raises, so a
silent divergence cannot happen.

R source: ``latentcor/R/latentcor.R`` + ``latentcor/R/internal.R``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bridge import r_switch22
from .kendall import kendall_tau_a_lower
from .nearpd import near_pd

__all__ = ["LatentCorResult", "latentcor"]


@dataclass
class LatentCorResult:
    """Result of :func:`latentcor` (mirrors the R list fields used by SPRING)."""

    K: np.ndarray
    R: np.ndarray
    Rpointwise: np.ndarray
    zratios: np.ndarray


def _check_types(types) -> None:
    if isinstance(types, (list, tuple, np.ndarray)):
        vals = set(np.asarray(types).tolist())
    else:
        vals = {types}
    if vals != {"tru"}:
        raise NotImplementedError(
            "This port implements only types='tru' (truncated), which is the "
            f"only type SPRING uses; got types={types!r}."
        )


def latentcor(
    X,
    types="tru",
    method: str = "approx",
    use_nearPD: bool = True,
    nu: float = 0.001,
    tol: float = 1e-6,
    ratio: float = 0.9,
) -> LatentCorResult:
    """Estimate the latent correlation matrix under a truncated Gaussian copula.

    Parameters
    ----------
    X : array-like of shape (n, p)
        Data matrix.
    types : {"tru"}, default "tru"
        Only ``"tru"`` is supported (SPRING's only mode).
    method : {"approx", "original"}, default "approx"
        ``"approx"`` uses the precomputed multilinear interpolation table;
        ``"original"`` inverts the bridge function numerically.
    use_nearPD : bool, default True
        Project the pointwise estimate onto the nearest correlation matrix.
    nu : float, default 0.001
        Shrinkage towards identity after ``nearPD``.
    tol : float, default 1e-6
        Accuracy of the bridge inversion when ``method="original"``.
    ratio : float, default 0.9
        Interpolation/optimisation boundary when ``method="approx"``.
    """
    if not (0.0 <= nu <= 1.0):
        raise ValueError("nu must be between 0 and 1.")
    if tol <= 0:
        raise ValueError("tol for optimization should be a positive value.")
    if not (0.0 <= ratio <= 1.0):
        raise ValueError("ratio must be between 0 and 1.")
    _check_types(types)

    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError("Input data matrix should be a 2-D numeric matrix.")
    n, p = X.shape

    # encodeX for type "tru": X <- X - colmin  (negate is all ones -> no sign flip)
    colmin = X.min(axis=0)
    Xe = X - colmin[None, :]

    # zratios for "tru": proportion of zeros per column (on the encoded matrix)
    zratio = (Xe == 0).mean(axis=0)

    rows, cols, k_lower = kendall_tau_a_lower(Xe)

    K = np.zeros((p, p))
    K[rows, cols] = k_lower
    K = K + K.T + np.eye(p)

    # all pairs are comb "22" -> r_switch on the lower triangle
    z1 = zratio[rows]
    z2 = zratio[cols]
    r_lower = r_switch22(k_lower, z1, z2, method=method, tol=tol, ratio=ratio)

    R = np.zeros((p, p))
    R[rows, cols] = r_lower
    R = R + R.T + np.eye(p)

    Rpointwise = R.copy()
    if use_nearPD:
        R = near_pd(R, corr=True)
        R = (1.0 - nu) * R + nu * np.eye(p)

    return LatentCorResult(K=K, R=R, Rpointwise=Rpointwise, zratios=zratio)
