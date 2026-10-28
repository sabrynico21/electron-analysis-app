"""``hugeKmb`` — rank-based MB graph estimation (port of ``R/helpers.R:44``)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .huge import huge_mb, symmetrize_path
from .latentcor import latentcor

__all__ = ["HugeKmbResult", "huge_kmb"]


@dataclass
class HugeKmbResult:
    """Mirrors the list returned by ``hugeKmb``."""

    beta: list[np.ndarray]
    path: list[np.ndarray]
    df: np.ndarray
    sparsity: np.ndarray
    lambda_: np.ndarray

    def __getitem__(self, key: str):
        """R-style ``fit$path`` / ``fit$beta`` access."""
        if key == "lambda":
            return self.lambda_
        return getattr(self, key)


def huge_kmb(
    data,
    lambda_=None,
    type: str = "tru",
    sym: str = "or",
    verbose: bool = True,
    Rmethod: str = "approx",
    tol: float = 1e-6,
    use_nearPD: bool = True,
    nu: float = 0.001,
    ratio: float = 0.9,
    **kwargs,
) -> HugeKmbResult:
    """Estimate the latent correlation and run MB neighborhood selection.

    The R argument is named ``lambda``; since that is a Python keyword it is
    accepted positionally or as the keyword ``lambda`` via ``**kwargs``.

    R source: ``SPRING/R/helpers.R:44``.
    """
    if lambda_ is None:
        lambda_ = kwargs.pop("lambda", None)
    if kwargs:
        raise TypeError(f"unused arguments: {', '.join(sorted(kwargs))}")
    if lambda_ is None:
        raise TypeError("huge_kmb() missing required argument: 'lambda'")
    data = np.asarray(data, dtype=float)
    lam = np.asarray(lambda_, dtype=float)

    lc = latentcor(
        data,
        types=type,
        method=Rmethod,
        use_nearPD=use_nearPD,
        nu=nu,
        tol=tol,
        ratio=ratio,
    )
    S = lc.R
    beta = huge_mb(S, lam)
    path = symmetrize_path(beta, sym)

    d = S.shape[0]
    df = np.asarray((np.abs(beta) != 0).sum(axis=1)).T  # (d, nlambda)
    if d <= 1:
        sparsity = np.zeros(lam.size)
    else:
        sparsity = path.sum(axis=(1, 2)) / d / (d - 1)

    return HugeKmbResult(
        beta=[beta[i] for i in range(lam.size)],
        path=[path[i] for i in range(lam.size)],
        df=df,
        sparsity=sparsity,
        lambda_=lam,
    )
