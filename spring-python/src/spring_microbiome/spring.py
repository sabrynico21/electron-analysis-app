"""SPRING main entry point (port of ``SPRING.default``, ``R/SPRING.R:98``).

R argument names contain dots (``lambda.min.ratio``); Python identifiers cannot,
so they are exposed in snake_case. The default values are unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .huge_kmb import HugeKmbResult, huge_kmb
from .latentcor import latentcor
from .mclr import mclr
from .pulsar import PulsarResult, RefitResult, pulsar
from .pulsar import refit as _refit
from .validation import validate_counts

__all__ = ["SpringResult", "spring", "default_lambdaseq"]


@dataclass
class SpringResult:
    """Return value of :func:`spring` (mirrors the R list)."""

    output: PulsarResult
    fit: RefitResult
    lambdaseq: np.ndarray

    @property
    def adjacency(self) -> np.ndarray:
        """Final adjacency matrix at the StARS-selected lambda."""
        return self.fit.refit["stars"]

    @property
    def beta(self) -> HugeKmbResult:
        """Full-data MB fit (``fit$est``)."""
        return self.fit.est


def default_lambdaseq(lambda_min_ratio: float, nlambda: int) -> np.ndarray:
    """R default: ``exp(seq(log(0.6), log(0.6 * lambda.min.ratio), length.out = nlambda))``."""
    return np.exp(
        np.linspace(np.log(0.6), np.log(0.6 * lambda_min_ratio), num=nlambda)
    )


def spring(
    data,
    quantitative: bool = False,
    method: str = "mb",
    lambda_min_ratio: float = 1e-2,
    nlambda: int = 20,
    lambdaseq=None,
    seed: int = 10010,
    ncores: int = 1,
    thresh: float = 0.1,
    subsample_ratio: float = 0.8,
    rep_num: int = 20,
    Rtol: float = 1e-6,
    verbose: bool = True,
    Rmethod: str = "approx",
    use_nearPD: bool = True,
    nu: float = 0.001,
    ratio: float = 0.9,
    *,
    subsample_indices=None,
    rng: np.random.Generator | None = None,
    validate_input: bool = True,
) -> SpringResult:
    """Estimate a sparse microbial association network with SPRING.

    Parameters mirror the R function of the same name. In addition:

    subsample_indices : sequence of int arrays, optional
        The StARS subsample indices. Providing them removes all randomness
        (used by the equivalence tests to reuse R's exact indices).
    rng : numpy.random.Generator, optional
        Used to draw subsample indices when ``subsample_indices`` is None.
    validate_input : bool, default True
        Apply the stricter input validation in :mod:`spring_microbiome.validation`.
    """
    if method != "mb":
        raise ValueError('method must be "mb".')
    if lambdaseq is not None and not (
        isinstance(lambdaseq, str) or np.ndim(lambdaseq) == 1
    ):
        raise ValueError("lambdaseq must be None, a sequence, or 'data-specific'.")

    if validate_input:
        data = validate_counts(data)
    else:
        data = np.asarray(data, dtype=float)
        if np.any(data < 0):
            raise ValueError(
                "Negative values are detected, but either quantitative or "
                "compositional counts are expected.\n"
            )

    p = data.shape[1]
    if quantitative:
        if np.max(data.sum(axis=1)) <= 1 or np.isclose(
            np.max(data.sum(axis=1)), 1.0
        ):
            import warnings

            warnings.warn(
                "The input data is normalized, but quantitative count data is "
                "expected.\n",
                stacklevel=2,
            )
        qdat = data
    else:
        qdat = mclr(data)

    if isinstance(lambdaseq, str):
        if lambdaseq != "data-specific":
            raise ValueError("The input for lambdaseq is not correct.")
        Kcor = latentcor(
            qdat,
            types="tru",
            method=Rmethod,
            use_nearPD=use_nearPD,
            nu=nu,
            tol=Rtol,
            ratio=ratio,
        ).R
        off = Kcor - np.eye(p)
        lambda_max = max(off.max(), -off.min())
        lambda_min = lambda_min_ratio * lambda_max
        lam = np.exp(np.linspace(np.log(lambda_max), np.log(lambda_min), num=nlambda))
    elif lambdaseq is None:
        lam = default_lambdaseq(lambda_min_ratio, nlambda)
    else:
        lam = np.asarray(lambdaseq, dtype=float)

    fargs = {
        "lambda": lam,
        "Rmethod": Rmethod,
        "tol": Rtol,
        "verbose": verbose,
        "use_nearPD": use_nearPD,
        "nu": nu,
        "ratio": ratio,
    }

    out = pulsar(
        qdat,
        fun=huge_kmb,
        fargs=fargs,
        criterion="stars",
        thresh=thresh,
        subsample_ratio=subsample_ratio,
        rep_num=rep_num,
        seed=seed,
        ncores=ncores,
        refit=True,
        subsample_indices=subsample_indices,
        rng=rng,
    )
    fit = _refit(out)
    return SpringResult(output=out, fit=fit, lambdaseq=lam)
