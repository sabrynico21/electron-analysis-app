"""Modified central log-ratio (mclr) transformation.

Faithful port of ``mclr`` in SPRING ``R/helpers.R:71``.
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = ["mclr"]


def mclr(
    dat,
    base: float | None = None,
    tol: float = 1e-16,
    eps: float | None = None,
    atleast: float = 1.0,
) -> np.ndarray:
    """Apply the modified central log-ratio transformation.

    Parameters
    ----------
    dat : array-like of shape (n, p)
        Raw count or compositional data. Non-negative.
    base : float, optional
        Logarithm base. Defaults to ``e`` (natural log), matching R.
    tol : float, default 1e-16
        A value is considered zero when ``dat < tol`` (R uses ``dat >= tol``).
    eps : float, optional
        If ``None`` (default) an automatic positive shift makes every non-zero
        entry strictly positive. ``eps = 0`` returns the unshifted CLR (may be
        negative). ``eps > 0`` adds a user-defined shift. ``eps < 0`` raises.
    atleast : float, default 1.0
        Constant added to ``abs(min(clrdat))`` when the automatic shift is used.

    Returns
    -------
    numpy.ndarray
        Matrix with the same shape as ``dat``; original zeros are exactly 0.

    Notes
    -----
    R source: ``SPRING/R/helpers.R:71``.
    """
    base = float(np.e) if base is None else float(base)
    x = np.asarray(dat, dtype=float)
    if x.ndim != 2:
        raise ValueError("mclr expects a 2-D array (n samples x p taxa).")

    nzero = x >= tol  # R: nzero <- (dat >= tol)
    with np.errstate(divide="ignore", invalid="ignore"):
        logdat = np.log(x) / np.log(base)  # R: log(dat, base)
    logdat = np.where(nzero, logdat, 0.0)  # R: ifelse(nzero, log(dat, base), 0)

    n_row = x.shape[0]
    if n_row > 1:
        # R: clrdat <- ifelse(nzero, LOG - rowMeans(LOG)/rowMeans(nzero), 0)
        row_mean_log = logdat.sum(axis=1) / x.shape[1]
        row_mean_nz = nzero.sum(axis=1) / x.shape[1]
        with np.errstate(divide="ignore", invalid="ignore"):
            const = row_mean_log / row_mean_nz
        clrdat = np.where(nzero, logdat - const[:, None], 0.0)
    elif n_row == 1:
        # R: clrdat <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            const = logdat.mean() / nzero.mean()
        clrdat = np.where(nzero, logdat - const, 0.0)
    else:
        raise ValueError("mclr expects at least one sample (row).")

    if eps is None:
        if atleast < 0:
            warnings.warn(
                "atleast should be positive. The functions uses default value 1 instead.",
                stacklevel=2,
            )
            atleast = 1.0
        if clrdat.min() < 0:
            positivecst = abs(clrdat.min()) + atleast
        else:
            positivecst = 0.0
        return np.where(nzero, clrdat + positivecst, 0.0)

    if eps == 0:
        # no shift. clr transform applied to non-zero proportions only.
        return clrdat
    if eps > 0:
        return np.where(nzero, clrdat + eps, 0.0)
    raise ValueError(
        "check your eps value for additional positive shift. Otherwise, leave it as NULL."
    )
