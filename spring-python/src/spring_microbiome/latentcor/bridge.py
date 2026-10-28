"""Bridge functions and inverse-bridge inversion for truncated variables.

Faithful port of the parts of ``latentcor/R/internal.R`` that SPRING can
actually reach. SPRING always calls ``latentcor(..., types = "tru")`` for every
column, so the only type combination that ever occurs is ``"22"``
(truncated/truncated). ``r_sol`` and ``bound_switch`` for ``"22"`` are ported
verbatim; the other combinations are intentionally *not* implemented and raise a
clear error, so a silent divergence is impossible.

R source: ``latentcor/R/internal.R`` (``r_sol`` -> "22" branch,
``bound_switch`` -> "22" branch, ``r_switch``, ``r_ml``).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import multivariate_normal, norm

from .interp import load_ipol22

__all__ = ["bound22", "r_sol22", "r_ml22", "r_switch22"]

_SQRT2 = float(np.sqrt(2.0))


def brent_fmin(ax: float, bx: float, f, tol: float) -> float:
    """Exact port of R's ``Brent_fmin`` (``stats/src/optimize.c``).

    This is the algorithm behind ``stats::optimize(f, lower, upper, tol)``, which
    ``latentcor`` uses to invert the bridge function. Using SciPy's Brent instead
    gives materially different minimizers on this flat objective (~1e-5), so the
    R algorithm is reproduced verbatim.
    """
    c = (3.0 - math.sqrt(5.0)) * 0.5
    eps = math.sqrt(np.finfo(float).eps)

    a, b = ax, bx
    v = a + c * (b - a)
    w = v
    x = v
    d = 0.0
    e = 0.0
    fx = f(x)
    fv = fx
    fw = fx
    tol3 = tol / 3.0

    while True:
        xm = (a + b) * 0.5
        tol1 = eps * abs(x) + tol3
        t2 = tol1 * 2.0

        if abs(x - xm) <= t2 - (b - a) * 0.5:
            break

        p = 0.0
        q = 0.0
        r = 0.0
        if abs(e) > tol1:  # fit parabola
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0
            if q > 0.0:
                p = -p
            else:
                q = -q
            r = e
            e = d

        if abs(p) >= abs(q * 0.5 * r) or p <= q * (a - x) or p >= q * (b - x):
            # golden-section step
            e = (b - x) if x < xm else (a - x)
            d = c * e
        else:
            # parabolic-interpolation step
            d = p / q
            u = x + d
            if u - a < t2 or b - u < t2:  # not too close to ax or bx
                d = tol1
                if x >= xm:
                    d = -d

        if abs(d) >= tol1:  # not too close to x
            u = x + d
        elif d > 0.0:
            u = x + tol1
        else:
            u = x - tol1

        fu = f(u)

        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu

    return x


def _pmnorm4(x: np.ndarray, cov: np.ndarray) -> float:
    """``mnormt::pmnorm(x, mean=0, varcov=cov)`` for 4 dimensions.

    ``scipy.stats.multivariate_normal.cdf`` uses randomized quasi-Monte-Carlo
    integration and is therefore NON-DETERMINISTIC by default (its result
    changes from call to call, which would make the whole pipeline depend on
    unrelated call history). We pass a fresh, fixed-seed generator on *every*
    call and tighten the tolerances, which makes each evaluation reproducible
    to ~1e-7 independently of ordering.
    """
    return float(
        multivariate_normal.cdf(
            x,
            mean=np.zeros(4),
            cov=cov,
            rng=np.random.default_rng(0),
            maxpts=1_000_000,
            abseps=1e-10,
            releps=1e-10,
        )
    )


def bridge22(r: float, z1: float, z2: float) -> float:
    """Bridge function F(r) for truncated/truncated (R ``r_sol`` comb "22")."""
    de1 = norm.ppf(z1)
    de2 = norm.ppf(z2)
    s = 1.0 / _SQRT2
    mat1 = np.array(
        [
            [1.0, 0.0, s, -r * s],
            [0.0, 1.0, -r * s, s],
            [s, -r * s, 1.0, -r],
            [-r * s, s, -r, 1.0],
        ]
    )
    mat2 = np.array(
        [
            [1.0, r, s, r * s],
            [r, 1.0, r * s, s],
            [s, r * s, 1.0, r],
            [r * s, s, r, 1.0],
        ]
    )
    x = np.array([-de1, -de2, 0.0, 0.0])
    return -2.0 * _pmnorm4(x, mat1) + 2.0 * _pmnorm4(x, mat2)


def bound22(z1, z2):
    """R ``bound_switch`` for comb "22": ``1 - pmax(z1, z2)**2``."""
    z1 = np.asarray(z1, dtype=float)
    z2 = np.asarray(z2, dtype=float)
    return 1.0 - np.maximum(z1, z2) ** 2


def r_sol22(K, z1, z2, tol: float = 1e-6):
    """Numerically invert the bridge function (R ``r_sol``, comb "22").

    Minimises ``(F(r) - K)**2`` over ``r in (-0.999, 0.999)`` using R's
    ``stats::optimize`` (see :func:`brent_fmin`).
    """
    K = np.atleast_1d(np.asarray(K, dtype=float))
    z1 = np.atleast_1d(np.asarray(z1, dtype=float))
    z2 = np.atleast_1d(np.asarray(z2, dtype=float))
    out = np.empty(K.size, dtype=float)
    for i in range(K.size):
        k = K[i]
        a = z1[i]
        b = z2[i]

        def f(r, a=a, b=b, k=k):
            return (bridge22(r, a, b) - k) ** 2

        out[i] = brent_fmin(-0.999, 0.999, f, tol)
    return out


def r_ml22(K, z1, z2):
    """Multilinear approximation (R ``r_ml`` / ``ipol_22``), divided by 1e7."""
    ipol = load_ipol22()
    K = np.atleast_1d(np.asarray(K, dtype=float))
    z1 = np.atleast_1d(np.asarray(z1, dtype=float))
    z2 = np.atleast_1d(np.asarray(z2, dtype=float))
    pts = np.vstack([K, z1, z2])
    return ipol(pts) / getattr(ipol, "scale", 1e7)


def r_switch22(K, z1, z2, method: str = "approx", tol: float = 1e-6, ratio: float = 0.9):
    """R ``r_switch`` for comb "22"."""
    if method == "original":
        return r_sol22(K, z1, z2, tol)
    if method != "approx":
        raise ValueError("method must be 'approx' or 'original'")

    K = np.atleast_1d(np.asarray(K, dtype=float))
    z1 = np.atleast_1d(np.asarray(z1, dtype=float))
    z2 = np.atleast_1d(np.asarray(z2, dtype=float))
    bound = bound22(z1, z2)
    cutoff = np.abs(K) > ratio * bound
    out = np.empty(K.size, dtype=float)
    if np.all(cutoff):
        return r_sol22(K, z1, z2, tol)
    if not np.any(cutoff):
        return r_ml22(K / bound, z1, z2)
    out[cutoff] = r_sol22(K[cutoff], z1[cutoff], z2[cutoff], tol)
    out[~cutoff] = r_ml22(K[~cutoff] / bound[~cutoff], z1[~cutoff], z2[~cutoff])
    return out
