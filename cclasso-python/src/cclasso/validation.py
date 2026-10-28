"""Input validation for :func:`cclasso.cclasso`.

R silently produces ``NaN``/``Inf`` for malformed inputs; the port fails loudly with
clear messages instead (the numerical algorithm itself is unchanged).
"""

from __future__ import annotations

import numpy as np

__all__ = ["validate_input", "CCLassoInputError"]


class CCLassoInputError(ValueError):
    """Raised for invalid input to :func:`cclasso.cclasso`."""


def validate_input(x, *, counts, pseudo, k_cv, lam_int, n_boot) -> np.ndarray:
    # --- matrix coercion ---
    try:
        x = np.asarray(x, dtype=float)
    except (TypeError, ValueError) as exc:  # non-numeric content
        raise CCLassoInputError("x must be a numeric matrix; could not convert to float.") from exc

    if x.ndim != 2:
        raise CCLassoInputError(f"x must be 2-dimensional (n x p); got ndim={x.ndim}.")
    n, p = x.shape
    if n < 2 or p < 2:
        raise CCLassoInputError(f"x must have at least 2 rows and 2 columns; got shape {x.shape}.")
    if not np.all(np.isfinite(x)):
        raise CCLassoInputError("x contains NaN or infinite values.")

    # --- counts / pseudo ---
    if pseudo < 0:
        raise CCLassoInputError(f"pseudo must be >= 0; got {pseudo}.")
    if counts:
        if np.any(x < 0):
            raise CCLassoInputError(
                "counts=True but x contains negative values; counts must be >= 0."
            )
        row_sums = x.sum(axis=1) + p * pseudo
        if np.any(row_sums <= 0):
            raise CCLassoInputError(
                "counts=True with pseudo=0 and one or more all-zero rows "
                "(row sum is 0); cannot normalise."
            )
    else:
        if np.any(x <= 0):
            raise CCLassoInputError(
                "counts=False but x contains non-positive values; log(x) is undefined. "
                "Use counts=True with a pseudo-count for count data containing zeros."
            )

    # --- cross-validation ---
    if int(k_cv) != k_cv or k_cv < 1:
        raise CCLassoInputError(f"k_cv must be a positive integer; got {k_cv}.")
    k_cv = int(k_cv)
    n_b = n // k_cv
    if n_b < 2 or (n - n_b) < 2:
        raise CCLassoInputError(
            f"k_cv={k_cv} is too large for n={n}: cross-validation needs at least "
            f"2 test rows and 2 training rows (floor(n/k_cv)={n_b})."
        )

    # --- tuning interval ---
    lam_int = np.asarray(lam_int, dtype=float).ravel()
    if lam_int.size != 2:
        raise CCLassoInputError(
            f"lam_int must be a length-2 interval; got {lam_int.size} value(s)."
        )
    if not np.all(np.isfinite(lam_int)) or np.any(lam_int <= 0):
        raise CCLassoInputError(f"lam_int values must be finite and > 0; got {lam_int}.")
    if lam_int[0] > lam_int[1]:
        raise CCLassoInputError(f"lam_int must be increasing (min, max); got {lam_int}.")

    # --- bootstrap ---
    if int(n_boot) != n_boot or n_boot < 1:
        raise CCLassoInputError(f"n_boot must be a positive integer; got {n_boot}.")

    return np.ascontiguousarray(x, dtype=float)
