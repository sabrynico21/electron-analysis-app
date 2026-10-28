"""Input validation for SPRING.

The R implementation only rejects negative values. This module adds clearer,
*stricter* checks requested for the Python/Electron interface. They are applied
by :func:`spring_microbiome.spring.spring` when ``validate_input=True`` (the
default); set it to ``False`` to reproduce R's permissive behaviour exactly.

R source: ``SPRING/R/SPRING.R:98`` (only the ``any(data < 0)`` check exists).
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = ["validate_counts"]


def validate_counts(data, *, allow_zero_rows: bool = False) -> np.ndarray:
    """Validate a count/relative-abundance matrix and return it as float64.

    Raises
    ------
    TypeError
        If ``data`` is not numeric.
    ValueError
        For negative values, non-finite values, all-zero columns, all-zero
        rows (unless ``allow_zero_rows``), or fewer than two samples.
    """
    arr = np.asarray(data)
    if arr.ndim != 2:
        raise ValueError("data must be a 2-D matrix (samples x taxa).")
    if not np.issubdtype(arr.dtype, np.number):
        raise TypeError("data must be numeric.")
    arr = arr.astype(float, copy=False)

    if np.any(~np.isfinite(arr)):
        raise ValueError("data contains NaN or infinite values.")
    if np.any(arr < 0):
        raise ValueError(
            "Negative values are detected, but either quantitative or "
            "compositional counts are expected."
        )
    n, p = arr.shape
    if n < 2:
        raise ValueError(f"at least 2 samples are required, got {n}.")
    if p < 2:
        raise ValueError(f"at least 2 taxa are required, got {p}.")

    zero_cols = np.flatnonzero(arr.sum(axis=0) == 0)
    if zero_cols.size:
        raise ValueError(
            f"columns {zero_cols.tolist()} are all zero; remove them before running SPRING."
        )
    if not allow_zero_rows:
        zero_rows = np.flatnonzero(arr.sum(axis=1) == 0)
        if zero_rows.size:
            raise ValueError(
                f"rows {zero_rows.tolist()} have zero total counts; remove them "
                "before running SPRING."
            )

    if p >= n:
        warnings.warn(
            f"number of taxa (p={p}) is >= number of samples (n={n}); "
            "the estimate may be unstable.",
            stacklevel=2,
        )
    return arr
