"""Compositional data transformations used by SPIEC-EASI."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

try:  # pragma: no cover - pandas is an optional convenience dependency
    import pandas as pd
except ImportError:  # pragma: no cover - exercised only without pandas
    pd = None

from ._types import ArrayLike2D
from ._validation import validate_count_matrix, validate_taxa_names

ZeroMethod = Literal["additive", "multiplicative"]


@dataclass(frozen=True)
class CLRResult:
    """Result of CLR transformation.

    Attributes
    ----------
    clr
        CLR-transformed matrix with shape ``(n_samples, n_taxa)``.
    taxa_names
        Taxa names aligned to CLR columns.
    """

    clr: np.ndarray
    taxa_names: list[str]


def _to_numpy_and_names(data: ArrayLike2D) -> tuple[np.ndarray, list[str]]:
    if pd is not None and isinstance(data, pd.DataFrame):
        matrix = data.to_numpy(dtype=float, copy=True)
        taxa_names = [str(c) for c in data.columns]
    else:
        matrix = np.asarray(data, dtype=float)
        taxa_names = [f"taxon_{i}" for i in range(matrix.shape[1])]

    validate_count_matrix(
        matrix,
        min_samples=1,
        require_samples_gt_taxa=False,
    )
    validate_taxa_names(taxa_names, matrix.shape[1])
    return matrix, taxa_names


def multiplicative_zero_replacement(
    proportions: np.ndarray,
    *,
    delta: float,
) -> np.ndarray:
    """Perform multiplicative replacement of zeros in a compositional matrix.

    Parameters
    ----------
    proportions
        Row-normalized compositional matrix where each row sums to one.
    delta
        Replacement value for zero entries in proportion space.

    Returns
    -------
    np.ndarray
        Matrix with zeros replaced and rows re-normalized to sum to one.
    """
    X = np.asarray(proportions, dtype=float).copy()
    n_taxa = X.shape[1]
    max_allowed_delta = (1.0 / n_taxa) * 0.99
    if not 0.0 < delta < max_allowed_delta:
        msg = (
            "For multiplicative replacement, pseudocount (delta) must satisfy "
            f"0 < delta < {max_allowed_delta:.4g} for n_taxa={n_taxa}."
        )
        raise ValueError(msg)

    for i in range(X.shape[0]):
        row = X[i, :]
        zero_mask = row == 0
        z = int(np.sum(zero_mask))
        if z == 0:
            continue

        nonzero_mask = ~zero_mask
        nonzero_sum = float(np.sum(row[nonzero_mask]))
        if nonzero_sum <= 0:
            msg = "A sample row is entirely zero after normalization."
            raise ValueError(msg)

        total_zero_mass = z * delta
        if total_zero_mass >= 1.0:
            msg = "Replacement mass for zeros is too large; decrease pseudocount."
            raise ValueError(msg)

        row[zero_mask] = delta
        row[nonzero_mask] = row[nonzero_mask] * ((1.0 - total_zero_mass) / nonzero_sum)
        X[i, :] = row

    X /= np.sum(X, axis=1, keepdims=True)
    return X


def clr_transform(
    data: ArrayLike2D,
    *,
    pseudocount: float = 1.0,
    zero_method: ZeroMethod = "additive",
) -> CLRResult:
    """Apply centered log-ratio transform to compositional count data.

    Parameters
    ----------
    data
        Count matrix with shape ``(n_samples, n_taxa)``.
    pseudocount
        Zero-handling parameter. With ``zero_method='additive'``, this value is
        added directly to each count before row normalization. With
        ``zero_method='multiplicative'``, this acts as replacement mass in
        compositional space.
    zero_method
        Strategy for handling zeros before the logarithm:
        ``'additive'`` or ``'multiplicative'``.

    Returns
    -------
    CLRResult
        Dataclass containing the CLR matrix and taxa names.
    """
    X, taxa_names = _to_numpy_and_names(data)
    if pseudocount <= 0:
        msg = "pseudocount must be strictly positive."
        raise ValueError(msg)

    if zero_method == "additive":
        adjusted = X + pseudocount
        proportions = adjusted / np.sum(adjusted, axis=1, keepdims=True)
    elif zero_method == "multiplicative":
        proportions = X / np.sum(X, axis=1, keepdims=True)
        proportions = multiplicative_zero_replacement(proportions, delta=pseudocount)
    else:
        msg = f"Unsupported zero_method '{zero_method}'."
        raise ValueError(msg)

    log_props = np.log(proportions)
    clr = log_props - np.mean(log_props, axis=1, keepdims=True)
    return CLRResult(clr=clr, taxa_names=taxa_names)


def coerce_taxa_names(taxa_names: Sequence[str] | None, n_taxa: int) -> list[str]:
    """Return validated taxa names with default fallback values."""
    if taxa_names is None:
        return [f"taxon_{i}" for i in range(n_taxa)]
    names = [str(name) for name in taxa_names]
    validate_taxa_names(names, n_taxa)
    return names
