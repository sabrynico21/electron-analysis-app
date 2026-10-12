"""Validation utilities for SPIEC-EASI MB inputs."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def validate_count_matrix(
    X: np.ndarray,
    *,
    min_samples: int = 1,
    require_samples_gt_taxa: bool = False,
) -> None:
    """Validate count matrix for downstream MB estimation.

    Parameters
    ----------
    X
        Input count matrix of shape ``(n_samples, n_taxa)``.
    min_samples
        Minimum required sample count.

    Raises
    ------
    ValueError
        If any constraints are violated.
    """
    if X.ndim != 2:
        msg = "Input data must be a 2D matrix with shape (samples, taxa)."
        raise ValueError(msg)

    n_samples, n_taxa = X.shape
    if n_taxa < 2:
        msg = "At least two taxa (columns) are required for network inference."
        raise ValueError(msg)

    if n_samples < min_samples:
        msg = f"At least {min_samples} samples are required, received {n_samples}."
        raise ValueError(msg)

    if require_samples_gt_taxa and n_samples <= n_taxa:
        msg = (
            "The sample count must exceed the number of taxa for stable MB "
            f"regressions (n_samples={n_samples}, n_taxa={n_taxa})."
        )
        raise ValueError(msg)

    if np.any(X < 0):
        msg = "Count data must be non-negative."
        raise ValueError(msg)

    if np.any(np.all(X == 0, axis=0)):
        msg = "At least one taxon has all-zero counts across all samples."
        raise ValueError(msg)

    if np.any(np.sum(X, axis=1) == 0):
        msg = "At least one sample has all-zero counts across all taxa."
        raise ValueError(msg)


def validate_taxa_names(taxa_names: Sequence[str], n_taxa: int) -> None:
    """Validate that taxa names match the number of taxa."""
    if len(taxa_names) != n_taxa:
        msg = (
            "Taxa name count does not match matrix columns "
            f"(len={len(taxa_names)}, n_taxa={n_taxa})."
        )
        raise ValueError(msg)
