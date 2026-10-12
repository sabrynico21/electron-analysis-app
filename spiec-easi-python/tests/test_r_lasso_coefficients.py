from __future__ import annotations

import numpy as np
import pytest

from spieceasi.neighborhood import _fit_single_neighborhood
from tests.validation_config import (
    LASSO_COEFFICIENT_TOLERANCE,
    MIN_LASSO_LAMBDAS_PER_TARGET,
    MIN_LASSO_TARGET_COMPARISONS,
    RELATIVE_ERROR_MIN_MAGNITUDE,
    SUPPORT_ZERO_THRESHOLD,
)


def _target_specific_alpha_max(X: np.ndarray, target_idx: int) -> float:
    n_samples, n_taxa = X.shape
    mask = np.ones(n_taxa, dtype=bool)
    mask[target_idx] = False
    X_other = X[:, mask]
    y = X[:, target_idx]

    X_center = X_other - X_other.mean(axis=0)
    X_scale = X_center.std(axis=0)
    X_scale[X_scale == 0.0] = 1.0
    X_std = X_center / X_scale
    y_center = y - y.mean()

    alpha_max = float(np.max(np.abs((X_std.T @ y_center) / n_samples)))
    return max(alpha_max, 1e-12)


def _relative_error_max(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    if not np.any(mask):
        return 0.0
    denom = np.maximum(np.abs(b[mask]), 1e-15)
    return float(np.max(np.abs(a[mask] - b[mask]) / denom))


def _make_gaussian_matrix(
    *,
    n_samples: int,
    n_taxa: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    latent = rng.normal(size=(n_samples, n_taxa))

    # Mild correlation structure while keeping a full-rank, non-compositional
    # design matrix suitable for strict coefficient-level solver parity.
    for j in range(1, n_taxa):
        latent[:, j] = 0.20 * latent[:, j - 1] + np.sqrt(1.0 - 0.20**2) * latent[:, j]
    return latent


@pytest.mark.r_validation
@pytest.mark.slow
def test_removed_sklearn_branch_raises(r_runtime) -> None:
    """The sklearn/glmnet-style coefficient branch was removed.

    R's ``method='mb'`` uses huge's coordinate descent, not glmnet, so the
    ``coefficient_scale='original'`` and ``standardize=False`` options now
    raise ``NotImplementedError``. ``r_runtime`` is used only so the test is
    skipped when the R stack is unavailable (parity is covered by
    ``test_r_components`` / ``test_r_fixed_lambda_network``).
    """
    X = _make_gaussian_matrix(n_samples=120, n_taxa=12, seed=5)
    with pytest.raises(NotImplementedError):
        _fit_single_neighborhood(
            X,
            0,
            alpha=0.05,
            standardize=True,
            max_iter=200_000,
            tol=1e-12,
            random_state=0,
            coefficient_scale="original",
        )
    with pytest.raises(NotImplementedError):
        _fit_single_neighborhood(
            X,
            0,
            alpha=0.05,
            standardize=False,
            max_iter=200_000,
            tol=1e-12,
            random_state=0,
            coefficient_scale="standardized",
        )
