"""Input validation tests."""

from __future__ import annotations

import numpy as np
import pytest

from cclasso import CCLassoInputError, cclasso


def _good(n=20, p=5, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((n, p))
    f = np.exp(z) / np.exp(z).sum(axis=1, keepdims=True)
    return f


def test_non_numeric_matrix():
    with pytest.raises(CCLassoInputError):
        cclasso([["a", "b"], ["c", "d"]])


def test_wrong_ndim():
    with pytest.raises(CCLassoInputError):
        cclasso(np.ones(5))


def test_too_small():
    with pytest.raises(CCLassoInputError):
        cclasso(np.ones((5, 1)))


def test_nan_rejected():
    x = _good()
    x[0, 0] = np.nan
    with pytest.raises(CCLassoInputError):
        cclasso(x)


def test_inf_rejected():
    x = _good()
    x[0, 0] = np.inf
    with pytest.raises(CCLassoInputError):
        cclasso(x)


def test_negative_counts_rejected():
    x = _good() * 100
    x[0, 0] = -3
    with pytest.raises(CCLassoInputError):
        cclasso(x, counts=True)


def test_zero_rows_with_zero_pseudo_rejected():
    x = np.zeros((20, 5))
    with pytest.raises(CCLassoInputError):
        cclasso(x, counts=True, pseudo=0.0)


def test_nonpositive_without_counts_rejected():
    x = _good()
    x[0, 0] = 0.0
    with pytest.raises(CCLassoInputError):
        cclasso(x, counts=False)


def test_bad_lam_int_order():
    with pytest.raises(CCLassoInputError):
        cclasso(_good(), lam_int=(1.0, 1e-4))


def test_bad_lam_int_length():
    with pytest.raises(CCLassoInputError):
        cclasso(_good(), lam_int=(1e-4, 1e-2, 1.0))


def test_kcv_too_large():
    with pytest.raises(CCLassoInputError):
        cclasso(_good(n=8, p=5), k_cv=8)


def test_n_boot_invalid():
    with pytest.raises(CCLassoInputError):
        cclasso(_good(), n_boot=0)
