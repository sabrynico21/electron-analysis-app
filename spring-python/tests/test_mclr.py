"""Unit tests for the mclr transformation (R: helpers.R:71)."""

from __future__ import annotations

import numpy as np
import pytest

from spring_microbiome.mclr import mclr


def test_preserves_zeros():
    dat = np.array([[0.0, 1.0, 0.0], [2.0, 3.0, 0.0]])
    out = mclr(dat)
    assert np.all(out[dat == 0] == 0.0)


def test_default_positive_shift():
    dat = np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    out = mclr(dat)
    assert np.all(out[dat > 0] > 0.0)


def test_eps_zero_is_unshifted_clr():
    dat = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    out = mclr(dat, eps=0)
    logd = np.log(dat)
    clr = logd - logd.mean(axis=1, keepdims=True)
    np.testing.assert_allclose(out, clr, rtol=1e-15, atol=1e-15)


def test_eps_positive_shift():
    dat = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    out0 = mclr(dat, eps=0)
    out2 = mclr(dat, eps=2)
    np.testing.assert_allclose(out2, out0 + 2.0, rtol=1e-15)


def test_single_row():
    dat = np.array([[0.0, 1.0, 2.0, 3.0]])
    out = mclr(dat)
    assert out.shape == dat.shape
    assert np.all(out[dat == 0] == 0)
    assert np.all(out[dat > 0] > 0)


def test_geometric_mean_of_nonzeros_is_removed():
    dat = np.array([[1.0, 0.0, 100.0]])
    out = mclr(dat, eps=0)
    # nonzero entries: log(1)-mean(log([1,100])) and log(100)-mean(log([1,100]))
    m = np.log(np.array([1.0, 100.0])).mean()
    np.testing.assert_allclose(out[0, 0], -m, rtol=1e-15)
    np.testing.assert_allclose(out[0, 2], np.log(100.0) - m, rtol=1e-15)


def test_negative_eps_raises():
    with pytest.raises(ValueError):
        mclr(np.array([[1.0, 2.0]]), eps=-1)


def test_negative_atleast_warns_and_uses_default():
    with pytest.warns(UserWarning, match="atleast should be positive"):
        mclr(np.array([[1.0, 2.0], [3.0, 4.0]]), atleast=-1)


def test_accepts_list_and_returns_2d():
    out = mclr([[1.0, 2.0], [3.0, 4.0]])
    assert out.shape == (2, 2)


def test_wrong_dimension_raises():
    with pytest.raises(ValueError):
        mclr(np.array([1.0, 2.0, 3.0]))
