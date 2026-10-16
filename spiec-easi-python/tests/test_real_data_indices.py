"""Unit tests for the lambda/sample index conversions and Wilson interval.

These pin the pre-registered rule: R lambda indices are 1-based, Python lambda
indices are 0-based, and the conversion is explicit and total (rejects invalid
inputs). No R is involved and no golden data is touched.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.real_data.reference_indices import (
    python_lambda_index_to_r,
    r_lambda_index_to_python,
    r_sample_indices_to_python,
    wilson_interval,
)


def test_r_to_python_lambda_index_roundtrip() -> None:
    for r_index in range(1, 21):
        py = r_lambda_index_to_python(r_index)
        assert py == r_index - 1
        assert python_lambda_index_to_r(py) == r_index


def test_python_to_r_lambda_index_roundtrip() -> None:
    for py_index in range(20):
        r = python_lambda_index_to_r(py_index)
        assert r == py_index + 1
        assert r_lambda_index_to_python(r) == py_index


def test_lambda_index_zero_based_boundary() -> None:
    # Python index 0 (lambda.max) is R index 1; this is the documented boundary.
    assert r_lambda_index_to_python(1) == 0
    assert python_lambda_index_to_r(0) == 1


@pytest.mark.parametrize("bad", [0, -1, -5])
def test_r_lambda_index_rejects_non_positive(bad: int) -> None:
    with pytest.raises(ValueError):
        r_lambda_index_to_python(bad)


@pytest.mark.parametrize("bad", [-1, -3])
def test_python_lambda_index_rejects_negative(bad: int) -> None:
    with pytest.raises(ValueError):
        python_lambda_index_to_r(bad)


def test_r_sample_indices_to_python() -> None:
    out = r_sample_indices_to_python([1, 5, 96])
    assert np.array_equal(out, np.array([0, 4, 95]))
    with pytest.raises(ValueError):
        r_sample_indices_to_python([0, 1, 2])


def test_wilson_interval_properties() -> None:
    lo, hi = wilson_interval(10, 10)
    assert 0.0 <= lo < hi <= 1.0
    assert lo > 0.5  # 10/10 is high but not degenerate
    # Point estimate lies within the interval.
    lo2, hi2 = wilson_interval(3, 10)
    assert lo2 < 0.3 < hi2
    # Symmetry: 0 successes and n successes mirror.
    a = wilson_interval(0, 10)
    b = wilson_interval(10, 10)
    assert abs(a[1] - (1.0 - b[0])) < 1e-12
    nan_lo, nan_hi = wilson_interval(0, 0)
    assert np.isnan(nan_lo) and np.isnan(nan_hi)
