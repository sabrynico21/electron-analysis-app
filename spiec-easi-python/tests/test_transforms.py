from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from spieceasi.transforms import clr_transform


def test_clr_additive_rows_sum_to_zero() -> None:
    data = pd.DataFrame([[10, 2, 0], [1, 4, 9]], columns=["a", "b", "c"])
    res = clr_transform(data, pseudocount=1.0, zero_method="additive")
    row_sums = np.sum(res.clr, axis=1)
    assert np.allclose(row_sums, 0.0, atol=1e-10)
    assert res.clr.shape == (2, 3)


def test_clr_multiplicative_no_infinite_values() -> None:
    X = np.array([[5, 0, 3], [1, 4, 8], [7, 2, 0]], dtype=float)
    res = clr_transform(X, pseudocount=0.01, zero_method="multiplicative")
    assert np.isfinite(res.clr).all()


def test_clr_raises_on_nonpositive_pseudocount() -> None:
    X = np.array([[1, 2], [3, 4]], dtype=float)
    with pytest.raises(ValueError, match="pseudocount"):
        clr_transform(X, pseudocount=0.0)
