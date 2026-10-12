from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from spieceasi.pipeline import spiec_easi_mb


def test_raises_when_single_taxon() -> None:
    counts = pd.DataFrame({"taxon_a": [1, 2, 3, 4, 5]})
    with pytest.raises(ValueError, match="At least two taxa"):
        spiec_easi_mb(counts)


def test_raises_when_taxon_all_zero() -> None:
    counts = pd.DataFrame(
        {
            "a": [10, 9, 11, 10, 12, 8],
            "b": [0, 0, 0, 0, 0, 0],
            "c": [2, 3, 1, 4, 2, 1],
        }
    )
    with pytest.raises(ValueError, match="all-zero"):
        spiec_easi_mb(counts)


def test_raises_when_samples_not_greater_than_taxa() -> None:
    X = np.array(
        [
            [10, 5, 3, 2],
            [9, 4, 2, 2],
            [11, 6, 4, 1],
            [8, 3, 2, 2],
        ],
        dtype=float,
    )
    with pytest.raises(ValueError, match="sample count must exceed"):
        spiec_easi_mb(X)
