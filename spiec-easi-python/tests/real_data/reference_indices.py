"""Shared helpers for the real-data R<->Python verification battery.

* ``r_lambda_index_to_python`` / ``python_lambda_index_to_r`` implement the
  explicit 1-based (R) <-> 0-based (Python) lambda index conversion. This is the
  single place where the conversion happens; a dedicated unit test
  (``tests/test_real_data_indices.py``) pins the behaviour.
* ``r_sample_indices_to_python`` converts R's 1-based sample indices to 0-based.
* ``wilson_interval`` returns the 95% Wilson score interval for a binomial
  proportion, used when reporting per-case frequencies.

No production code lives here; this is test/analysis scaffolding only.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np


def r_lambda_index_to_python(r_index: int) -> int:
    """Convert an R 1-based lambda index (``opt.index``) to a Python 0-based index."""
    if r_index < 1:
        msg = f"R lambda indices are 1-based (>=1); got {r_index}."
        raise ValueError(msg)
    return int(r_index) - 1


def python_lambda_index_to_r(py_index: int) -> int:
    """Convert a Python 0-based lambda index to an R 1-based index."""
    if py_index < 0:
        msg = f"Python lambda indices are 0-based (>=0); got {py_index}."
        raise ValueError(msg)
    return int(py_index) + 1


def r_sample_indices_to_python(indices: Sequence[int]) -> np.ndarray:
    """Convert R 1-based sample indices to 0-based (as used by the pipeline)."""
    arr = np.asarray(indices, dtype=int)
    if arr.size and int(arr.min()) < 1:
        msg = f"R sample indices are 1-based (>=1); got min {int(arr.min())}."
        raise ValueError(msg)
    return arr - 1


def wilson_interval(
    successes: int, trials: int, z: float = 1.959963984540054
) -> tuple[float, float]:
    """95% Wilson score interval for ``successes / trials`` (default z = 1.96)."""
    if trials <= 0:
        return (float("nan"), float("nan"))
    if not 0 <= successes <= trials:
        msg = f"successes must be in [0, trials]; got {successes}/{trials}."
        raise ValueError(msg)
    phat = successes / trials
    denom = 1.0 + z * z / trials
    center = (phat + z * z / (2.0 * trials)) / denom
    half = (z / denom) * math.sqrt(
        phat * (1.0 - phat) / trials + z * z / (4.0 * trials * trials)
    )
    return (max(0.0, center - half), min(1.0, center + half))
