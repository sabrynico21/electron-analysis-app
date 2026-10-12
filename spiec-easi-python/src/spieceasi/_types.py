"""Internal type aliases."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd

    ArrayLike2D = np.ndarray | pd.DataFrame
else:  # pragma: no cover - runtime value is unused (annotations are strings)
    # ``pandas`` is an optional convenience dependency: the runtime fallback
    # keeps the package importable without it.
    ArrayLike2D = np.ndarray
