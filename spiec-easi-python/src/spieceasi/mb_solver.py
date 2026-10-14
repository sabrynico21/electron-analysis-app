"""MB neighborhood selection — public entry point.

Two interchangeable implementations of the ``huge::mb`` coordinate descent:

* :mod:`spieceasi._mb_core` — compiled C++ (near-literal translation of
  ``huge_core.cpp`` lines 640-810). Used when importable. This is the default.
* :mod:`spieceasi.mb_solver_reference` — the pure-Python reference, always
  available and used as an automatic fallback (with a ``RuntimeWarning``).

Only the solver differs between the two. The correlation matrix, the lambda path
and everything downstream stay in Python and are shared, so switching
implementations cannot change any other stage of the pipeline.

See ``docs/AUDIT.md`` for the declared deviations and ``setup.py`` for the
compile flags.
"""

from __future__ import annotations

import warnings

import numpy as np

from . import mb_solver_reference as _reference

# Re-exported unchanged: these are shared by both implementations.
from .mb_solver_reference import (
    HUGE_MAX_ITER,
    HUGE_THOL,
    huge_correlation,
    soft_threshold,
    symmetrize_path,
)

__all__ = [
    "HUGE_MAX_ITER",
    "HUGE_THOL",
    "HAVE_CPP_SOLVER",
    "active_solver",
    "huge_correlation",
    "soft_threshold",
    "solve_mb_path_huge_scale",
    "symmetrize_path",
]

_IMPORT_ERROR: str | None = None
try:  # pragma: no cover - depends on whether the extension was built
    from . import _mb_core as _core
except Exception as _exc:  # noqa: BLE001 - any import failure means "no C++"
    _core = None
    _IMPORT_ERROR = f"{type(_exc).__name__}: {_exc}"

HAVE_CPP_SOLVER: bool = _core is not None
"""True when the compiled MB core is importable and will be used."""

_warned = False


def active_solver() -> str:
    """Name of the solver in use: ``"cpp"`` or ``"python"``."""
    return "cpp" if HAVE_CPP_SOLVER else "python"


def _warn_fallback() -> None:
    global _warned
    if _warned:
        return
    _warned = True
    warnings.warn(
        "spieceasi._mb_core is not available "
        f"({_IMPORT_ERROR}); falling back to the pure-Python MB solver, which is "
        "roughly two orders of magnitude slower. Build the extension with "
        "`pip install --no-build-isolation -e .` (needs pybind11 and a C++17 "
        "compiler).",
        RuntimeWarning,
        stacklevel=3,
    )


def solve_mb_path_huge_scale(
    X: np.ndarray,
    alpha_path: np.ndarray,
    *,
    max_iter: int = HUGE_MAX_ITER,
    tol: float = HUGE_THOL,
) -> np.ndarray:
    """Solve the MB nodewise Lasso path on the correlation scale, warm-started.

    Returns a ``(n_lambda, d, d)`` array where ``directed_path[i, m, k]`` is the
    coefficient of predictor ``k`` in the regression for response ``m`` at
    ``alpha_path[i]`` (the Python ``[node, predictor]`` convention; ``huge``
    stores the transpose).
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] < 2:
        msg = "At least two features are required for MB neighborhood selection."
        raise ValueError(msg)

    alpha_path = np.asarray(alpha_path, dtype=float)
    if alpha_path.ndim != 1 or alpha_path.size < 1:
        msg = "alpha_path must be a non-empty 1-D array."
        raise ValueError(msg)

    if not HAVE_CPP_SOLVER:
        _warn_fallback()
        return _reference.solve_mb_path_huge_scale(
            X, alpha_path, max_iter=max_iter, tol=tol
        )

    # Resolve thol / max_iter exactly as the reference does, so both
    # implementations receive identical values.
    max_cd_iter = min(int(max_iter), HUGE_MAX_ITER)
    if max_cd_iter < 1:
        max_cd_iter = 1
    thol = max(float(tol), HUGE_THOL)

    S = huge_correlation(X)
    coefficients, _hit_max_iter = _core.solve_mb_path(
        S, alpha_path, thol, max_cd_iter
    )
    return np.ascontiguousarray(coefficients)
