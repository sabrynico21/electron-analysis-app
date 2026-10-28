"""cclasso — Python port of CCLasso (R package, v2.0).

Faithful translation of ``CCLasso-master/R/cclasso.R`` (Fang et al., 2015/2016).
See ``docs/PORTING_NOTES.md`` for the function-by-function mapping and
``docs/EQUIVALENCE.md`` for what is guaranteed identical vs only numerically close.
"""

from __future__ import annotations

from ._core import boot_cclasso, cclasso, cclasso_sub, cv_loss_cclasso
from .validation import CCLassoInputError

__version__ = "0.1.0"

__all__ = [
    "cclasso",
    "cclasso_sub",
    "cv_loss_cclasso",
    "boot_cclasso",
    "CCLassoInputError",
    "__version__",
]
