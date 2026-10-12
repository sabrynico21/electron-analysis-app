from __future__ import annotations

import pytest


@pytest.fixture(scope="session")
def r_runtime():
    """Provide an initialized R runtime or skip R validation tests.

    This fixture intentionally performs availability checks at runtime instead of
    collection time, so environments without R can still collect the suite.
    """
    try:
        from rpy2 import robjects
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects.packages import importr
    except Exception as exc:  # pragma: no cover - depends on external runtime
        pytest.skip(f"rpy2 is not available: {exc}")

    try:
        robjects.r('.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))')
        for pkg in ("SpiecEasi", "pulsar", "glmnet"):
            importr(pkg)
    except Exception as exc:  # pragma: no cover - depends on external runtime
        pytest.skip(f"R package stack unavailable: {exc}")

    robjects.r("set.seed(1)")
    return {
        "robjects": robjects,
        "numpy2ri": numpy2ri,
        "pandas2ri": pandas2ri,
        "localconverter": localconverter,
    }
