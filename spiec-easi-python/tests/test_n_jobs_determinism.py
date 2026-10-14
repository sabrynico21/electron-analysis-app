"""``n_jobs`` must not change the result: fixed-order, deterministic aggregation.

The StARS subsamples are independent and each is solved by the compiled MB core
(which releases the GIL), so a thread pool gives real parallelism. The results
must be *bit-identical* for any worker count and for either joblib backend: the
adjacency paths are stacked in submission order and summed in that order, and no
global mutable state is involved.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spieceasi.pipeline import spiec_easi_mb  # noqa: E402
from spieceasi.runtime import default_n_jobs, physical_core_count  # noqa: E402

CASE = Path(__file__).parent / "golden_extended" / "rand_p10_n30_s201"


def _run(n_jobs: int, prefer: str):
    counts = np.loadtxt(CASE / "counts.tsv", delimiter="\t")
    return spiec_easi_mb(
        counts,
        n_lambdas=8,
        lambda_min_ratio=1e-2,
        n_subsamples=4,
        stars_threshold=0.05,
        sym_rule="or",
        random_state=17,
        max_iter=10_000,
        n_jobs=n_jobs,
        prefer=prefer,
        require_samples_gt_taxa=False,
    )


@pytest.mark.parametrize("n_jobs", [2, 4, 8])
def test_n_jobs_does_not_change_results(n_jobs: int) -> None:
    """Threaded runs are bit-identical to the serial run."""
    base = _run(1, "threads")
    parallel = _run(n_jobs, "threads")

    assert np.array_equal(parallel.adjacency, base.adjacency)
    assert np.array_equal(parallel.weights, base.weights)
    assert np.array_equal(parallel.directed_coefficients, base.directed_coefficients)
    assert np.array_equal(parallel.instability, base.instability)
    assert np.array_equal(parallel.monotone_instability, base.monotone_instability)
    assert np.array_equal(parallel.alpha_path, base.alpha_path)
    assert parallel.selected_lambda == base.selected_lambda


def test_process_backend_still_matches_threads() -> None:
    """The legacy process pool remains available and equally deterministic."""
    threaded = _run(2, "threads")
    processed = _run(2, "processes")

    assert np.array_equal(processed.adjacency, threaded.adjacency)
    assert np.array_equal(processed.instability, threaded.instability)
    assert processed.selected_lambda == threaded.selected_lambda


def test_default_n_jobs_uses_physical_cores() -> None:
    """The default worker count is never more than the physical core count."""
    default = default_n_jobs()
    assert default >= 1

    physical = physical_core_count()
    if physical is not None:
        assert default == physical
        logical = os.cpu_count()
        if logical:
            # sanity: physical cores can never exceed the logical ones
            assert physical <= logical
    else:
        logical = os.cpu_count()
        assert default == (max(1, logical // 2) if logical else 1)


def test_default_n_jobs_stays_serial_when_nothing_is_detectable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Conservative fallback: 1 worker when neither count can be read."""
    from spieceasi import runtime

    monkeypatch.setattr(runtime, "physical_core_count", lambda: None)
    monkeypatch.setattr(runtime.os, "cpu_count", lambda: None)
    assert runtime.default_n_jobs() == 1

    monkeypatch.setattr(runtime.os, "cpu_count", lambda: 8)
    assert runtime.default_n_jobs() == 4  # assume SMT, stay on the safe side
