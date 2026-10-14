"""Phase D — application integration tests (`python/analysis.py`) against R.

These tests drive the *application* entry point the Electron app calls,
``python/analysis.py::_compute_spieceasi_for_matrix``, and compare its output
with R goldens in ``tests/golden_app/`` (generate with
``Rscript tests/r_reference/generate_app_golden.R``).

Design notes
------------
* R's StARS subsample indices are injected through the optional ``subsamples``
  argument added to the application function (default ``None`` keeps the app
  behaviour unchanged);
* a case with an extra constant feature checks that the app drops it (R would
  error on it) and still matches the golden computed on the reduced matrix;
* the "tiny + high-dimensional" regime (``n <= 10 & p >= 200``) is covered at
  the parameter level: such a fit is far too slow for a test, so a stub records
  the arguments the application would forward.

No golden file is modified and no tolerance is widened.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
GOLDEN_APP = Path(__file__).parent / "golden_app"
MANIFEST = GOLDEN_APP / "manifest.json"

if str(REPO_ROOT / "python") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "python"))

import analysis  # noqa: E402  (path set above)

# same a-priori tolerance as the other equivalence tests (no widening)
TOL_OPT_LAMBDA = dict(rtol=1e-9, atol=1e-12)

pytestmark = pytest.mark.skipif(
    not MANIFEST.is_file(),
    reason=(
        "App golden not found. Generate with "
        "`Rscript tests/r_reference/generate_app_golden.R`."
    ),
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _cases() -> list[str]:
    return list(json.loads(MANIFEST.read_text()))


def _read_meta(path: Path) -> dict:
    meta: dict[str, str] = {}
    for line in path.read_text().split():
        if "=" in line:
            key, value = line.split("=", 1)
            meta[key] = value
    return meta


def _read_counts(path: Path) -> tuple[list[str], list[str], dict]:
    """Read a samples x features TSV with header/row names into app structures."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    header = lines[0].split("\t")
    features = header[1:]
    sample_ids: list[str] = []
    matrix_by_sample: dict[str, dict[str, float]] = {}
    for line in lines[1:]:
        cells = line.split("\t")
        sample_id = cells[0]
        sample_ids.append(sample_id)
        matrix_by_sample[sample_id] = {
            feat: float(val) for feat, val in zip(features, cells[1:])
        }
    return features, sample_ids, matrix_by_sample


def _app_params(meta: dict) -> dict:
    return {
        "spieceasiNLambdas": int(meta["nlambda"]),
        "spieceasiLambdaMinRatio": float(meta["lambda_min_ratio"]),
        "spieceasiNSubsamples": int(meta["rep_num"]),
        "spieceasiStarsThreshold": float(meta["threshold"]),
        "spieceasiJobs": 1,
    }


def _edge_pairs_from_pairs(features: list[str], is_edge) -> set[tuple[str, str]]:
    return {
        (features[i], features[j])
        for i in range(len(features))
        for j in range(i + 1, len(features))
        if is_edge(i, j)
    }


def _argmax_abs_corr_pair(clr: np.ndarray) -> tuple[int, int]:
    from spieceasi.mb_solver import huge_correlation

    corr = np.abs(huge_correlation(clr))
    np.fill_diagonal(corr, 0.0)
    i, j = np.unravel_index(int(np.argmax(corr)), corr.shape)
    return (min(int(i), int(j)), max(int(i), int(j)))


# ---------------------------------------------------------------------------
# Slow: real application run vs R golden
# ---------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("case", _cases())
def test_app_matches_r_golden(case: str) -> None:
    cdir = GOLDEN_APP / case
    meta = _read_meta(cdir / "meta.txt")
    kept = (cdir / "expected_kept_features.txt").read_text().split()

    features, sample_ids, matrix_by_sample = _read_counts(cdir / "counts_full.tsv")
    indices = np.loadtxt(cdir / "subsample_indices.tsv", delimiter="\t").astype(int)
    subsamples = [indices[i] for i in range(indices.shape[0])]

    corr, _pvals, diag = analysis._compute_spieceasi_for_matrix(
        matrix_by_sample,
        features,
        sample_ids,
        _app_params(meta),
        subsamples=subsamples,
    )

    # the app must have dropped exactly the constant feature(s)
    assert sorted(corr) == sorted(kept)

    # selected lambda
    r_lambda = float((cdir / "opt_lambda.txt").read_text().strip())
    assert np.isclose(diag["selected_lambda"], r_lambda, **TOL_OPT_LAMBDA), (
        f"{case}: app lambda {diag['selected_lambda']!r} != R {r_lambda!r}"
    )

    # network
    r_adj = np.loadtxt(cdir / "final_adj.tsv", delimiter="\t").astype(int)
    r_edges = _edge_pairs_from_pairs(
        kept, lambda i, j: r_adj[i, j] != 0 or r_adj[j, i] != 0
    )
    py_edges = _edge_pairs_from_pairs(kept, lambda i, j: corr[kept[i]][kept[j]] != 0)

    grid = np.loadtxt(cdir / "lambda_path.txt")
    py_index = int(np.argmin(np.abs(grid - diag["selected_lambda"])))

    differing = py_edges ^ r_edges
    if py_index == 0:
        # lambda = lambda.max: numerically indeterminate even in R (see
        # docs/AUDIT.md §2.E). At most one differing edge, and only the
        # argmax|cor| pair, is tolerated here.
        clr = np.loadtxt(cdir / "clr.tsv", delimiter="\t")
        am = _argmax_abs_corr_pair(clr)
        argmax_pair = (kept[am[0]], kept[am[1]])
        assert differing <= {argmax_pair}, (
            f"{case}: lambda=lambda.max, differing={sorted(differing)}, "
            f"allowed only {{{argmax_pair}}}"
        )
    else:
        assert not differing, (
            f"{case}: {len(differing)} edge(s) differ: {sorted(differing)}"
        )

    assert diag["edge_count"] == len(py_edges)


# ---------------------------------------------------------------------------
# Fast: application parameter handling
# ---------------------------------------------------------------------------
class _Recorder:
    """Stand-in for spiec_easi_mb that records the forwarded kwargs."""

    def __init__(self) -> None:
        self.kwargs: dict | None = None

    def __call__(self, data, **kwargs):  # noqa: ANN001, ANN003
        self.kwargs = kwargs
        n_taxa = data.shape[1]

        class _Result:
            adjacency = np.zeros((n_taxa, n_taxa), dtype=int)
            weights = np.zeros((n_taxa, n_taxa))
            selected_lambda = 1.0
            instability = np.zeros(1)
            monotone_instability = np.zeros(1)

        return _Result()


def _synthetic(n_samples: int, n_taxa: int, seed: int = 0):
    rng = np.random.default_rng(seed)
    counts = rng.integers(1, 200, size=(n_samples, n_taxa))
    features = [f"f{i:03d}" for i in range(n_taxa)]
    sample_ids = [f"s{i:02d}" for i in range(n_samples)]
    matrix_by_sample = {
        sample_ids[i]: {features[j]: float(counts[i, j]) for j in range(n_taxa)}
        for i in range(n_samples)
    }
    return features, sample_ids, matrix_by_sample


def _capture(monkeypatch, n_samples: int, n_taxa: int, params: dict) -> dict:
    import spieceasi

    recorder = _Recorder()
    monkeypatch.setattr(spieceasi, "spiec_easi_mb", recorder)
    features, sample_ids, matrix_by_sample = _synthetic(n_samples, n_taxa)
    analysis._compute_spieceasi_for_matrix(
        matrix_by_sample,
        features,
        sample_ids,
        params,
        subsamples=[np.arange(min(4, n_samples)), np.arange(min(4, n_samples))],
    )
    assert recorder.kwargs is not None
    return recorder.kwargs


def test_app_forwards_pulsar_subsample_rule_by_default(monkeypatch) -> None:
    kwargs = _capture(monkeypatch, n_samples=40, n_taxa=6, params={})
    # None -> the library applies pulsar::.ratcheck (0.8 for n <= 144)
    assert kwargs["subsample_fraction"] is None

    kwargs = _capture(
        monkeypatch, n_samples=40, n_taxa=6, params={"spieceasiSubsampleFraction": 0.5}
    )
    assert kwargs["subsample_fraction"] == 0.5


def test_app_defaults_match_r_reference(monkeypatch) -> None:
    kwargs = _capture(monkeypatch, n_samples=40, n_taxa=6, params={})
    assert kwargs["n_lambdas"] == 20
    assert kwargs["lambda_min_ratio"] == 1e-3
    assert kwargs["n_subsamples"] == 20
    assert kwargs["stars_threshold"] == 0.05
    assert kwargs["sym_rule"] == "or"
    assert kwargs["max_iter"] == 10_000
    assert kwargs["require_samples_gt_taxa"] is False


def test_app_tiny_highdim_regime_caps_defaults(monkeypatch) -> None:
    kwargs = _capture(monkeypatch, n_samples=8, n_taxa=220, params={})
    assert kwargs["max_iter"] == 3000
    assert kwargs["n_lambdas"] == 10
    assert kwargs["n_subsamples"] == 8


def test_app_tiny_highdim_regime_respects_explicit_overrides(monkeypatch) -> None:
    kwargs = _capture(
        monkeypatch,
        n_samples=8,
        n_taxa=220,
        params={
            "spieceasiNLambdas": 7,
            "spieceasiNSubsamples": 3,
            "spieceasiMaxIter": 123,
        },
    )
    assert kwargs["max_iter"] == 123
    assert kwargs["n_lambdas"] == 7
    assert kwargs["n_subsamples"] == 3


# ---------------------------------------------------------------------------
# Slow: the frozen (PyInstaller) backend must run SPIEC-EASI at its defaults
# ---------------------------------------------------------------------------
FROZEN_BACKEND = REPO_ROOT / "python-dist" / "analysis-backend" / "analysis-backend"
TINY_DIR = REPO_ROOT / "benchmarks" / "tiny"


def _tiny_params(job_id: str, **extra) -> dict:
    return {
        "jobId": job_id,
        "files": {
            "bacteria": [str(TINY_DIR / "counts.txt")],
            "bacteria_metadata": [str(TINY_DIR / "metadata.txt")],
        },
        "params": {
            "datasetType": "bacteria",
            "correlationMethod": "spieceasi",
            "conditionColumn": "Combination_treat2",
            "conditionA": "C1",
            "conditionB": "C2",
            "groupingMode": "all",
            "pValueThreshold": 0.05,
            "spieceasiNLambdas": 20,
            "spieceasiLambdaMinRatio": 0.001,
            "spieceasiNSubsamples": 20,
            "spieceasiStarsThreshold": 0.05,
            **extra,
        },
    }


def _run_backend(binary: list[str], params: dict, tmp_path: Path, name: str) -> tuple[int, str, dict]:
    params_path = tmp_path / f"{name}.params.json"
    result_path = tmp_path / f"{name}.result.json"
    params_path.write_text(json.dumps(params))
    proc = subprocess.run(
        [*binary, "analyze", "--params", str(params_path), "--output", str(result_path)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result_path.read_text()) if result_path.is_file() else {}
    return proc.returncode, proc.stdout + proc.stderr, payload


@pytest.mark.slow
@pytest.mark.skipif(
    not FROZEN_BACKEND.is_file(),
    reason="frozen backend not built (run `npm run build:backend`)",
)
def test_frozen_backend_runs_spieceasi_with_app_defaults(tmp_path: Path) -> None:
    """The app never sends `spieceasiJobs`; that default must not crash frozen.

    A process pool would re-execute the frozen binary. Regression test for
    "unknown command: -B" (joblib/loky in a PyInstaller bundle).
    """
    binary = [str(FROZEN_BACKEND)]

    rc_default, log_default, payload_default = _run_backend(
        binary, _tiny_params("frozen-default"), tmp_path, "default"
    )
    assert rc_default == 0, log_default
    # the frozen build must select the thread backend (never a process pool)
    assert "parallel_backend=threads" in log_default, log_default
    assert "unknown command" not in log_default

    rc_one, log_one, payload_one = _run_backend(
        binary, _tiny_params("frozen-default", spieceasiJobs=1), tmp_path, "jobs1"
    )
    assert rc_one == 0, log_one

    # identical payload, not merely identical shape
    assert payload_default == payload_one
    assert payload_default["graph_pipeline"] == payload_one["graph_pipeline"]
    assert payload_default["summary"] == payload_one["summary"]


@pytest.mark.slow
@pytest.mark.skipif(
    not FROZEN_BACKEND.is_file(),
    reason="frozen backend not built (run `npm run build:backend`)",
)
def test_frozen_backend_uses_the_compiled_mb_core(tmp_path: Path) -> None:
    """The packaged app must ship *and use* the compiled MB core.

    Without it the backend silently falls back to the pure-Python reference
    solver, which is roughly two orders of magnitude slower; that would be a
    silent performance regression in the shipped application.
    """
    rc, log, _payload = _run_backend(
        [str(FROZEN_BACKEND)], _tiny_params("frozen-solver"), tmp_path, "solver"
    )
    assert rc == 0, log
    assert "SPIEC-EASI solver: cpp" in log, log
    assert "compiled MB core (spieceasi._mb_core) is not available" not in log, log
