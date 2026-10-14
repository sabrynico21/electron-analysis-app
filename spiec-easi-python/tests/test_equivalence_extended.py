"""Phase B — extended R/Python equivalence tests on a larger case matrix.

Golden files live in ``tests/golden_extended/`` (regenerate with
``Rscript tests/r_reference/generate_golden_extended.R``); the case list and
metadata are in ``tests/golden_extended/manifest.json``.

Design mirrors ``tests/test_equivalence.py``:
* R's StARS subsample indices are injected into Python (R RNG is not
  reproducible with numpy);
* the same a-priori tolerances are used (imported from ``test_equivalence``);
* heavy end-to-end cases are marked ``slow`` (excluded by default), with the
  per-case wall time printed as ``CASEREPORT`` for the coverage table;
* the first diverging step is reported (with subsample/lambda where relevant).

One narrow, explicitly scoped exception is allowed at step 9: when the selected
lambda is ``lambda.max`` (grid index 0), the result is numerically indeterminate
even in R (see ``docs/AUDIT.md`` §2.E). Only there is a difference of at most ONE
edge tolerated, and only if that edge is exactly the ``argmax|cor|`` pair. In
every other case the comparison stays exact (zero differing edges).

Existing golden files and tolerances are untouched.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pytest

from tests.test_equivalence import (  # reuse the SAME a-priori tolerances
    TOL_COEFF,
    TOL_LAMBDA,
    TOL_OPT_LAMBDA,
    TOL_STARS,
    TOL_TRANSFORM,
    _check,
    _load_matrix,
    _load_params,
    _load_scalar,
    _load_vector,
    _max_abs_diff,
    FirstDivergenceError,
)

GOLDEN_EXT_DIR = Path(__file__).parent / "golden_extended"
MANIFEST = GOLDEN_EXT_DIR / "manifest.json"

pytestmark = pytest.mark.skipif(
    not MANIFEST.is_file(),
    reason=(
        "Extended golden not found. Generate with "
        "`Rscript tests/r_reference/generate_golden_extended.R`."
    ),
)


def _symmetric_edge_pairs(adj: np.ndarray) -> set[tuple[int, int]]:
    """Unordered edge pairs (i < j) present in a possibly asymmetric adjacency."""
    d = adj.shape[0]
    return {
        (i, j)
        for i in range(d)
        for j in range(i + 1, d)
        if adj[i, j] != 0 or adj[j, i] != 0
    }


def _argmax_abs_corr_pair(X: np.ndarray) -> tuple[int, int]:
    """Index pair (i < j) attaining the largest absolute correlation in ``X``."""
    from spieceasi.mb_solver import huge_correlation

    corr = np.abs(huge_correlation(X))
    np.fill_diagonal(corr, 0.0)
    i, j = np.unravel_index(int(np.argmax(corr)), corr.shape)
    return (min(int(i), int(j)), max(int(i), int(j)))


def _manifest() -> list[dict]:
    return json.loads(MANIFEST.read_text())


def _case_ids() -> list[str]:
    return [c["id"] for c in _manifest()]


# ---------------------------------------------------------------------------
# Fast consistency checks (run by default, no R required)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("case", _case_ids())
def test_extended_golden_consistent(case: str) -> None:
    cdir = GOLDEN_EXT_DIR / case
    params = _load_params(cdir / "solver_dims.txt")
    d, n_lambda, n_sub = int(params["d"]), int(params["n_lambda"]), int(params["n_sub"])
    grid = _load_vector(cdir / "lambda_path.txt")
    assert grid.size == n_lambda
    assert np.all(np.diff(grid) < 0), f"{case}: lambda path not decreasing"
    indices = _load_matrix(cdir / "subsample_indices.tsv")
    assert indices.shape[0] == n_sub
    adj = _load_matrix(cdir / "solver_adj_path.tsv")
    assert adj.shape == (n_sub * n_lambda, d * d)
    final = _load_matrix(cdir / "final_adj.tsv")
    assert final.shape == (d, d)


def _pulsar_subsample_width(n: int) -> int:
    """pulsar::.ratcheck size: ``floor(ratio * n)``, ratio 0.8 (n<=144) else 10/sqrt(n)."""
    if n > 144:
        return int(np.floor(10 * np.sqrt(n)))
    return int(np.floor(0.8 * n))


@pytest.mark.parametrize("case", _case_ids())
def test_extended_subsample_size_follows_pulsar_rule(case: str) -> None:
    """R's StARS subsample width must follow ``pulsar::.ratcheck`` (incl. n > 144)."""
    cdir = GOLDEN_EXT_DIR / case
    params = _load_params(cdir / "solver_dims.txt")
    n = int(params["n"])
    meta = next(c for c in _manifest() if c["id"] == case)
    indices = _load_matrix(cdir / "subsample_indices.tsv").astype(int)
    expected = _pulsar_subsample_width(n)
    assert indices.shape[0] == int(meta["rep_num"])
    assert indices.shape[1] == expected, (
        f"{case}: R subsample width {indices.shape[1]} != pulsar::.ratcheck "
        f"expectation {expected} for n={n}"
    )
    assert indices.shape[1] > 1
    assert indices.min() >= 0 and indices.max() < n


def test_infer_subsample_size_matches_pulsar() -> None:
    from spieceasi.stars import infer_subsample_size

    # n <= 144 -> ratio 0.8
    assert infer_subsample_size(100) == 80
    assert infer_subsample_size(144) == int(np.floor(0.8 * 144))
    # n > 144 -> floor(10*sqrt(n))
    assert infer_subsample_size(145) == int(np.floor(10 * np.sqrt(145)))
    assert infer_subsample_size(200) == 141
    assert infer_subsample_size(300) == 173
    # explicit overrides still win
    assert infer_subsample_size(200, subsample_fraction=0.8) == 160
    assert infer_subsample_size(200, subsample_size=77) == 77
    with pytest.raises(ValueError):
        infer_subsample_size(200, subsample_fraction=1.5)


# ---------------------------------------------------------------------------
# Slow: full parity per case
# ---------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("case", _case_ids())
def test_extended_case_matches_r(case: str) -> None:
    from spieceasi.mb_solver import solve_mb_path_huge_scale, symmetrize_path
    from spieceasi.neighborhood import (
        alpha_max_for_mb,
        make_alpha_path,
    )
    from spieceasi.pipeline import spiec_easi_mb
    from spieceasi.transforms import clr_transform

    cdir = GOLDEN_EXT_DIR / case
    t0 = time.perf_counter()
    diffs: dict[str, float] = {}

    dims = _load_params(cdir / "solver_dims.txt")
    d, n_lambda, n_sub = int(dims["d"]), int(dims["n_lambda"]), int(dims["n_sub"])
    case_meta = next(c for c in _manifest() if c["id"] == case)
    nl = int(case_meta["nlambda"])
    ratio = float(case_meta["lambda_min_ratio"])
    rep_num = int(case_meta["rep_num"])

    counts = _load_matrix(cdir / "counts.tsv")

    try:
        # 1. CLR
        r_clr = _load_matrix(cdir / "clr.tsv")
        py_clr = clr_transform(counts).clr
        diffs["clr"] = _max_abs_diff(py_clr, r_clr)
        _check("1_transform_clr",
               py_clr.shape == r_clr.shape and np.allclose(py_clr, r_clr, **TOL_TRANSFORM),
               f"max|diff|={diffs['clr']:.3e}")

        # 2. lambda.max
        r_lmax = _load_scalar(cdir / "lambda_max.txt")
        py_lmax = alpha_max_for_mb(py_clr)
        diffs["lambda_max"] = abs(py_lmax - r_lmax)
        _check("2_lambda_max", np.isclose(py_lmax, r_lmax, **TOL_LAMBDA),
               f"py={py_lmax!r} r={r_lmax!r}")

        # 3. lambda grid
        r_grid = _load_vector(cdir / "lambda_path.txt")
        py_grid = make_alpha_path(py_clr, n_lambdas=nl, lambda_min_ratio=ratio)
        diffs["grid"] = _max_abs_diff(py_grid, r_grid)
        _check("3_lambda_path",
               py_grid.shape == r_grid.shape and np.allclose(py_grid, r_grid, **TOL_LAMBDA),
               f"max|diff|={diffs['grid']:.3e}")

        indices = _load_matrix(cdir / "subsample_indices.tsv").astype(int)
        subsamples = [indices[i] for i in range(indices.shape[0])]

        # 4. solver per subsample (adjacency per lambda) + coef path (subsample 0)
        r_adj = _load_matrix(cdir / "solver_adj_path.tsv").astype(np.int8)
        r_adj = r_adj.reshape(n_sub, n_lambda, d, d)
        solver_max_edges = 0
        first_bad = None
        for s in range(n_sub):
            directed = solve_mb_path_huge_scale(py_clr[indices[s]], py_grid,
                                                max_iter=10_000, tol=1e-4)
            adj = symmetrize_path(directed, "or").astype(np.int8)
            bad = int(np.sum(adj != r_adj[s]) // 2)
            if bad:
                solver_max_edges = max(solver_max_edges, bad)
                if first_bad is None:
                    k, i, j = np.argwhere(adj != r_adj[s])[0]
                    first_bad = (s, int(k), int(i), int(j), bad)
        diffs["solver_edges"] = float(solver_max_edges)
        _check("4_solver_path",
               first_bad is None,
               f"subsample {first_bad[0]}, lambda idx {first_bad[1]}, "
               f"edge ({first_bad[2]},{first_bad[3]}): {first_bad[4]} edge(s) differ"
               if first_bad else "ok")

        r_coef = _load_matrix(cdir / "solver_coef_path.tsv").reshape(n_lambda, d, d)
        directed0 = solve_mb_path_huge_scale(py_clr[indices[0]], py_grid,
                                             max_iter=10_000, tol=1e-4)
        diffs["coef"] = _max_abs_diff(directed0, r_coef)
        _check("4b_solver_coef", np.allclose(directed0, r_coef, **TOL_COEFF),
               f"max|diff|={diffs['coef']:.3e}")

        # 5-9. full pipeline with injected subsamples
        result = spiec_easi_mb(
            counts,
            n_lambdas=nl,
            lambda_min_ratio=ratio,
            n_subsamples=rep_num,
            alpha_path=py_grid,
            subsamples=subsamples,
            stars_threshold=0.05,
            sym_rule="or",
            random_state=0,
            max_iter=10_000,
            tol=1e-4,
            n_jobs=-1,
            require_samples_gt_taxa=False,
        )

        r_raw = _load_vector(cdir / "stars_raw.txt")
        diffs["stars_raw"] = _max_abs_diff(result.instability, r_raw)
        _check("5_stars_raw", np.allclose(result.instability, r_raw, **TOL_STARS),
               f"max|diff|={diffs['stars_raw']:.3e}")

        r_mono = _load_vector(cdir / "stars_monotone.txt")
        diffs["stars_mono"] = _max_abs_diff(result.monotone_instability, r_mono)
        _check("6_stars_monotone",
               np.allclose(result.monotone_instability, r_mono, **TOL_STARS),
               f"max|diff|={diffs['stars_mono']:.3e}")

        r_opt_index = int(_load_scalar(cdir / "opt_index.txt"))
        r_opt_lambda = _load_scalar(cdir / "opt_lambda.txt")
        py_idx = int(np.argmax(np.isclose(py_grid, result.selected_lambda, **TOL_OPT_LAMBDA)))
        diffs["opt_lambda"] = abs(result.selected_lambda - r_opt_lambda)
        _check("7_stars_opt_index", py_idx == r_opt_index - 1,
               f"py={py_idx + 1} r={r_opt_index}")
        _check("8_stars_opt_lambda",
               np.isclose(result.selected_lambda, r_opt_lambda, **TOL_OPT_LAMBDA),
               f"py={result.selected_lambda:.17g} r={r_opt_lambda:.17g}")

        r_final = _load_matrix(cdir / "final_adj.tsv").astype(int)
        py_final = result.adjacency.astype(int)
        differing = _symmetric_edge_pairs(py_final) ^ _symmetric_edge_pairs(r_final)
        diffs["final_edges"] = float(len(differing))

        # Numerical indeterminacy at lambda = lambda.max (docs/AUDIT.md §2.E):
        # huge's MB coordinate descent keeps the argmax|cor| edge for some
        # datasets and drops it for others at exactly this lambda. Only here,
        # and only for that single edge, a difference is tolerated.
        if py_idx == 0:
            argmax_pair = _argmax_abs_corr_pair(py_clr)
            allowed = differing <= {argmax_pair}
            detail = (
                f"{len(differing)} edge(s) differ; differing={sorted(differing)}; "
                f"argmax|cor|={argmax_pair} (lambda=lambda.max: indeterminate in R)"
            )
        else:
            allowed = not differing
            detail = (
                f"{len(differing)} edge(s) differ; differing={sorted(differing)}"
            )
        _check("9_final_adjacency", allowed, detail)

    except FirstDivergenceError as exc:
        print(
            f"CASEREPORT {case} status=FAIL first_divergence=[{exc}] "
            f"t={time.perf_counter() - t0:.1f}s " + " ".join(f"{k}={v:.3e}" for k, v in diffs.items())
        )
        pytest.fail(f"{case}: first divergence {exc}", pytrace=False)

    elapsed = time.perf_counter() - t0
    print(
        f"CASEREPORT {case} status=OK t={elapsed:.1f}s " + " ".join(f"{k}={v:.3e}" for k, v in diffs.items())
    )
