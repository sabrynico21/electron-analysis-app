"""Phase 2 — R/Python equivalence tests for SPIEC-EASI MB.

Compares the Python implementation (``spiec_easi_mb`` and its components) against
golden files produced by the reference R stack via
``tests/r_reference/generate_golden.R``.

Key design decisions:

* R's StARS subsample indices are **injected** into Python (the R RNG is not
  reproducible with numpy). This isolates every deterministic step and the
  deterministic consequence of identical subsampling.
* Tolerances are declared a priori (see constants below) and are **not** widened
  to force passing tests.
* If a step diverges, the test fails at the **first** diverging step and names it.

Run the generator first::

    Rscript tests/r_reference/generate_golden.R

If the golden directory is missing the whole module is skipped.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

GOLDEN_DIR = Path(__file__).parent / "golden"
SELECTION_DIR = GOLDEN_DIR / "selection"

# joblib worker count for StARS (parallel over subsamples). Parallel execution
# does not change any per-subsample result (the MB solver is deterministic).
N_JOBS = int(os.environ.get("SPIECEASI_TEST_JOBS", "-1"))

pytestmark = pytest.mark.skipif(
    not GOLDEN_DIR.is_dir(),
    reason=(
        "Golden files not found. Generate them with "
        "`Rscript tests/r_reference/generate_golden.R`."
    ),
)

# ---------------------------------------------------------------------------
# A-priori tolerances (declared before running; not tuned to outcomes).
# ---------------------------------------------------------------------------
# Well-conditioned deterministic transforms.
TOL_TRANSFORM = dict(rtol=1e-12, atol=1e-12)
# lambda.max and the geometric lambda grid.
TOL_LAMBDA = dict(rtol=1e-12, atol=1e-15)
# MB lasso coefficients at a fixed lambda: huge C++ CD vs Python CD replica.
# Both solve the same strictly-convex problem to KKT certification, so agreement
# is expected to be near machine precision; the tolerance only absorbs BLAS
# accumulation-order differences.
TOL_COEFF = dict(rtol=1e-6, atol=1e-8)
# StARS instability summaries with identical subsample indices + grid.
TOL_STARS = dict(rtol=1e-6, atol=1e-8)
# Selected optimal lambda.
TOL_OPT_LAMBDA = dict(rtol=1e-9, atol=0.0)


class FirstDivergenceError(Exception):
    """Raised with the name of the first step where Python diverges from R."""


def _check(step: str, ok: bool, detail: str) -> None:
    if not ok:
        raise FirstDivergenceError(f"[{step}] {detail}")


def _resolve_golden(path: Path) -> Path:
    """Return ``path``, or its ``.gz`` companion when only that exists.

    The large golden matrices (``solver_adj_path.tsv``) are stored gzipped to
    keep the fixture tree manageable; ``numpy.loadtxt`` decompresses a ``.gz``
    path transparently.  The decompressed bytes are identical to the original
    files, see ``golden_extended/CHECKSUMS.sha256``.
    """
    if path.is_file():
        return path
    compressed = path.with_name(path.name + ".gz")
    return compressed if compressed.is_file() else path


def _load_matrix(path: Path) -> np.ndarray:
    arr = np.loadtxt(_resolve_golden(path), delimiter="\t")
    return np.atleast_2d(arr)


def _load_vector(path: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(_resolve_golden(path), delimiter="\t"))


def _load_scalar(path: Path) -> float:
    return float(np.loadtxt(_resolve_golden(path)))


def _load_params(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def _cases() -> list[str]:
    params = sorted(GOLDEN_DIR.glob("*.params.txt"))
    return [p.name[: -len(".params.txt")] for p in params]


def _max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b))) if a.size else 0.0


SOLVER_CASES = [c for c in _cases() if c.startswith(("ds3_", "ds4_"))]


@pytest.mark.slow
@pytest.mark.parametrize("case", SOLVER_CASES)
def test_solver_path_matches_r(case: str) -> None:
    """Compare the warm-started MB path solver with huge **before** StARS.

    One warm-started solve of the full lambda path per pulsar subsample, exactly
    as ``pulsar``'s ``fun`` does. If this fails, the divergence originates in the
    solver; if it passes but the StARS test fails, the divergence is in StARS.
    """
    from spieceasi.mb_solver import solve_mb_path_huge_scale, symmetrize_path
    from spieceasi.transforms import clr_transform

    dims = _load_params(GOLDEN_DIR / f"{case}.solver_dims.txt")
    d = int(dims["d"])
    n_lambda = int(dims["n_lambda"])
    n_sub = int(dims["n_sub"])

    counts = _load_matrix(GOLDEN_DIR / f"{case}.counts.tsv")
    r_grid = _load_vector(GOLDEN_DIR / f"{case}.lambda_path.txt")
    indices = _load_matrix(GOLDEN_DIR / f"{case}.subsample_indices.tsv").astype(int)
    X_clr = clr_transform(counts).clr

    r_adj = _load_matrix(GOLDEN_DIR / f"{case}.solver_adj_path.tsv").astype(int)
    r_adj = r_adj.reshape(n_sub, n_lambda, d, d)

    try:
        for s in range(n_sub):
            directed = solve_mb_path_huge_scale(
                X_clr[indices[s]], r_grid, max_iter=10_000, tol=1e-4
            )
            adj = symmetrize_path(directed, "or").astype(int)
            mismatch = adj != r_adj[s]
            n_bad = int(mismatch.sum() // 2)
            if n_bad:
                k, i, j = np.argwhere(mismatch)[0]
                raise FirstDivergenceError(
                    f"solver adjacency path: subsample {s}, lambda index {int(k)} "
                    f"(alpha={r_grid[int(k)]:.10g}), edge ({int(i)},{int(j)}); "
                    f"{n_bad} differing edge(s) at that lambda"
                )

        r_coef = _load_matrix(GOLDEN_DIR / f"{case}.solver_coef_path.tsv")
        r_coef = r_coef.reshape(n_lambda, d, d)
        directed0 = solve_mb_path_huge_scale(
            X_clr[indices[0]], r_grid, max_iter=10_000, tol=1e-4
        )
        _check(
            "solver coefficient path (subsample 0)",
            np.allclose(directed0, r_coef, **TOL_COEFF),
            f"max|diff|={_max_abs_diff(directed0, r_coef):.3e}",
        )
    except FirstDivergenceError as exc:
        pytest.fail(f"{case}: first divergence [{exc}]", pytrace=False)


@pytest.mark.slow
@pytest.mark.parametrize("case", _cases())
def test_python_matches_r_golden(case: str) -> None:
    from spieceasi.neighborhood import (
        alpha_max_for_mb,
        make_alpha_path,
        mb_neighborhood_selection,
    )
    from spieceasi.pipeline import spiec_easi_mb
    from spieceasi.transforms import clr_transform

    params = _load_params(GOLDEN_DIR / f"{case}.params.txt")
    nl = int(params["nlambda"])
    ratio = float(params["lambda_min_ratio"])
    rep_num = int(params["rep.num"])
    threshold = float(params["threshold"])
    fixed_index = int(params["fixed_index"])  # 1-based in R

    counts = _load_matrix(GOLDEN_DIR / f"{case}.counts.tsv")

    try:
        # -- Step 1: CLR transform -----------------------------------------
        r_clr = _load_matrix(GOLDEN_DIR / f"{case}.clr.tsv")
        py_clr = clr_transform(counts, pseudocount=1.0, zero_method="additive").clr
        _check(
            "1_transform_clr",
            py_clr.shape == r_clr.shape
            and np.allclose(py_clr, r_clr, **TOL_TRANSFORM),
            f"CLR max|diff|={_max_abs_diff(py_clr, r_clr):.3e} "
            f"(shape py={py_clr.shape} r={r_clr.shape})",
        )

        # -- Step 2: lambda.max --------------------------------------------
        r_lmax = _load_scalar(GOLDEN_DIR / f"{case}.lambda_max.txt")
        py_lmax = alpha_max_for_mb(py_clr)
        _check(
            "2_lambda_max",
            np.isclose(py_lmax, r_lmax, **TOL_LAMBDA),
            f"lambda.max py={py_lmax!r} r={r_lmax!r} "
            f"(rel={abs(py_lmax - r_lmax) / max(abs(r_lmax), 1e-300):.3e})",
        )

        # -- Step 3: lambda grid -------------------------------------------
        r_grid = _load_vector(GOLDEN_DIR / f"{case}.lambda_path.txt")
        py_grid = make_alpha_path(py_clr, n_lambdas=nl, lambda_min_ratio=ratio)
        _check(
            "3_lambda_path",
            py_grid.shape == r_grid.shape
            and np.allclose(py_grid, r_grid, **TOL_LAMBDA),
            f"grid max|diff|={_max_abs_diff(py_grid, r_grid):.3e}",
        )

        # -- Step 4: MB coefficients at fixed lambda (isolate the lasso) ----
        fixed_lambda = _load_scalar(GOLDEN_DIR / f"{case}.fixed_lambda.txt")
        r_coef = _load_matrix(GOLDEN_DIR / f"{case}.coef_fixed.tsv")
        iso = mb_neighborhood_selection(
            py_clr,
            alpha=fixed_lambda,
            sym_rule="or",
            standardize=True,
            max_iter=10_000,
            tol=1e-4,
        )
        py_coef = iso.directed_coefficients
        _check(
            "4_mb_coefficients_fixed_lambda",
            py_coef.shape == r_coef.shape
            and np.allclose(py_coef, r_coef, **TOL_COEFF),
            f"coef max|diff|={_max_abs_diff(py_coef, r_coef):.3e} "
            f"(lambda={fixed_lambda:.6g})",
        )

        # -- Step 5: symmetrized adjacency at fixed lambda ------------------
        r_adj_fixed = _load_matrix(GOLDEN_DIR / f"{case}.adj_fixed.tsv").astype(int)
        py_adj_fixed = iso.adjacency.astype(int)
        n_diff_edges = int(np.sum(py_adj_fixed != r_adj_fixed) // 2)
        _check(
            "5_adjacency_fixed_lambda",
            np.array_equal(py_adj_fixed, r_adj_fixed),
            f"{n_diff_edges} edge(s) differ (lambda={fixed_lambda:.6g})",
        )

        # -- Steps 6-10: full pipeline with injected R subsample indices ----
        r_indices = _load_matrix(
            GOLDEN_DIR / f"{case}.subsample_indices.tsv"
        ).astype(int)  # 0-based, shape (rep.num, size)
        subsamples = [r_indices[i] for i in range(r_indices.shape[0])]

        result = spiec_easi_mb(
            counts,
            n_lambdas=nl,
            lambda_min_ratio=ratio,
            n_subsamples=rep_num,
            alpha_path=py_grid,
            subsamples=subsamples,
            stars_threshold=threshold,
            sym_rule="or",
            random_state=0,
            max_iter=10_000,
            tol=1e-4,
            n_jobs=N_JOBS,
            require_samples_gt_taxa=False,
        )

        # -- Step 6: StARS instability (raw) --------------------------------
        r_raw = _load_vector(GOLDEN_DIR / f"{case}.stars_raw.txt")
        _check(
            "6_stars_instability_raw",
            result.instability.shape == r_raw.shape
            and np.allclose(result.instability, r_raw, **TOL_STARS),
            f"instability max|diff|={_max_abs_diff(result.instability, r_raw):.3e}",
        )

        # -- Step 7: monotone instability -----------------------------------
        r_mono = _load_vector(GOLDEN_DIR / f"{case}.stars_monotone.txt")
        _check(
            "7_stars_instability_monotone",
            result.monotone_instability.shape == r_mono.shape
            and np.allclose(result.monotone_instability, r_mono, **TOL_STARS),
            f"monotone max|diff|="
            f"{_max_abs_diff(result.monotone_instability, r_mono):.3e}",
        )

        # -- Step 8/9: optimal index & lambda -------------------------------
        r_opt_index = int(_load_scalar(GOLDEN_DIR / f"{case}.opt_index.txt"))  # 1-based
        r_opt_lambda = _load_scalar(GOLDEN_DIR / f"{case}.opt_lambda.txt")
        py_opt_matches = np.isclose(
            result.selected_lambda, r_opt_lambda, **TOL_OPT_LAMBDA
        )
        py_index_matches = np.isclose(
            py_grid, result.selected_lambda, **TOL_OPT_LAMBDA
        )
        py_opt_index = int(np.argmax(py_index_matches))
        _check(
            "8_stars_opt_index",
            py_index_matches[py_opt_index] and py_opt_index == r_opt_index - 1,
            f"opt index py={py_opt_index + 1} r={r_opt_index} "
            f"(lambda py={result.selected_lambda:.10g} r={r_opt_lambda:.10g})",
        )
        _check(
            "9_stars_opt_lambda",
            py_opt_matches,
            f"opt lambda py={result.selected_lambda:.17g} r={r_opt_lambda:.17g}",
        )

        # -- Step 10: final network (full-data refit at optimal lambda) -----
        r_final = _load_matrix(GOLDEN_DIR / f"{case}.final_adj.tsv").astype(int)
        py_final = result.adjacency.astype(int)
        n_diff_final = int(np.sum(py_final != r_final) // 2)
        _check(
            "10_final_adjacency",
            np.array_equal(py_final, r_final),
            f"{n_diff_final} edge(s) differ "
            f"(selected lambda py={result.selected_lambda:.10g})",
        )

    except FirstDivergenceError as exc:
        pytest.fail(f"{case}: first divergence {exc}", pytrace=False)


# ---------------------------------------------------------------------------
# Selection-rule boundary test (pulsar:::.starsind vs Python selection rule).
# Mirrors the current Python selection logic in spieceasi/stars.py:
#   monotone = np.maximum.accumulate(instability)
#   feasible = np.where(monotone <= threshold)[0]
#   selected_idx = feasible[-1] if feasible.size else 0
# ---------------------------------------------------------------------------
def _python_select_index(monotone: np.ndarray, threshold: float) -> int:
    # Exercise the production selection rule (pulsar:::.starsind replica).
    from spieceasi.stars import _stars_index

    return _stars_index(monotone, threshold)


def _selection_cases() -> list[str]:
    if not SELECTION_DIR.is_dir():
        return []
    return sorted(p.stem for p in SELECTION_DIR.glob("*.txt"))


@pytest.mark.parametrize("name", _selection_cases())
def test_selection_index_matches_r_starsind(name: str) -> None:
    lines = (SELECTION_DIR / f"{name}.txt").read_text().splitlines()
    threshold = float(lines[0])
    r_index_1b = int(lines[1])
    summary = np.array([float(v) for v in lines[2:]], dtype=float)

    py_index_0b = _python_select_index(summary, threshold)
    r_index_0b = r_index_1b - 1

    assert py_index_0b == r_index_0b, (
        f"selection boundary case '{name}': Python index (0-based)={py_index_0b} "
        f"vs R index (0-based)={r_index_0b} for summary={summary.tolist()}, "
        f"threshold={threshold}"
    )
