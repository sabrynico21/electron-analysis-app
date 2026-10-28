"""Equivalence tests: Python port vs. R (CCLasso v2.0) golden fixtures.

Run with:  pytest tests/test_equivalence.py -m slow -s

Design:
- The R fixtures in ``tests/golden/`` include the *injected* bootstrap indices and
  artificial-null draws, so the random parts are bit-comparable across languages.
- Intermediate quantities use tight tolerances; the final ADMM-driven outputs use a
  looser (documented) tolerance.  Chosen lambda and ADMM iteration counts must match
  exactly.  Any failing comparison reports the FIRST step that diverges.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from cclasso import cclasso

GOLDEN = Path(os.environ.get("CCLASSO_GOLDEN_DIR", Path(__file__).parent / "golden"))
CASES = sorted(p for p in GOLDEN.iterdir() if p.is_dir() and (p / "params.json").is_file())
FAILURE_CASES = (
    sorted(p for p in (GOLDEN / "failures").iterdir() if p.is_dir())
    if (GOLDEN / "failures").is_dir()
    else []
)

# Tolerances (see docs/EQUIVALENCE.md for the justification of each).
TOL_INTERMEDIATE = dict(rtol=1e-10, atol=1e-12)  # transform, vx2, wd
TOL_LAMS = dict(rtol=1e-12, atol=1e-15)  # pure arithmetic (golden-section points)
TOL_FVALS = dict(rtol=1e-6, atol=1e-8)  # ADMM-dependent CV loss
TOL_FINAL = dict(rtol=1e-6, atol=1e-8)  # var_w / cor_w / p_vals
TOL_LAMBDA = dict(rtol=1e-12, atol=1e-15)

KNOWN_DATASETS = {
    "QMP",
    "American Gut Project Round 1",
    "GlobalPatterns / Caporaso et al. 2011",
    "Synthetic logistic-normal fractions",
}
PUBLIC_DATASETS = {
    "American Gut Project Round 1",
    "GlobalPatterns / Caporaso et al. 2011",
}
PROVENANCE_FIELDS = (
    "source",
    "source_package",
    "source_package_version",
    "source_license",
    "source_citation",
    "retrieval_date",
    "redistribution_cleared",
)


def _vec(path: Path) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",", missing_values="NA", filling_values=np.nan).reshape(
        -1
    )


def _mat(path: Path) -> np.ndarray:
    return np.genfromtxt(
        path,
        delimiter=",",
        missing_values="NA",
        filling_values=np.nan,
        ndmin=2,
    )


def _load(case: Path) -> dict:
    params = json.loads((case / "params.json").read_text())
    return {
        "params": params,
        "input": _mat(case / "input.csv"),
        "boot": _mat(case / "boot_indices.csv").astype(np.intp) - 1,  # R 1-based -> 0-based
        "null": _mat(case / "null_norm.csv"),
        "x_log": _mat(case / "x_log.csv"),
        "vx2": _mat(case / "vx2.csv"),
        "wd": _vec(case / "wd.csv"),
        "lams": _vec(case / "lams.csv"),
        "fvals": _vec(case / "fvals.csv"),
        "info_cv": json.loads((case / "info_cv.json").read_text()),
        "lambda": json.loads((case / "lambda.json").read_text()),
        "var_w": _vec(case / "var_w.csv"),
        "cor_w": _mat(case / "cor_w.csv"),
        "p_vals": _mat(case / "p_vals.csv"),
        "admm_iters": _vec(case / "admm_iters.csv").astype(np.intp),
        "admm_relative_errors": _vec(case / "admm_relative_errors.csv"),
    }


def _first_divergence(steps: list[tuple[str, np.ndarray, np.ndarray, dict]]) -> str | None:
    """Return a description of the first step that fails its tolerance, else None."""
    for name, got, exp, tol in steps:
        got = np.asarray(got, dtype=float)
        exp = np.asarray(exp, dtype=float)
        if got.shape != exp.shape:
            return f"{name}: shape {got.shape} != {exp.shape}"
        if not np.allclose(got, exp, equal_nan=True, **tol):
            diff = np.nanmax(np.abs(got - exp)) if got.size else 0.0
            return f"{name}: max|diff|={diff:.3e} (tol={tol})"
    return None


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_equivalence(case: Path):
    g = _load(case)
    p = g["params"]
    assert p["dataset"] in KNOWN_DATASETS

    # ---- Step 1: data transform (pseudo-count + log) ----
    x = g["input"]
    xt = x + p["pseudo"] if p["counts"] else x
    if p["counts"]:
        xt = xt / xt.sum(axis=1, keepdims=True)
    xt = np.log(xt)
    assert np.allclose(xt, g["x_log"], **TOL_INTERMEDIATE), "step 1: log-transform"

    # ---- Step 2: run the full pipeline with injected random draws ----
    res = cclasso(
        x,
        counts=p["counts"],
        pseudo=p["pseudo"],
        k_cv=p["k_cv"],
        lam_int=tuple(p["lam_int"]),
        k_max=p["k_max"],
        n_boot=p["n_boot"],
        boot_indices=g["boot"],
        null_norm=g["null"],
    )
    # Non-finite outputs are permitted only when R has the same non-finite
    # status at the same positions. Such estimates remain limitations, not
    # valid scientific results.
    for key in ("var_w", "cor_w", "p_vals"):
        assert np.array_equal(~np.isfinite(res[key]), ~np.isfinite(g[key])), (
            f"{case.name}: {key} non-finite mask differs from R"
        )

    # ---- Step 3: golden-section path (must be exact-ish; pure arithmetic) ----
    lams = res["info_cv"]["lams"]
    fvals = res["info_cv"]["fvals"]

    assert res["info_cv"]["k"] == g["info_cv"]["k"], (
        f"step 3: number of CV evaluations {res['info_cv']['k']} != "
        f"{g['info_cv']['k']} (golden-section path diverged)"
    )
    assert lams.shape == g["lams"].shape, (
        f"step 3: lambda path length {lams.shape} != {g['lams'].shape}"
    )

    steps = [
        ("step 3 lams", lams, g["lams"], TOL_LAMS),
        ("step 4 fvals", fvals, g["fvals"], TOL_FVALS),
        ("step 5 lambda", res["lambda"], np.array(g["lambda"]["lambda"]), TOL_LAMBDA),
        ("step 6 var_w", res["var_w"], g["var_w"], TOL_FINAL),
        ("step 7 cor_w", res["cor_w"], g["cor_w"], TOL_FINAL),
        ("step 8 p_vals", res["p_vals"], g["p_vals"], TOL_FINAL),
    ]
    msg = _first_divergence(steps)
    assert msg is None, f"{case.name}: first divergence at {msg}"

    # ---- Step 9: ADMM iteration counts (exact) ----
    iters = np.asarray(res["_admm_iters"], dtype=np.intp)
    assert iters.shape == g["admm_iters"].shape, (
        f"{case.name}: step 9 ADMM call count {iters.shape} != {g['admm_iters'].shape}"
    )
    if not np.array_equal(iters, g["admm_iters"]):
        idx = int(np.argmax(iters != g["admm_iters"]))
        raise AssertionError(
            f"{case.name}: step 9 ADMM iterations differ first at call {idx}: "
            f"got {iters[idx]}, expected {g['admm_iters'][idx]}"
        )

    # Compare the terminal relative residual from every ADMM solve, including
    # solves that hit k_max. This is diagnostic-only data; it does not change
    # the solver's control flow or relax output parity.
    residuals = np.asarray(res["_admm_relative_errors"], dtype=float)
    expected_residuals = g["admm_relative_errors"]
    assert residuals.shape == expected_residuals.shape, (
        f"step 10: ADMM residual count {residuals.shape} != {expected_residuals.shape}"
    )
    np.testing.assert_allclose(
        residuals,
        expected_residuals,
        rtol=1e-10,
        atol=1e-13,
        err_msg=f"{case.name}: terminal ADMM relative residuals",
    )


def test_golden_matrix_covers_counts_false_and_parameter_variation():
    """The committed R goldens must exercise the full parameter surface.

    The original fixture set used counts=True with default k_cv/n_boot/lam_int
    everywhere, so several branches were only self-tested. This asserts the
    coverage gap stays closed.
    """
    params = [
        json.loads((case / "params.json").read_text())
        for case in CASES
        if (case / "params.json").is_file()
    ]
    assert params, "no golden fixtures found"

    assert any(not p["counts"] for p in params), "no counts=False golden case"
    assert any(p["counts"] for p in params), "no counts=True golden case"

    assert len({p["k_cv"] for p in params}) >= 2, "k_cv is not varied"
    assert len({p["n_boot"] for p in params}) >= 2, "n_boot is not varied"
    assert len({tuple(p["lam_int"]) for p in params}) >= 2, "lam_int is not varied"

    defaults = (3, 3, (0.0001, 1.0))
    assert any((p["k_cv"], p["n_boot"], tuple(p["lam_int"])) != defaults for p in params), (
        "no case departs from the default k_cv/n_boot/lam_int triple"
    )


def test_public_fixtures_record_provenance():
    """Public-data cases must carry source, package/version, license and a
    redistribution-clearance flag. Terms are audited in
    ``docs/DATASET_PROVENANCE.md``; fixtures with ``redistribution_cleared``
    false are not part of a distributable artifact until cleared.
    """
    seen = set()
    for case in CASES:
        params = json.loads((case / "params.json").read_text())
        dataset = params["dataset"]
        if dataset not in PUBLIC_DATASETS:
            continue
        seen.add(dataset)
        for field in PROVENANCE_FIELDS:
            assert field in params, f"{case.name}: missing provenance field {field}"
            assert params[field] not in (None, ""), f"{case.name}: empty provenance field {field}"
        assert isinstance(params["redistribution_cleared"], bool)

    assert seen == PUBLIC_DATASETS, f"public fixtures missing for: {PUBLIC_DATASETS - seen}"


@pytest.mark.parametrize("case", FAILURE_CASES, ids=lambda path: path.name)
def test_r_failure_status_matches_python(case: Path):
    """A reference-side convergence error must not become a Python estimate."""
    params = json.loads((case / "params.json").read_text())
    assert params["expected_status"] == "error"
    x = _mat(case / "input.csv")
    boot = _mat(case / "boot_indices.csv").astype(np.intp) - 1
    null = _mat(case / "null_norm.csv")

    with pytest.raises(ValueError, match="missing value where TRUE/FALSE needed"):
        cclasso(
            x,
            counts=params["counts"],
            pseudo=params["pseudo"],
            k_cv=params["k_cv"],
            lam_int=tuple(params["lam_int"]),
            k_max=params["k_max"],
            n_boot=params["n_boot"],
            boot_indices=boot,
            null_norm=null,
        )
