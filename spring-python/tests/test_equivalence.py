"""R-equivalence tests against the committed golden fixtures.

These tests replay the exact R computation step by step and fail with a pointer
to the *first divergent step*. Tolerances are deliberately tight; they must not
be widened to make a test pass (see docs/EQUIVALENCE.md).

Run with:  pytest tests/test_equivalence.py -m slow
"""

from __future__ import annotations

import numpy as np
import pytest

from spring_microbiome import spring
from spring_microbiome.huge import huge_mb, symmetrize_path
from spring_microbiome.latentcor import latentcor
from spring_microbiome.latentcor.bridge import bound22
from spring_microbiome.mclr import mclr

from .conftest import CASES, LATENTCOR_CASES, load_case, load_latentcor_case

pytestmark = pytest.mark.slow

TOL_STRICT = dict(rtol=1e-10, atol=1e-12)
TOL_FINAL = dict(rtol=1e-6, atol=1e-8)

# ---------------------------------------------------------------------------
# Why a second, looser tolerance for the latent correlation exists
# ---------------------------------------------------------------------------
# For pairs with |K| > ratio * bound, `latentcor` does NOT interpolate: it
# inverts the bridge function numerically, which needs the 4-variate normal CDF.
# R uses `mnormt::pmnorm`, Python uses `scipy.stats.multivariate_normal.cdf`.
# On the (mildly ill-conditioned) 4x4 bridge matrices that occur here the two
# CDFs disagree by up to ~8e-5 -- mnormt's value is stable from abseps=1e-6 to
# 1e-12, so this is an algorithmic difference, not a tolerance setting we can
# match. Because the bridge objective is very flat near its root, that CDF
# difference shifts the inverted correlation by up to ~3e-5 (up to ~2e-4
# relative). We therefore:
#   * require the STRICT tolerance on every pair that is interpolated
#     (i.e. |K| <= ratio*bound) -- this is a sharp, cause-based check;
#   * allow the documented tolerance below only on the pairs that go through the
#     numerical inversion;
#   * `R` (after nearPD) is compared with the documented tolerance globally,
#     because nearPD mixes all entries of the matrix.
# Crucially, the *discrete* outputs (opt.index, edge sets, final adjacency) are
# still matched EXACTLY in every case. See docs/EQUIVALENCE.md section 2.
TOL_TRANSFORM = dict(rtol=1e-3, atol=1e-4)
# MB coefficients for Rmethod="original" inherit the lambda-grid difference.
TOL_FINAL_ORIGINAL = dict(rtol=1e-4, atol=1e-6)
# A data-specific lambda grid inherits the documented bridge-inversion error in
# latentcor. The QMP p>n fixture measures a 4.9e-8 relative grid deviation.
TOL_DATA_SPECIFIC_GRID = dict(rtol=1e-7, atol=1e-9)
# When the *full-data* latentcor contains a bridge-inversion pair, lambda-max is
# one of its entries, so the documented ~2e-4 bridge-CDF difference propagates
# into the grid. GlobalPatterns measures 1.27e-6 relative; keep a small margin.
TOL_DATA_SPECIFIC_GRID_BRIDGE = dict(rtol=5e-6, atol=1e-8)

_PIPELINE_RESULTS = {}


def _beta_tol(params):
    return TOL_FINAL_ORIGINAL if params["Rmethod"] == "original" else TOL_FINAL


def _max_kkt_violation(S, beta, lam):
    """Largest active/inactive KKT residual for all nodewise Lasso fits."""
    active_max = 0.0
    inactive_max = 0.0
    p = S.shape[0]
    for m in range(p):
        w = np.asarray(beta[:, m], dtype=float).copy()
        w[m] = 0.0
        residual = S @ w - S[:, m]
        active = np.abs(w) > 0.0
        active[m] = False
        inactive = ~active
        inactive[m] = False  # the response coordinate is constrained to zero
        if np.any(active):
            active_max = max(
                active_max,
                float(np.max(np.abs(residual[active] + lam * np.sign(w[active])))),
            )
        if np.any(inactive):
            inactive_max = max(
                inactive_max,
                float(np.max(np.maximum(np.abs(residual[inactive]) - lam, 0.0))),
            )
    return active_max, inactive_max


def _max_scaled_objective_gap(S_py, beta_py, lam_py, S_r, beta_r, lam_r):
    """Largest nodewise objective difference, scaled by max(1, |R objective|)."""
    gaps = []
    for node in range(S_py.shape[0]):
        wp = np.asarray(beta_py[:, node], dtype=float).copy()
        wr = np.asarray(beta_r[:, node], dtype=float).copy()
        wp[node] = wr[node] = 0.0
        obj_py = 0.5 * wp @ S_py @ wp - S_py[:, node] @ wp + lam_py * np.abs(wp).sum()
        obj_r = 0.5 * wr @ S_r @ wr - S_r[:, node] @ wr + lam_r * np.abs(wr).sum()
        gaps.append(abs(float(obj_py - obj_r)) / max(1.0, abs(float(obj_r))))
    return max(gaps, default=0.0)


def _zratios(data):
    """zratios for `types="tru"`: proportion of entries equal to the column min."""
    data = np.asarray(data, dtype=float)
    return (data == data.min(axis=0)).mean(axis=0)


def _bridge_solver_mask(data, params, K):
    """Off-diagonal pairs that go through the numerical bridge inversion."""
    p = K.shape[0]
    off = ~np.eye(p, dtype=bool)
    if params["Rmethod"] == "original":
        return off
    zr = _zratios(data)
    rows, cols = np.tril_indices(p, -1)
    bound = bound22(zr[rows], zr[cols])
    cut = np.abs(K[rows, cols]) > params["ratio"] * bound
    mask = np.zeros((p, p), dtype=bool)
    mask[rows[cut], cols[cut]] = True
    return (mask | mask.T) & off


def _assert_transform(res, gold, data, params, ctx):
    """Compare latentcor outputs, splitting by the bridge-solver region."""
    np.testing.assert_allclose(res.K, gold["K"], err_msg=f"{ctx} K", **TOL_STRICT)
    mask = _bridge_solver_mask(data, params, res.K)
    strict = ~mask & ~np.eye(res.K.shape[0], dtype=bool)
    if strict.any():
        np.testing.assert_allclose(
            res.Rpointwise[strict],
            gold["Rpointwise"][strict],
            err_msg=f"{ctx} Rpointwise (interpolated pairs)",
            **TOL_STRICT,
        )
    if mask.any():
        np.testing.assert_allclose(
            res.Rpointwise[mask],
            gold["Rpointwise"][mask],
            err_msg=f"{ctx} Rpointwise (bridge-inversion pairs)",
            **TOL_TRANSFORM,
        )
    np.testing.assert_allclose(res.R, gold["R"], err_msg=f"{ctx} R (after nearPD)", **TOL_TRANSFORM)


def _p(params, key):
    return params[key]


def _pipeline_kwargs(params, indices):
    lam_spec = params["lambdaseq_spec"]
    lam = (
        np.asarray(params["lambdaseq"], dtype=float)
        if lam_spec == "injected_explicit"
        else None
        if lam_spec == "explicit"
        else "data-specific"
    )
    return dict(
        quantitative=params["quantitative"],
        nlambda=params["nlambda"],
        lambdaseq=lam,
        seed=params["seed"],
        thresh=params["thresh"],
        subsample_ratio=params["subsample.ratio"],
        rep_num=params["rep.num"],
        Rtol=params["Rtol"],
        verbose=False,
        Rmethod=params["Rmethod"],
        use_nearPD=params["use.nearPD"],
        nu=params["nu"],
        ratio=params["ratio"],
        subsample_indices=indices,
        validate_input=False,
    )


def _spring_result(name, case):
    """Reuse the expensive full StARS result across its two assertions."""
    if name not in _PIPELINE_RESULTS:
        params = case["params"]
        indices = [np.asarray(ix) for ix in case["subsample_indices"]]
        _PIPELINE_RESULTS[name] = spring(
            case["input"], **_pipeline_kwargs(params, indices)
        )
    return _PIPELINE_RESULTS[name]


@pytest.mark.parametrize("name", CASES)
def test_step1_mclr(name):
    case = load_case(name)
    got = mclr(case["input"]) if not case["params"]["quantitative"] else case["input"]
    assert case["params"]["dataset"] in {
        "QMP",
        "American Gut Project Round 1",
        "GlobalPatterns / Caporaso et al. 2011",
        "synthetic",
    }
    np.testing.assert_allclose(got, case["qdat"], **TOL_STRICT)


@pytest.mark.parametrize("name", CASES)
def test_step2_kendall_and_latent_correlation(name):
    case = load_case(name)
    p = case["params"]
    res = latentcor(
        case["qdat"],
        types="tru",
        method=p["Rmethod"],
        use_nearPD=p["use.nearPD"],
        nu=p["nu"],
        tol=p["Rtol"],
        ratio=p["ratio"],
    )
    _assert_transform(res, case, case["qdat"], p, f"{name} full")


@pytest.mark.parametrize("name", LATENTCOR_CASES)
def test_step2b_latentcor_original_matches_r(name):
    """Direct R comparison of ``Rmethod="original"`` (numeric bridge inversion).

    A full StARS pipeline is impractical here: every pair is inverted with Brent
    over the 4-D normal CDF (~20 s/pair). This dedicated fixture compares the
    latent-correlation stage itself, which is exactly the divergent branch.
    """
    case = load_latentcor_case(name)
    p = case["params"]
    res = latentcor(
        case["qdat"],
        types="tru",
        method=p["Rmethod"],
        use_nearPD=p["use.nearPD"],
        nu=p["nu"],
        tol=p["Rtol"],
        ratio=p["ratio"],
    )
    np.testing.assert_allclose(res.K, case["K"], err_msg=f"{name} K", **TOL_STRICT)
    np.testing.assert_allclose(
        res.Rpointwise,
        case["Rpointwise"],
        err_msg=f"{name} Rpointwise (original)",
        **TOL_TRANSFORM,
    )
    np.testing.assert_allclose(
        res.R, case["R"], err_msg=f"{name} R (original, after nearPD)", **TOL_TRANSFORM
    )


@pytest.mark.parametrize("name", CASES)
def test_step3_subsample_latent_correlation(name):
    case = load_case(name)
    p = case["params"]
    for k, rows in enumerate(case["subsample_indices"]):
        sub = case["qdat"][rows]
        res = latentcor(
            sub,
            method=p["Rmethod"],
            use_nearPD=p["use.nearPD"],
            nu=p["nu"],
            tol=p["Rtol"],
            ratio=p["ratio"],
        )
        np.testing.assert_allclose(
            res.R,
            case["subsample_R"][k],
            err_msg=f"{name} subsample {k} R",
            **TOL_TRANSFORM,
        )


@pytest.mark.parametrize("name", CASES)
def test_step4_mb_paths_per_subsample(name):
    """Exact edge-set agreement and tight coefficient agreement per subsample."""
    case = load_case(name)
    lam = case["lambdaseq"]
    d = case["qdat"].shape[1]

    # rebuild the golden long-format structures
    beta_gold = _long_to_set(case["dir"] / "subsample_beta_long.csv", 5)
    path_gold = _long_to_set(case["dir"] / "subsample_paths_long.csv", 4)

    for k in range(1, len(case["subsample_indices"]) + 1):
        sub_r_matrix = case["subsample_R"][k - 1]
        beta = huge_mb(sub_r_matrix, lam)
        path = symmetrize_path(beta, "or")

        # edge sets (unsigned adjacency) must be identical
        py_edges = {
            (ell, i, j)
            for ell in range(len(lam))
            for i in range(d)
            for j in range(d)
            if path[ell, i, j] != 0
        }
        g_edges = {(ell - 1, i, j) for (kk, ell, i, j) in path_gold if kk == k}
        assert py_edges == g_edges, (
            f"{name} subsample {k}: edge sets differ "
            f"(py-only {sorted(py_edges - g_edges)[:5]}, "
            f"gold-only {sorted(g_edges - py_edges)[:5]})"
        )

        # coefficients
        for kk, ell, i, j, v in beta_gold:
            if kk != k:
                continue
            np.testing.assert_allclose(
                beta[ell - 1, i, j], v, err_msg=f"{name} beta subsample {k}", **TOL_FINAL
            )


@pytest.mark.parametrize("name", CASES)
def test_step5_stars_and_final_adjacency(name):
    case = load_case(name)
    p = case["params"]
    result = _spring_result(name, case)

    if p["Rmethod"] == "original":
        lam_tol = dict(rtol=1e-4, atol=1e-9)
    elif p["lambdaseq_spec"] == "data-specific":
        # lambda-max is an entry of the latent correlation; if that matrix needs
        # numeric bridge inversion, the grid inherits the bridge-CDF difference.
        bridge = _bridge_solver_mask(case["qdat"], p, case["K"])
        lam_tol = (
            TOL_DATA_SPECIFIC_GRID_BRIDGE if bridge.any() else TOL_DATA_SPECIFIC_GRID
        )
    else:
        lam_tol = TOL_STRICT
    np.testing.assert_allclose(result.lambdaseq, case["lambdaseq"], err_msg="lambdaseq", **lam_tol)

    stars = result.output.stars
    assert stars.opt_index == case["opt_index"], (
        f"{name}: opt.index {stars.opt_index} != R {case['opt_index']}"
    )
    np.testing.assert_allclose(
        stars.summary, case["summary"], err_msg="instability curve", **TOL_STRICT
    )
    for i, merge in enumerate(stars.merge):
        np.testing.assert_allclose(
            merge,
            case["merge"][i],
            err_msg=f"selection frequency lambda {i + 1}",
            **TOL_STRICT,
        )

    # final adjacency must be exactly equal
    diff = np.argwhere(np.asarray(result.adjacency) != case["refit"])
    assert diff.size == 0, (
        f"{name}: final adjacency differs at {diff.tolist()[:5]} "
        f"(py={np.asarray(result.adjacency)[tuple(diff[0])] if diff.size else None}, "
        f"R={case['refit'][tuple(diff[0])] if diff.size else None})"
    )


@pytest.mark.parametrize("name", CASES)
def test_step6_beta_weights_at_opt(name):
    case = load_case(name)
    p = case["params"]
    result = _spring_result(name, case)
    opt = case["opt_index"]
    got = np.asarray(result.beta.beta[opt - 1])
    inversion_mask = _bridge_solver_mask(case["input"], p, case["K"])
    if p["n"] < p["p"] or inversion_mask.any():
        # Rank-deficient p>n designs can have non-unique coefficients. On any
        # case that also uses numeric bridge inversion, the R/Python covariance
        # inputs differ at the measured bridge-CDF level, which can amplify into
        # beta differences even when both solve their own objective correctly.
        # Compare exact support, each solver's KKT residual, and objective gap.
        ref = case["beta_full"][opt - 1]
        np.testing.assert_array_equal(got != 0.0, ref != 0.0)
        S_py = latentcor(
            case["qdat"],
            method=p["Rmethod"],
            use_nearPD=p["use.nearPD"],
            nu=p["nu"],
            tol=p["Rtol"],
            ratio=p["ratio"],
        ).R
        lam_py = float(result.lambdaseq[opt - 1])
        lam_r = float(case["lambdaseq"][opt - 1])
        py_kkt = _max_kkt_violation(S_py, got, lam_py)
        r_kkt = _max_kkt_violation(case["R"], ref, lam_r)
        assert max(py_kkt) <= 2e-4, f"Python KKT residuals {py_kkt} exceed huge's thol"
        assert max(r_kkt) <= 2e-4, f"R KKT residuals {r_kkt} exceed huge's thol"
        objective_gap = _max_scaled_objective_gap(
            S_py, got, lam_py, case["R"], ref, lam_r
        )
        assert objective_gap <= 1e-6, f"relative objective gap {objective_gap:.3e}"
        return
    np.testing.assert_allclose(got, case["beta_full"][opt - 1], err_msg="beta", **_beta_tol(p))


def _long_to_set(path, ncols):
    if not path.exists():
        return []
    raw = np.genfromtxt(path, delimiter=",", skip_header=1, ndmin=2)
    if raw.size == 0:
        return []
    if ncols == 4:
        return [tuple(int(v) for v in row) for row in raw]
    return [(*[int(v) for v in row[:4]], float(row[4])) for row in raw]
