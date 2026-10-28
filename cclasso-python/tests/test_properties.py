"""Mathematical-property and invariance tests for the cclasso port."""

from __future__ import annotations

import numpy as np
import pytest

from cclasso import cclasso


def _lognormal_fractions(n, p, seed, block=0.5):
    """Logistic-normal fractions whose *latent* (absolute) correlation is positive on
    the true edges (0,1),(1,2),(0,2),(3,4)."""
    rng = np.random.default_rng(seed)
    corr = np.eye(p)
    for i, j in [(0, 1), (1, 2), (0, 2), (3, 4)]:
        if j < p:
            corr[i, j] = corr[j, i] = block
    z = rng.standard_normal((n, p)) @ np.linalg.cholesky(corr).T
    return np.exp(z) / np.exp(z).sum(axis=1, keepdims=True)


def _counts_from_fractions(frac, seed, totals=(2000, 5000)):
    rng = np.random.default_rng(seed)
    tot = rng.integers(totals[0], totals[1] + 1, size=frac.shape[0])
    out = np.empty_like(frac)
    for i in range(frac.shape[0]):
        out[i] = rng.multinomial(tot[i], frac[i])
    return out


@pytest.fixture(scope="module")
def frac_result():
    x = _lognormal_fractions(40, 8, seed=1)
    return x, cclasso(x, k_cv=3, n_boot=8, seed=7)


def test_correlation_is_symmetric_unit_diagonal(frac_result):
    _, res = frac_result
    c = res["cor_w"]
    assert np.allclose(c, c.T, atol=1e-12)
    assert np.allclose(np.diag(c), 1.0, atol=1e-12)
    assert np.all(np.abs(c) <= 1.0 + 1e-9)


def test_correlation_is_psd(frac_result):
    _, res = frac_result
    eig = np.linalg.eigvalsh(res["cor_w"])
    assert eig.min() >= -1e-8


def test_pvalues_properties(frac_result):
    _, res = frac_result
    pv = res["p_vals"]
    assert np.allclose(pv, pv.T, atol=1e-12)
    assert np.all((pv >= 0.0) & (pv <= 1.0))
    assert np.allclose(np.diag(pv), 1.0)


def test_variances_positive(frac_result):
    _, res = frac_result
    assert np.all(res["var_w"] > 0)


def test_lambda_inside_or_bounded_interval(frac_result):
    _, res = frac_result
    lo, hi = res["info_cv"]["lam_int"]
    # final lambda is the midpoint of the inner points, within the final bracket
    assert lo - 1e-9 <= res["lambda"] <= hi + 1e-9


def test_scale_invariance_for_fractions():
    """log(c*x) = log(x) + const: variance-based pipeline is exactly scale-invariant."""
    x = _lognormal_fractions(30, 6, seed=2)
    r1 = cclasso(x, k_cv=3, n_boot=5, seed=3)
    r2 = cclasso(x * 1234.5, k_cv=3, n_boot=5, seed=3)
    assert np.allclose(r1["cor_w"], r2["cor_w"], rtol=1e-12, atol=1e-14)
    assert np.allclose(r1["var_w"], r2["var_w"], rtol=1e-10, atol=1e-12)


def test_counts_pipeline_equals_pre_normalised_fractions():
    """counts=True with pseudo=0 == feeding the row-normalised fractions directly."""
    rng = np.random.default_rng(5)
    counts = rng.integers(1, 50, size=(35, 6)).astype(float)  # strictly positive
    frac = counts / counts.sum(axis=1, keepdims=True)
    r_counts = cclasso(counts, counts=True, pseudo=0.0, k_cv=3, n_boot=5, seed=6)
    r_frac = cclasso(frac, counts=False, k_cv=3, n_boot=5, seed=6)
    assert np.allclose(r_counts["cor_w"], r_frac["cor_w"], rtol=1e-12, atol=1e-14)


def test_known_nonfinite_bootstrap_is_deterministic_and_reproduced_in_r():
    """Pin a known shared R/Python failure instead of treating equal NaNs as success.

    The matching R replay is available via tests/r_reference/replay_numpy_nan.R;
    its report shows the same negative bootstrap variance and non-finite outputs.
    """
    x = _lognormal_fractions(25, 5, seed=8)
    with pytest.warns(RuntimeWarning, match="invalid value encountered in sqrt"):
        a = cclasso(x, k_cv=3, n_boot=6, seed=42)
    with pytest.warns(RuntimeWarning, match="invalid value encountered in sqrt"):
        b = cclasso(x, k_cv=3, n_boot=6, seed=42)

    assert np.isfinite(a["var_w"]).all()
    assert np.isnan(a["cor_w"]).sum() == 8
    assert np.isnan(a["p_vals"]).sum() == 8
    negative_variances = np.argwhere(a["_bootstrap_variances"] < 0.0)
    np.testing.assert_array_equal(negative_variances, np.array([[0, 2]]))
    np.testing.assert_array_equal(a["_bootstrap_variances"], b["_bootstrap_variances"])
    assert np.array_equal(a["cor_w"], b["cor_w"], equal_nan=True)
    assert np.array_equal(a["p_vals"], b["p_vals"], equal_nan=True)
    assert a["lambda"] == b["lambda"]


def test_recovers_known_block_structure():
    x = _lognormal_fractions(120, 6, seed=11, block=0.6)
    res = cclasso(x, k_cv=3, n_boot=10, seed=12)
    c = res["cor_w"]
    # true edges (0,1),(1,2),(0,2),(3,4) should be positive; unrelated pairs smaller
    assert c[0, 1] > 0 and c[1, 2] > 0 and c[0, 2] > 0
    assert c[3, 4] > 0
    edge_mean = np.mean([c[0, 1], c[1, 2], c[0, 2], c[3, 4]])
    off_mean = np.mean([c[0, 5], c[1, 5], c[2, 5], c[3, 5]])
    assert edge_mean > off_mean


def test_minimal_two_part_composition_is_finite():
    """The minimum supported feature dimension must survive the full pipeline."""
    rng = np.random.default_rng(71)
    z = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.35], [0.35, 1.0]], size=36)
    x = np.exp(z)
    x /= x.sum(axis=1, keepdims=True)

    result = cclasso(x, k_cv=3, n_boot=3, seed=72)
    for key in ("var_w", "cor_w", "p_vals"):
        assert np.isfinite(result[key]).all(), f"non-finite {key} for p=2"
    assert result["cor_w"].shape == (2, 2)
    assert np.allclose(result["cor_w"], result["cor_w"].T, atol=1e-12)
    assert np.allclose(np.diag(result["cor_w"]), 1.0, atol=1e-12)


def test_extreme_common_count_scale_is_numerically_invariant():
    """Large but finite counts must not overflow the row normalization."""
    rng = np.random.default_rng(5)
    counts = rng.integers(1, 50, size=(35, 6)).astype(float)
    kwargs = dict(counts=True, pseudo=0.0, k_cv=3, n_boot=5, seed=6)
    base = cclasso(counts, **kwargs)
    scaled = cclasso(counts * 1e100, **kwargs)

    for key in ("var_w", "cor_w", "p_vals"):
        assert np.isfinite(scaled[key]).all(), f"non-finite {key} after scaling"
        np.testing.assert_allclose(scaled[key], base[key], rtol=1e-8, atol=1e-10)


def test_student_t_cdf_is_finite_in_both_tails():
    from cclasso._core import _student_t_cdf

    cdf = _student_t_cdf(np.array([-1e200, -10.0, 0.0, 10.0, 1e200]), df=3)
    assert np.isfinite(cdf).all()
    assert np.all((cdf >= 0.0) & (cdf <= 1.0))
    np.testing.assert_allclose(cdf + cdf[::-1], 1.0, atol=1e-14)
