"""Unit tests for the MB solver, nearPD, pulsar and validation."""

from __future__ import annotations

import numpy as np
import pytest

from spring_microbiome.huge import huge_mb, symmetrize_path
from spring_microbiome.latentcor.nearpd import near_pd
from spring_microbiome.pulsar import pulsar, ratcheck, stars_stability, starsind
from spring_microbiome.validation import validate_counts


# ---------------------------------------------------------------- nearPD -----
def test_near_pd_returns_valid_correlation():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(6, 6))
    approx = (a + a.T) / 2
    d = np.sqrt(np.abs(np.diag(approx)) + 1.0)
    approx = approx / np.outer(d, d)
    out = near_pd(approx, corr=True)
    np.testing.assert_allclose(np.diag(out), 1.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(out, out.T, rtol=0, atol=1e-14)
    assert np.linalg.eigvalsh(out).min() > -1e-8


def test_near_pd_identity_is_fixed_point():
    out = near_pd(np.eye(4), corr=True)
    np.testing.assert_allclose(out, np.eye(4), rtol=0, atol=1e-12)


# ------------------------------------------------------------------- MB ------
def test_huge_mb_path_is_symmetric_and_sparser_with_lambda():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(8, 8))
    s = (a @ a.T) / 8
    d = np.sqrt(np.diag(s))
    s = s / np.outer(d, d)
    lam = np.array([0.6, 0.4, 0.2])
    beta = huge_mb(s, lam)
    path = symmetrize_path(beta, "or")
    assert beta.shape == (3, 8, 8)
    for i in range(3):
        np.testing.assert_allclose(path[i], path[i].T, rtol=0, atol=0)
    assert path[0].sum() <= path[1].sum() <= path[2].sum()


def test_huge_mb_rejects_increasing_lambda():
    with pytest.raises(ValueError, match="non-increasing"):
        huge_mb(np.eye(3), np.array([0.1, 0.2]))


def test_huge_mb_rejects_non_positive_lambda():
    with pytest.raises(ValueError, match="strictly positive"):
        huge_mb(np.eye(3), np.array([0.1, 0.0]))


def test_huge_mb_empty_is_all_zero():
    s = np.eye(4)
    s[0, 1] = s[1, 0] = 0.3
    beta = huge_mb(s, np.array([10.0]))
    assert np.all(beta == 0)


# --------------------------------------------------------------- pulsar ------
def test_ratcheck_defaults():
    assert ratcheck(None, 100) == 0.8
    assert ratcheck(None, 200) == pytest.approx(10 * np.sqrt(200) / 200)
    assert ratcheck(0.5, 100) == 0.5


def test_starsind_first_threshold_minus_one():
    assert starsind(np.array([0.0, 0.0, 0.2, 0.3]), 0.1) == 2
    assert starsind(np.array([0.2, 0.3]), 0.1) == 1


def test_starsind_warns_and_returns_length():
    with pytest.warns(UserWarning):
        assert starsind(np.array([0.0, 0.0]), 0.1) == 2


def test_stars_stability_summary_and_cummax():
    p = 3
    empty = np.zeros((p, p))
    one = np.zeros((p, p))
    one[0, 1] = one[1, 0] = 1
    premerge = [[empty, empty], [one, empty]]
    out = stars_stability(premerge, 0.1, 2, p)
    # lambda 1: p=0 everywhere -> instability 0.
    # lambda 2: p=0.5 on the symmetric edge (0,1) and (1,0):
    #   sum p(1-p) = 2 * 0.25 = 0.5  ->  4*0.5/(p*(p-1)) = 1/3
    np.testing.assert_allclose(out.summary[0], 0.0, atol=1e-15)
    np.testing.assert_allclose(out.summary[1], 4 * 0.5 / (p * (p - 1)), atol=1e-15)
    assert out.opt_index == 1


def test_pulsar_injects_indices_deterministically():
    rng = np.random.default_rng(0)
    data = rng.normal(size=(10, 5))

    def fake(d, **kw):
        return {"path": [np.ones((d.shape[1], d.shape[1])), np.zeros((d.shape[1], d.shape[1]))]}

    idx = [np.array([0, 1, 2, 3, 4, 5, 6, 7]) for _ in range(4)]
    res1 = pulsar(data, fake, {"lambda": np.array([0.5, 0.1])}, rep_num=4, subsample_indices=idx)
    res2 = pulsar(data, fake, {"lambda": np.array([0.5, 0.1])}, rep_num=4, subsample_indices=idx)
    for a, b in zip(res1.stars.merge, res2.stars.merge):
        np.testing.assert_array_equal(a, b)


# --------------------------------------------------------------- bridge ------
def test_bridge22_is_deterministic_and_order_independent():
    """SciPy's normal CDF is randomized by default; the port must not leak that.

    Regression test for a non-determinism bug: ``multivariate_normal.cdf`` uses
    quasi-Monte-Carlo integration driven by an RNG, so without an explicit fixed
    generator the bridge value changes from call to call and the whole pipeline
    depends on unrelated call history.
    """
    from spring_microbiome.latentcor.bridge import bridge22, r_sol22

    a = bridge22(0.6, 0.125, 0.125)
    b = bridge22(0.6, 0.125, 0.125)
    # interleave unrelated evaluations, then repeat
    bridge22(0.91, 0.02, 0.8)
    c = bridge22(0.6, 0.125, 0.125)
    assert a == b == c

    v1 = r_sol22(np.array([0.8916666666666666]), 0.125, 0.125, 1e-6)[0]
    r_sol22(np.array([0.3]), 0.2, 0.4, 1e-6)
    v2 = r_sol22(np.array([0.8916666666666666]), 0.125, 0.125, 1e-6)[0]
    assert v1 == v2


def test_latentcor_original_is_deterministic():
    from spring_microbiome.latentcor import latentcor

    rng = np.random.default_rng(0)
    x = rng.normal(size=(25, 4))
    x[x < -0.5] = 0.0
    a = latentcor(x, method="original", use_nearPD=False).Rpointwise
    latentcor(rng.normal(size=(25, 4)), method="original", use_nearPD=False)
    b = latentcor(x, method="original", use_nearPD=False).Rpointwise
    np.testing.assert_array_equal(a, b)


def test_latentcor_without_nearpd_returns_pointwise_matrix():
    from spring_microbiome.latentcor import latentcor

    rng = np.random.default_rng(8)
    x = rng.poisson(3.0, size=(32, 6)).astype(float)
    result = latentcor(x, method="approx", use_nearPD=False)
    np.testing.assert_array_equal(result.R, result.Rpointwise)
    assert np.isfinite(result.R).all()


def test_latentcor_rejects_unimplemented_variable_types():
    from spring_microbiome.latentcor import latentcor

    with pytest.raises(NotImplementedError, match="only types='tru'"):
        latentcor(np.ones((5, 3)), types="con")


def test_brent_fmin_matches_known_minimum():
    from spring_microbiome.latentcor.bridge import brent_fmin

    x = brent_fmin(-2.0, 5.0, lambda t: (t - 1.25) ** 2, 1e-8)
    assert x == pytest.approx(1.25, abs=1e-6)


# ------------------------------------------------------------ validation -----
def test_validate_counts_rejects_negative():
    with pytest.raises(ValueError, match="Negative"):
        validate_counts(np.array([[1.0, -1.0], [2.0, 3.0]]))


def test_validate_counts_rejects_nan():
    with pytest.raises(ValueError, match="NaN"):
        validate_counts(np.array([[1.0, np.nan], [2.0, 3.0]]))


def test_validate_counts_rejects_all_zero_column():
    with pytest.raises(ValueError, match="all zero"):
        validate_counts(np.array([[0.0, 1.0], [0.0, 2.0]]))


def test_validate_counts_rejects_all_zero_row():
    with pytest.raises(ValueError, match="zero total"):
        validate_counts(np.array([[0.0, 0.0], [1.0, 2.0]]))


def test_validate_counts_rejects_non_numeric():
    with pytest.raises(TypeError):
        validate_counts(np.array([["a", "b"], ["c", "d"]]))


def test_validate_counts_rejects_too_few_samples():
    with pytest.raises(ValueError, match="at least 2 samples"):
        validate_counts(np.array([[1.0, 2.0]]))


def test_validate_counts_warns_when_p_ge_n():
    with pytest.warns(UserWarning, match="taxa"):
        validate_counts(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
