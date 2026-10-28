"""Unit tests for Kendall tau and the interpolation table."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from spring_microbiome.latentcor.interp import load_ipol22
from spring_microbiome.latentcor.kendall import kendall_tau_b, tie_count


def test_tie_count():
    assert tie_count(np.array([1.0, 2.0, 3.0])) == 0
    assert tie_count(np.array([1.0, 1.0, 2.0])) == 1
    assert tie_count(np.array([1.0, 1.0, 1.0, 2.0])) == 3


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_tau_b_matches_scipy_no_ties(seed):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=40)
    y = rng.normal(size=40)
    assert kendall_tau_b(x, y) == pytest.approx(
        stats.kendalltau(x, y).statistic, rel=1e-12, abs=1e-14
    )


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_tau_b_matches_scipy_with_ties(seed):
    rng = np.random.default_rng(seed)
    x = np.round(rng.normal(size=40), 1)
    y = np.round(rng.normal(size=40), 1)
    assert kendall_tau_b(x, y) == pytest.approx(
        stats.kendalltau(x, y).statistic, rel=1e-12, abs=1e-14
    )


def test_tau_b_is_symmetric():
    rng = np.random.default_rng(3)
    x = np.round(rng.normal(size=30), 1)
    y = np.round(rng.normal(size=30), 1)
    assert kendall_tau_b(x, y) == pytest.approx(kendall_tau_b(y, x), rel=1e-15)


def test_ipol22_exact_on_grid_nodes():
    ipol = load_ipol22()
    gk, gz1, gz2 = ipol.grids
    values = ipol.values
    for i in (0, 5, 20, 40):
        for j in (0, 10, 48):
            for k in (0, 24, 48):
                pt = np.array([[gk[i]], [gz1[j]], [gz2[k]]])
                got = ipol(pt)[0] / ipol.scale
                want = values[i, j, k] / ipol.scale
                assert got == pytest.approx(want, rel=1e-12, abs=1e-12)


def test_ipol22_linear_along_first_axis():
    ipol = load_ipol22()
    gk, gz1, gz2 = ipol.grids
    # midpoint between two grid nodes on the K axis: linear interpolation
    a, b = 5, 6
    mid = 0.5 * (gk[a] + gk[b])
    j, k = 20, 20
    got = ipol(np.array([[mid], [gz1[j]], [gz2[k]]]))[0]
    want = 0.5 * (ipol.values[a, j, k] + ipol.values[b, j, k])
    assert got == pytest.approx(want, rel=1e-12, abs=1e-9)


def test_ipol22_interpolates_inside_a_cell():
    ipol = load_ipol22()
    gk, gz1, gz2 = ipol.grids
    a, b = 5, 6
    for s in (0.25, 0.5, 0.75):
        x0 = gk[a] + s * (gk[b] - gk[a])
        got = ipol(np.array([[x0], [gz1[24]], [gz2[24]]]))[0]
        want = (1 - s) * ipol.values[a, 24, 24] + s * ipol.values[b, 24, 24]
        assert got == pytest.approx(want, rel=1e-10, abs=1e-6)


def test_ipol22_extrapolates_linearly_outside():
    ipol = load_ipol22()
    gk, gz1, gz2 = ipol.grids
    # above the last node: uses (n-2, n-1) with weight > 1
    x0 = gk[-1] + 0.5 * (gk[-1] - gk[-2])
    got = ipol(np.array([[x0], [gz1[24]], [gz2[24]]]))[0]
    w = (x0 - gk[-2]) / (gk[-1] - gk[-2])
    want = (1 - w) * ipol.values[-2, 24, 24] + w * ipol.values[-1, 24, 24]
    assert got == pytest.approx(want, rel=1e-10, abs=1e-6)

    # below the first node: binsearch clamps to index 1, weight < 0
    x1 = gk[0] - 0.5 * (gk[1] - gk[0])
    got = ipol(np.array([[x1], [gz1[24]], [gz2[24]]]))[0]
    w = (x1 - gk[0]) / (gk[1] - gk[0])
    want = (1 - w) * ipol.values[0, 24, 24] + w * ipol.values[1, 24, 24]
    assert got == pytest.approx(want, rel=1e-10, abs=1e-6)
