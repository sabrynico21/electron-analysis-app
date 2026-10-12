from __future__ import annotations

import numpy as np
import pytest

import spieceasi.stars as stars_module
from spieceasi.neighborhood import make_alpha_path
from spieceasi.stars import _instability_from_adjacencies, _stars_index, stars_select
from spieceasi.transforms import clr_transform
from tests.synthetic_data import generate_compositional_counts


def test_stars_returns_alpha_from_path() -> None:
    counts, _ = generate_compositional_counts(n_samples=70, n_taxa=7)
    X = clr_transform(counts).clr
    alpha_path = make_alpha_path(X, n_lambdas=10)
    result = stars_select(
        X,
        alpha_path=alpha_path,
        n_subsamples=8,
        subsample_fraction=0.8,
        random_state=123,
        n_jobs=1,
    )
    assert result.selected_alpha in alpha_path
    assert result.instability.shape == alpha_path.shape
    assert np.all(np.diff(result.monotone_instability) >= -1e-12)


def test_stars_is_deterministic_with_fixed_seed() -> None:
    counts, _ = generate_compositional_counts(n_samples=90, n_taxa=9, seed=5)
    X = clr_transform(counts).clr
    alpha_path = make_alpha_path(X, n_lambdas=12)

    kwargs = dict(
        alpha_path=alpha_path,
        n_subsamples=12,
        subsample_fraction=0.8,
        random_state=77,
        n_jobs=1,
    )
    res1 = stars_select(X, **kwargs)
    res2 = stars_select(X, **kwargs)

    assert np.array_equal(res1.alpha_path, res2.alpha_path)
    assert np.allclose(res1.instability, res2.instability)
    assert np.allclose(res1.monotone_instability, res2.monotone_instability)
    assert res1.selected_alpha == res2.selected_alpha


def test_stars_default_for_tiny_n_matches_ratio_rule() -> None:
    # With pulsar-compatible defaults and n<=144, subsample size is floor(0.8*n).
    X = np.array(
        [
            [0.1, -0.1],
            [0.0, 0.2],
            [-0.2, 0.0],
            [0.3, -0.1],
            [0.2, 0.4],
            [-0.1, -0.3],
        ],
        dtype=float,
    )
    alpha_path = np.array([0.5, 0.1], dtype=float)
    res = stars_select(
        X, alpha_path=alpha_path, n_subsamples=5, random_state=3, n_jobs=1
    )
    assert res.selected_alpha in alpha_path


def test_instability_formula_matches_hand_computation() -> None:
    # Two-node graph: one possible undirected edge, selected in half subsamples.
    adjs = [
        np.array([[0, 1], [1, 0]], dtype=int),
        np.array([[0, 1], [1, 0]], dtype=int),
        np.array([[0, 0], [0, 0]], dtype=int),
        np.array([[0, 0], [0, 0]], dtype=int),
    ]
    instability = _instability_from_adjacencies(adjs)

    # pulsar::stars.stability uses 4*p*(1-p) per undirected edge;
    # at p=0.5 this equals 1.0.
    assert np.isclose(instability, 1.0, atol=1e-12)


def test_selection_rule_matches_pulsar_starsind() -> None:
    import warnings

    # Case 1: normal interior crossing (first summary >= thresh at index 2).
    m1 = np.array([0.01, 0.03, 0.07, 0.11])
    assert _stars_index(m1, 0.05) == 1

    # Case 2: all lambdas above threshold -> first (most regularized) lambda.
    m2 = np.array([0.06, 0.07, 0.09])
    assert _stars_index(m2, 0.05) == 0

    # Case 3: no lambda reaches the threshold -> last (least regularized),
    # with pulsar's warning.
    m3 = np.array([0.01, 0.02, 0.03])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert _stars_index(m3, 0.05) == 2
    assert any("Optimal lambda" in str(w.message) for w in caught)

    # Case 4: exact-equality boundary. pulsar uses `>= thresh`, so it returns
    # the index just before the crossing (the old `<= thresh` rule returned 2).
    m4 = np.array([0.01, 0.02, 0.05, 0.07, 0.10])
    assert _stars_index(m4, 0.05) == 1


def test_default_subsample_size_matches_pulsar_rule(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recorded_sizes: list[int] = []

    class DummyRng:
        def choice(
            self, n_samples: int, size: int, replace: bool = False
        ) -> np.ndarray:
            assert replace is False
            recorded_sizes.append(int(size))
            return np.arange(size, dtype=int)

    def fake_default_rng(_seed: int):
        return DummyRng()

    def fake_infer(*_args, **_kwargs) -> np.ndarray:
        return np.zeros((2, 4, 4), dtype=np.uint8)  # (n_lambda, d, d)

    monkeypatch.setattr(np.random, "default_rng", fake_default_rng)
    monkeypatch.setattr(
        stars_module, "_infer_subsample_adjacency_path", fake_infer
    )

    alpha_path = np.array([1.0, 0.5], dtype=float)

    X_small = np.zeros((100, 4), dtype=float)
    stars_select(
        X_small, alpha_path=alpha_path, n_subsamples=3, random_state=1, n_jobs=1
    )
    # pulsar::.ratcheck: n<=144 => subsample.ratio=0.8
    assert all(size == 80 for size in recorded_sizes[:3])

    X_large = np.zeros((225, 4), dtype=float)
    stars_select(
        X_large, alpha_path=alpha_path, n_subsamples=3, random_state=1, n_jobs=1
    )
    # pulsar::.ratcheck: n>144 => subsample size floor(10*sqrt(n)).
    assert all(size == 150 for size in recorded_sizes[3:6])
