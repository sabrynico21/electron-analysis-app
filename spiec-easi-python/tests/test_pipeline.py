from __future__ import annotations

import numpy as np

from spieceasi.pipeline import spiec_easi_mb
from tests.synthetic_data import (
    generate_compositional_counts,
    generate_strong_signal_counts,
)


def _upper_triangle_edges(A: np.ndarray) -> set[tuple[int, int]]:
    idx = np.triu_indices(A.shape[0], k=1)
    return {(int(i), int(j)) for i, j in zip(*idx, strict=True) if A[i, j] != 0}


def test_pipeline_recovers_reasonable_edges() -> None:
    # Edge recovery is only meaningful with a strong signal: on weak, high-variance
    # compositional data even the reference R SpiecEasi implementation recovers no
    # ground-truth edges, so this test uses the strong-signal generator.
    counts, ground_truth = generate_strong_signal_counts(
        n_samples=220, n_taxa=24, depth=6000, rho=0.45, seed=31
    )
    result = spiec_easi_mb(
        counts,
        n_lambdas=20,
        lambda_min_ratio=1e-2,
        n_subsamples=20,
        subsample_fraction=0.8,
        stars_threshold=0.05,
        sym_rule="or",
        random_state=31,
        n_jobs=1,
    )

    pred_edges = _upper_triangle_edges(result.adjacency)
    true_edges = _upper_triangle_edges(ground_truth)
    intersection = pred_edges.intersection(true_edges)

    precision = len(intersection) / len(pred_edges)
    recall = len(intersection) / len(true_edges)

    assert len(pred_edges) > 0
    assert precision >= 0.6
    assert recall >= 0.7


def test_pipeline_reproducible_with_seed() -> None:
    counts, _ = generate_compositional_counts(n_samples=80, n_taxa=8, seed=15)
    kwargs = dict(
        n_lambdas=10,
        n_subsamples=8,
        subsample_fraction=0.8,
        stars_threshold=0.1,
        random_state=11,
        n_jobs=1,
    )
    r1 = spiec_easi_mb(counts, **kwargs)
    r2 = spiec_easi_mb(counts, **kwargs)

    assert np.array_equal(r1.adjacency, r2.adjacency)
    assert np.allclose(r1.weights, r2.weights)
    assert r1.selected_lambda == r2.selected_lambda


def test_pipeline_invariant_to_n_jobs_with_fixed_seed() -> None:
    counts, _ = generate_compositional_counts(n_samples=120, n_taxa=10, seed=19)
    kwargs = dict(
        n_lambdas=12,
        n_subsamples=12,
        subsample_fraction=0.8,
        stars_threshold=0.1,
        random_state=7,
    )
    r1 = spiec_easi_mb(counts, n_jobs=1, **kwargs)
    r2 = spiec_easi_mb(counts, n_jobs=2, **kwargs)

    assert np.array_equal(r1.adjacency, r2.adjacency)
    assert np.allclose(r1.weights, r2.weights)
    assert np.isclose(r1.selected_lambda, r2.selected_lambda)


def test_pipeline_taxa_permutation_invariance() -> None:
    # Permutation invariance is asserted with strict equality, so the scenario must
    # be non-degenerate. On weak, near-empty networks a coefficient can sit exactly
    # on the Lasso soft-threshold boundary, which makes the coordinate-descent refit
    # order-dependent (a 0- vs 1-edge flip at the *same* selected alpha). With a
    # strong signal the recovered network is non-empty and the property holds
    # exactly. Weak compositional data is a documented degenerate regime.
    counts, _ = generate_strong_signal_counts(
        n_samples=200, n_taxa=16, depth=8000, rho=0.50, seed=21
    )
    rng = np.random.default_rng(21)
    perm = rng.permutation(counts.shape[1])
    inv = np.argsort(perm)

    r_ref = spiec_easi_mb(
        counts,
        n_lambdas=12,
        n_subsamples=12,
        subsample_fraction=0.8,
        stars_threshold=0.1,
        random_state=11,
        n_jobs=1,
    )
    r_perm = spiec_easi_mb(
        counts.iloc[:, perm],
        n_lambdas=12,
        n_subsamples=12,
        subsample_fraction=0.8,
        stars_threshold=0.1,
        random_state=11,
        n_jobs=1,
    )

    A_back = r_perm.adjacency[np.ix_(inv, inv)]
    W_back = r_perm.weights[np.ix_(inv, inv)]

    # Adjacency (the structural claim) must be exactly permutation-invariant.
    assert np.array_equal(r_ref.adjacency, A_back)
    # Symmetrized weights are only determined to the coordinate-descent solver
    # tolerance (``thol = max(tol, 1e-4)`` in the huge-scale CD), so compare with a
    # matching absolute tolerance rather than exact equality.
    assert np.allclose(r_ref.weights, W_back, atol=1e-4)
