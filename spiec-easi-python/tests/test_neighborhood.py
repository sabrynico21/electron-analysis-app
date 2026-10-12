from __future__ import annotations

import numpy as np

from spieceasi.neighborhood import make_alpha_path, mb_neighborhood_selection
from spieceasi.transforms import clr_transform
from tests.synthetic_data import generate_compositional_counts


def test_mb_outputs_symmetric_adjacency() -> None:
    counts, _ = generate_compositional_counts(n_samples=60, n_taxa=6)
    X = clr_transform(counts).clr
    alpha = make_alpha_path(X, n_lambdas=8)[3]
    result = mb_neighborhood_selection(X, alpha=float(alpha), sym_rule="or")

    assert np.array_equal(result.adjacency, result.adjacency.T)
    assert np.all(np.diag(result.adjacency) == 0)
    assert result.weight_matrix.shape == result.adjacency.shape


def test_sym_rule_and_not_denser_than_or() -> None:
    counts, _ = generate_compositional_counts(n_samples=60, n_taxa=6, seed=12)
    X = clr_transform(counts).clr
    alpha = make_alpha_path(X, n_lambdas=8)[4]

    res_or = mb_neighborhood_selection(X, alpha=float(alpha), sym_rule="or")
    res_and = mb_neighborhood_selection(X, alpha=float(alpha), sym_rule="and")

    assert int(np.sum(res_and.adjacency)) <= int(np.sum(res_or.adjacency))
