"""CLI and end-to-end property tests."""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pytest

from spring_microbiome import spring
from spring_microbiome.cli import run


def _counts(n=30, p=6, seed=0):
    rng = np.random.default_rng(seed)
    return np.round(np.exp(rng.normal(2.0, 0.5, size=(n, p)))).astype(int)


def test_cli_run_returns_expected_keys():
    out = run({"data": _counts().tolist(), "nlambda": 5, "rep_num": 4, "seed": 1})
    for key in ("lambdaseq", "opt_index", "adjacency", "weights", "instability"):
        assert key in out
    assert len(out["lambdaseq"]) == 5
    assert len(out["adjacency"]) == 6
    assert 1 <= out["opt_index"] <= 5


def test_cli_run_accepts_counts_alias():
    out = run({"counts": _counts().tolist(), "nlambda": 4, "rep_num": 4})
    assert out["n"] == 30


def test_cli_run_rejects_unknown_param():
    with pytest.raises(ValueError, match="unknown parameter"):
        run({"data": _counts().tolist(), "nope": 1})


def test_cli_run_rejects_missing_data():
    with pytest.raises(ValueError, match="data"):
        run({"nlambda": 5})


def test_cli_subprocess_stdin_stdout():
    payload = json.dumps({"data": _counts().tolist(), "nlambda": 5, "rep_num": 4})
    proc = subprocess.run(
        [sys.executable, "-m", "spring_microbiome"],
        input=payload,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert out["p"] == 6
    # stdout must be pure JSON (warnings, if any, go to stderr)
    assert proc.stdout.strip().startswith("{")


def test_cli_subprocess_error_exit_code():
    proc = subprocess.run(
        [sys.executable, "-m", "spring_microbiome"],
        input=json.dumps({"data": [[1, -1], [2, 3]]}),
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert proc.stdout == ""
    assert "error" in proc.stderr.lower()


def test_adjacency_is_symmetric_binary():
    res = spring(_counts(), nlambda=6, rep_num=5, seed=3, verbose=False)
    adj = np.asarray(res.adjacency)
    np.testing.assert_array_equal(adj, adj.T)
    assert set(np.unique(adj)).issubset({0.0, 1.0})
    assert np.all(np.diag(adj) == 0)


def test_determinism_with_injected_indices():
    data = _counts()
    idx = [np.random.default_rng(0).choice(30, 24, replace=False) for _ in range(5)]
    a = spring(data, nlambda=5, rep_num=5, subsample_indices=idx, verbose=False)
    b = spring(data, nlambda=5, rep_num=5, subsample_indices=idx, verbose=False)
    np.testing.assert_array_equal(np.asarray(a.adjacency), np.asarray(b.adjacency))


def test_reproducible_with_seed_and_rng():
    data = _counts()
    a = spring(data, nlambda=5, rep_num=5, seed=42, verbose=False)
    b = spring(data, nlambda=5, rep_num=5, seed=42, verbose=False)
    assert a.output.stars.opt_index == b.output.stars.opt_index
    np.testing.assert_array_equal(np.asarray(a.adjacency), np.asarray(b.adjacency))


def test_scale_invariance_of_compositional_input():
    """mclr is scale invariant, so scaling counts by a constant changes nothing."""
    data = _counts()
    idx = [np.random.default_rng(0).choice(30, 24, replace=False) for _ in range(5)]
    a = spring(data, nlambda=5, rep_num=5, subsample_indices=idx, verbose=False)
    b = spring(data * 7, nlambda=5, rep_num=5, subsample_indices=idx, verbose=False)
    np.testing.assert_allclose(a.lambdaseq, b.lambdaseq, rtol=1e-9, atol=1e-12)
    np.testing.assert_array_equal(np.asarray(a.adjacency), np.asarray(b.adjacency))


def test_synthetic_recovery_known_edges():
    """SPRING recovers a sparse planted network from synthetic counts."""
    n, p = 300, 8

    def build(seed, strength=0.7, scale=1.0):
        rng = np.random.default_rng(seed)
        omega = np.eye(p)
        for i, j in ((0, 1), (2, 3)):
            omega[i, j] = omega[j, i] = -strength
        assert np.linalg.eigvalsh(omega).min() > 0
        cov = np.linalg.inv(omega)
        d = np.sqrt(np.diag(cov))
        corr = cov / np.outer(d, d)
        z = rng.multivariate_normal(np.zeros(p), corr, size=n)
        return np.round(np.exp(2.0 + scale * z)).astype(int)

    for seed in (0, 1, 2):
        counts = build(seed)
        res = spring(
            counts,
            lambdaseq="data-specific",
            nlambda=25,
            rep_num=20,
            seed=1,
            verbose=False,
        )
        adj = np.asarray(res.adjacency)
        # the two planted edges must be recovered
        assert adj[0, 1] == 1 and adj[2, 3] == 1, f"seed {seed} missed a planted edge"
        # and the network must stay sparse (few false positives)
        assert adj.sum() / 2 <= 4, f"seed {seed} too dense: {adj.sum() / 2} edges"
