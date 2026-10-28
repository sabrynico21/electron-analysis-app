"""CLI contract tests (JSON in / JSON out; clean stdout)."""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np

from cclasso.cli import run_request


def _request():
    rng = np.random.default_rng(0)
    z = rng.standard_normal((25, 5))
    frac = np.exp(z) / np.exp(z).sum(axis=1, keepdims=True)
    return {
        "x": frac.tolist(),
        "counts": False,
        "k_cv": 3,
        "n_boot": 5,
        "seed": 1,
    }


def test_run_request_keys():
    res = run_request(_request())
    for key in ("var_w", "cor_w", "p_vals", "lambda", "info_cv"):
        assert key in res
    assert len(res["cor_w"]) == 5
    assert set(res["info_cv"]) == {"lams", "fvals", "k", "lam_int"}


def test_cli_stdin_stdout_clean():
    proc = subprocess.run(
        [sys.executable, "-m", "cclasso"],
        input=json.dumps(_request()),
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    # stdout must be exactly one JSON object (no warnings/logs)
    payload = json.loads(proc.stdout)
    assert "cor_w" in payload and "lambda" in payload


def test_cli_error_goes_to_stderr_nonzero_exit():
    bad = {"x": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]], "counts": False}
    proc = subprocess.run(
        [sys.executable, "-m", "cclasso"],
        input=json.dumps(bad),
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert proc.stdout.strip() == ""
    assert "error" in proc.stderr.lower()


def test_cli_missing_x():
    proc = subprocess.run(
        [sys.executable, "-m", "cclasso"],
        input=json.dumps({"counts": False}),
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert proc.stdout.strip() == ""
