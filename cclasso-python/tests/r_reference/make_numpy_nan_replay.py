#!/usr/bin/env python3
"""Recreate the known seeded CCLasso non-finite bootstrap case for R replay.

The generated directory is diagnostic input, not a golden output. It contains
the exact Python input and NumPy-generated random draws; invoke
``replay_numpy_nan.R`` on it to determine whether R shares the failure.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import warnings
from pathlib import Path

import numpy as np

from cclasso import cclasso
from cclasso._core import _resolve_random_inputs


def _input_matrix() -> np.ndarray:
    n, p, seed = 25, 5, 8
    rng = np.random.default_rng(seed)
    corr = np.eye(p)
    for i, j in ((0, 1), (1, 2), (0, 2), (3, 4)):
        if j < p:
            corr[i, j] = corr[j, i] = 0.5
    latent = rng.standard_normal((n, p)) @ np.linalg.cholesky(corr).T
    values = np.exp(latent)
    return values / values.sum(axis=1, keepdims=True)


def _write_matrix(path: Path, values: np.ndarray) -> None:
    np.savetxt(path, np.asarray(values), delimiter=",", fmt="%.17g")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    x = _input_matrix()
    n, p = x.shape
    params = {
        "case": "numpy_seeded_nonfinite_replay",
        "input_seed": 8,
        "counts": False,
        "pseudo": 0.5,
        "k_cv": 3,
        "n_boot": 6,
        "seed": 42,
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "platform": platform.platform(),
    }
    boot_indices, null_norm = _resolve_random_inputs(
        n, p, params["n_boot"], None, None, params["seed"]
    )
    _write_matrix(out_dir / "input.csv", x)
    # R's sample() uses one-based indices.
    _write_matrix(out_dir / "boot_indices.csv", boot_indices + 1)
    _write_matrix(out_dir / "null_norm.csv", null_norm)
    (out_dir / "params.json").write_text(json.dumps(params, indent=2) + "\n", encoding="utf-8")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = cclasso(
            x,
            counts=params["counts"],
            pseudo=params["pseudo"],
            k_cv=params["k_cv"],
            n_boot=params["n_boot"],
            seed=params["seed"],
        )
    _write_matrix(out_dir / "python_bootstrap_variances.csv", result["_bootstrap_variances"])
    _write_matrix(out_dir / "python_var_w.csv", np.asarray(result["var_w"])[:, None])
    _write_matrix(out_dir / "python_cor_w.csv", result["cor_w"])
    _write_matrix(out_dir / "python_p_vals.csv", result["p_vals"])
    summary = {
        key: {
            "finite": bool(np.isfinite(result[key]).all()),
            "nan_count": int(np.isnan(result[key]).sum()),
        }
        for key in ("var_w", "cor_w", "p_vals")
    }
    summary["negative_bootstrap_variance_count"] = int(
        np.count_nonzero(result["_bootstrap_variances"] < 0.0)
    )
    (out_dir / "python_status.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
