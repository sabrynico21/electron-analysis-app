#!/usr/bin/env python
"""Point F (p=1135) - Python side: reproduce R's huge MB path on the SAME S and
SAME lambda path (no timing comparison) and compare the adjacency.

Reads the artifacts written by ``point_f.R`` (S, lambda path, R adjacency),
verifies the lambda path is the real one (equal to ``make_alpha_path`` on the
full CLR matrix), runs the compiled core on R's S, and reports per-lambda
differing edges. Analysis tooling only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PKG = HERE.parents[1]
sys.path.insert(0, str(PKG))

from spieceasi import _mb_core  # noqa: E402
from spieceasi.mb_solver import HUGE_MAX_ITER, HUGE_THOL  # noqa: E402
from spieceasi.neighborhood import make_alpha_path  # noqa: E402
from spieceasi.transforms import clr_transform  # noqa: E402

OUT = HERE / "results" / "point_f"


def lm(p: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(p, delimiter="\t"))


def lv(p: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(p, delimiter="\t"))


def sym_or(coef: np.ndarray) -> np.ndarray:
    A = coef != 0
    return A | A.T


def main() -> int:
    S = lm(OUT / "S.tsv")
    grid = lv(OUT / "lambda_path.txt")
    lmax = float(lv(OUT / "lambda_max.txt")[0])
    adj_r = lm(OUT / "spmb_adj_r.tsv")  # (n_lambda, d*d)
    d = S.shape[0]
    n_lambda = grid.size

    counts = pd.read_csv(
        HERE / "data" / "genere_table.txt", sep="\t", index_col=0
    ).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    X_clr = clr_transform(counts.to_numpy().T).clr  # samples x genera

    # verify the lambda path is the real one
    py_grid = make_alpha_path(X_clr, n_lambdas=n_lambda, lambda_min_ratio=1e-3)
    path_ok = np.allclose(py_grid, grid, rtol=1e-12, atol=1e-15)
    lmax_ok = np.isclose(make_alpha_path(X_clr, n_lambdas=2, lambda_min_ratio=1e-3)[0], lmax,
                         rtol=1e-12, atol=1e-15)

    thol = max(1e-4, HUGE_THOL)
    coef, _ = _mb_core.solve_mb_path(S, grid, thol, min(10_000, HUGE_MAX_ITER))
    coef = np.ascontiguousarray(coef)

    total_diff = 0
    per_lambda = []
    for k in range(n_lambda):
        A_py = sym_or(coef[k])
        A_r = adj_r[k].reshape(d, d).astype(bool)
        diff = int(np.sum(A_py != A_r)) // 2
        per_lambda.append(diff)
        total_diff += diff

    print(f"POINT_F(python) d={d} n_lambda={n_lambda} lambda_max_R={lmax:.10g}")
    print(f"  real_lambda_path_match={path_ok} lambda_max_match={lmax_ok}")
    print(f"  total_differing_edges_over_path={total_diff}")
    print(f"  per_lambda_diff={per_lambda}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
