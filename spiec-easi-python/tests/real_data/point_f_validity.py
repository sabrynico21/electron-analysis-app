#!/usr/bin/env python
"""p=1135 validity (KKT/objective) + "S di R" diagnostic at lambda idx 1..5.

Reads the artifacts written by ``point_f_validity.R`` (R's S, the real lambda
path, R's huge MB coefficients at idx 1..5) and:

* runs the compiled core on R's S and on OUR OWN S (X_clr subsample) --- the
  latter is the "S di R" comparison;
* computes the Lasso KKT violation and objective for R's solution and ours
  (declared criteria: KKT <= 1e-3, objective gap <= 1e-5).

Analysis tooling only; no production code touched.
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
from spieceasi.mb_solver import HUGE_MAX_ITER, HUGE_THOL, huge_correlation  # noqa: E402
from spieceasi.transforms import clr_transform  # noqa: E402

OUT = HERE / "results" / "point_f_validity"
IDX = [1, 2, 3, 4, 5]  # 0-based lambda indices


def lm(p: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(p, delimiter="\t"))


def lv(p: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(p, delimiter="\t"))


def kkt_and_objective(S: np.ndarray, beta: np.ndarray, lam: float) -> tuple[float, float]:
    d = S.shape[0]
    viol = 0.0
    obj = 0.0
    for i in range(d):
        b = beta[i].copy()
        b[i] = 0.0
        g = S @ b - S[:, i]
        mask = np.ones(d, dtype=bool)
        mask[i] = False
        gj = g[mask]
        bj = b[mask]
        nz = np.abs(bj) > 0.0
        if nz.any():
            viol = max(viol, float(np.max(np.abs(gj[nz] + lam * np.sign(bj[nz])))))
        if (~nz).any():
            viol = max(viol, float(np.max(np.maximum(0.0, np.abs(gj[~nz]) - lam))))
        obj += float(0.5 * b @ (S @ b) - S[:, i] @ b + lam * np.sum(np.abs(b)))
    return viol, obj


def main() -> int:
    S_R = lm(OUT / "S.tsv")
    grid = lv(OUT / "lambda_path.txt")
    beta_R = lm(OUT / "beta_idx1_5.tsv").reshape(len(IDX), S_R.shape[0], S_R.shape[0])

    counts = pd.read_csv(
        HERE / "data" / "genere_table.txt", sep="\t", index_col=0
    ).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    X_clr = clr_transform(counts.to_numpy().T).clr  # samples x genera
    idx = lm(HERE / "results" / "point_f" / "subsample_indices.tsv").ravel().astype(int) - 1
    own_S = huge_correlation(X_clr[idx])
    sdiff = float(np.max(np.abs(own_S - S_R)))

    thol = max(1e-4, HUGE_THOL)
    mi = min(10_000, HUGE_MAX_ITER)
    coef_rs, _ = _mb_core.solve_mb_path(S_R, grid, thol, mi)
    coef_own, _ = _mb_core.solve_mb_path(own_S, grid, thol, mi)
    coef_rs = np.ascontiguousarray(coef_rs)
    coef_own = np.ascontiguousarray(coef_own)

    print(f"p=1135 S_di_R: max|own_S - R_S| = {sdiff:.3e}  (our S vs R's S)")
    print(f"{'idx':>3} | {'lam':>10} | {'kkt_R':>9} {'kkt_PyS':>9} {'kkt_ownS':>9} | "
          f"{'objgap_R_vs_PyS':>15} {'both_valid':>10} | {'edgediff(R,ourS)':>16}")
    for k, li in enumerate(IDX):
        lam = float(grid[li])
        vR, oR = kkt_and_objective(S_R, beta_R[k], lam)
        vP, oP = kkt_and_objective(S_R, coef_rs[li], lam)
        vO, oO = kkt_and_objective(own_S, coef_own[li], lam)
        gap = abs(oR - oP)
        both = (vR <= 1e-3) and (vP <= 1e-3) and (gap <= 1e-5)
        eR = beta_R[k] != 0
        eR = eR | eR.T
        eP = coef_rs[li] != 0
        eP = eP | eP.T
        ed = int(np.sum(eR != eP)) // 2
        print(f"{li:>3} | {lam:>10.6f} | {vR:>9.2e} {vP:>9.2e} {vO:>9.2e} | "
              f"{gap:>15.3e} {str(both):>10} | {ed:>16}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
