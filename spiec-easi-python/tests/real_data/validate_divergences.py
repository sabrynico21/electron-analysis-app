#!/usr/bin/env python
"""KKT / objective-gap validity check + "S di R" diagnostic for solver divergences.

For a case (L2 or L3) and the seeds present in ``results/battery_results.csv``,
it locates every (subsample, lambda) where the Python solver adjacency differs
from R's, then, using R's own correlation matrix S (from
``validate_solver_helper.R``):

* computes the Lasso KKT violation and the objective for BOTH R's and Python's
  coefficients at that (sub, lambda);
* re-solves with the C++ core fed R's S ("S di R") and checks whether the
  divergence disappears;
* reports, separately, divergences at lambda index 0 (lambda.max, covered by the
  documented exception) and at internal lambdas.

Declared validity criteria (docs/VERIFICA_DATI_REALI.md): KKT <= 1e-3, objective
gap <= 1e-5.  This is analysis tooling only; no production code is touched.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PKG = HERE.parents[1]
sys.path.insert(0, str(PKG))

from spieceasi.mb_solver import (  # noqa: E402
    HUGE_MAX_ITER,
    HUGE_THOL,
    solve_mb_path_huge_scale,
    symmetrize_path,
)
from spieceasi.transforms import clr_transform  # noqa: E402

GOLDEN = HERE / "golden"
RESULTS = HERE / "results"
HELPER = HERE / "validate_solver_helper.R"


def load_matrix(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path, delimiter="\t"))


def load_vector(path: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(path, delimiter="\t"))


def kkt_and_objective(S: np.ndarray, beta: np.ndarray, lam: float) -> tuple[float, float]:
    """Max KKT violation and total objective for MB coefficients ``beta`` (node x predictor)."""
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
            v = float(np.max(np.abs(gj[nz] + lam * np.sign(bj[nz]))))
            viol = max(viol, v)
        if (~nz).any():
            v = float(np.max(np.maximum(0.0, np.abs(gj[~nz]) - lam)))
            viol = max(viol, v)
        obj += float(0.5 * b @ (S @ b) - S[:, i] @ b + lam * np.sum(np.abs(b)))
    return viol, obj


def symmetric_pairs(adj: np.ndarray) -> set[tuple[int, int]]:
    d = adj.shape[0]
    return {(i, j) for i in range(d) for j in range(i + 1, d) if adj[i, j] or adj[j, i]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", default="L2")
    ap.add_argument("--out", default="")
    ap.add_argument("--seed", default="")
    args = ap.parse_args()
    case = args.case
    out_path = Path(args.out) if args.out else RESULTS / f"validity_{case}.csv"

    seeds = []
    with (RESULTS / "battery_results.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["case"] == case:
                seeds.append(int(r["seed"]))
    seeds = sorted(set(seeds))
    if args.seed:
        wanted = {int(x) for x in args.seed.split(",")}
        seeds = [s for s in seeds if s in wanted]

    fields = ["case", "seed", "sub", "lam_idx", "lam_1based", "at_lambda_max",
              "edges_R", "edges_Py", "edges_diff", "kkt_R", "kkt_Py",
              "obj_R", "obj_Py", "obj_gap", "both_valid", "sdiR_edges_diff"]
    rows_out = []
    for seed in seeds:
        d = GOLDEN / case / f"seed_{seed}"
        grid = load_vector(d / "lambda_path.txt")
        n_lambda = grid.size
        counts = load_matrix(d / "counts.tsv")
        X = clr_transform(counts).clr
        indices = load_matrix(d / "subsample_indices.tsv").astype(int)
        r_adj = load_matrix(d / "solver_adj_path.tsv").astype(np.int8).reshape(
            indices.shape[0], n_lambda, X.shape[1], X.shape[1])

        # find divergences
        diverging: dict[int, list[int]] = {}
        py_coef: dict[int, np.ndarray] = {}
        py_adj: dict[int, np.ndarray] = {}
        for s in range(indices.shape[0]):
            directed = solve_mb_path_huge_scale(X[indices[s]], grid, max_iter=10_000, tol=1e-4)
            adj = symmetrize_path(directed, "or").astype(np.int8)
            py_coef[s] = directed
            py_adj[s] = adj
            for k in range(n_lambda):
                if np.any(adj[k] != r_adj[s, k]):
                    diverging.setdefault(s, []).append(k)
        if not diverging:
            print(f"{case} seed={seed}: no solver divergence")
            continue

        subs = ",".join(str(s) for s in sorted(diverging))
        subprocess.run(["Rscript", str(HELPER), "--case", case, "--seed", str(seed),
                        "--subs", subs, "--outdir", str(d)], check=True,
                       capture_output=True, text=True)

        for s, ks in sorted(diverging.items()):
            S_R = load_matrix(d / f"s_{s}.tsv")
            beta_R = load_matrix(d / f"beta_{s}.tsv").reshape(n_lambda, X.shape[1], X.shape[1])
            beta_Py = py_coef[s]
            # S di R: feed the C++ core R's S
            from spieceasi import _mb_core
            thol = max(1e-4, HUGE_THOL)
            diag_coef, _ = _mb_core.solve_mb_path(S_R, grid, thol, min(10_000, HUGE_MAX_ITER))
            diag_adj = symmetrize_path(np.ascontiguousarray(diag_coef), "or").astype(np.int8)
            for k in ks:
                lam = float(grid[k])
                eR = symmetric_pairs(r_adj[s, k].astype(int))
                eP = symmetric_pairs(py_adj[s][k].astype(int))
                kkt_R, obj_R = kkt_and_objective(S_R, beta_R[k], lam)
                kkt_Py, obj_Py = kkt_and_objective(S_R, beta_Py[k], lam)
                gap = abs(obj_R - obj_Py)
                both = (kkt_R <= 1e-3) and (kkt_Py <= 1e-3) and (gap <= 1e-5)
                ediag = len(symmetric_pairs(diag_adj[k].astype(int)) ^ eR)
                rows_out.append({
                    "case": case, "seed": seed, "sub": s, "lam_idx": k,
                    "lam_1based": k + 1, "at_lambda_max": int(k == 0),
                    "edges_R": len(eR), "edges_Py": len(eP),
                    "edges_diff": len(eR ^ eP), "kkt_R": f"{kkt_R:.3e}",
                    "kkt_Py": f"{kkt_Py:.3e}", "obj_R": f"{obj_R:.10g}",
                    "obj_Py": f"{obj_Py:.10g}", "obj_gap": f"{gap:.3e}",
                    "both_valid": int(both), "sdiR_edges_diff": ediag,
                })
        print(f"{case} seed={seed}: divergences={sum(len(v) for v in diverging.values())} "
              f"subs={sorted(diverging)}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows_out)
    print(f"wrote {out_path} ({len(rows_out)} divergences)")

    # summary
    if rows_out:
        n = len(rows_out)
        n_lmax = sum(r["at_lambda_max"] for r in rows_out)
        n_int = n - n_lmax
        n_both = sum(r["both_valid"] for r in rows_out)
        n_diag_fixed = sum(r["sdiR_edges_diff"] == 0 for r in rows_out)
        print(f"SUMMARY {case}: total={n} at_lambda_max={n_lmax} internal_lambda={n_int}")
        print(f"  both_valid={n_both}/{n}; S_di_R_makes_it_match={n_diag_fixed}/{n}")
        print(f"  internal-lambda divergences NOT both valid: "
              f"{sum(1 for r in rows_out if not r['at_lambda_max'] and not r['both_valid'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
