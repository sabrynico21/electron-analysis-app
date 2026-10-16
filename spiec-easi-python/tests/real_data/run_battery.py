#!/usr/bin/env python
"""Driver + comparator for the real-data verification battery.

For every (case, seed) in the declared order it:

1. ensures the R golden exists (runs ``generate_r_golden.R`` if missing);
2. runs the Python pipeline with R's subsample indices injected;
3. compares the 9 steps against R (same tolerances as the extended parity
   suite) and classifies the outcome (section D of docs/VERIFICA_DATI_REALI.md);
4. appends one CSV row and flushes ``results/battery_results.csv`` after EVERY
   (case, seed).

The script is idempotent: (case, seed) pairs already present in the CSV are
skipped, so the run can be interrupted (Ctrl-C / killed) and resumed.

Usage (from ``spiec-easi-python/``):

    .venv/bin/python tests/real_data/run_battery.py            # block (a)
    .venv/bin/python tests/real_data/run_battery.py --cases R6,L3
    .venv/bin/python tests/real_data/run_battery.py --status
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PKG = HERE.parents[1]  # spiec-easi-python/
sys.path.insert(0, str(PKG))

from tests.test_equivalence import (  # noqa: E402
    TOL_COEFF,
    TOL_LAMBDA,
    TOL_OPT_LAMBDA,
    TOL_STARS,
    TOL_TRANSFORM,
)

from spieceasi.mb_solver import (  # noqa: E402
    huge_correlation,
    solve_mb_path_huge_scale,
    symmetrize_path,
)
from spieceasi.neighborhood import alpha_max_for_mb, make_alpha_path  # noqa: E402
from spieceasi.pipeline import spiec_easi_mb  # noqa: E402
from spieceasi.transforms import clr_transform  # noqa: E402

GOLDEN = HERE / "golden"
RESULTS = HERE / "results"
CSV_PATH = RESULTS / "battery_results.csv"
R_SCRIPT = HERE / "generate_r_golden.R"
BASE_SEEDS = list(range(1001, 1011))

KIND = {
    "R1": "real", "R2": "real", "R3": "real", "R4": "real", "R5": "real",
    "R6": "real", "L1": "limit", "L2": "limit", "L3": "real", "L4": "limit",
    "S1": "synthetic", "S2": "synthetic", "S3": "synthetic", "S4": "synthetic",
    "P1": "pathological", "P2": "pathological", "P3": "pathological",
}

# Block (a) order; R6 and L3 are run separately afterwards.
BLOCK_A = ["R1", "R2", "R3", "R4", "R5", "L1", "L2", "L4", "S1", "S2", "S3", "S4", "P1", "P2", "P3"]
ALL_ORDER = [*BLOCK_A, "R6", "L3"]

FIELDS = [
    "case", "seed", "kind", "status", "classification", "first_div_step",
    "p", "n", "opt_index_r", "opt_index_py", "opt_lambda_r", "opt_lambda_py",
    "n_edges_r", "n_edges_py", "final_diff_edges", "clr_maxdiff", "lmax_diff",
    "grid_maxdiff", "solver_diff_edges", "coef_maxdiff", "stars_raw_maxdiff",
    "stars_mono_maxdiff", "opt_index_match", "opt_lambda_diff", "elapsed_s",
    "note",
]


def load_matrix(path: Path) -> np.ndarray:
    arr = np.loadtxt(path, delimiter="\t")
    return np.atleast_2d(arr)


def load_vector(path: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(path, delimiter="\t"))


def load_scalar(path: Path) -> float:
    return float(np.loadtxt(path))


def load_params(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def symmetric_edge_pairs(adj: np.ndarray) -> set[tuple[int, int]]:
    d = adj.shape[0]
    return {(i, j) for i in range(d) for j in range(i + 1, d)
            if adj[i, j] != 0 or adj[j, i] != 0}


def argmax_abs_corr_pair(X: np.ndarray) -> tuple[int, int]:
    corr = np.abs(huge_correlation(X))
    np.fill_diagonal(corr, 0.0)
    i, j = np.unravel_index(int(np.argmax(corr)), corr.shape)
    return (min(int(i), int(j)), max(int(i), int(j)))


def max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def compare(case: str, seed: int, d: Path) -> dict:
    dims = load_params(d / "solver_dims.txt")
    n_lambda = int(dims["n_lambda"])
    n_sub = int(dims["n_sub"])
    p = int(dims["p"])
    n = int(dims["n"])
    counts = load_matrix(d / "counts.tsv")

    row: dict = {"case": case, "seed": seed, "p": p, "n": n,
                 "status": "ok", "first_div_step": ""}
    first_div = None

    # 1. CLR
    r_clr = load_matrix(d / "clr.tsv")
    py_clr = clr_transform(counts).clr
    row["clr_maxdiff"] = max_abs_diff(py_clr, r_clr)
    if not (py_clr.shape == r_clr.shape and np.allclose(py_clr, r_clr, **TOL_TRANSFORM)):
        first_div = first_div or "1_clr"

    # 2. lambda.max
    r_lmax = load_scalar(d / "lambda_max.txt")
    py_lmax = alpha_max_for_mb(py_clr)
    row["lmax_diff"] = abs(py_lmax - r_lmax)
    if not np.isclose(py_lmax, r_lmax, **TOL_LAMBDA):
        first_div = first_div or "2_lambda_max"

    # 3. lambda grid
    r_grid = load_vector(d / "lambda_path.txt")
    py_grid = make_alpha_path(py_clr, n_lambdas=n_lambda, lambda_min_ratio=1e-3)
    row["grid_maxdiff"] = max_abs_diff(py_grid, r_grid)
    if not (py_grid.shape == r_grid.shape and np.allclose(py_grid, r_grid, **TOL_LAMBDA)):
        first_div = first_div or "3_lambda_path"

    indices = load_matrix(d / "subsample_indices.tsv").astype(int)
    subsamples = [indices[i] for i in range(indices.shape[0])]

    # 4. solver per subsample (adjacency per lambda)
    r_adj = load_matrix(d / "solver_adj_path.tsv").astype(np.int8).reshape(n_sub, n_lambda, p, p)
    solver_max = 0
    first_bad = None
    for s in range(n_sub):
        directed = solve_mb_path_huge_scale(py_clr[indices[s]], py_grid, max_iter=10_000, tol=1e-4)
        adj = symmetrize_path(directed, "or").astype(np.int8)
        bad = int(np.sum(adj != r_adj[s]) // 2)
        if bad:
            solver_max = max(solver_max, bad)
            if first_bad is None:
                k, i, j = np.argwhere(adj != r_adj[s])[0]
                first_bad = (s, int(k), int(i), int(j))
    row["solver_diff_edges"] = solver_max
    if first_bad is not None:
        first_div = first_div or f"4_solver(sub={first_bad[0]},lam={first_bad[1]})"

    # 4b. coefficients, subsample 0
    r_coef = load_matrix(d / "solver_coef_path.tsv").reshape(n_lambda, p, p)
    directed0 = solve_mb_path_huge_scale(py_clr[indices[0]], py_grid, max_iter=10_000, tol=1e-4)
    row["coef_maxdiff"] = max_abs_diff(directed0, r_coef)
    if not np.allclose(directed0, r_coef, **TOL_COEFF):
        first_div = first_div or "4b_coef"

    # 5-9. full pipeline with injected subsamples
    result = spiec_easi_mb(
        counts, n_lambdas=n_lambda, lambda_min_ratio=1e-3, n_subsamples=n_sub,
        alpha_path=py_grid, subsamples=subsamples, stars_threshold=0.05,
        sym_rule="or", random_state=0, max_iter=10_000, tol=1e-4, n_jobs=-1,
        require_samples_gt_taxa=False,
    )
    r_raw = load_vector(d / "stars_raw.txt")
    row["stars_raw_maxdiff"] = max_abs_diff(result.instability, r_raw)
    if not np.allclose(result.instability, r_raw, **TOL_STARS):
        first_div = first_div or "5_stars_raw"

    r_mono = load_vector(d / "stars_monotone.txt")
    row["stars_mono_maxdiff"] = max_abs_diff(result.monotone_instability, r_mono)
    if not np.allclose(result.monotone_instability, r_mono, **TOL_STARS):
        first_div = first_div or "6_stars_monotone"

    r_opt_index = int(load_scalar(d / "opt_index.txt"))
    r_opt_lambda = load_scalar(d / "opt_lambda.txt")
    py_idx = int(np.argmax(np.isclose(py_grid, result.selected_lambda, **TOL_OPT_LAMBDA)))
    row["opt_index_r"] = r_opt_index
    row["opt_index_py"] = py_idx + 1
    row["opt_index_match"] = int(py_idx == r_opt_index - 1)
    row["opt_lambda_r"] = r_opt_lambda
    row["opt_lambda_py"] = float(result.selected_lambda)
    row["opt_lambda_diff"] = abs(result.selected_lambda - r_opt_lambda)

    r_final = load_matrix(d / "final_adj.tsv").astype(int)
    py_final = result.adjacency.astype(int)
    differing = symmetric_edge_pairs(py_final) ^ symmetric_edge_pairs(r_final)
    row["final_diff_edges"] = len(differing)
    row["n_edges_r"] = len(symmetric_edge_pairs(r_final))
    row["n_edges_py"] = len(symmetric_edge_pairs(py_final))

    if py_idx == 0:
        allowed = differing <= {argmax_abs_corr_pair(py_clr)}
    else:
        allowed = not differing

    # classification (section D)
    lambda_max_exc = (py_idx == 0)
    if py_idx != r_opt_index - 1 or not allowed:
        if not allowed:
            first_div = first_div or "9_final_adjacency"
        if py_idx != r_opt_index - 1:
            first_div = first_div or "7_opt_index"
        classification = "diverso"
    else:
        identical = (
            first_div is None
            and row["stars_raw_maxdiff"] <= 1e-10
            and row["stars_mono_maxdiff"] <= 1e-10
            and row["solver_diff_edges"] == 0
        )
        classification = "identico" if identical else "equivalente"

    row["classification"] = classification
    if first_div:
        row["first_div_step"] = first_div
    if lambda_max_exc:
        row["note"] = "lambda.max (indeterminate exception)"
    return row


def r_list_cases() -> dict[str, int]:
    out = subprocess.run(["Rscript", str(R_SCRIPT), "--list"], capture_output=True, text=True, cwd=str(PKG))
    cases: dict[str, int] = {}
    for line in out.stdout.strip().splitlines():
        parts = line.split()
        if len(parts) == 2:
            cases[parts[0]] = int(parts[1])
    return cases


def ensure_golden(case: str, seed: int) -> Path:
    d = GOLDEN / case / f"seed_{seed}"
    if (d / "done.txt").is_file() or (d / "error.txt").is_file():
        return d
    subprocess.run(["Rscript", str(R_SCRIPT), "--case", case, "--seed", str(seed)],
                   cwd=str(PKG), check=False)
    return d


def read_done_keys() -> set[tuple[str, int]]:
    if not CSV_PATH.is_file():
        return set()
    keys = set()
    with CSV_PATH.open() as fh:
        for r in csv.DictReader(fh):
            keys.add((r["case"], int(r["seed"])))
    return keys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", default=",".join(BLOCK_A))
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)
    done = read_done_keys()
    if args.status:
        case_counts: dict[str, int] = {}
        for c, _ in sorted(done):
            case_counts[c] = case_counts.get(c, 0) + 1
        print("completed (case, seeds):")
        for c, k in sorted(case_counts.items()):
            print(f"  {c}: {k}")
        print(f"total rows: {len(done)}")
        return 0

    # only include cases that actually exist in the R table
    avail = r_list_cases()
    cases = [c.strip() for c in args.cases.split(",") if c.strip()]
    cases = [c for c in cases if c in avail]

    write_header = not CSV_PATH.is_file()
    with CSV_PATH.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if write_header:
            w.writeheader()
            fh.flush()
        for case in cases:
            for seed in BASE_SEEDS[: avail[case]]:
                if (case, seed) in done:
                    continue
                t0 = time.time()
                d = ensure_golden(case, seed)
                errfile = d / "error.txt"
                if errfile.is_file():
                    msg = errfile.read_text().strip()
                    row = {k: "" for k in FIELDS}
                    row.update({"case": case, "seed": seed, "status": "r_error",
                                "note": msg[:200], "first_div_step": "r_error",
                                "classification": "diverso", "elapsed_s": round(time.time() - t0, 2)})
                else:
                    try:
                        row = compare(case, seed, d)
                        row["elapsed_s"] = round(time.time() - t0, 2)
                    except Exception as exc:  # noqa: BLE001
                        row = {k: "" for k in FIELDS}
                        row.update({"case": case, "seed": seed, "status": "py_error",
                                    "note": f"{type(exc).__name__}: {exc}"[:200],
                                    "first_div_step": "py_error", "classification": "diverso",
                                    "elapsed_s": round(time.time() - t0, 2)})
                row["kind"] = KIND.get(case, "")
                w.writerow({k: row.get(k, "") for k in FIELDS})
                fh.flush()
                os.fsync(fh.fileno())
                print(f"ROW {case} seed={seed} class={row['classification']} "
                      f"status={row['status']} t={row['elapsed_s']}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
