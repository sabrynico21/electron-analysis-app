#!/usr/bin/env python3
"""Seeded R/Python status and discrete-output stress sweep for SPIEC-EASI/SPRING.

Small count tables are generated once, R subsample indices are saved/injected,
and both implementations are classified on identical data. Run with
``--cases 1000`` for the registered sweep; smaller counts smoke-test the runner.
The complete case and per-language status archive is retained at ``out_dir``.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np


def _counts(rng: np.random.Generator, index: int) -> np.ndarray:
    n = int(rng.integers(8, 21))
    p = int(rng.integers(2, min(11, n + 5)))
    mode = index % 5
    latent = rng.normal(0.0, 0.75, size=(n, p))
    composition = np.exp(latent - latent.max(axis=1, keepdims=True))
    composition /= composition.sum(axis=1, keepdims=True)
    totals = rng.integers(100, 2501, size=n)
    x = np.vstack([rng.multinomial(int(totals[i]), composition[i]) for i in range(n)])
    if mode in (1, 3):
        x[rng.random(x.shape) < (0.35 if mode == 1 else 0.7)] = 0
    if mode == 2 and p > 2:
        x[:, -1] = x[:, 0]  # duplicate column, deliberately singular design
    if mode == 4 and p > 2:
        x[:, -1] = 1  # constant feature stress case
        x[rng.random(x.shape) < 0.2] = 0
    for row in range(n):
        if x[row].sum() == 0:
            x[row, int(rng.integers(0, p))] = 1
    return x.astype(float)


def _finite_result(method: str, result: dict) -> bool:
    if method == "spieceasi":
        arrays = [result["beta"], result["stability"], result["refit"]]
    else:
        arrays = [result.beta.beta, result.output.stars.summary, result.adjacency]
    return all(np.isfinite(np.asarray(array, dtype=float)).all() for array in arrays)


def _run_python(method: str, case_dir: Path, params: dict) -> dict:
    x = np.loadtxt(case_dir / "input.csv", delimiter=",", ndmin=2)
    idx = np.loadtxt(case_dir / "r_subsample_indices.csv", delimiter=",", dtype=int, ndmin=2) - 1
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if method == "spieceasi":
                from spieceasi import spiec_easi

                result = spiec_easi(
                    x,
                    method="mb",
                    sel_criterion="stars",
                    nlambda=3,
                    lambda_min_ratio=0.1,
                    lambda_log=True,
                    pulsar_params={
                        "rep_num": params["rep_num"],
                        "subsample_ratio": 0.8,
                        "thresh": 0.2,
                        "seed": params["seed"],
                    },
                    subsample_indices=[row for row in idx],
                    n_jobs=1,
                )
                status = "finite" if _finite_result(method, result) else "nonfinite"
                return {
                    "status": status,
                    "opt_index": int(result["opt_index"]),
                    "adjacency": np.asarray(result["refit"], dtype=int),
                    "beta": np.asarray(result["beta"][int(result["opt_index"]) - 1]),
                    "stability": np.asarray(result["stability"], dtype=float),
                    "lambda_path": np.asarray(result["lambda"], dtype=float),
                }

            from spring_microbiome import spring

            lambda_path_file = case_dir / "r_lambda_path.csv"
            lambdaseq = (
                np.loadtxt(lambda_path_file, delimiter=",").reshape(-1)
                if lambda_path_file.exists()
                else "data-specific"
            )
            result = spring(
                x,
                quantitative=False,
                lambda_min_ratio=0.1,
                nlambda=3,
                lambdaseq=lambdaseq,
                seed=params["seed"],
                thresh=0.2,
                subsample_ratio=0.8,
                rep_num=params["rep_num"],
                Rtol=1e-6,
                verbose=False,
                Rmethod="approx",
                use_nearPD=True,
                nu=0.001,
                ratio=0.9,
                subsample_indices=[row for row in idx],
            )
            status = "finite" if _finite_result(method, result) else "nonfinite"
            return {
                "status": status,
                "opt_index": int(result.output.stars.opt_index),
                "adjacency": np.asarray(result.adjacency, dtype=int),
                "beta": np.asarray(result.beta.beta[int(result.output.stars.opt_index) - 1]),
                "stability": np.asarray(result.output.stars.summary, dtype=float),
                "lambda_path": np.asarray(result.lambdaseq, dtype=float),
            }
    except Exception as exc:
        return {
            "status": "error",
            "error_class": type(exc).__name__,
            "error": str(exc),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=("spieceasi", "spring"))
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--cases", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error("--cases must be positive")
    out_dir = args.out_dir
    if out_dir.exists() and any(out_dir.iterdir()):
        parser.error(f"refusing to overwrite non-empty output directory: {out_dir}")
    cases_dir = out_dir / "cases"
    cases_dir.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]

    source_dir = root / ("spiec-easi-python/src" if args.method == "spieceasi" else "spring-python/src")
    sys.path.insert(0, str(source_dir))

    rng = np.random.default_rng(args.seed)
    for i in range(args.cases):
        case_rng = np.random.default_rng(int(rng.integers(0, 2**32 - 1)))
        x = _counts(case_rng, i)
        case_dir = cases_dir / f"case-{i:05d}"
        case_dir.mkdir()
        np.savetxt(case_dir / "input.csv", x, delimiter=",", fmt="%.17g")
        params = {"case_index": i, "seed": int(rng.integers(0, 2**31 - 1)), "rep_num": 2}
        (case_dir / "params.json").write_text(
            json.dumps(params, indent=2) + "\n", encoding="utf-8"
        )

    r_status_path = out_dir / "r_status.csv"
    r_script = root / "scripts/network_status_stress.R"
    r_run = subprocess.run(
        ["Rscript", str(r_script), args.method, str(cases_dir), str(out_dir)],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if r_run.returncode:
        raise RuntimeError(
            "R network stress runner failed:\n"
            f"stdout: {r_run.stdout[-4000:]}\nstderr: {r_run.stderr[-4000:]}"
        )

    with r_status_path.open(newline="", encoding="utf-8") as stream:
        r_rows = {row["case"]: row for row in csv.DictReader(stream)}
    python_rows = []
    mismatches = []
    boundary_ties = []
    for case_dir in sorted(cases_dir.iterdir()):
        params = json.loads((case_dir / "params.json").read_text())
        py = _run_python(args.method, case_dir, params)
        ref = r_rows.get(case_dir.name)
        record = {
            "case": case_dir.name,
            "status": py["status"],
            "error_class": py.get("error_class", ""),
            "error": py.get("error", ""),
            "opt_index": py.get("opt_index", ""),
        }
        if "adjacency" in py:
            np.savetxt(case_dir / "python_adjacency.csv", py["adjacency"], delimiter=",", fmt="%d")
            np.savetxt(case_dir / "python_beta_selected.csv", py["beta"], delimiter=",", fmt="%.17g")
            np.savetxt(case_dir / "python_stability.csv", py["stability"], delimiter=",", fmt="%.17g")
            np.savetxt(case_dir / "python_lambda_path.csv", py["lambda_path"], delimiter=",", fmt="%.17g")
        python_rows.append(record)
        if ref is None or ref["status"] != py["status"]:
            mismatches.append(
                {
                    "case": case_dir.name,
                    "python": py["status"],
                    "python_error": py.get("error"),
                    "R": ref,
                }
            )
            continue
        if py["status"] == "finite":
            if int(ref["opt_index"]) != py["opt_index"]:
                mismatches.append(
                    {
                        "case": case_dir.name,
                        "field": "opt_index",
                        "python": py["opt_index"],
                        "R": ref["opt_index"],
                        "python_stability": py["stability"].tolist(),
                        "R_stability": np.loadtxt(
                            case_dir / "r_stability.csv", delimiter=",", ndmin=1
                        ).tolist(),
                    }
                )
            r_adj = np.loadtxt(case_dir / "r_adjacency.csv", delimiter=",", ndmin=2).astype(int)
            if not np.array_equal(r_adj, py["adjacency"]):
                r_beta_path = case_dir / "r_beta_selected.csv"
                edge_details = []
                if r_beta_path.exists() and "beta" in py:
                    r_beta = np.loadtxt(r_beta_path, delimiter=",", ndmin=2)
                    for i, j in zip(*np.triu_indices(r_adj.shape[0], k=1)):
                        if r_adj[i, j] != py["adjacency"][i, j]:
                            edge_details.append(
                                {
                                    "edge": [int(i), int(j)],
                                    "r_abs_beta": float(max(abs(r_beta[i, j]), abs(r_beta[j, i]))),
                                    "python_abs_beta": float(max(abs(py["beta"][i, j]), abs(py["beta"][j, i]))),
                                }
                            )
                is_boundary = bool(edge_details) and all(
                    detail["r_abs_beta"] <= 1e-8
                    and detail["python_abs_beta"] <= 1e-8
                    for detail in edge_details
                )
                if is_boundary:
                    boundary_ties.append(
                        {"case": case_dir.name, "edges": edge_details}
                    )
                    continue
                mismatches.append(
                    {
                        "case": case_dir.name,
                        "field": "adjacency",
                        "symmetric_difference": int(np.count_nonzero(r_adj != py["adjacency"])),
                        "edge_details": edge_details[:10],
                    }
                )

    with (out_dir / "python_status.csv").open("w", newline="", encoding="utf-8") as stream:
        columns = sorted({key for row in python_rows for key in row})
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(python_rows)
    summary = {
        "method": args.method,
        "seed": args.seed,
        "cases": args.cases,
        "python_status_counts": {
            status: sum(row["status"] == status for row in python_rows)
            for status in ("finite", "nonfinite", "error")
        },
        "r_status_counts": {
            status: sum(row["status"] == status for row in r_rows.values())
            for status in ("finite", "nonfinite", "error")
        },
        "mismatch_count": len(mismatches),
        "mismatch_field_counts": {
            field: sum(item.get("field", "status") == field for item in mismatches)
            for field in sorted({item.get("field", "status") for item in mismatches})
        },
        "boundary_tie_count": len(boundary_ties),
        "boundary_ties_preview": boundary_ties[:20],
        "mismatches_preview": mismatches[:20],
    }
    (out_dir / "mismatches.json").write_text(json.dumps(mismatches, indent=2) + "\n")
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
