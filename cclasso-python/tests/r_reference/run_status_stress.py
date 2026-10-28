#!/usr/bin/env python3
"""Compare R/Python CCLasso failure categories on reproducible random cases.

This is the CCLasso P3 runner. It creates a self-contained input/draw archive,
classifies each case in Python, asks the pinned R source to classify the same
arrays, and exits nonzero on any status mismatch. Run with ``--cases 1000`` for
the registered sweep; a smaller count is useful as a smoke test.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import warnings
from pathlib import Path

import numpy as np

from cclasso import cclasso


def _write_matrix(path: Path, values: np.ndarray) -> None:
    np.savetxt(path, np.asarray(values), delimiter=",", fmt="%.17g")


def _make_counts(rng: np.random.Generator, case_index: int) -> np.ndarray:
    n = int(rng.integers(8, 19))
    p = int(rng.integers(2, min(10, n + 4)))
    mode = case_index % 4
    latent = rng.normal(0.0, 0.8, size=(n, p))
    library = rng.integers(12, 2001, size=(n, 1))
    composition = np.exp(latent - latent.max(axis=1, keepdims=True))
    composition /= composition.sum(axis=1, keepdims=True)
    counts = np.vstack(
        [rng.multinomial(int(library[i, 0]), composition[i]) for i in range(n)]
    ).astype(float)

    if mode == 0:
        counts[rng.random(counts.shape) < 0.15] = 0.0
    elif mode == 1:
        counts[rng.random(counts.shape) < 0.65] = 0.0
    elif mode == 2:
        counts[:, -1] = counts[0, -1]  # constant feature / rank deficiency
        counts[rng.random(counts.shape) < 0.35] = 0.0
    else:
        counts = np.maximum(0.0, np.rint(counts / 250.0))  # low library sizes
        counts[rng.random(counts.shape) < 0.45] = 0.0

    return counts


def _python_status(x, boot, null_norm, params) -> dict[str, object]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            result = cclasso(
                x,
                counts=True,
                pseudo=params["pseudo"],
                k_cv=params["k_cv"],
                lam_int=tuple(params["lam_int"]),
                k_max=params["k_max"],
                n_boot=params["n_boot"],
                boot_indices=boot,
                null_norm=null_norm,
            )
        except Exception as exc:  # status is compared, not the language-specific text
            return {"status": "error", "error_class": type(exc).__name__}

    nonfinite = {
        key: int((~np.isfinite(result[key])).sum()) for key in ("var_w", "cor_w", "p_vals")
    }
    return {
        "status": "nonfinite" if any(nonfinite.values()) else "finite",
        "error_class": "",
        "nonfinite_var": nonfinite["var_w"],
        "nonfinite_cor": nonfinite["cor_w"],
        "nonfinite_p": nonfinite["p_vals"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--cases", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261019)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error("--cases must be positive")
    out_dir = args.out_dir
    if out_dir.exists() and any(out_dir.iterdir()):
        parser.error(f"refusing to overwrite non-empty output directory: {out_dir}")
    cases_dir = out_dir / "cases"
    cases_dir.mkdir(parents=True, exist_ok=True)

    master_rng = np.random.default_rng(args.seed)
    python_rows = []
    for i in range(args.cases):
        case_rng = np.random.default_rng(int(master_rng.integers(0, 2**32 - 1)))
        x = _make_counts(case_rng, i)
        n, p = x.shape
        n_boot = 3
        draw_rng = np.random.default_rng(int(master_rng.integers(0, 2**32 - 1)))
        boot = draw_rng.integers(0, n, size=(n, n_boot), dtype=np.intp)
        null_norm = draw_rng.standard_normal((n, p))
        params = {
            "counts": True,
            "pseudo": 0.5,
            "k_cv": 3,
            "lam_int": [1e-4, 1.0],
            "k_max": 3,
            "n_boot": n_boot,
            "case_seed": args.seed,
            "case_index": i,
        }
        case_dir = cases_dir / f"case-{i:05d}"
        case_dir.mkdir()
        _write_matrix(case_dir / "input.csv", x)
        _write_matrix(case_dir / "boot_indices.csv", boot + 1)  # R indices
        _write_matrix(case_dir / "null_norm.csv", null_norm)
        (case_dir / "params.json").write_text(json.dumps(params, indent=2) + "\n", encoding="utf-8")
        python_rows.append({"case": case_dir.name, **_python_status(x, boot, null_norm, params)})

    python_status = out_dir / "python_status.csv"
    with python_status.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(python_rows[0]))
        writer.writeheader()
        writer.writerows(python_rows)

    r_status = out_dir / "r_status.csv"
    r_runner = Path(__file__).with_name("stress_status_R.R")
    r_run = subprocess.run(
        ["Rscript", str(r_runner), str(cases_dir), str(r_status)],
        cwd=Path(__file__).resolve().parents[2],
        capture_output=True,
        text=True,
    )
    if r_run.returncode:
        raise RuntimeError(
            "R stress classification failed:\n"
            f"stdout: {r_run.stdout[-4000:]}\nstderr: {r_run.stderr[-4000:]}"
        )
    with r_status.open(newline="", encoding="utf-8") as stream:
        r_rows = {row["case"]: row for row in csv.DictReader(stream)}
    mismatches = []
    for row in python_rows:
        ref = r_rows.get(row["case"])
        if ref is None or ref["status"] != row["status"]:
            mismatches.append({"case": row["case"], "python": row["status"], "R": ref})
        elif row["status"] == "nonfinite":
            for py_key, r_key in (
                ("nonfinite_var", "nonfinite_var"),
                ("nonfinite_cor", "nonfinite_cor"),
                ("nonfinite_p", "nonfinite_p"),
            ):
                if int(row[py_key]) != int(ref[r_key]):
                    mismatches.append(
                        {
                            "case": row["case"],
                            "python": row["status"],
                            "R": ref["status"],
                            "field": py_key,
                            "python_count": row[py_key],
                            "R_count": ref[r_key],
                        }
                    )

    summary = {
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
        "mismatches_preview": mismatches[:20],
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
