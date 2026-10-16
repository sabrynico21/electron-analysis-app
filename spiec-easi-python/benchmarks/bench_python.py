#!/usr/bin/env python3
"""FASE 1 baseline — timing of the *complete* MB+StARS pipeline (Python, serial).

Runs the production Python pipeline on fixed golden cases, injecting the R
subsample indices from the goldens so that R and Python perform the same work.
``n_jobs=1`` keeps it serial, to compare against serial R (no pulsar parallelism).

Usage:
    ./.venv/bin/python benchmarks/bench_python.py [--reps 3] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "spiec-easi-python" / "src"))

GOLD = REPO / "spiec-easi-python" / "tests" / "golden"
GOLD_EXT = REPO / "spiec-easi-python" / "tests" / "golden_extended"

# Same pipeline settings everywhere: nlambda=20, lambda.min.ratio=1e-3, rep.num=20.
CASES = [
    {
        "id": "ds2_p30_n100",
        "counts": GOLD / "ds2_r.counts.tsv",
        "subs": GOLD / "ds2_r.subsample_indices.tsv",
        "nlambda": 20,
        "ratio": 1e-3,
        "rep": 20,
    },
    {
        "id": "rand_p50_n50_s211",
        "counts": GOLD_EXT / "rand_p50_n50_s211" / "counts.tsv",
        "subs": GOLD_EXT / "rand_p50_n50_s211" / "subsample_indices.tsv",
        "nlambda": 20,
        "ratio": 1e-3,
        "rep": 20,
    },
    {
        "id": "largep_p100_n80_s215",
        "counts": GOLD_EXT / "largep_p100_n80_s215" / "counts.tsv",
        "subs": GOLD_EXT / "largep_p100_n80_s215" / "subsample_indices.tsv",
        "nlambda": 20,
        "ratio": 1e-3,
        "rep": 20,
    },
]


def environment() -> dict:
    import numpy

    blas = "unknown"
    try:
        info = numpy.show_config(mode="dicts")
        blas = info.get("Build Dependencies", {}).get("blas", {}).get("name", "unknown")
    except Exception:  # pragma: no cover
        pass
    model = "unknown"
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    except Exception:  # pragma: no cover
        pass
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_model": model,
        "cpu_count": os.cpu_count(),
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "blas": blas,
    }


def run_case(case: dict, n_jobs: int) -> tuple[float, dict]:
    from spieceasi import spiec_easi_mb

    counts = np.loadtxt(case["counts"], delimiter="\t")
    idx = np.loadtxt(case["subs"], delimiter="\t").astype(int)
    subsamples = [idx[i] for i in range(idx.shape[0])]

    t0 = time.perf_counter()
    result = spiec_easi_mb(
        counts,
        n_lambdas=case["nlambda"],
        lambda_min_ratio=case["ratio"],
        n_subsamples=case["rep"],
        subsamples=subsamples,
        stars_threshold=0.05,
        sym_rule="or",
        random_state=17,
        max_iter=10_000,
        tol=1e-4,
        n_jobs=n_jobs,
        require_samples_gt_taxa=False,
    )
    elapsed = time.perf_counter() - t0
    n_edges = int(
        (np.count_nonzero(result.adjacency) - np.count_nonzero(np.diag(result.adjacency))) // 2
    )
    return elapsed, {
        "selected_lambda": float(result.selected_lambda),
        "edge_count": n_edges,
        "n_samples": int(counts.shape[0]),
        "n_taxa": int(counts.shape[1]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--only", default=None, help="substring filter on case id")
    args = ap.parse_args()

    env = environment()
    print("== environment ==")
    for key, value in env.items():
        print(f"  {key}: {value}")
    print(f"== python pipeline, n_jobs={args.n_jobs}, reps={args.reps} ==")

    out: dict = {"kind": "python", "environment": env, "n_jobs": args.n_jobs, "cases": {}}
    for case in CASES:
        if args.only and args.only not in case["id"]:
            continue
        times: list[float] = []
        info: dict = {}
        for r in range(args.reps):
            elapsed, info = run_case(case, args.n_jobs)
            times.append(elapsed)
            print(f"  {case['id']} rep {r + 1}/{args.reps}: {elapsed:.1f} s", flush=True)
        entry = {
            "min": min(times),
            "median": statistics.median(times),
            "times": times,
            **info,
        }
        out["cases"][case["id"]] = entry
        print(
            f"  => {case['id']}: min={entry['min']:.1f}s median={entry['median']:.1f}s "
            f"edges={info['edge_count']} lambda={info['selected_lambda']:.6g}",
            flush=True,
        )

    if args.json:
        args.json.write_text(json.dumps(out, indent=2))
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
