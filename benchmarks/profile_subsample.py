#!/usr/bin/env python3
"""FASE 1 — profile a single MB path solve (one StARS subsample).

Profiles ``solve_mb_path_huge_scale`` on one subsample of a large-p golden case
with cProfile (always) and line_profiler (when installed).

Usage:
    ./.venv/bin/python benchmarks/profile_subsample.py [--case largep_p100_n80_s215]
                                                       [--subsample 0] [--line]
"""

from __future__ import annotations

import argparse
import cProfile
import io
import pstats
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "spiec-easi-python" / "src"))

GOLD_EXT = REPO / "spiec-easi-python" / "tests" / "golden_extended"
GOLD = REPO / "spiec-easi-python" / "tests" / "golden"


def load_case(case: str):
    if case == "ds2_p30_n100":
        clr = np.loadtxt(GOLD / "ds2_r.clr.tsv", delimiter="\t")
        grid = np.loadtxt(GOLD / "ds2_r.lambda_path.txt")
        idx = np.loadtxt(GOLD / "ds2_r.subsample_indices.tsv", delimiter="\t").astype(int)
    else:
        cdir = GOLD_EXT / case
        clr = np.loadtxt(cdir / "clr.tsv", delimiter="\t")
        grid = np.loadtxt(cdir / "lambda_path.txt")
        idx = np.loadtxt(cdir / "subsample_indices.tsv", delimiter="\t").astype(int)
    return clr, grid, idx


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="largep_p100_n80_s215")
    ap.add_argument("--subsample", type=int, default=0)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--line", action="store_true")
    args = ap.parse_args()

    clr, grid, idx = load_case(args.case)
    Xs = clr[idx[args.subsample]]

    from spieceasi.mb_solver import solve_mb_path_huge_scale

    print(
        f"case={args.case} subsample={args.subsample} "
        f"n={Xs.shape[0]} p={Xs.shape[1]} lambdas={grid.size}"
    )

    t0 = time.perf_counter()
    directed = solve_mb_path_huge_scale(Xs, grid, max_iter=10_000, tol=1e-4)
    t1 = time.perf_counter()
    print(f"single subsample solve: {t1 - t0:.3f} s   (output {directed.shape})")

    profiler = cProfile.Profile()
    profiler.enable()
    solve_mb_path_huge_scale(Xs, grid, max_iter=10_000, tol=1e-4)
    profiler.disable()

    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).sort_stats("tottime")
    stats.print_stats(args.top)
    print("== cProfile (sorted by tottime) ==")
    print(stream.getvalue())

    cum = io.StringIO()
    pstats.Stats(profiler, stream=cum).sort_stats("cumulative").print_stats(15)
    print("== cProfile (sorted by cumulative) ==")
    print(cum.getvalue())

    if args.line:
        try:
            from line_profiler import LineProfiler
        except ImportError:
            print("line_profiler not installed -> skipping line profile")
        else:
            lp = LineProfiler()
            lp.add_function(solve_mb_path_huge_scale)
            lp.runcall(solve_mb_path_huge_scale, Xs, grid, max_iter=10_000, tol=1e-4)
            print("== line_profiler ==")
            lp.print_stats()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
