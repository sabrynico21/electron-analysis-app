#!/usr/bin/env python3
"""FASE 1 — cold/warm start-up timing of the analysis backend.

Measures, on a tiny dataset (p=5, n=20) in ``benchmarks/tiny/``:

* interpreter + PyInstaller bootstrap only (backend invoked with no arguments);
* a full ``analyze`` run (bootstrap + imports + computation).

Both are done for the *frozen* backend (``python-dist/…/analysis-backend``) and
for the *source* backend (``python/backend.py``), so the PyInstaller overhead can
be separated from the Python/numpy import cost. Exit codes are recorded so a
crashing configuration is not silently reported as a fast one.

Run it **alone** (no other load) for meaningful numbers.

Usage:
    ./.venv/bin/python benchmarks/bench_startup.py [--params FILE] [--repeats 5]
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FROZEN = REPO / "python-dist" / "analysis-backend" / "analysis-backend"
SOURCE = REPO / "python" / "backend.py"


def run(cmd: list[str], repeats: int, label: str) -> dict:
    times: list[float] = []
    codes: list[int] = []
    for i in range(repeats):
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
        elapsed = time.perf_counter() - t0
        times.append(elapsed)
        codes.append(proc.returncode)
        if proc.returncode not in (0, 2):
            tail = (proc.stderr or "").strip().splitlines()[-1:] or ["(no stderr)"]
            print(f"  !! {label} rep {i + 1}: rc={proc.returncode}  {tail[0][:110]}")
    warm = [t for t, c in zip(times[1:], codes[1:]) if c in (0, 2)] or times[1:]
    return {
        "cold": times[0],
        "warm_median": statistics.median(warm) if warm else times[0],
        "min": min(times),
        "times": times,
        "exit_codes": codes,
        "all_ok": all(c in (0, 2) for c in codes),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", type=Path, default=REPO / "benchmarks" / "tiny" / "params_jobs1.json")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    out_json = REPO / "benchmarks" / "tiny" / "result.json"
    params_path = args.params if args.params.is_absolute() else (REPO / args.params)
    analyze = ["analyze", "--params", str(args.params), "--output", str(out_json)]
    try:
        label = str(params_path.resolve().relative_to(REPO))
    except ValueError:
        label = str(params_path)
    results: dict = {"params": label}

    print(f"== params: {results['params']} ==")
    if FROZEN.is_file():
        print("== frozen (PyInstaller) ==")
        results["frozen_startup_only"] = run([str(FROZEN)], args.repeats, "frozen startup")
        results["frozen_analyze_tiny"] = run([str(FROZEN), *analyze], args.repeats, "frozen analyze")
    else:
        print(f"!! frozen backend not found at {FROZEN}")

    print("== source (python/backend.py) ==")
    results["source_startup_only"] = run([sys.executable, str(SOURCE)], args.repeats, "source startup")
    results["source_analyze_tiny"] = run(
        [sys.executable, str(SOURCE), *analyze], args.repeats, "source analyze"
    )

    print("== summary (warm median) ==")
    for name, entry in results.items():
        if not isinstance(entry, dict):
            continue
        flag = "" if entry["all_ok"] else "   <-- NOT ALL RUNS OK"
        print(f"  {name}: {entry['warm_median']:.3f} s{flag}")
    if "frozen_analyze_tiny" in results and "source_analyze_tiny" in results:
        delta = results["frozen_analyze_tiny"]["warm_median"] - results["source_analyze_tiny"]["warm_median"]
        print(f"  frozen - source (analyze): {delta:+.3f} s")
    if "source_analyze_tiny" in results:
        compute = results["source_analyze_tiny"]["warm_median"] - results["source_startup_only"]["warm_median"]
        print(f"  source: analyze - startup (imports + compute): {compute:.3f} s")

    if args.json:
        args.json.write_text(json.dumps(results, indent=2))
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
