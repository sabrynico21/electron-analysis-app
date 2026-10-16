#!/usr/bin/env python
"""Part B timing calibration: measure Python wall-clock time for one seed.

Runs the full default Python protocol on a real-data sub-case (a
prevalence-threshold genus subset of the real matrix) and prints the elapsed
seconds, so that the per-case Python time of the battery can be estimated.

Protocol matches the R reference defaults: ``n_lambdas=20``,
``lambda_min_ratio=1e-3``, ``n_subsamples=20``, ``stars_threshold=0.05``,
``sym_rule="or"``, ``max_iter=10000``, ``tol=1e-4``. ``random_state`` is the
seed (the battery injects R's subsample indices, which yields the same amount
of work).

Usage (from ``spiec-easi-python/``):

    .venv/bin/python tests/real_data/calibrate_python_time.py \
        --threshold 0.70 [--group all] [--seed 17]
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_GENERA = DATA_DIR / "genere_table.txt"
DEFAULT_METADATA = DATA_DIR / "sample_metadatab.txt"


def load_subcase(genera_path: Path, metadata_path: Path, group: str, threshold: float):
    counts = (
        pd.read_csv(genera_path, sep="\t", index_col=0)
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
    )
    if group == "all":
        samples = list(counts.columns)
    else:
        meta = pd.read_csv(metadata_path, sep="\t")
        meta.columns = [c.strip() for c in meta.columns]
        ids = meta.loc[meta["Condition"] == group, meta.columns[0]].astype(str)
        samples = [s for s in counts.columns if s in set(ids)]
    sub = counts[samples]
    prevalence = (sub.to_numpy() > 0).mean(axis=1)
    kept = sub.loc[prevalence >= threshold]
    # samples x genera
    return kept.to_numpy().T.astype(float), samples


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--group", default="all")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--genera", type=Path, default=DEFAULT_GENERA)
    ap.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    args = ap.parse_args()

    from spieceasi.pipeline import spiec_easi_mb

    Y, samples = load_subcase(args.genera, args.metadata, args.group, args.threshold)
    n, p = Y.shape
    zero_frac = float((Y == 0).mean())

    t0 = time.perf_counter()
    result = spiec_easi_mb(
        Y,
        n_lambdas=20,
        lambda_min_ratio=1e-3,
        n_subsamples=20,
        stars_threshold=0.05,
        sym_rule="or",
        max_iter=10_000,
        tol=1e-4,
        random_state=args.seed,
        n_jobs=-1,
        require_samples_gt_taxa=False,
    )
    dt = time.perf_counter() - t0

    opt_index = int(np.argmin(np.abs(result.alpha_path - result.selected_lambda)))
    n_edges = int(result.adjacency.sum() // 2)
    print(
        f"PYCALIB threshold={args.threshold:.2f} group={args.group} seed={args.seed} "
        f"n={n} p={p} zero_frac={zero_frac:.4f} opt_index={opt_index} "
        f"n_edges={n_edges} elapsed_s={dt:.2f}"
    )


if __name__ == "__main__":
    main()
