#!/usr/bin/env python
"""Part A of the real-data R<->Python verification battery: prevalence/dimensions.

Reads the real count matrix (``genere_table.txt``, genera x samples) and the
sample metadata (``sample_metadatab.txt``) and, for a list of prevalence
thresholds, reports for each group (all / Treated / Not_Treated):

* ``p`` — number of genera present (count > 0) in at least ``threshold``
  fraction of the group's samples;
* ``zero_fraction`` — fraction of zeros in the resulting sub-matrix;
* ``n`` — number of samples in the group.

Prevalence of a genus = fraction of the group's samples with count > 0.

Outputs a JSON file (with a SHA-256 checksum of itself and of the inputs) and
prints a Markdown table to stdout. No production code is touched.

Usage (from ``spiec-easi-python/``):

    .venv/bin/python tests/real_data/estimate_dimensions.py \
        --out tests/real_data/dimensions.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_GENERA = DATA_DIR / "genere_table.txt"
DEFAULT_METADATA = DATA_DIR / "sample_metadatab.txt"
DEFAULT_THRESHOLDS = "0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def load_counts(genera_path: Path) -> pd.DataFrame:
    """Load the genera x samples count matrix (genera are rows)."""
    df = pd.read_csv(genera_path, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return df


def load_groups(metadata_path: Path, sample_order: list[str]) -> dict[str, list[str]]:
    meta = pd.read_csv(metadata_path, sep="\t")
    id_col = meta.columns[0]
    cond_col = "Condition"
    meta[id_col] = meta[id_col].astype(str)
    by_cond: dict[str, list[str]] = {}
    for cond, sub in meta.groupby(cond_col):
        by_cond[str(cond)] = [s for s in sub[id_col].tolist() if s in set(sample_order)]
    groups = {"all": list(sample_order)}
    groups.update(by_cond)
    return groups


def dimensions_for_group(
    counts: pd.DataFrame, samples: list[str], thresholds: list[float]
) -> list[dict]:
    sub = counts[samples]
    present = (sub.to_numpy() > 0).astype(float)  # genera x samples
    prevalence = present.mean(axis=1)  # fraction of samples where genus is present
    rows = []
    for thr in thresholds:
        mask = prevalence >= thr
        kept = sub.loc[mask]
        rows.append(
            {
                "threshold": thr,
                "n": len(samples),
                "p": int(mask.sum()),
                "zero_fraction": (
                    float((kept.to_numpy() == 0).mean()) if kept.size else float("nan")
                ),
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genera", type=Path, default=DEFAULT_GENERA)
    ap.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    ap.add_argument("--thresholds", default=DEFAULT_THRESHOLDS)
    ap.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "spiec-easi-python/tests/real_data/dimensions.json",
    )
    args = ap.parse_args()

    thresholds = [float(t) for t in args.thresholds.split(",") if t.strip()]

    counts = load_counts(args.genera)
    groups = load_groups(args.metadata, list(counts.columns))

    result: dict = {
        "genera_file": str(args.genera),
        "metadata_file": str(args.metadata),
        "genera_sha256": sha256_file(args.genera),
        "metadata_sha256": sha256_file(args.metadata),
        "n_samples": int(counts.shape[1]),
        "p_genera_total": int(counts.shape[0]),
        "group_sizes": {g: len(s) for g, s in groups.items()},
        "thresholds": thresholds,
        "rows": [],
    }

    for group in ["all", *[g for g in groups if g != "all"]]:
        samples = groups[group]
        for row in dimensions_for_group(counts, samples, thresholds):
            result["rows"].append({"group": group, **row})

    # Markdown table
    lines = [
        "| group | threshold | n | p | zero_fraction |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in result["rows"]:
        lines.append(
            f"| {row['group']} | {row['threshold']:.2f} | {row['n']} | {row['p']} | "
            f"{row['zero_fraction']:.4f} |"
        )
    table = "\n".join(lines)
    print(table)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    result["out_sha256"] = sha256_file(args.out)
    print(f"\nSHA-256 {args.genera.name}: {result['genera_sha256']}", flush=True)
    print(f"SHA-256 {args.metadata.name}: {result['metadata_sha256']}")
    print(f"SHA-256 {args.out}: {result['out_sha256']}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
