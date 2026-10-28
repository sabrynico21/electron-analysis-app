"""Load the two static test networks used by the clustering experiments.

The files in ``data/`` are the exported ``*_full_graph_edges.csv`` edge
lists (one row per pair: condition, source, target, correlation, ...).  They
are loaded into the same ``BacteriaGraph`` wrapper the Electron app passes to
``run_clustering.py``, so the experiments exercise the algorithm and only the
algorithm.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Optional

TEST_DIR = Path(__file__).resolve().parent
CLUSTERING_DIR = TEST_DIR.parent
if str(CLUSTERING_DIR) not in sys.path:
    sys.path.insert(0, str(CLUSTERING_DIR))

from bacteria_graph import BacteriaGraph  # noqa: E402

DATA_DIR = TEST_DIR / "data"
NETWORKS = {
    "rs4-101": DATA_DIR / "RS4_101_Not_Treated_full_graph_edges.csv",
    "101": DATA_DIR / "101_Not_Treated_full_graph_edges.csv",
}

REQUIRED_COLUMNS = ("source", "target", "correlation")


def resolve_network(name_or_path: str) -> Path:
    """Accept a network alias (``101``, ``rs4-101``) or a direct file path."""

    key = name_or_path.strip().lower()
    if key in NETWORKS:
        return NETWORKS[key]
    path = Path(name_or_path).expanduser()
    if not path.is_absolute():
        candidate = TEST_DIR / path
        if candidate.exists():
            return candidate
    if path.exists():
        return path
    aliases = ", ".join(sorted(NETWORKS))
    raise FileNotFoundError(
        f"Unknown network {name_or_path!r}. Use one of: {aliases}, or a CSV path."
    )


def load_network(name_or_path: str) -> tuple[BacteriaGraph, dict]:
    """Build a ``BacteriaGraph`` from an exported edge CSV.

    Returns ``(graph, info)`` where ``info`` records the file and the row
    counts.  The edge weight follows the app semantics: absolute correlation
    (``BacteriaGraph.add_edge`` uses ``use_abs_weight=True``); unweighted runs
    simply ignore the weight.
    """

    path = resolve_network(name_or_path)
    nodes: set[str] = set()
    edges: list[dict] = []
    skipped_self_loops = 0
    skipped_blank = 0

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
        if missing:
            raise ValueError(
                f"{path.name}: missing required columns {missing}; found {fieldnames}"
            )

        for row in reader:
            source = (row.get("source") or "").strip()
            target = (row.get("target") or "").strip()
            correlation_raw = (row.get("correlation") or "").strip()
            if not source or not target or not correlation_raw:
                skipped_blank += 1
                continue
            if source == target:
                skipped_self_loops += 1
                continue
            p_raw = (row.get("p_value") or "").strip()
            p_value: Optional[float] = float(p_raw) if p_raw else None
            nodes.add(source)
            nodes.add(target)
            edges.append(
                {
                    "source": source,
                    "target": target,
                    "correlation": float(correlation_raw),
                    "p_value": p_value,
                }
            )

    graph = BacteriaGraph.from_nodes_edges(sorted(nodes), edges)
    info = {
        "file": str(path),
        "rows_used": len(edges),
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "skipped_self_loops": skipped_self_loops,
        "skipped_blank": skipped_blank,
    }
    return graph, info
