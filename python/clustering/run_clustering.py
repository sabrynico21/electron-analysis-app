#!/usr/bin/env python3
"""Run adaptive PageRank clustering over a cached condition graph payload."""

from __future__ import annotations

import argparse
import json
import sys

from bacteria_graph import BacteriaGraph
from bacteria_pagerank import pagerank_nibble_adaptive


def build_graph(payload: dict) -> BacteriaGraph:
    graph_payload = payload.get("graph", {})
    nodes_payload = graph_payload.get("nodes", [])
    edges_payload = graph_payload.get("edges", [])

    node_ids = [str(node.get("id")) for node in nodes_payload if node.get("id")]
    graph = BacteriaGraph.from_nodes_edges(node_ids, edges_payload)
    return graph


def extract_subgraph(graph_payload: dict, cluster_nodes: set[str], ranked_nodes: list[str] | None = None) -> dict:
    rank_order = ranked_nodes or []
    rank_index = {node_id: idx for idx, node_id in enumerate(rank_order)}
    nodes = [node for node in graph_payload.get("nodes", []) if node.get("id") in cluster_nodes]
    if rank_order:
        nodes.sort(key=lambda node: rank_index.get(node.get("id"), len(rank_index)))
    edges = [
        edge
        for edge in graph_payload.get("edges", [])
        if edge.get("source") in cluster_nodes and edge.get("target") in cluster_nodes
    ]
    return {
        "nodes": nodes,
        "edges": edges,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to clustering input JSON")
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    seed = payload.get("seedNode")
    if not seed:
        raise ValueError("seedNode is required")

    weighted_mode = bool(payload.get("weighted", False))

    graph_payload = payload.get("graphData", {})
    graph = build_graph(graph_payload)

    cluster_nodes, metadata = pagerank_nibble_adaptive(
        graph,
        seed,
        target_conductance=float(payload.get("targetConductance", 0.2)),
        b=(int(payload["b"]) if payload.get("b") is not None else None),
        teleport_alpha=(float(payload["teleportAlpha"]) if payload.get("teleportAlpha") is not None else None),
        weight="weight" if weighted_mode else None,
        verbose=False,
    )

    metadata = dict(metadata)
    metadata["weighted_mode"] = weighted_mode

    ranked_nodes = metadata.get("ranked_nodes") if isinstance(metadata, dict) else None
    subgraph = extract_subgraph(graph_payload.get("graph", {}), cluster_nodes, ranked_nodes)

    output = {
        "success": True,
        "seedNode": seed,
        "clusterNodes": sorted(list(cluster_nodes)),
        "subgraph": subgraph,
        "metadata": metadata,
    }

    print(json.dumps(output), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(json.dumps({"success": False, "error": str(exc)}), flush=True)
        sys.exit(1)
