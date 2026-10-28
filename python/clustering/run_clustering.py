#!/usr/bin/env python3
"""Run PageRank-Nibble clustering over a cached condition graph payload.

Payload keys
------------
``seedNode`` (required), ``graphData``, ``weighted``,
``targetConductance``, ``b``, ``teleportAlpha``, ``signed``,
``signLambda`` and the optional ``parameterOverrides``.

``parameterProfile`` accepts only ``"app"`` (the default): the former
``recommended`` profile was removed; and ``selection`` is ignored: the
five-step ladder of ``selection_ladder.py`` is the only selection rule.
"""

from __future__ import annotations

import argparse
import json
import math
import sys

from bacteria_graph import BacteriaGraph
from bacteria_pagerank import pagerank_nibble_adaptive
from clustering_advanced import pagerank_cluster
from clustering_profiles import DEFAULT_TARGET_CONDUCTANCE, resolve_parameters


def _json_safe(value):
    """Recursively replace non-finite floats with ``None``.

    JSON has no ``Infinity``/``NaN`` (RFC 8259) and the Electron side parses
    this stream with ``JSON.parse``: emitting ``Infinity`` (the honest Python
    representation of an undefined conductance, e.g. the seed-singleton
    fallback) would make the whole response unreadable for the app.
    """
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


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
    profile = str(payload.get("parameterProfile", "app"))
    signed = bool(payload.get("signed", False))
    sign_lambda = float(payload.get("signLambda", 1.0))

    graph_payload = payload.get("graphData", {})
    graph = build_graph(graph_payload)

    overrides = {}
    if payload.get("targetConductance") is not None:
        overrides["target_conductance"] = float(payload["targetConductance"])
    if payload.get("b") is not None:
        overrides["b"] = int(payload["b"])
    if payload.get("teleportAlpha") is not None:
        overrides["teleport_alpha"] = float(payload["teleportAlpha"])
    for key, value in (payload.get("parameterOverrides") or {}).items():
        if value is not None:
            overrides[key] = value
    params = resolve_parameters(
        graph, profile, weighted=weighted_mode, overrides=overrides
    )

    # Boundary hygiene: an empty conductance field in the UI reaches this
    # process as 0 (``Number('') === 0`` in JS), and every ladder step
    # rejects a non-positive bar - only the singleton would survive.  Fall
    # back to the app default and record the raw value for transparency.
    phi_clamped_from = None
    if not (params["target_conductance"] > 0):
        phi_clamped_from = params["target_conductance"]
        params["target_conductance"] = DEFAULT_TARGET_CONDUCTANCE
    alpha_clamped_from = None
    if params["teleport_alpha"] is not None and params["teleport_alpha"] < 0:
        alpha_clamped_from = params["teleport_alpha"]
        params["teleport_alpha"] = None

    if signed:
        cluster_nodes, metadata = pagerank_cluster(
            graph,
            seed,
            target_conductance=params["target_conductance"],
            b=params["b"],
            teleport_alpha=params["teleport_alpha"],
            weighted=weighted_mode,
            signed=True,
            sign_lambda=sign_lambda,
        )
    else:
        cluster_nodes, metadata = pagerank_nibble_adaptive(
            graph,
            seed,
            target_conductance=params["target_conductance"],
            b=params["b"],
            teleport_alpha=params["teleport_alpha"],
            weight=params["weight"],
            verbose=False,
        )

    metadata = dict(metadata)
    metadata["weighted_mode"] = weighted_mode
    metadata["parameters"] = params
    metadata["signed_mode"] = signed
    if phi_clamped_from is not None:
        metadata["phi_clamped_from"] = phi_clamped_from
    if alpha_clamped_from is not None:
        metadata["alpha_clamped_from"] = alpha_clamped_from

    ranked_nodes = metadata.get("ranked_nodes") if isinstance(metadata, dict) else None
    subgraph = extract_subgraph(graph_payload.get("graph", {}), cluster_nodes, ranked_nodes)

    output = {
        "success": True,
        "seedNode": seed,
        "clusterNodes": sorted(list(cluster_nodes)),
        "subgraph": subgraph,
        "metadata": metadata,
    }

    print(json.dumps(_json_safe(output), allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(json.dumps({"success": False, "error": str(exc)}), flush=True)
        sys.exit(1)
