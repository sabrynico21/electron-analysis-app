"""Signed scoring for PageRank-Nibble and the advanced clustering entry.

Two independent concerns live in this module:

* ``signed=True`` - positive correlations drive the random walk
  (``A+ = max(A, 0)``) while negative correlations are used only as
  violations in the cluster score:
  ``(cut_pos + lambda * neg_inside) / min(vol_pos(S), vol_pos(V-S))``.
* :func:`pagerank_cluster` - the clustering entry used for signed runs
  (unsigned runs go through ``bacteria_pagerank.pagerank_nibble_adaptive``).
  Both entries execute the *same* five-step selection ladder of
  ``selection_ladder.py``; the ladder replaced the former
  ``selection="first" / "min_conductance"`` rules, and the profile
  ``recommended`` was removed together with them.

The plain ``|correlation|`` weighting and the parameter defaults of the
Electron app are unchanged: see ``clustering_profiles.py`` and
``tests/docs/ALGORITHM.md``; the ladder itself is documented in
``tests/docs/SELECTION_LADDER.md``.
"""

from __future__ import annotations

from typing import Optional

try:
    from .bacteria_graph import BacteriaGraph
    from .selection_ladder import (
        DEFAULT_MIN_CUT_SIZE,
        SignedScorer,
        cluster_page_rank_nibble,
    )
except ImportError:
    from bacteria_graph import BacteriaGraph
    from selection_ladder import (
        DEFAULT_MIN_CUT_SIZE,
        SignedScorer,
        cluster_page_rank_nibble,
    )


# ---------------------------------------------------------------------------
# Signed graph helpers
# ---------------------------------------------------------------------------
def build_signed_index(graph: BacteriaGraph, weighted: bool = True) -> dict:
    """Pre-compute positive/negative adjacency lists and degrees.

    ``weighted=False`` treats every edge as weight 1, both in the walk and in
    the signed score.
    """

    pos_adj: dict[str, list[tuple[str, float]]] = {node: [] for node in graph.available_nodes()}
    neg_adj: dict[str, list[tuple[str, float]]] = {node: [] for node in graph.available_nodes()}
    pos_deg = {node: 0.0 for node in pos_adj}
    neg_deg = {node: 0.0 for node in pos_adj}
    positive_edges = 0
    negative_edges = 0

    for source, target, data in graph.graph.edges(data=True):
        correlation = float(data.get("correlation", 0.0))
        if correlation > 0:
            weight = correlation if weighted else 1.0
            pos_adj[source].append((target, weight))
            pos_adj[target].append((source, weight))
            pos_deg[source] += weight
            pos_deg[target] += weight
            positive_edges += 1
        elif correlation < 0:
            weight = -correlation if weighted else 1.0
            neg_adj[source].append((target, weight))
            neg_adj[target].append((source, weight))
            neg_deg[source] += weight
            neg_deg[target] += weight
            negative_edges += 1

    return {
        "pos_adj": pos_adj,
        "neg_adj": neg_adj,
        "pos_deg": pos_deg,
        "neg_deg": neg_deg,
        "total_pos_volume": sum(pos_deg.values()),
        "total_neg_volume": sum(neg_deg.values()),
        "positive_edges": positive_edges,
        "negative_edges": negative_edges,
    }


def signed_metrics(index: dict, nodes: set[str], sign_lambda: float = 1.0) -> dict:
    """Signed conductance and sign composition of a candidate cluster.

    Thin wrapper over :meth:`selection_ladder.SignedScorer.composition`,
    which is the single implementation of the formulas (it is also used
    by the sweep itself and by the run metadata).
    """

    return SignedScorer(index, sign_lambda).composition(nodes)


def compute_signed_metrics(
    graph: BacteriaGraph,
    nodes: set[str],
    sign_lambda: float = 1.0,
    weighted: bool = True,
) -> dict:
    return signed_metrics(
        build_signed_index(graph, weighted=weighted), nodes, sign_lambda
    )


def positive_walk_graph(graph: BacteriaGraph) -> BacteriaGraph:
    """Graph used by the random walk when signs are enabled: positive edges only."""

    nodes = graph.available_nodes()
    edges = []
    for source, target, data in graph.graph.edges(data=True):
        correlation = float(data.get("correlation", 0.0))
        if correlation > 0:
            edges.append(
                {
                    "source": source,
                    "target": target,
                    "correlation": correlation,
                    "p_value": data.get("p_value"),
                }
            )
    return BacteriaGraph.from_nodes_edges(nodes, edges)


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------
def pagerank_cluster(
    graph: BacteriaGraph,
    seed: str,
    *,
    target_conductance: float,
    b: int,
    teleport_alpha: float,
    weighted: bool = True,
    signed: bool = False,
    sign_lambda: float = 1.0,
    random_seed: Optional[int] = None,
    max_iterations: int = 100_000,
    min_cut_size: int = DEFAULT_MIN_CUT_SIZE,
    max_cut_size: Optional[int] = None,
) -> tuple[set[str], dict]:
    """Cluster around ``seed`` with the five-step ladder (signed or not).

    ``signed=True`` runs the walk on the positive sub-graph and replaces
    the conductance by the signed score everywhere in the ladder
    (acceptance, relaxation rungs, penalty).  All other behaviour -
    the sweep scan, the five steps and the metadata schema - is shared
    with ``bacteria_pagerank.pagerank_nibble_adaptive`` through
    ``selection_ladder.cluster_page_rank_nibble``.
    """

    if not graph.has_node(seed):
        raise KeyError(f"Unknown seed node: {seed}")

    weight = "weight" if weighted else None
    walk_graph = None
    scorer = None
    if signed:
        walk_graph = positive_walk_graph(graph)
        scorer = SignedScorer(
            build_signed_index(graph, weighted=weighted), sign_lambda
        )

    return cluster_page_rank_nibble(
        graph,
        seed,
        phi=target_conductance,
        b=b,
        alpha=teleport_alpha,
        weight=weight,
        walk_graph=walk_graph,
        scorer=scorer,
        signed=signed,
        sign_lambda=sign_lambda,
        min_cut_size=min_cut_size,
        max_cut_size=max_cut_size,
        max_iterations=max_iterations,
        random_seed=random_seed,
    )
