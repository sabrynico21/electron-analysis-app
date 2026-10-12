"""Graph conversion utilities for inferred microbial networks."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    import igraph as ig


def to_networkx_graph(
    adjacency: np.ndarray,
    weights: np.ndarray,
    taxa_names: Sequence[str],
) -> nx.Graph:
    """Convert adjacency and weight matrices into a networkx graph."""
    n_taxa = adjacency.shape[0]
    if adjacency.shape != (n_taxa, n_taxa) or weights.shape != (n_taxa, n_taxa):
        msg = "adjacency and weights must be square matrices of identical shape."
        raise ValueError(msg)
    if len(taxa_names) != n_taxa:
        msg = "taxa_names length must match matrix size."
        raise ValueError(msg)

    G: nx.Graph = nx.Graph()
    for i, name in enumerate(taxa_names):
        G.add_node(i, label=str(name))

    upper = np.triu_indices(n_taxa, k=1)
    for row, col in zip(*upper, strict=True):
        if adjacency[row, col] != 0:
            w = float(weights[row, col])
            G.add_edge(int(row), int(col), weight=w, abs_weight=abs(w))
    return G


def to_igraph_graph(
    adjacency: np.ndarray,
    weights: np.ndarray,
    taxa_names: Sequence[str],
) -> ig.Graph:
    """Convert adjacency and weights into an igraph graph.

    Raises
    ------
    ImportError
        If python-igraph is not installed.
    """
    try:
        import igraph as ig
    except ImportError as exc:
        msg = "python-igraph is required for to_igraph_graph()."
        raise ImportError(msg) from exc

    n_taxa = adjacency.shape[0]
    g = ig.Graph(n=n_taxa, directed=False)
    g.vs["label"] = list(taxa_names)

    edges: list[tuple[int, int]] = []
    edge_weights: list[float] = []
    upper = np.triu_indices(n_taxa, k=1)
    for i, j in zip(*upper, strict=True):
        if adjacency[i, j] != 0:
            edges.append((int(i), int(j)))
            edge_weights.append(float(weights[i, j]))

    g.add_edges(edges)
    g.es["weight"] = edge_weights
    g.es["abs_weight"] = [abs(w) for w in edge_weights]
    return g
