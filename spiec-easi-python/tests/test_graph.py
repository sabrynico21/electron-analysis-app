from __future__ import annotations

import numpy as np

from spieceasi.graph import to_networkx_graph


def test_to_networkx_graph_basic() -> None:
    adjacency = np.array(
        [
            [0, 1, 0],
            [1, 0, 1],
            [0, 1, 0],
        ],
        dtype=int,
    )
    weights = np.array(
        [
            [0.0, 0.4, 0.0],
            [0.4, 0.0, -0.2],
            [0.0, -0.2, 0.0],
        ]
    )
    G = to_networkx_graph(adjacency, weights, ["a", "b", "c"])
    assert G.number_of_nodes() == 3
    assert G.number_of_edges() == 2
    assert G.nodes[0]["label"] == "a"
