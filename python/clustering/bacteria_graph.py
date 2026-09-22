"""Utility classes for building bacteria co-occurrence graphs.

The new dataset expresses weighted relationships between bacteria species.
Each row of the input text file is expected to follow the schema::

    row_id    source_name    target_name    correlation    p_value

Values can be separated via spaces or tabs. During ingestion we keep the
original correlation (possibly negative) and derive a non-negative weight
based on its absolute value. This makes the graph compatible with the
PageRank Nibble routine that assumes positive transition probabilities.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import networkx as nx


@dataclass(frozen=True)
class EdgeRecord:
    """Structured representation of a single bacteria relationship."""

    row_id: str
    source: str
    target: str
    correlation: float
    p_value: float


class BacteriaGraph:
    """Graph wrapper tailored for bacteria-only metadata.

    The class exposes just the subset of helpers that the PageRank Nibble
    implementation and the PyQt interface rely on (neighbors, weights, degree
    lookups, etc.).
    """

    def __init__(self) -> None:
        self.graph: nx.Graph = nx.Graph()
        self._degree = defaultdict(int)
        self._weighted_degree = defaultdict(float)

    @classmethod
    def from_edge_file(
        cls,
        path: Path | str,
        *,
        min_corr: float = 0.0,
        max_p_value: Optional[float] = None,
        use_abs_weight: bool = True,
    ) -> "BacteriaGraph":
        graph = cls()
        graph.load_from_file(
            path,
            min_corr=min_corr,
            max_p_value=max_p_value,
            use_abs_weight=use_abs_weight,
        )
        return graph

    @classmethod
    def from_nodes_edges(cls, nodes: Iterable[str], edges: Iterable[dict]) -> "BacteriaGraph":
        graph = cls()
        for node in nodes:
            graph.graph.add_node(str(node))
        for edge in edges:
            source = str(edge.get("source"))
            target = str(edge.get("target"))
            if not source or not target:
                continue
            correlation = edge.get("correlation", edge.get("weight", 0.0))
            graph.add_edge(
                row_id=str(edge.get("row_id", "")),
                source=source,
                target=target,
                correlation=float(correlation),
                p_value=float(edge.get("p_value", 0.0)),
                use_abs_weight=True,
            )
        return graph

    def load_from_file(
        self,
        path: Path | str,
        *,
        min_corr: float = 0.0,
        max_p_value: Optional[float] = None,
        use_abs_weight: bool = True,
    ) -> None:
        """Populate the graph from a whitespace-delimited text file."""

        file_path = Path(path)
        if not file_path.exists():
            raise FileNotFoundError(f"Edge file not found: {file_path}")

        with file_path.open("r", encoding="utf-8") as handle:
            header = handle.readline().strip()
            if not header:
                raise ValueError("The edge file appears to be empty.")

            for line in handle:
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) < 5:
                    raise ValueError(
                        "Malformed row. Expected: row_id source target correlation p-value"
                    )
                row_id, source, target, corr_raw, p_raw = parts[:5]
                correlation = float(corr_raw)
                p_value = float(p_raw)

                if abs(correlation) < min_corr:
                    continue
                if max_p_value is not None and p_value > max_p_value:
                    continue
                self.add_edge(
                    row_id=row_id,
                    source=source,
                    target=target,
                    correlation=correlation,
                    p_value=p_value,
                    use_abs_weight=use_abs_weight,
                )

    def add_edge(
        self,
        *,
        row_id: str,
        source: str,
        target: str,
        correlation: float,
        p_value: float,
        use_abs_weight: bool = True,
    ) -> None:
        weight = abs(correlation) if use_abs_weight else correlation
        if weight <= 0:
            weight = abs(correlation)

        existing = self.graph.get_edge_data(source, target, default=None)
        if existing:
            # Keep the strongest relationship observed so far.
            if weight <= existing.get("weight", 0.0):
                return

        self.graph.add_edge(
            source,
            target,
            row_id=row_id,
            correlation=correlation,
            p_value=p_value,
            weight=weight,
        )

        self._degree[source] += 1
        self._degree[target] += 1
        self._weighted_degree[source] += weight
        self._weighted_degree[target] += weight

    # ------------------------------------------------------------------
    # Graph helpers consumed by PageRank Nibble and the GUI
    # ------------------------------------------------------------------
    def get_neighbors(self, node: str) -> List[str]:
        return list(self.graph.neighbors(node)) if self.graph.has_node(node) else []

    def get_weight(self, source: str, target: str) -> float:
        data = self.graph.get_edge_data(source, target)
        return float(data.get("weight", 0.0)) if data else 0.0

    def get_correlation(self, source: str, target: str) -> Optional[float]:
        """Return the stored correlation between two nodes, if available."""

        data = self.graph.get_edge_data(source, target)
        if not data:
            return None
        try:
            return float(data.get("correlation"))
        except (TypeError, ValueError):
            return None

    def get_deg(self, nodes: Optional[str] = None, weight: Optional[str] = None):
        store = self._weighted_degree if weight == "weight" else self._degree
        if nodes is None:
            return store
        return store[nodes]

    def has_node(self, node: str) -> bool:
        return self.graph.has_node(node)

    def number_of_nodes(self) -> int:
        return self.graph.number_of_nodes()

    def number_of_edges(self) -> int:
        return self.graph.number_of_edges()

    def get_subgraph(self, nodes: Sequence[str]) -> nx.Graph:
        return self.graph.subgraph(nodes).copy()

    # ------------------------------------------------------------------
    # Metadata helpers for the GUI layer
    # ------------------------------------------------------------------
    def available_nodes(self) -> List[str]:
        return sorted(self.graph.nodes())

    def node_stats(self, node: str) -> Dict[str, float | int | List[Tuple[str, float]]]:
        neighbors = self.get_neighbors(node)
        weighted_neighbors = [(nb, self.get_weight(node, nb)) for nb in neighbors]
        weighted_neighbors.sort(key=lambda item: item[1], reverse=True)
        return {
            "degree": self._degree[node],
            "weighted_degree": self._weighted_degree[node],
            "top_neighbors": weighted_neighbors[:10],
        }

    def edge_records(self, nodes: Iterable[str]) -> Iterator[EdgeRecord]:
        seen = set()
        for u in nodes:
            for v, data in self.graph[u].items():
                pair = tuple(sorted((u, v)))
                if pair in seen:
                    continue
                seen.add(pair)
                yield EdgeRecord(
                    row_id=str(data.get("row_id", "")),
                    source=u,
                    target=v,
                    correlation=float(data.get("correlation", 0.0)),
                    p_value=float(data.get("p_value", 0.0)),
                )
