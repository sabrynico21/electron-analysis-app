"""PageRank-Nibble implementation based on the original paper.

Reference:
    Andersen, R., Chung, F., & Lang, K. (2006).
    Local Graph Partitioning using PageRank Vectors.
    https://www.cs.cmu.edu/~15859n/RelatedWork/local_partitioning_full.pdf
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
import random
from typing import Dict, Optional, Set, Tuple

try:
    from .bacteria_graph import BacteriaGraph
except ImportError:
    from bacteria_graph import BacteriaGraph


@dataclass(frozen=True)
class PageRankResult:
    """Result from computing an approximate PageRank vector."""
    seed: str
    scores: Dict[str, float]
    residuals: Dict[str, float]
    epsilon: float
    alpha: float
    iterations: int

    def top_nodes(self, k: int) -> list[tuple[str, float]]:
        ordered = sorted(self.scores.items(), key=lambda item: item[1], reverse=True)
        return ordered[:k]

    def get_sweep_sets(self, graph: BacteriaGraph) -> list[tuple[str, float, Set[str]]]:
        deg_store = graph.get_deg(weight=None)
        ratios = []
        for node, pr_value in self.scores.items():
            degree = max(deg_store.get(node, 1), 1)
            ratio = pr_value / degree
            ratios.append((node, ratio))
        ratios.sort(key=lambda x: x[1], reverse=True)

        sweep_sets = []
        current_set = set()
        for node, ratio in ratios:
            current_set.add(node)
            sweep_sets.append((node, ratio, current_set.copy()))

        return sweep_sets


def compute_epsilon_paper(b: int, m: int) -> float:
    log_m_ceil = math.ceil(math.log2(max(m, 2)))
    return 1.0 / (2**b * 48 * log_m_ceil)


def compute_alpha_paper(phi: float, m: int) -> float:
    return phi**2 / (225 * math.log(100 * math.sqrt(m)))


def approximate_pagerank(
    graph: BacteriaGraph,
    seed: str,
    *,
    alpha: float = 0.15,
    epsilon: Optional[float] = None,
    max_iterations: int = 100_000,
    weight: Optional[str] = None,
    random_seed: Optional[int] = None,
) -> PageRankResult:
    rng = random.Random(random_seed) if random_seed is not None else None

    if not graph.has_node(seed):
        raise KeyError(f"Unknown seed node: {seed}")

    deg_store = graph.get_deg(weight=weight)

    if epsilon is None:
        m = max(graph.number_of_edges(), 1)
        degree = max(deg_store.get(seed, 1), 1)
        b = max(1, int(math.log2(degree)))
        epsilon = compute_epsilon_paper(b, m)

    p: Dict[str, float] = {}
    r: Dict[str, float] = {seed: 1.0}

    queue: deque[str] = deque([seed])
    in_queue: Set[str] = {seed}

    iterations = 0

    while queue and iterations < max_iterations:
        if rng is None:
            u = queue.popleft()
        else:
            idx = rng.randrange(len(queue))
            queue.rotate(-idx)
            u = queue.popleft()
        in_queue.discard(u)

        r_u = r.get(u, 0.0)
        d_u = max(deg_store.get(u, 1), 1)

        if r_u < epsilon * d_u:
            continue

        iterations += 1

        p[u] = p.get(u, 0.0) + alpha * r_u

        new_r_u = (1 - alpha) * r_u / 2
        r[u] = new_r_u

        if new_r_u >= epsilon * d_u and u not in in_queue:
            queue.append(u)
            in_queue.add(u)

        neighbors = graph.get_neighbors(u)
        if neighbors and rng is not None:
            neighbors = list(neighbors)
            rng.shuffle(neighbors)
        if neighbors:
            if weight is None:
                mass_per_neighbor = (1 - alpha) * r_u / (2 * d_u)

                for v in neighbors:
                    r[v] = r.get(v, 0.0) + mass_per_neighbor

                    d_v = max(deg_store.get(v, 1), 1)
                    if r[v] >= epsilon * d_v and v not in in_queue:
                        queue.append(v)
                        in_queue.add(v)
            else:
                total_mass = (1 - alpha) * r_u / 2

                for v in neighbors:
                    edge_weight = graph.get_weight(u, v)
                    mass_to_v = total_mass * edge_weight / d_u

                    r[v] = r.get(v, 0.0) + mass_to_v

                    d_v = max(deg_store.get(v, 1), 1)
                    if r[v] >= epsilon * d_v and v not in in_queue:
                        queue.append(v)
                        in_queue.add(v)

    return PageRankResult(
        seed=seed,
        scores=p,
        residuals=r,
        epsilon=epsilon,
        alpha=alpha,
        iterations=iterations,
    )


def pagerank_nibble(
    graph: BacteriaGraph,
    seed: str,
    *,
    phi: float = 0.3,
    b: Optional[int] = None,
    alpha_override: Optional[float] = None,
    max_iterations: int = 100_000,
    weight: Optional[str] = None,
    random_seed: Optional[int] = None,
) -> Tuple[Optional[Set[str]], PageRankResult]:
    if not graph.has_node(seed):
        raise KeyError(f"Unknown seed node: {seed}")

    m = max(graph.number_of_edges(), 1)

    alpha = alpha_override if alpha_override is not None else compute_alpha_paper(phi, m)

    if b is None:
        B = max(1, int(math.ceil(math.log2(m))))
        b = max(1, B // 2)

    epsilon = compute_epsilon_paper(b, m)

    pr_result = approximate_pagerank(
        graph,
        seed,
        alpha=alpha,
        epsilon=epsilon,
        max_iterations=max_iterations,
        weight=weight,
        random_seed=random_seed,
    )

    deg_store = graph.get_deg(weight=weight)
    total_volume = sum(deg_store.values())
    sweep_sets = pr_result.get_sweep_sets(graph)

    best_cut = None
    best_conductance = float("inf")

    log_m_ceil = math.ceil(math.log2(max(m, 2)))
    volume_min = 2**(b - 1)
    volume_max = (2 * total_volume) // 3

    for _node_id, _ratio, sweep_set in sweep_sets:
        vol_S = sum(deg_store.get(v, 1) for v in sweep_set)

        if vol_S <= volume_min or vol_S >= volume_max:
            continue

        boundary_weight = 0
        for v in sweep_set:
            for neighbor in graph.get_neighbors(v):
                if neighbor not in sweep_set:
                    if weight is None:
                        boundary_weight += 1
                    else:
                        boundary_weight += graph.get_weight(v, neighbor)

        vol_complement = total_volume - vol_S
        min_vol = min(vol_S, vol_complement)
        if min_vol == 0:
            continue

        conductance = boundary_weight / min_vol

        if conductance >= phi:
            continue

        p_at_vol = sum(pr_result.scores.get(v, 0.0) for v in sweep_set)

        if p_at_vol > 1.0 / (48 * log_m_ceil):
            if conductance < best_conductance:
                best_conductance = conductance
                best_cut = sweep_set

    return best_cut, pr_result


def compute_conductance(
    graph: BacteriaGraph,
    node_set: Set[str],
    weight: Optional[str] = None,
) -> float:
    deg_store = graph.get_deg(weight=weight)

    vol_S = sum(deg_store.get(v, 1) for v in node_set)

    boundary = 0
    for v in node_set:
        for neighbor in graph.get_neighbors(v):
            if neighbor not in node_set:
                if weight is None:
                    boundary += 1
                else:
                    boundary += graph.get_weight(v, neighbor)

    total_volume = sum(deg_store.values())
    vol_complement = total_volume - vol_S

    min_vol = min(vol_S, vol_complement)
    if min_vol == 0:
        return float("inf")

    return boundary / min_vol
