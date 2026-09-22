"""
Adaptive PageRank-Nibble
================================================
This provides adaptive parameter selection that scales with graph size
and implements a fallback strategy to ALWAYS return a cut.
"""

import math
from typing import Dict, Set, Tuple, Optional
from dataclasses import dataclass

try:
    from .bacteria_graph import BacteriaGraph
except ImportError:
    from bacteria_graph import BacteriaGraph


@dataclass
class AdaptiveParams:
    """Auto-computed parameters based on graph size."""
    n: int
    m: int
    alpha: float
    b_min: int
    b_max: int
    b_default: int
    epsilon_at_default: float

    def __str__(self):
        return f"""Adaptive Parameters for n={self.n}, m={self.m}:
  - Alpha (alpha): {self.alpha:.6f}
  - Scale range: b in [{self.b_min}, {self.b_max}]
  - Default scale: b = {self.b_default} (targets volume ~{2**self.b_default})
      - Epsilon at b={self.b_default}: {self.epsilon_at_default:.8f}"""


def compute_adaptive_params(
    graph: BacteriaGraph,
    target_conductance: float = 0.2,
    conservative: bool = False,
) -> AdaptiveParams:
    n = graph.number_of_nodes()
    m = max(graph.number_of_edges(), 1)

    alpha = target_conductance**2 / (225 * math.log(100 * math.sqrt(m)))

    if conservative:
        alpha *= 3

    b_max = max(1, int(math.ceil(math.log2(m))))
    b_min = max(1, int(math.log2(math.sqrt(max(n, 1)))))

    avg_degree = (2 * m) / n if n > 0 else 1

    if avg_degree < 10:
        b_default = max(b_min, int(math.log2(max(n / 10, 1))))
    elif avg_degree > 50:
        b_default = max(b_min, int(math.log2(max(n / 3, 1))))
    else:
        b_default = max(b_min, int(math.log2(max(math.sqrt(n), 1))))

    b_default = min(b_default, b_max)

    log_m_ceil = math.ceil(math.log2(max(m, 2)))
    epsilon_at_default = 1.0 / (2**b_default * 48 * log_m_ceil)

    return AdaptiveParams(
        n=n,
        m=m,
        alpha=alpha,
        b_min=b_min,
        b_max=b_max,
        b_default=b_default,
        epsilon_at_default=epsilon_at_default,
    )


def pagerank_nibble_adaptive(
    graph: BacteriaGraph,
    seed: str,
    *,
    target_conductance: float = 0.2,
    b: Optional[int] = None,
    teleport_alpha: Optional[float] = None,
    min_cut_size: int = 5,
    max_cut_size: Optional[int] = None,
    conservative: bool = False,
    max_iterations: int = 100_000,
    verbose: bool = True,
    random_seed: Optional[int] = None,
    weight: Optional[str] = None,
) -> Tuple[Set[str], Dict]:
    del min_cut_size, max_cut_size

    if not graph.has_node(seed):
        raise KeyError(f"Unknown seed node: {seed}")

    deg_store = graph.get_deg(weight=weight)
    total_volume = sum(deg_store.values())

    params = compute_adaptive_params(graph, target_conductance, conservative)

    if verbose:
        print(params)
        print()

    b_values = [int(b)] if b is not None else [params.b_default]
    if verbose:
        print(f"Using single scale b={b_values[0]}")

    def _annotate_and_print(cut_set, metadata):
        try:
            from .pagerank_nibble import compute_conductance
        except ImportError:
            from pagerank_nibble import compute_conductance

        b_local = metadata.get("b", params.b_default)
        log_m_ceil = math.ceil(math.log2(max(max(graph.number_of_edges(), 1), 2)))
        volume_min = 2 ** (b_local - 1)
        volume_max = (2 * total_volume) // 3

        vol_S = sum(deg_store.get(v, 1) for v in cut_set)
        volume_ok = (vol_S > volume_min) and (vol_S < volume_max)

        phi_target = metadata.get("phi_target", target_conductance)
        prob_threshold = 1.0 / (48 * log_m_ceil)
        p_at_vol = metadata.get("p_at_vol", None)
        prob_ok = (p_at_vol is not None) and (p_at_vol > prob_threshold)

        conductance = compute_conductance(graph, cut_set, weight=weight)
        conductance_ok = conductance < phi_target
        quality_ok = conductance_ok and prob_ok

        all_ok = volume_ok and quality_ok

        metadata = dict(metadata)
        metadata["constraints"] = {
            "volume_ok": volume_ok,
            "quality_ok": quality_ok,
        }
        metadata["constraints_satisfied"] = all_ok
        return cut_set, metadata

    if verbose:
        print(f"Stage 1: Target phi={target_conductance:.3f}")

    fallback_candidate = None
    fallback_rank = None
    last_pr_scores = {}
    evaluated_sweeps = 0

    def _attach_ranked_nodes(cut_set, metadata, score_map):
        ranked_nodes = sorted(
            list(cut_set),
            key=lambda node_id: score_map.get(node_id, 0.0),
            reverse=True,
        )
        ranked_scores = [
            {
                "id": node_id,
                "score": float(score_map.get(node_id, 0.0)),
            }
            for node_id in ranked_nodes
        ]
        enriched = dict(metadata)
        enriched["ranked_nodes"] = ranked_nodes
        enriched["ranked_scores"] = ranked_scores
        return cut_set, enriched

    for b_local in b_values:
        try:
            from .pagerank_nibble import pagerank_nibble
        except ImportError:
            from pagerank_nibble import pagerank_nibble

        _, pr_result = pagerank_nibble(
            graph,
            seed,
            phi=target_conductance,
            b=b_local,
            alpha_override=teleport_alpha,
            max_iterations=max_iterations,
            weight=weight,
            random_seed=random_seed,
        )

        sweep_sets = pr_result.get_sweep_sets(graph)
        last_pr_scores = pr_result.scores
        log_m_ceil = math.ceil(math.log2(max(max(graph.number_of_edges(), 1), 2)))
        volume_min = 2**(b_local - 1)
        volume_max = (2 * total_volume) // 3
        prob_threshold = 1.0 / (48 * log_m_ceil)

        for sweep_idx, (_, _, sweep_set) in enumerate(sweep_sets):
            if not sweep_set:
                continue

            evaluated_sweeps += 1

            vol_S = sum(deg_store.get(v, 1) for v in sweep_set)
            volume_ok = (vol_S > volume_min) and (vol_S < volume_max)

            boundary_edges = 0.0
            for v in sweep_set:
                for neighbor in graph.get_neighbors(v):
                    if neighbor not in sweep_set:
                        if weight is None:
                            boundary_edges += 1.0
                        else:
                            boundary_edges += graph.get_weight(v, neighbor)
            vol_complement = total_volume - vol_S
            min_vol = min(vol_S, vol_complement)
            if min_vol == 0:
                continue
            conductance = boundary_edges / min_vol

            p_at_vol = sum(pr_result.scores.get(v, 0.0) for v in sweep_set)
            prob_ok = p_at_vol > prob_threshold
            phi_ok = conductance < target_conductance
            conditions_met_count = int(volume_ok) + int(prob_ok) + int(phi_ok)

            rank = (
                conditions_met_count,
                -conductance,
                p_at_vol,
                -len(sweep_set),
                -sweep_idx,
            )

            if fallback_rank is None or rank > fallback_rank:
                failed_conditions = []
                if not volume_ok:
                    failed_conditions.append("volume_ok")
                if not prob_ok:
                    failed_conditions.append("prob_ok")
                if not phi_ok:
                    failed_conditions.append("phi_ok")

                fallback_rank = rank
                fallback_candidate = (
                    sweep_set,
                    {
                        "strategy": "target_conductance_fallback",
                        "selection": "max_condition_coverage_fallback",
                        "b": b_local,
                        "conductance": conductance,
                        "size": len(sweep_set),
                        "phi_target": target_conductance,
                        "p_at_vol": p_at_vol,
                        "conditions_met_count": conditions_met_count,
                        "failed_conditions": failed_conditions,
                        "sweep_index": sweep_idx,
                        "evaluated_sweeps": evaluated_sweeps,
                    },
                )

            if volume_ok and prob_ok and phi_ok:
                cut_set, strict_meta = _annotate_and_print(
                    sweep_set,
                    {
                        "strategy": "target_conductance",
                        "selection": "strict_first_match",
                        "b": b_local,
                        "conductance": conductance,
                        "size": len(sweep_set),
                        "phi_target": target_conductance,
                        "p_at_vol": p_at_vol,
                        "conditions_met_count": 3,
                        "failed_conditions": [],
                        "sweep_index": sweep_idx,
                    },
                )
                return _attach_ranked_nodes(cut_set, strict_meta, pr_result.scores)

    if fallback_candidate is not None:
        candidate_set, candidate_meta = fallback_candidate
        candidate_meta["satisfied"] = False
        candidate_meta["evaluated_sweeps"] = evaluated_sweeps
        cut_set, fallback_meta = _annotate_and_print(candidate_set, candidate_meta)
        return _attach_ranked_nodes(cut_set, fallback_meta, last_pr_scores)

    seed_p = last_pr_scores.get(seed, 0.0)
    cut_set, seed_meta = _annotate_and_print(
        {seed},
        {
            "strategy": "seed_singleton_fallback",
            "selection": "no_sweep_candidate",
            "b": params.b_default,
            "conductance": float("inf"),
            "size": 1,
            "phi_target": target_conductance,
            "p_at_vol": seed_p,
            "conditions_met_count": 0,
            "failed_conditions": ["volume_ok", "prob_ok", "phi_ok"],
            "satisfied": False,
            "evaluated_sweeps": evaluated_sweeps,
        },
    )
    return _attach_ranked_nodes(cut_set, seed_meta, last_pr_scores)
