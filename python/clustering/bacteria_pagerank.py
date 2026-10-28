"""
Adaptive PageRank-Nibble
===============================================
This provides adaptive parameter selection that scales with graph size
and always returns a cut through the five-step selection ladder of
``selection_ladder.py`` (strict -> quality -> multi-scale -> relaxed
phi -> best-effort penalty -> honest singleton).

The three admission criteria (volume window, PageRank mass, quality
bar) and the ladder's parameters are documented in
``tests/docs/SELECTION_LADDER.md``.
"""

import math
from typing import Dict, Set, Tuple, Optional
from dataclasses import dataclass

try:
    from .bacteria_graph import BacteriaGraph
    from .clustering_profiles import (
        DEFAULT_TELEPORT_ALPHA,
        DEFAULT_TARGET_CONDUCTANCE,
        compute_default_b,
    )
    from .selection_ladder import (
        DEFAULT_MIN_CUT_SIZE,
        cluster_page_rank_nibble,
    )
except ImportError:
    from bacteria_graph import BacteriaGraph
    from clustering_profiles import (
        DEFAULT_TELEPORT_ALPHA,
        DEFAULT_TARGET_CONDUCTANCE,
        compute_default_b,
    )
    from selection_ladder import (
        DEFAULT_MIN_CUT_SIZE,
        cluster_page_rank_nibble,
    )


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
    target_conductance: float = DEFAULT_TARGET_CONDUCTANCE,
    conservative: bool = False,
) -> AdaptiveParams:
    """Default parameters of the reference set, with graph-size bookkeeping.

    ``alpha`` is the experimentally validated constant (0.05, x3 when
    ``conservative``) and ``b_default`` the dimension-based scale of
    ``clustering_profiles.compute_default_b`` (10 on both bundled networks);
    ``b_min``/``b_max`` are only informational bounds.
    """
    n = graph.number_of_nodes()
    m = max(graph.number_of_edges(), 1)

    alpha = DEFAULT_TELEPORT_ALPHA

    if conservative:
        alpha *= 3

    b_max = max(1, int(math.ceil(math.log2(m))))
    b_min = max(1, int(math.log2(math.sqrt(max(n, 1)))))
    b_default = compute_default_b(n, m)

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
    target_conductance: float = DEFAULT_TARGET_CONDUCTANCE,
    b: Optional[int] = None,
    teleport_alpha: Optional[float] = None,
    min_cut_size: int = DEFAULT_MIN_CUT_SIZE,
    max_cut_size: Optional[int] = None,
    conservative: bool = False,
    max_iterations: int = 100_000,
    verbose: bool = True,
    random_seed: Optional[int] = None,
    weight: Optional[str] = None,
) -> Tuple[Set[str], Dict]:
    """Cluster around ``seed`` with a single push and the 5-step ladder.

    The approximate PageRank vector is computed once (teleport
    probability ``teleport_alpha``, falling back to the reference default
    0.05 when ``None``), the sweep is scanned once incrementally, and the
    five-step ladder selects the cluster.  ``min_cut_size`` is a universal
    admissibility floor: no step can return a cluster with fewer taxa than
    that; ``None``-like disabling is achieved with ``min_cut_size=1``.
    """
    if not graph.has_node(seed):
        raise KeyError(f"Unknown seed node: {seed}")

    params = compute_adaptive_params(graph, target_conductance, conservative)

    b_local = int(b) if b is not None else params.b_default
    alpha = teleport_alpha if teleport_alpha is not None else params.alpha

    if verbose:
        print(params)
        print()
        print(
            f"Stage 1: Target phi={target_conductance:.3f}, "
            f"b={b_local}, alpha={alpha:g}"
        )

    cluster_nodes, metadata = cluster_page_rank_nibble(
        graph,
        seed,
        phi=target_conductance,
        b=b_local,
        alpha=alpha,
        weight=weight,
        min_cut_size=min_cut_size,
        max_cut_size=max_cut_size,
        max_iterations=max_iterations,
        random_seed=random_seed,
    )

    if verbose:
        conductance = metadata.get("conductance")
        conductance_text = (
            f"{conductance:.4f}" if math.isfinite(conductance) else "inf"
        )
        print(
            f"Step {metadata['ladder_step']} ({metadata['ladder_rule']}): "
            f"size={metadata['size']}, conductance={conductance_text}, "
            f"conditions={metadata['conditions_met_count']}/3, "
            f"failed={metadata['failed_conditions']}"
        )

    return cluster_nodes, metadata
