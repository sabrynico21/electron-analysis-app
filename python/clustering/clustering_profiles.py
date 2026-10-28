"""Parameter resolution for the PageRank-Nibble clustering.

A single parameter set is supported: the ``app`` profile, i.e. the
reference defaults validated by the published experiments
(``phi = 0.4``, ``b = compute_default_b(n, m)`` (10 on both bundled
networks), ``alpha = 0.05`` independent of ``phi``, ``|correlation|``
weights).  The legacy app values (``phi = 0.2``, ``b = 4``,
``alpha ~ 2e-5`` from the paper's bound) were retired together with the
former ``recommended`` profile; both remain reproducible through explicit
overrides.

Cluster selection is exclusively the five-step ladder documented in
``tests/docs/SELECTION_LADDER.md`` (the ``first`` / ``min_conductance``
selection rules were removed with the profile).  The experiment harness
passes its base values as explicit overrides
(``tests/stability_experiment.py``), so the published experiments stay
reproducible.

Every value can still be overridden through
``resolve_parameters(..., overrides={...})``.
"""

from __future__ import annotations

import math
from typing import Optional

from bacteria_graph import BacteriaGraph

# Reference defaults: the parameter set of the published tuning/comparison
# experiments (phi = 0.4, b = 10, alpha = 0.05 on both bundled networks).
# The legacy app values (phi = 0.2, b = 4, alpha ~ 2e-5 from the paper's
# bound) were retired: on network 101 they never pass strictly and the tiny
# teleport consumes the whole iteration budget.  Callers can still request
# them through explicit overrides.
DEFAULT_TARGET_CONDUCTANCE = 0.4
DEFAULT_TELEPORT_ALPHA = 0.05
DEFAULT_TARGET_CLUSTER_FRACTION = 0.10
DEFAULT_WEIGHTED = True


def compute_default_b(node_count: int, edge_count: int) -> int:
    """Scale of the default cluster, derived from the graph dimensions.

    Targets a cluster of ~10% of the taxa (clamped to ``[10, n/3]``), turns
    that size into a target volume with the average degree, and returns
    ``b = round(log2(target_volume))`` capped so the paper's volume window
    stays inside the graph.  **10 on both bundled networks**; keep in sync
    with ``computeDefaultB`` in ``AnalysisWorkspace.jsx``.
    """

    n = max(1, int(node_count or 1))
    m = max(1, int(edge_count or 1))

    average_degree = (2 * m) / n
    target_size = max(10, min(max(n // 3, 10), round(DEFAULT_TARGET_CLUSTER_FRACTION * n)))
    target_volume = max(2.0, average_degree * target_size)

    b_max = max(1, math.floor(math.log2(max(2 * m, 2))))
    return max(1, min(b_max, int(round(math.log2(target_volume)))))


def compute_default_teleport_alpha(
    target_conductance: float = DEFAULT_TARGET_CONDUCTANCE,
    edge_count: int = 1,
) -> float:
    """Default teleport probability: a constant, independent of ``phi``.

    The paper's ``phi^2 / (225 ln(100 sqrt(m)))`` is a worst-case bound for
    the theorem, not an operating point: it yields ~2e-5, which barely
    localises the walk and burns the iteration budget.  The reference
    experiments select **0.05**; the parameters are kept for signature
    compatibility with the UI port.
    """

    del target_conductance, edge_count
    return DEFAULT_TELEPORT_ALPHA


def app_parameters(
    graph: BacteriaGraph,
    target_conductance: float = DEFAULT_TARGET_CONDUCTANCE,
    weighted: bool = DEFAULT_WEIGHTED,
) -> dict:
    n = graph.number_of_nodes()
    m = graph.number_of_edges()
    return {
        "profile": "app",
        "target_conductance": float(target_conductance),
        "b": compute_default_b(n, m),
        "teleport_alpha": compute_default_teleport_alpha(target_conductance, m),
        "weighted": bool(weighted),
        "weight": "weight" if weighted else None,
        "node_count": n,
        "edge_count": m,
    }


def resolve_parameters(
    graph: BacteriaGraph,
    profile: str = "app",
    *,
    weighted: bool = DEFAULT_WEIGHTED,
    overrides: Optional[dict] = None,
) -> dict:
    """Build the parameter dict for ``profile`` and apply explicit overrides."""

    if profile != "app":
        raise ValueError(
            f"Unknown parameter profile: {profile!r}. The 'recommended' "
            "profile was removed; pass the values as explicit overrides "
            "(target_conductance, b, teleport_alpha) instead."
        )
    params = app_parameters(graph, weighted=weighted)
    for key, value in (overrides or {}).items():
        if value is not None:
            params[key] = value
    return params
