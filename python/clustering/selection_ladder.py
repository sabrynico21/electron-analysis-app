"""Five-step cluster selection for the PageRank-Nibble sweep.

This module implements the *selection stage* of the PageRank-Nibble
pipeline.  The stage receives the approximate PageRank vector produced by
the push algorithm of Andersen, Chung and Lang and must answer a single
question: **which sweep prefix should be returned as the cluster?**

Pipeline
--------

1. **Push** - :func:`pagerank_nibble.approximate_pagerank` computes the
   vector ``p`` with teleport probability ``alpha`` and residual
   threshold ``epsilon``.
2. **Sweep** - nodes are ordered by ``p(v) / d(v)`` (descending; ties
   are broken lexicographically by node identifier, so the order is a
   deterministic function of ``p``) and the nested prefixes
   ``S_1, ..., S_K`` of that order are scanned **incrementally**: one
   pass maintains the volume, the cut and the PageRank mass of every
   prefix in ``O(sum_v d(v))`` time, instead of the ``O(K * sum_v d(v))``
   cost of recomputing each cut from scratch.
3. **Select** - the five-step ladder implemented here returns the first
   cluster it can justify and records in the metadata exactly how strict
   the returned cut is with respect to the three admission criteria.

Admission criteria
------------------

Evaluated for every prefix ``S`` given the quality bar ``phi``, the
scale ``b`` and a graph with ``m`` edges:

* **volume**  ``vol(S) in (2^(b-1), vol(V)/2]``
* **mass**    ``p(S) > tau`` with ``tau = 1 / (48 * ceil(log2 m))``
* **quality** ``score(S) < phi`` where the score is the plain
  conductance, or the signed score when signed clustering is enabled.

The volume ceiling is one half of the total volume (not two thirds):
conductance is invariant under complementation
(``phi(S) == phi(V \\ S)``), so every prefix beyond half of the volume
duplicates the value of a smaller set, and the cap also bounds the
plain conductance of any admissible prefix to ``[0, 1]``.

The five steps
--------------

=============  ==========================================================
step 1         **strict** - volume, mass and quality all hold; the
               candidate with the lowest score is returned, with ties
               resolved towards the smallest prefix (see
               ``epsilon_tie``).
step 2         **quality** - the mass certificate is dropped: ``tau`` is
               a device of the paper's proof and is nearly vacuous at
               this graph scale; a prefix with the right size and a
               conductance below ``phi`` is a valid cluster even if its
               mass is marginally small.
step 3         **multi-scale** - the volume floor is halved up to
               ``scale_relaxations`` times (``b-1``, ``b-2``).  This is
               equivalent to re-running the algorithm at a coarser
               scale: a finer ``epsilon`` (the one used for the original
               ``b``) is always admissible for a coarser scale, so the
               *same* PageRank vector and the *same* sweep scan are
               reused - no second push is executed.
step 4         **relaxed quality** - the bar is raised geometrically
               (``x1.25, x1.5, x1.75, x2.0``, capped at ``2 * phi``)
               inside the original volume window; the rung that first
               admits a prefix is recorded as ``phi_effective``.
step 5         **best effort** - every admissible prefix is graded by a
               normalised penalty and the minimum is returned with its
               breakdown (see :func:`penalty_terms`).
=============  ==========================================================

If no admissible prefix exists (for example because the graph has fewer
than ``min_cut_size`` taxa), the seed itself is returned with
``ladder_step = 0`` and ``strategy = "seed_singleton_fallback"``: an
honest "no cluster found" instead of a degenerate micro-cluster.

Reference
---------

Andersen, R., Chung, F., & Lang, K. (2006). Local Graph Partitioning
using PageRank Vectors.
https://www.cs.cmu.edu/~15859n/RelatedWork/local_partitioning_full.pdf
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from .bacteria_graph import BacteriaGraph
    from .pagerank_nibble import (
        PageRankResult,
        approximate_pagerank,
        compute_conductance,
        compute_epsilon_paper,
    )
except ImportError:  # executed when run as a flat module (scripts/tests)
    from bacteria_graph import BacteriaGraph
    from pagerank_nibble import (
        PageRankResult,
        approximate_pagerank,
        compute_conductance,
        compute_epsilon_paper,
    )


# ---------------------------------------------------------------------------
# Ladder parameters (documented in tests/docs/SELECTION_LADDER.md)
# ---------------------------------------------------------------------------
DEFAULT_MIN_CUT_SIZE = 5
DEFAULT_SCALE_RELAXATIONS = 2
DEFAULT_PHI_RELAX_FACTORS: Tuple[float, ...] = (1.25, 1.5, 1.75, 2.0)
DEFAULT_EPSILON_TIE = 0.10
DEFAULT_PENALTY_VOLUME_WEIGHT = 0.5
DEFAULT_PENALTY_MASS_WEIGHT = 0.10


@dataclass(frozen=True)
class LadderParameters:
    """Configuration of the five-step selection ladder.

    ``phi`` and ``b`` are the algorithm's quality bar and scale; the
    remaining fields are the ladder's own controls and are identical for
    every graph unless a caller overrides them explicitly.
    """

    phi: float
    b: int
    min_cut_size: int = DEFAULT_MIN_CUT_SIZE
    max_cut_size: Optional[int] = None
    scale_relaxations: int = DEFAULT_SCALE_RELAXATIONS
    phi_relax_factors: Tuple[float, ...] = DEFAULT_PHI_RELAX_FACTORS
    epsilon_tie: float = DEFAULT_EPSILON_TIE
    penalty_volume_weight: float = DEFAULT_PENALTY_VOLUME_WEIGHT
    penalty_mass_weight: float = DEFAULT_PENALTY_MASS_WEIGHT
    weight: Optional[str] = None

    @property
    def volume_min(self) -> float:
        """Strict volume floor ``2^(b-1)`` of the paper's window."""
        return float(2 ** (self.b - 1))

    def relaxed_volume_min(self, steps: int) -> float:
        """Volume floor after ``steps`` halvings (``steps >= 1``)."""
        return float(2 ** (self.b - 1 - steps))


# ---------------------------------------------------------------------------
# Incremental prefix scores
# ---------------------------------------------------------------------------
class ConductanceScorer:
    """Streaming conductance of nested prefixes (plain clustering).

    ``add(node)`` folds one node into the current prefix: the volume
    grows by ``d(v)`` and the cut changes by ``d(v) - 2 * w(v, S)``,
    which makes a full sweep cost ``O(sum_v d(v)``) instead of
    ``O(K * sum_v d(v))``.
    """

    is_signed = False

    def __init__(self, graph: BacteriaGraph, weight: Optional[str] = None) -> None:
        self._graph = graph
        self._weight = weight
        self._deg = graph.get_deg(weight=weight)
        self._inside: set[str] = set()
        self._volume = 0.0
        self._cut = 0.0
        self.total_volume = float(sum(self._deg.values()))

    def degree(self, node: str) -> float:
        return float(self._deg.get(node, 0))

    def add(self, node: str) -> None:
        d_v = self.degree(node)
        inside_weight = 0.0
        for neighbor in self._graph.get_neighbors(node):
            if neighbor in self._inside:
                if self._weight is None:
                    inside_weight += 1.0
                else:
                    inside_weight += self._graph.get_weight(node, neighbor)
        self._volume += d_v
        self._cut += d_v - 2.0 * inside_weight
        self._inside.add(node)

    @property
    def volume(self) -> float:
        return self._volume

    @property
    def cut(self) -> float:
        return self._cut

    @property
    def score(self) -> float:
        denominator = min(self._volume, self.total_volume - self._volume)
        if denominator <= 0.0:
            return float("inf")
        return self._cut / denominator


class SignedScorer:
    """Streaming signed score of nested prefixes (signed clustering).

    The score is ``(cut_pos + lambda * neg_inside) / min(vol_pos(S),
    vol_pos(V \\ S))`` over the positive-edge volume; see
    ``tests/docs/ALGORITHM.md`` for the rationale.  ``composition()``
    reports the sign breakdown of a finished prefix and is the single
    implementation behind ``clustering_advanced.signed_metrics``.
    """

    is_signed = True

    def __init__(self, index: dict, sign_lambda: float = 1.0) -> None:
        self.index = index
        self.sign_lambda = float(sign_lambda)
        self._pos_adj: dict = index["pos_adj"]
        self._neg_adj: dict = index["neg_adj"]
        self._pos_deg: dict = index["pos_deg"]
        self.total_volume = float(index["total_pos_volume"])
        self._inside: set[str] = set()
        self._volume = 0.0
        self._cut_pos = 0.0
        self._neg_inside = 0.0

    def degree(self, node: str) -> float:
        return float(self._pos_deg.get(node, 0))

    def add(self, node: str) -> None:
        pos_deg_v = self.degree(node)
        inside_pos = 0.0
        for neighbor, weight in self._pos_adj.get(node, ()):
            if neighbor in self._inside:
                inside_pos += weight
        inside_neg = 0.0
        for neighbor, weight in self._neg_adj.get(node, ()):
            if neighbor in self._inside:
                inside_neg += weight
        self._volume += pos_deg_v
        self._cut_pos += pos_deg_v - 2.0 * inside_pos
        self._neg_inside += inside_neg
        self._inside.add(node)

    @property
    def volume(self) -> float:
        return self._volume

    @property
    def cut(self) -> float:
        """Positive cut of the current prefix (kept for symmetry)."""
        return self._cut_pos

    @property
    def score(self) -> float:
        denominator = min(self._volume, self.total_volume - self._volume)
        if denominator <= 0.0:
            return float("inf")
        violations = self._cut_pos + self.sign_lambda * self._neg_inside
        return violations / denominator

    def composition(self, nodes: set[str]) -> dict:
        """Sign breakdown of an arbitrary node set (metadata support)."""
        pos_adj = self._pos_adj
        neg_adj = self._neg_adj
        pos_deg = self._pos_deg

        vol_pos = sum(pos_deg.get(node, 0.0) for node in nodes)
        cut_pos = 0.0
        positive_inside = 0.0
        for node in nodes:
            for neighbor, weight in pos_adj.get(node, ()):
                if neighbor in nodes:
                    positive_inside += weight
                else:
                    cut_pos += weight
        positive_inside /= 2.0

        neg_inside = 0.0
        neg_cut = 0.0
        for node in nodes:
            for neighbor, weight in neg_adj.get(node, ()):
                if neighbor in nodes:
                    if node < neighbor:
                        neg_inside += weight
                else:
                    neg_cut += weight

        denominator = min(vol_pos, self.total_volume - vol_pos)
        violations = cut_pos + self.sign_lambda * neg_inside
        score = violations / denominator if denominator > 0 else float("inf")
        return {
            "score": score,
            "cut_pos": cut_pos,
            "neg_inside": neg_inside,
            "neg_cut": neg_cut,
            "positive_inside": positive_inside,
            "violations": violations,
            "vol_pos": vol_pos,
            "total_pos_volume": self.total_volume,
            "sign_lambda": self.sign_lambda,
        }


# ---------------------------------------------------------------------------
# Sweep scan
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PrefixStats:
    """Pre-computed statistics of one sweep prefix."""

    size: int          # number of nodes |S| (the prefix ends at order[|S|-1])
    volume: float      # vol(S) in the scorer's degree measure
    cut: float         # boundary weight of S
    score: float       # quality metric (conductance or signed score)
    mass: float        # PageRank mass p(S)


@dataclass(frozen=True)
class SweepScan:
    """Result of the single incremental sweep pass."""

    order: List[str]
    prefixes: List[PrefixStats]
    total_volume: float
    volume_max: float
    mass_threshold: float
    seed_index: Optional[int]
    min_admissible_size: int

    @property
    def evaluated_sweeps(self) -> int:
        """Number of admissible prefixes (seed present, size floor)."""
        return sum(
            1
            for prefix in self.prefixes
            if prefix.size >= self.min_admissible_size
        )

    def nodes_of(self, prefix: PrefixStats) -> set[str]:
        return set(self.order[: prefix.size])


def mass_threshold_for(graph: BacteriaGraph) -> float:
    """``tau = 1 / (48 * ceil(log2 m))`` of the paper's mass certificate."""
    m = max(graph.number_of_edges(), 1)
    return 1.0 / (48.0 * math.ceil(math.log2(max(m, 2))))


def evaluate_sweep(
    graph: BacteriaGraph,
    pr_result: PageRankResult,
    scorer,
    *,
    seed: str,
    ladder: LadderParameters,
) -> SweepScan:
    """Scan every sweep prefix once and cache volume, cut, score, mass.

    The sweep stops as soon as a prefix exceeds ``vol(V)/2``: all later
    prefixes exceed it too (the prefixes are nested), so none of them
    can be admissible in any step of the ladder.
    """
    scores = pr_result.scores
    total_volume = scorer.total_volume
    volume_max = float(int(total_volume // 2))

    order: List[str] = sorted(
        scores,
        key=lambda node: (-scores[node] / max(scorer.degree(node), 1.0), node),
    )
    seed_index = order.index(seed) if seed in scores else None
    if seed_index is None:
        min_admissible = max(ladder.min_cut_size, 1)
    else:
        # Admissible prefixes must contain the seed and respect the size
        # floor; ``max_cut_size`` only adds an upper bound (applied in
        # ``select_prefix``), so the admissible range can legitimately be
        # empty and the ladder falls through to the singleton.
        min_admissible = max(ladder.min_cut_size, seed_index + 1)

    mass_threshold = mass_threshold_for(graph)
    prefixes: List[PrefixStats] = []
    mass = 0.0
    for index, node in enumerate(order):
        mass += scores.get(node, 0.0)
        scorer.add(node)
        volume = scorer.volume
        if volume > volume_max:
            # Nested prefixes only grow: the rest of the sweep is out of
            # the volume window for every step of the ladder.
            break
        prefixes.append(
            PrefixStats(
                size=index + 1,
                volume=volume,
                cut=scorer.cut,
                score=scorer.score,
                mass=mass,
            )
        )

    return SweepScan(
        order=order,
        prefixes=prefixes,
        total_volume=total_volume,
        volume_max=volume_max,
        mass_threshold=mass_threshold,
        seed_index=seed_index,
        min_admissible_size=min_admissible,
    )


# ---------------------------------------------------------------------------
# Step 5: graded penalty
# ---------------------------------------------------------------------------
def penalty_terms(
    prefix: PrefixStats, ladder: LadderParameters, scan: SweepScan
) -> dict:
    """Normalised penalty of a prefix (all terms dimensionless).

    * ``quality``   - ``score / phi``: 1.0 sits exactly on the bar.
    * ``volume``    - octaves (``log2``) below the strict floor; 0 when
      the prefix is inside the window.  One octave (a factor two in
      volume) is priced at ``penalty_volume_weight`` units of quality
      ratio.
    * ``mass``      - ``-weight * mass``: a mild bonus that only breaks
      near-ties, because the mass certificate is the weakest criterion.
    """
    quality = prefix.score / ladder.phi if ladder.phi > 0 else float("inf")
    floor = ladder.volume_min
    octaves = 0.0
    if prefix.volume <= 0.0:
        octaves = float("inf")
    elif prefix.volume < floor:
        octaves = math.log2(floor / prefix.volume)
    mass_bonus = -ladder.penalty_mass_weight * min(prefix.mass, 1.0)
    total = (
        quality
        + ladder.penalty_volume_weight * octaves
        + mass_bonus
    )
    return {
        "quality": quality,
        "volume_octaves": octaves,
        "mass_bonus": mass_bonus,
        "penalty": total,
    }


# ---------------------------------------------------------------------------
# The ladder
# ---------------------------------------------------------------------------
@dataclass
class SelectionOutcome:
    """Which step produced the answer, and with which relaxations."""

    step: int                       # 0 = none (singleton), 1..5
    rule: str                       # human-readable rule name
    prefix: Optional[PrefixStats]
    volume_min_used: float          # floor applied at this step
    phi_effective: float            # quality bar applied at this step
    phi_relax_factor: Optional[float] = None   # step 4 rung
    scale_relaxation: Optional[int] = None     # step 3 halvings
    penalty: Optional[dict] = None             # step 5 breakdown
    evaluated: int = 0


def _pick_smallest_within(
    candidates: Sequence[PrefixStats], epsilon_tie: float
) -> Optional[PrefixStats]:
    """Lowest-score candidate, ties resolved towards the smallest prefix.

    The winner is the first candidate (in sweep order) whose score is
    within a relative ``epsilon_tie`` of the best score observed.  With
    a flat conductance curve this returns the smallest of the
    equivalent prefixes; a larger prefix is chosen only when it improves
    the score by more than ``epsilon_tie`` (10% by default), which is
    what prevents the "argument minimum drifts to a huge set" failure
    mode of a raw ``argmin``.
    """
    if not candidates:
        return None
    best = min(candidate.score for candidate in candidates)
    limit = best * (1.0 + epsilon_tie)
    for candidate in candidates:  # candidates are in increasing sweep order
        if candidate.score <= limit:
            return candidate
    return candidates[0]  # numerical safety: best itself satisfies best <= limit


def select_prefix(scan: SweepScan, ladder: LadderParameters) -> SelectionOutcome:
    """Run the five-step ladder over a finished sweep scan."""
    phi = ladder.phi
    tau = scan.mass_threshold
    floor_strict = ladder.volume_min

    admissible = [
        prefix
        for prefix in scan.prefixes
        if prefix.size >= scan.min_admissible_size
        and (
            ladder.max_cut_size is None
            or prefix.size <= ladder.max_cut_size
        )
    ]
    evaluated = len(admissible)

    def outcome(
        step: int,
        rule: str,
        prefix: PrefixStats,
        *,
        volume_min_used: float,
        phi_effective: float,
        **extra,
    ) -> SelectionOutcome:
        return SelectionOutcome(
            step=step,
            rule=rule,
            prefix=prefix,
            volume_min_used=volume_min_used,
            phi_effective=phi_effective,
            evaluated=evaluated,
            **extra,
        )

    # Step 1 - strict: volume, mass and quality.
    strict = [
        prefix
        for prefix in admissible
        if prefix.volume > floor_strict
        and prefix.mass > tau
        and prefix.score < phi
    ]
    chosen = _pick_smallest_within(strict, ladder.epsilon_tie)
    if chosen is not None:
        return outcome(1, "strict", chosen,
                       volume_min_used=floor_strict, phi_effective=phi)

    # Step 2 - quality without the (near-vacuous) mass certificate.
    quality = [
        prefix
        for prefix in admissible
        if prefix.volume > floor_strict and prefix.score < phi
    ]
    chosen = _pick_smallest_within(quality, ladder.epsilon_tie)
    if chosen is not None:
        return outcome(2, "quality_without_mass", chosen,
                       volume_min_used=floor_strict, phi_effective=phi)

    # Step 3 - multi-scale: halve the volume floor, reuse the same vector.
    for steps in range(1, max(ladder.scale_relaxations, 0) + 1):
        floor = ladder.relaxed_volume_min(steps)
        relaxed = [
            prefix
            for prefix in admissible
            if prefix.volume > floor and prefix.score < phi
        ]
        chosen = _pick_smallest_within(relaxed, ladder.epsilon_tie)
        if chosen is not None:
            return outcome(3, "multi_scale", chosen,
                           volume_min_used=floor, phi_effective=phi,
                           scale_relaxation=steps)

    # Step 4 - relax the quality bar geometrically inside the window.
    for factor in ladder.phi_relax_factors:
        rung = phi * factor
        relaxed = [
            prefix
            for prefix in admissible
            if prefix.volume > floor_strict and prefix.score < rung
        ]
        chosen = _pick_smallest_within(relaxed, ladder.epsilon_tie)
        if chosen is not None:
            return outcome(4, "relaxed_phi", chosen,
                           volume_min_used=floor_strict, phi_effective=rung,
                           phi_relax_factor=factor)

    # Step 5 - best effort: minimum normalised penalty.
    best_prefix: Optional[PrefixStats] = None
    best_penalty: Optional[dict] = None
    for prefix in admissible:
        terms = penalty_terms(prefix, ladder, scan)
        if not math.isfinite(terms["penalty"]):
            continue
        if best_penalty is None or terms["penalty"] < best_penalty["penalty"]:
            best_prefix = prefix
            best_penalty = terms
    if best_prefix is not None:
        return outcome(5, "penalty_best_effort", best_prefix,
                       volume_min_used=floor_strict, phi_effective=phi,
                       penalty=best_penalty)

    # Nothing admissible: honest "no cluster found".
    return SelectionOutcome(
        step=0,
        rule="no_admissible_cluster",
        prefix=None,
        volume_min_used=floor_strict,
        phi_effective=phi,
        evaluated=evaluated,
    )


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
def _strict_conditions(
    prefix: Optional[PrefixStats], ladder: LadderParameters, scan: SweepScan
) -> Tuple[bool, bool, bool]:
    if prefix is None:
        return False, False, False
    volume_ok = ladder.volume_min < prefix.volume <= scan.volume_max
    prob_ok = prefix.mass > scan.mass_threshold
    quality_ok = prefix.score < ladder.phi
    return volume_ok, prob_ok, quality_ok


def build_metadata(
    outcome: SelectionOutcome,
    scan: SweepScan,
    ladder: LadderParameters,
    *,
    graph: BacteriaGraph,
    seed: str,
    scores: Dict[str, float],
    scorer,
    alpha: float,
    epsilon: float,
    weight: Optional[str],
    signed: bool,
    sign_lambda: float,
) -> dict:
    """Assemble the public metadata of the returned cluster.

    The keys of the pre-ladder implementation (``strategy``,
    ``conditions_met_count``, ``failed_conditions``, ``constraints``,
    ``ranked_nodes``, ...) are preserved so that the reports and the UI
    keep working; the ladder adds ``ladder_step``, ``ladder_rule``,
    ``phi_effective``, ``volume_min_used`` and, on step 5, the
    ``penalty_breakdown``.
    """
    metadata: dict = {
        "b": ladder.b,
        "epsilon": epsilon,
        "alpha": alpha,
        "signed": signed,
        "phi_target": ladder.phi,
        "phi_effective": outcome.phi_effective,
        "ladder_step": outcome.step,
        "ladder_rule": outcome.rule,
        "volume_min": ladder.volume_min,
        "volume_min_used": outcome.volume_min_used,
        "volume_max": scan.volume_max,
        "mass_threshold": scan.mass_threshold,
        "min_cut_size": ladder.min_cut_size,
        "evaluated_sweeps": outcome.evaluated,
        "selection": outcome.rule,
    }
    if outcome.phi_relax_factor is not None:
        metadata["phi_relax_factor"] = outcome.phi_relax_factor
    if outcome.scale_relaxation is not None:
        metadata["scale_relaxation"] = outcome.scale_relaxation
    if outcome.penalty is not None:
        metadata["penalty_breakdown"] = dict(outcome.penalty)

    if outcome.prefix is None:
        # Honest singleton: no admissible prefix existed.
        metadata.update(
            {
                "strategy": "seed_singleton_fallback",
                "conductance": float("inf"),
                "score": float("inf"),
                "size": 1,
                "volume": scorer.degree(seed),
                "p_at_vol": scores.get(seed, 0.0),
                "conditions_met_count": 0,
                "failed_conditions": ["volume_ok", "prob_ok", "phi_ok"],
                "satisfied": False,
                "constraints": {
                    "volume_ok": False,
                    "prob_ok": False,
                    "quality_ok": False,
                },
                "constraints_satisfied": False,
                "ranked_nodes": [seed],
                "ranked_scores": [{"id": seed, "score": float(scores.get(seed, 0.0))}],
            }
        )
        if signed:
            metadata["sign_lambda"] = sign_lambda
        return metadata

    prefix = outcome.prefix
    nodes = scan.nodes_of(prefix)
    volume_ok, prob_ok, quality_ok = _strict_conditions(prefix, ladder, scan)
    failed = [
        name
        for name, ok in (
            ("volume_ok", volume_ok),
            ("prob_ok", prob_ok),
            ("phi_ok", quality_ok),
        )
        if not ok
    ]
    conditions_met = 3 - len(failed)
    satisfied = conditions_met == 3

    metadata.update(
        {
            # ``target_conductance`` marks a strict (all-three) success,
            # which is what the stability reports count; every relaxed
            # step reports ``target_conductance_fallback``.
            "strategy": (
                "target_conductance" if outcome.step == 1
                else "target_conductance_fallback"
            ),
            "conductance": compute_conductance(graph, nodes, weight=weight),
            "score": prefix.score,
            "size": prefix.size,
            "volume": prefix.volume,
            "p_at_vol": prefix.mass,
            "conditions_met_count": conditions_met,
            "failed_conditions": failed,
            "satisfied": satisfied,
            "sweep_index": prefix.size - 1,
            "constraints": {
                "volume_ok": volume_ok,
                "prob_ok": prob_ok,
                "quality_ok": quality_ok,
            },
            "constraints_satisfied": satisfied,
            "ranked_nodes": sorted(
                nodes, key=lambda node: scores.get(node, 0.0), reverse=True
            ),
            "ranked_scores": [
                {"id": node, "score": float(scores.get(node, 0.0))}
                for node in sorted(
                    nodes, key=lambda node: scores.get(node, 0.0), reverse=True
                )
            ],
        }
    )
    if signed:
        metadata["sign_lambda"] = sign_lambda
        composition = scorer.composition(nodes)
        metadata.update(
            {
                key: value
                for key, value in composition.items()
                if key not in {"score", "sign_lambda"}
            }
        )
        metadata["score"] = composition["score"]
    return metadata


# ---------------------------------------------------------------------------
# End-to-end entry point
# ---------------------------------------------------------------------------
def cluster_page_rank_nibble(
    graph: BacteriaGraph,
    seed: str,
    *,
    phi: float,
    b: int,
    alpha: float,
    weight: Optional[str] = None,
    walk_graph: Optional[BacteriaGraph] = None,
    scorer=None,
    signed: bool = False,
    sign_lambda: float = 1.0,
    min_cut_size: int = DEFAULT_MIN_CUT_SIZE,
    max_cut_size: Optional[int] = None,
    max_iterations: int = 100_000,
    random_seed: Optional[int] = None,
) -> Tuple[set, dict]:
    """Push + single sweep scan + five-step ladder.

    ``walk_graph`` is the graph the random walk runs on (the positive
    sub-graph when ``signed``); ``graph`` is always the full graph and
    provides the reported conductance.  ``scorer`` may be supplied by
    the caller (a :class:`SignedScorer` for signed clustering);
    otherwise a plain :class:`ConductanceScorer` over ``graph`` and
    ``weight`` is built.
    """
    if not graph.has_node(seed):
        raise KeyError(f"Unknown seed node: {seed}")

    walk = walk_graph if walk_graph is not None else graph
    if scorer is None:
        if signed:
            raise ValueError(
                "signed clustering requires a SignedScorer built from the "
                "signed edge index"
            )
        scorer = ConductanceScorer(graph, weight)

    ladder = LadderParameters(
        phi=phi,
        b=b,
        min_cut_size=min_cut_size,
        max_cut_size=max_cut_size,
        weight=weight,
    )

    m = max(graph.number_of_edges(), 1)
    epsilon = compute_epsilon_paper(b, m)
    pr_result = approximate_pagerank(
        walk,
        seed,
        alpha=alpha,
        epsilon=epsilon,
        max_iterations=max_iterations,
        weight=weight,
        random_seed=random_seed,
    )

    scan = evaluate_sweep(graph, pr_result, scorer, seed=seed, ladder=ladder)
    outcome = select_prefix(scan, ladder)
    metadata = build_metadata(
        outcome,
        scan,
        ladder,
        graph=graph,
        seed=seed,
        scores=pr_result.scores,
        scorer=scorer,
        alpha=alpha,
        epsilon=epsilon,
        weight=weight,
        signed=signed,
        sign_lambda=sign_lambda,
    )
    # Push diagnostics: with a very small teleport probability the vector
    # can hit the iteration budget before the residual mass is absorbed;
    # reporting it keeps truncated runs visible in the metadata.
    metadata["push_iterations"] = pr_result.iterations
    metadata["push_converged"] = pr_result.iterations < max_iterations

    if outcome.prefix is None:
        return {seed}, metadata
    return scan.nodes_of(outcome.prefix), metadata
