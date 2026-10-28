"""Core stability experiments shared by the tuning and comparison scripts.

Each configuration is evaluated with:

1. Repeatability - N runs from the same seed taxon, varying only the
   residual-queue order; reports intersection/union, inclusion frequency,
   pairwise Jaccard, chance-corrected Jaccard, conductance and signed score.
2. Nested consistency - take the deterministic cluster from the seed, re-run
   the algorithm from **every** taxon of that cluster and compare the
   resulting clusters.  There is no run cap here: the number of runs is the
   size of the base cluster.  (``max_nodes`` exists only to cap the cost of
   smoke tests; the published experiments never set it.)

Both ``weighted=True`` (``|correlation|`` edges) and ``weighted=False``
(unweighted edges, signed score counts every edge as 1) are supported.

The report generators are ``tuning_experiment.py`` and
``comparison_experiment.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import statistics
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

TEST_DIR = Path(__file__).resolve().parent
CLUSTERING_DIR = TEST_DIR.parent
if str(TEST_DIR) not in sys.path:
    sys.path.insert(0, str(TEST_DIR))
if str(CLUSTERING_DIR) not in sys.path:
    sys.path.insert(0, str(CLUSTERING_DIR))

from bacteria_pagerank import pagerank_nibble_adaptive  # noqa: E402
from clustering_advanced import compute_signed_metrics, pagerank_cluster  # noqa: E402
from clustering_profiles import resolve_parameters  # noqa: E402
from network_loader import NETWORKS, load_network, resolve_network  # noqa: E402
from pagerank_nibble import compute_conductance, compute_epsilon_paper  # noqa: E402
from selection_ladder import DEFAULT_MIN_CUT_SIZE  # noqa: E402

DEFAULT_RUNS = 100
DEFAULT_RUN_SEED = 42

# Base values of the experiment configurations.  They used to come from the
# removed ``recommended`` profile; the harness now passes them as explicit
# overrides so the published experiments stay reproducible.
EXPERIMENT_DEFAULTS = {
    "teleport_alpha": 0.05,
    "b": 10,
    "target_conductance": 0.4,
}

# Every configuration runs the same five-step selection ladder
# (``selection_ladder.py``); only the parameter values and the sign
# handling differ.
CONFIGS = {
    "baseline": {
        "label": "Baseline (|r|, app parameters)",
        "params": None,
        "signed": False,
    },
    "new_parameters": {
        "label": "Default parameters (|r|)",
        "params": EXPERIMENT_DEFAULTS,
        "signed": False,
    },
    "signed_weights": {
        "label": "Signed weights (app parameters)",
        "params": None,
        "signed": True,
    },
    "both": {
        "label": "Both (signed + default parameters)",
        "params": EXPERIMENT_DEFAULTS,
        "signed": True,
    },
}

_META_KEYS = (
    "strategy",
    "selection",
    "ladder_step",
    "ladder_rule",
    "b",
    "epsilon",
    "alpha",
    "conductance",
    "score",
    "signed",
    "sign_lambda",
    "size",
    "volume",
    "phi_target",
    "phi_effective",
    "phi_relax_factor",
    "scale_relaxation",
    "volume_min",
    "volume_min_used",
    "volume_max",
    "mass_threshold",
    "min_cut_size",
    "p_at_vol",
    "conditions_met_count",
    "failed_conditions",
    "sweep_index",
    "satisfied",
    "evaluated_sweeps",
    "penalty_breakdown",
    "push_iterations",
    "push_converged",
    "constraints",
    "constraints_satisfied",
    "cut_pos",
    "neg_inside",
    "neg_cut",
    "positive_inside",
    "violations",
    "vol_pos",
)


# ---------------------------------------------------------------------------
# Clustering dispatcher
# ---------------------------------------------------------------------------
def run_cluster(
    graph,
    seed: str,
    config: str,
    random_seed: Optional[int] = None,
    parameter_overrides: Optional[dict] = None,
    weighted: bool = True,
):
    """Run one clustering with a configuration's parameters.

    The configuration contributes its base values (``params``), the call
    contributes the sweep overrides, and both are applied on top of the
    ``app`` profile - the only profile left.  Signed configurations are
    dispatched to ``pagerank_cluster``, unsigned ones to
    ``pagerank_nibble_adaptive``; both run the same five-step ladder.
    """
    settings = CONFIGS[config]
    overrides = dict(settings.get("params") or {})
    if parameter_overrides:
        overrides.update(parameter_overrides)
    params = resolve_parameters(
        graph, "app", weighted=weighted, overrides=overrides or None
    )

    if settings["signed"]:
        nodes, metadata = pagerank_cluster(
            graph,
            seed,
            target_conductance=params["target_conductance"],
            b=params["b"],
            teleport_alpha=params["teleport_alpha"],
            weighted=weighted,
            signed=True,
            sign_lambda=1.0,
            random_seed=random_seed,
        )
    else:
        nodes, metadata = pagerank_nibble_adaptive(
            graph,
            seed,
            target_conductance=params["target_conductance"],
            b=params["b"],
            teleport_alpha=params["teleport_alpha"],
            weight=params["weight"],
            random_seed=random_seed,
            verbose=False,
        )

    signed_stats = compute_signed_metrics(
        graph, set(nodes), 1.0, weighted=weighted
    )
    return set(nodes), dict(metadata), params, signed_stats


# ---------------------------------------------------------------------------
# Parallel execution helpers
# ---------------------------------------------------------------------------
_WORKER_GRAPH = None


def _worker_init(network_path: str) -> None:
    global _WORKER_GRAPH
    _WORKER_GRAPH, _ = load_network(network_path)


def _worker_run(
    task: tuple[str, Optional[int], str, Optional[dict], bool]
) -> dict:
    seed, random_seed, config, overrides, weighted = task
    nodes, metadata, params, signed_stats = run_cluster(
        _WORKER_GRAPH, seed, config, random_seed, overrides, weighted=weighted
    )
    deg_store = _WORKER_GRAPH.get_deg(weight=params["weight"])
    volume = sum(deg_store.get(node, 1) for node in nodes)
    return {
        "seed": seed,
        "random_seed": random_seed,
        "config": config,
        "weighted": weighted,
        "cluster": sorted(nodes),
        "volume": volume,
        "conductance": compute_conductance(
            _WORKER_GRAPH, nodes, weight=params["weight"]
        ),
        "signed_score": signed_stats["score"],
        "cut_pos": signed_stats["cut_pos"],
        "neg_inside": signed_stats["neg_inside"],
        "violations": signed_stats["violations"],
        "metadata": {key: metadata.get(key) for key in _META_KEYS if key in metadata},
    }


def _run_tasks(
    network_path: Path,
    tasks: list[tuple[str, Optional[int], str, Optional[dict], bool]],
    jobs: int,
    label: str,
) -> list[dict]:
    if not tasks:
        return []

    results: list[Optional[dict]] = [None] * len(tasks)
    started = time.perf_counter()
    report_every = max(1, len(tasks) // 10)

    if jobs <= 1:
        _worker_init(str(network_path))
        for index, task in enumerate(tasks, start=1):
            results[index - 1] = _worker_run(task)
            if index % report_every == 0 or index == len(tasks):
                elapsed = time.perf_counter() - started
                print(f"  {label}: {index}/{len(tasks)} runs ({elapsed:.1f}s)", flush=True)
    else:
        with ProcessPoolExecutor(
            max_workers=jobs,
            initializer=_worker_init,
            initargs=(str(network_path),),
        ) as executor:
            future_to_index = {
                executor.submit(_worker_run, task): index
                for index, task in enumerate(tasks)
            }
            for completed, future in enumerate(as_completed(future_to_index), start=1):
                results[future_to_index[future]] = future.result()
                if completed % report_every == 0 or completed == len(tasks):
                    elapsed = time.perf_counter() - started
                    print(
                        f"  {label}: {completed}/{len(tasks)} runs "
                        f"({elapsed:.1f}s, {elapsed / len(tasks):.2f}s/run wall)",
                        flush=True,
                    )

    return [item for item in results if item is not None]


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------
def jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    left_set, right_set = set(left), set(right)
    union = left_set | right_set
    if not union:
        return 1.0
    return len(left_set & right_set) / len(union)


def _percentile(sorted_values: list[float], quantile: float) -> float:
    if not sorted_values:
        return float("nan")
    position = (len(sorted_values) - 1) * quantile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight


def describe(values: Iterable[float]) -> dict:
    data = [float(v) for v in values if v is not None and math.isfinite(v)]
    if not data:
        return {"count": 0, "mean": None, "median": None, "min": None, "max": None,
                "p05": None, "p95": None, "stdev": None}
    ordered = sorted(data)
    return {
        "count": len(data),
        "mean": statistics.fmean(data),
        "median": statistics.median(data),
        "min": ordered[0],
        "max": ordered[-1],
        "p05": _percentile(ordered, 0.05),
        "p95": _percentile(ordered, 0.95),
        "stdev": statistics.pstdev(data) if len(data) > 1 else 0.0,
    }


def pairwise_jaccard(clusters: list[set[str]]) -> dict:
    values = [
        jaccard(left, right)
        for index, left in enumerate(clusters)
        for right in clusters[index + 1:]
    ]
    return describe(values)


def cluster_hash(nodes: Iterable[str]) -> str:
    return hashlib.sha1("\x1f".join(sorted(nodes)).encode("utf-8")).hexdigest()[:12]


def _intersection(clusters: list[set[str]]) -> set[str]:
    if not clusters:
        return set()
    result = set(clusters[0])
    for nodes in clusters[1:]:
        result &= nodes
    return result


def _union(clusters: list[set[str]]) -> set[str]:
    result: set[str] = set()
    for nodes in clusters:
        result |= nodes
    return result


def _freq_table(frequencies: dict[str, int], total: int) -> dict:
    buckets = Counter()
    for count in frequencies.values():
        buckets[round(count / total, 4)] += 1
    return {str(ratio): buckets[ratio] for ratio in sorted(buckets, reverse=True)}


def graph_thresholds(graph, params: dict) -> dict:
    m = max(graph.number_of_edges(), 1)
    b = params["b"]
    log_m_ceil = math.ceil(math.log2(max(m, 2)))
    deg_store = graph.get_deg(weight=params["weight"])
    total_volume = sum(deg_store.values())
    return {
        "phi_target": params["target_conductance"],
        "b": b,
        "alpha": params["teleport_alpha"],
        "epsilon": compute_epsilon_paper(b, m),
        "total_volume": total_volume,
        "volume_min": 2 ** (b - 1),
        # the ladder caps the window at half of the volume: beyond it a
        # prefix duplicates the conductance of its complement.
        "volume_max": int(total_volume // 2),
        "prob_threshold": 1 / (48 * log_m_ceil),
        "log_m_ceil": log_m_ceil,
        "min_cut_size": DEFAULT_MIN_CUT_SIZE,
    }


def condition_summary(results: list[dict]) -> dict:
    total = len(results)
    strategies = Counter((result["metadata"] or {}).get("strategy") for result in results)
    met_histogram = Counter()
    failed = Counter()
    met_values: list[int] = []
    for result in results:
        metadata = result["metadata"] or {}
        count = metadata.get("conditions_met_count")
        met_histogram[count] += 1
        if count is not None:
            met_values.append(int(count))
        for name in metadata.get("failed_conditions") or []:
            failed[name] += 1

    return {
        "runs": total,
        "strict_success": strategies.get("target_conductance", 0),
        "strict_success_rate": (
            strategies.get("target_conductance", 0) / total if total else 0.0
        ),
        "fallback_runs": total - strategies.get("target_conductance", 0),
        "seed_singleton_runs": strategies.get("seed_singleton_fallback", 0),
        "strategy_counts": {
            str(strategy): count for strategy, count in strategies.items()
        },
        "conditions_met_mean": statistics.fmean(met_values) if met_values else None,
        "conditions_met_histogram": {
            str(key): value
            for key, value in sorted(
                met_histogram.items(), key=lambda item: (item[0] is None, item[0])
            )
        },
        "failed_condition_counts": dict(failed),
        "pass_rate": {
            name: (1 - failed[name] / total) if total else 0.0
            for name in ("volume_ok", "prob_ok", "phi_ok")
        },
    }


def edge_sign_stats(graph) -> dict:
    positive = negative = zero = 0
    positive_weight = negative_weight = 0.0
    for _source, _target, data in graph.graph.edges(data=True):
        correlation = float(data.get("correlation", 0.0))
        if correlation > 0:
            positive += 1
            positive_weight += correlation
        elif correlation < 0:
            negative += 1
            negative_weight += -correlation
        else:
            zero += 1
    total = positive + negative + zero
    return {
        "positive_edges": positive,
        "negative_edges": negative,
        "zero_edges": zero,
        "positive_fraction": positive / total if total else 0.0,
        "negative_fraction": negative / total if total else 0.0,
        "mean_abs_positive_corr": positive_weight / positive if positive else 0.0,
        "mean_abs_negative_corr": negative_weight / negative if negative else 0.0,
    }


# ---------------------------------------------------------------------------
# Experiment 1: repeatability
# ---------------------------------------------------------------------------
def repeatability_experiment(
    graph,
    seed: str,
    config: str,
    *,
    runs: int,
    run_seed: int,
    jobs: int,
    network_path: Path,
    overrides: Optional[dict] = None,
    weighted: bool = True,
) -> dict:
    order_seeds = [run_seed + offset for offset in range(runs)]
    results = _run_tasks(
        network_path,
        [(seed, order_seed, config, overrides, weighted)
         for order_seed in order_seeds],
        jobs,
        label=f"repeatability[{config}:{seed}]",
    )

    clusters = [set(result["cluster"]) for result in results]
    sizes = [len(nodes) for nodes in clusters]
    intersection = _intersection(clusters)
    union = _union(clusters)
    frequencies = Counter()
    for nodes in clusters:
        frequencies.update(nodes)

    deterministic_nodes, deterministic_metadata, params, deterministic_signed = run_cluster(
        graph, seed, config, random_seed=None,
        parameter_overrides=overrides, weighted=weighted
    )
    deterministic_conductance = compute_conductance(
        graph, deterministic_nodes, weight=params["weight"]
    )
    distinct = {frozenset(nodes) for nodes in clusters}
    thresholds = graph_thresholds(graph, params)
    n_nodes = graph.number_of_nodes()
    mean_size = sum(sizes) / len(sizes)
    mean_volume = sum(result["volume"] for result in results) / len(results)
    jaccard_mean = pairwise_jaccard(clusters)
    jaccard_null = mean_size / (2 * n_nodes - mean_size)
    jaccard_obs = jaccard_mean["mean"] or 0.0
    jaccard_above_chance = (
        (jaccard_obs - jaccard_null) / (1 - jaccard_null)
        if jaccard_null < 1
        else 0.0
    )

    return {
        "experiment": "repeatability",
        "config": config,
        "config_label": CONFIGS[config]["label"],
        "seed_taxon": seed,
        "runs": runs,
        "ordering_seeds": order_seeds,
        "parameters": params,
        "thresholds": thresholds,
        "cluster_size": describe(sizes),
        "size_fraction": mean_size / n_nodes,
        "volume_fraction": mean_volume / thresholds["total_volume"],
        "volume_cap_fraction": mean_volume / thresholds["volume_max"],
        "jaccard_null": jaccard_null,
        "jaccard_above_chance": jaccard_above_chance,
        "intersection_size": len(intersection),
        "union_size": len(union),
        "intersection": sorted(intersection),
        "nodes_in_all_runs": len(intersection),
        "nodes_in_at_least_90pct": sum(
            1 for count in frequencies.values() if count / runs >= 0.9
        ),
        "nodes_in_at_least_80pct": sum(
            1 for count in frequencies.values() if count / runs >= 0.8
        ),
        "consensus_90pct": sorted(
            node for node, count in frequencies.items() if count / runs >= 0.9
        ),
        "inclusion_frequency": dict(sorted(frequencies.items())),
        "inclusion_frequency_histogram": _freq_table(frequencies, runs),
        "pairwise_jaccard": jaccard_mean,
        "jaccard_vs_intersection": describe(
            jaccard(nodes, intersection) for nodes in clusters
        ),
        "distinct_clusters": len(distinct),
        "identical_runs": max(Counter(frozenset(nodes) for nodes in clusters).values()),
        "seed_present_in_all_runs": bool(seed in intersection),
        "seed_presence_count": frequencies.get(seed, 0),
        "conductance": describe(result["conductance"] for result in results),
        "signed_score": describe(result["signed_score"] for result in results),
        "neg_inside": describe(result["neg_inside"] for result in results),
        "violations": describe(result["violations"] for result in results),
        "conditions": condition_summary(results),
        "deterministic_run": {
            "size": len(deterministic_nodes),
            "conductance": deterministic_conductance,
            "signed_score": deterministic_signed["score"],
            "cluster": sorted(deterministic_nodes),
            "matches_any_randomized_run": any(
                deterministic_nodes == nodes for nodes in clusters
            ),
            "jaccard_vs_intersection": jaccard(deterministic_nodes, intersection),
            "metadata": {
                key: deterministic_metadata.get(key)
                for key in _META_KEYS
                if key in deterministic_metadata
            },
        },
        "per_run": [
            {
                "random_seed": result["random_seed"],
                "size": len(nodes),
                "volume": result["volume"],
                "conductance": result["conductance"],
                "signed_score": result["signed_score"],
                "neg_inside": result["neg_inside"],
                "cluster_hash": cluster_hash(nodes),
                "strategy": (result["metadata"] or {}).get("strategy"),
                "conditions_met_count": (result["metadata"] or {}).get("conditions_met_count"),
                "failed_conditions": (result["metadata"] or {}).get("failed_conditions"),
            }
            for result, nodes in zip(results, clusters)
        ],
    }


# ---------------------------------------------------------------------------
# Experiment 2: nested consistency
# ---------------------------------------------------------------------------
def nested_consistency_experiment(
    graph,
    seed: str,
    config: str,
    *,
    jobs: int,
    network_path: Path,
    max_nodes: Optional[int] = None,
    overrides: Optional[dict] = None,
    weighted: bool = True,
) -> dict:
    """Experiment 2: re-seed the clustering from the members of its own cluster.

    The deterministic cluster of ``seed`` is computed first (the "base"
    cluster), then the clustering is re-run once from **every** taxon of that
    cluster.  The fixed run count of Experiment 1 does not apply here: the
    number of runs equals the size of the base cluster, because every member
    must be given the chance to recover the community.

    ``max_nodes`` optionally caps the number of re-seeded members by random
    subsampling (fixed RNG seed, for reproducible smoke tests only); when it
    is ``None`` all members are tested and ``all_nodes_tested`` is ``True``.
    """
    base_nodes, base_metadata, params, base_signed = run_cluster(
        graph, seed, config, random_seed=None,
        parameter_overrides=overrides, weighted=weighted
    )
    base_conductance = compute_conductance(graph, base_nodes, weight=params["weight"])

    test_nodes = sorted(base_nodes)
    sampled = False
    if max_nodes is not None and len(test_nodes) > max_nodes:
        test_nodes = sorted(
            random.Random(20241023).sample(test_nodes, max_nodes)
        )
        sampled = True

    results = _run_tasks(
        network_path,
        [(node, None, config, overrides, weighted) for node in test_nodes],
        jobs,
        label=f"nested[{config}:{seed}]",
    )

    nested = {result["seed"]: set(result["cluster"]) for result in results}
    clusters = [nested[node] for node in test_nodes if node in nested]
    total = len(clusters)
    intersection = _intersection(clusters)
    union = _union(clusters)

    support = {
        node: sum(1 for nodes in clusters if node in nodes) / total if total else 0.0
        for node in sorted(base_nodes)
    }
    core_nodes = {node for node, value in support.items() if value == 1.0}
    seed_recovery = sum(1 for nodes in clusters if seed in nodes) / total if total else 0.0
    self_recovery = (
        sum(1 for node, nodes in zip(test_nodes, clusters) if node in nodes) / total
        if total else 0.0
    )
    jaccard_vs_base = [jaccard(nodes, base_nodes) for nodes in clusters]
    escape_counts = [len(nodes - base_nodes) for nodes in clusters]

    both = only_left = only_right = neither = 0
    for i, left in enumerate(test_nodes):
        for right in test_nodes[i + 1:]:
            left_finds_right = right in nested.get(left, set())
            right_finds_left = left in nested.get(right, set())
            if left_finds_right and right_finds_left:
                both += 1
            elif left_finds_right:
                only_left += 1
            elif right_finds_left:
                only_right += 1
            else:
                neither += 1
    pairs = both + only_left + only_right + neither
    thresholds = graph_thresholds(graph, params)
    n_nodes = graph.number_of_nodes()
    mean_nested_size = sum(len(nodes) for nodes in clusters) / len(clusters)
    support_mean = sum(support.values()) / len(support) if support else 0.0
    support_null = mean_nested_size / n_nodes
    support_above_chance = (
        (support_mean - support_null) / (1 - support_null)
        if support_null < 1
        else 0.0
    )
    jaccard_mean = pairwise_jaccard(clusters)
    jaccard_null = mean_nested_size / (2 * n_nodes - mean_nested_size)
    jaccard_obs = jaccard_mean["mean"] or 0.0
    jaccard_above_chance = (
        (jaccard_obs - jaccard_null) / (1 - jaccard_null)
        if jaccard_null < 1
        else 0.0
    )

    return {
        "experiment": "nested_consistency",
        "config": config,
        "config_label": CONFIGS[config]["label"],
        "seed_taxon": seed,
        "nodes_tested": total,
        "nodes_in_base_cluster": len(base_nodes),
        "all_nodes_tested": not sampled,
        "parameters": params,
        "thresholds": thresholds,
        "base_size_fraction": len(base_nodes) / n_nodes,
        "mean_nested_size": mean_nested_size,
        "mean_nested_size_fraction": mean_nested_size / n_nodes,
        "base_cluster": {
            "seed": seed,
            "size": len(base_nodes),
            "conductance": base_conductance,
            "signed_score": base_signed["score"],
            "nodes": sorted(base_nodes),
            "metadata": {
                key: base_metadata.get(key) for key in _META_KEYS if key in base_metadata
            },
        },
        "intersection_size": len(intersection),
        "intersection": sorted(intersection),
        "union_size": len(union),
        "coverage_of_base_cluster": (
            len(intersection & base_nodes) / len(base_nodes) if base_nodes else 0.0
        ),
        "core_nodes_always_found": len(core_nodes),
        "nodes_support_90pct": sum(
            1 for value in support.values() if value >= 0.9
        ),
        "nodes_support_80pct": sum(
            1 for value in support.values() if value >= 0.8
        ),
        "consensus_80pct": sorted(
            node for node, value in support.items() if value >= 0.8
        ),
        "node_support": support,
        "node_support_stats": describe(support.values()),
        "support_null": support_null,
        "support_above_chance": support_above_chance,
        "seed_recovery_rate": seed_recovery,
        "self_recovery_rate": self_recovery,
        "cluster_size": describe(len(nodes) for nodes in clusters),
        "pairwise_jaccard": jaccard_mean,
        "jaccard_null": jaccard_null,
        "jaccard_above_chance": jaccard_above_chance,
        "identical_runs": max(
            Counter(frozenset(nodes) for nodes in clusters).values()
        ),
        "jaccard_vs_base_cluster": describe(jaccard_vs_base),
        "escape_nodes": describe(escape_counts),
        "conductance": describe(result["conductance"] for result in results),
        "signed_score": describe(result["signed_score"] for result in results),
        "neg_inside": describe(result["neg_inside"] for result in results),
        "conditions": condition_summary(results),
        "asymmetry": {
            "pairs": pairs,
            "both_directions": both,
            "only_left_finds_right": only_left,
            "only_right_finds_left": only_right,
            "neither_direction": neither,
            "asymmetry_rate": (only_left + only_right) / pairs if pairs else 0.0,
        },
        "per_run": [
            {
                "node": result["seed"],
                "size": len(nested[result["seed"]]),
                "volume": result["volume"],
                "conductance": result["conductance"],
                "signed_score": result["signed_score"],
                "contains_base_seed": seed in nested[result["seed"]],
                "jaccard_vs_base": jaccard(nested[result["seed"]], base_nodes),
                "cluster_hash": cluster_hash(nested[result["seed"]]),
                "strategy": (result["metadata"] or {}).get("strategy"),
                "conditions_met_count": (result["metadata"] or {}).get("conditions_met_count"),
                "failed_conditions": (result["metadata"] or {}).get("failed_conditions"),
            }
            for result in results
        ],
    }
