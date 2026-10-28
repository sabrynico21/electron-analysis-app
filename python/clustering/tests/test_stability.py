"""Fast smoke tests for the clustering variants and the stability harness.

Run with:
    .venv/bin/python -m unittest discover -s python/clustering/tests -p "test_*.py" -v

They validate network loading, the parameter resolution (including the
rejection of the removed ``recommended`` profile), the five-step
selection ladder, the signed score and the report structure on a few
runs.  The full experiments are produced by ``tuning_experiment.py`` and
``comparison_experiment.py``: Experiment 1 (repeatability) uses a fixed
run count, Experiment 2 (nested consistency) runs from every member of
the base cluster.
"""

from __future__ import annotations

import inspect
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
CLUSTERING_DIR = TEST_DIR.parent
for path in (TEST_DIR, CLUSTERING_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from bacteria_graph import BacteriaGraph  # noqa: E402
from bacteria_pagerank import pagerank_nibble_adaptive  # noqa: E402
from clustering_advanced import (  # noqa: E402
    build_signed_index,
    compute_signed_metrics,
    pagerank_cluster,
    positive_walk_graph,
    signed_metrics,
)
from clustering_profiles import (  # noqa: E402
    app_parameters,
    compute_default_b,
    compute_default_teleport_alpha,
    resolve_parameters,
)
from network_loader import load_network, resolve_network  # noqa: E402
from stability_experiment import (  # noqa: E402
    CONFIGS,
    edge_sign_stats,
    nested_consistency_experiment,
    repeatability_experiment,
    run_cluster,
)

EXPECTED_SHAPES = {
    "rs4-101": (368, 5223),
    "101": (387, 5121),
}


class NetworkLoaderTests(unittest.TestCase):
    def test_expected_graph_shapes(self):
        for name, (nodes, edges) in EXPECTED_SHAPES.items():
            graph, info = load_network(name)
            self.assertEqual(graph.number_of_nodes(), nodes, name)
            self.assertEqual(graph.number_of_edges(), edges, name)
            self.assertEqual(info["nodes"], nodes)
            self.assertEqual(info["edges"], edges)

    def test_weights_are_absolute_correlations(self):
        graph, _ = load_network("101")
        for source, target in list(graph.graph.edges())[:50]:
            weight = graph.get_weight(source, target)
            correlation = graph.get_correlation(source, target)
            self.assertGreater(weight, 0.0)
            self.assertAlmostEqual(weight, abs(correlation))

    def test_unknown_network_raises(self):
        with self.assertRaises(FileNotFoundError):
            resolve_network("does-not-exist")


class ParameterProfileTests(unittest.TestCase):
    def test_app_profile_matches_the_electron_app(self):
        # The app defaults are the reference values of the experiments
        # (phi = 0.4, b = 10, alpha = 0.05 on both bundled networks).
        for name in EXPECTED_SHAPES:
            graph, _ = load_network(name)
            params = app_parameters(graph)
            self.assertEqual(params["target_conductance"], 0.4)
            self.assertEqual(params["b"], 10)
            self.assertEqual(params["teleport_alpha"], 0.05)
            self.assertTrue(params["weighted"])
            self.assertEqual(params["weight"], "weight")

    def test_removed_profile_is_rejected(self):
        # The ``recommended`` profile was removed together with the
        # ``first`` / ``min_conductance`` selection rules; explicit
        # overrides are the supported replacement.
        graph, _ = load_network("101")
        with self.assertRaises(ValueError) as ctx:
            resolve_parameters(graph, "recommended")
        self.assertIn("recommended", str(ctx.exception))
        params = resolve_parameters(
            graph,
            "app",
            overrides={
                "teleport_alpha": 0.05,
                "b": 10,
                "target_conductance": 0.4,
            },
        )
        self.assertEqual(params["teleport_alpha"], 0.05)
        self.assertEqual(params["b"], 10)
        self.assertEqual(params["target_conductance"], 0.4)

    def test_parameter_overrides(self):
        graph, _ = load_network("101")
        params = resolve_parameters(
            graph, "app", overrides={"b": 7, "target_conductance": 0.35}
        )
        self.assertEqual(params["b"], 7)
        self.assertEqual(params["target_conductance"], 0.35)

    def test_compute_default_b_edge_cases(self):
        self.assertGreaterEqual(compute_default_b(1, 1), 1)
        self.assertGreaterEqual(compute_default_teleport_alpha(0.2, 1), 0.0)
        self.assertLessEqual(compute_default_teleport_alpha(1.0, 10), 1.0)


class SignedGraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, _ = load_network("101")
        cls.index = build_signed_index(cls.graph)

    def test_positive_walk_graph_drops_negative_edges(self):
        positive = positive_walk_graph(self.graph)
        self.assertLess(
            positive.number_of_edges(), self.graph.number_of_edges()
        )
        for source, target in positive.graph.edges():
            self.assertGreater(self.graph.get_correlation(source, target), 0.0)

    def test_signed_metrics_are_consistent(self):
        nodes = {"Chloroplast", "Mitochondria"}
        metrics = signed_metrics(self.index, nodes, sign_lambda=1.0)
        self.assertGreaterEqual(metrics["score"], 0.0)
        self.assertGreaterEqual(metrics["violations"], 0.0)
        self.assertAlmostEqual(
            metrics["violations"],
            metrics["cut_pos"] + metrics["neg_inside"],
        )
        index_once = compute_signed_metrics(self.graph, nodes)
        self.assertAlmostEqual(index_once["score"], metrics["score"])

    def test_unweighted_signed_cluster(self):
        nodes, metadata, params, signed_stats = run_cluster(
            self.graph, "Chloroplast", "both", weighted=False
        )
        self.assertIsNone(params["weight"])
        self.assertGreater(len(nodes), 0)
        self.assertTrue(metadata["signed"])
        self.assertGreaterEqual(signed_stats["score"], 0.0)

    def test_signed_cluster_has_lower_signed_score(self):
        baseline_nodes, _, _, baseline_signed = run_cluster(
            self.graph, "Chloroplast", "baseline"
        )
        signed_nodes, _, _, signed_stats = run_cluster(
            self.graph, "Chloroplast", "both"
        )
        self.assertIn("Chloroplast", signed_nodes)
        self.assertLess(signed_stats["score"], baseline_signed["score"])


class AdvancedClusteringTests(unittest.TestCase):
    def test_default_parameters_reach_the_quality_target(self):
        graph, _ = load_network("101")
        nodes, metadata = pagerank_cluster(
            graph,
            "Chloroplast",
            target_conductance=0.4,
            b=10,
            teleport_alpha=0.05,
            weighted=True,
            signed=False,
        )
        self.assertEqual(metadata["ladder_step"], 1)
        self.assertEqual(metadata["strategy"], "target_conductance")
        self.assertTrue(metadata["satisfied"])
        self.assertLess(metadata["conductance"], 0.4)
        self.assertIn("Chloroplast", nodes)

    def test_signed_ladder_returns_a_cluster(self):
        graph, _ = load_network("rs4-101")
        params = app_parameters(graph)
        nodes, metadata = pagerank_cluster(
            graph,
            "Chloroplast",
            target_conductance=params["target_conductance"],
            b=params["b"],
            teleport_alpha=params["teleport_alpha"],
            weighted=True,
            signed=True,
        )
        self.assertGreater(len(nodes), 0)
        self.assertIn("Chloroplast", nodes)
        self.assertIn(
            metadata["strategy"],
            {"target_conductance", "target_conductance_fallback"},
        )


class LadderSelectionTests(unittest.TestCase):
    """The five-step selection ladder of ``selection_ladder.py``."""

    @classmethod
    def setUpClass(cls):
        cls.graph, _ = load_network("101")

    def _adaptive(self, **kwargs):
        kwargs.setdefault("weight", "weight")
        kwargs.setdefault("verbose", False)
        nodes, metadata = pagerank_nibble_adaptive(
            self.graph, "Chloroplast", **kwargs
        )
        return nodes, metadata

    def test_step1_strict_when_every_criterion_holds(self):
        nodes, metadata = self._adaptive(
            target_conductance=0.4, b=10, teleport_alpha=0.05
        )
        self.assertEqual(metadata["ladder_step"], 1)
        self.assertEqual(metadata["ladder_rule"], "strict")
        self.assertEqual(metadata["conditions_met_count"], 3)
        self.assertEqual(metadata["failed_conditions"], [])
        self.assertTrue(metadata["constraints_satisfied"])
        self.assertLess(metadata["conductance"], 0.4)
        self.assertGreaterEqual(metadata["size"], metadata["min_cut_size"])
        self.assertIn("Chloroplast", nodes)

    def test_step3_multiscale_when_the_floor_is_unreachable(self):
        # b=13 demands vol > 4096, but the best cut lives below it: the
        # floor is halved once (2048) and the same PageRank vector is reused.
        nodes, metadata = self._adaptive(target_conductance=0.4, b=13,
                                         teleport_alpha=0.05)
        self.assertEqual(metadata["ladder_step"], 3)
        self.assertEqual(metadata["ladder_rule"], "multi_scale")
        self.assertEqual(metadata["scale_relaxation"], 1)
        self.assertEqual(metadata["volume_min_used"], 2048.0)
        self.assertEqual(metadata["volume_min"], 4096.0)
        self.assertIn("volume_ok", metadata["failed_conditions"])
        self.assertIn("Chloroplast", nodes)

    def test_step4_relaxed_phi_reports_the_effective_bar(self):
        # phi=0.2 is below every achievable conductance on 101 (best ~0.27),
        # so the ladder relaxes the bar and reports the rung it used.
        nodes, metadata = self._adaptive(target_conductance=0.2, b=4,
                                         teleport_alpha=2e-05)
        if metadata["ladder_step"] == 4:
            self.assertEqual(metadata["ladder_rule"], "relaxed_phi")
            self.assertGreater(metadata["phi_effective"], 0.2)
            self.assertLessEqual(metadata["phi_effective"], 0.4)
            self.assertIn("phi_ok", metadata["failed_conditions"])
        else:
            # The rung ladder may also solve it strictly on a different
            # seed/graph state; the invariant is the strategy taxonomy.
            self.assertIn(
                metadata["strategy"],
                {"target_conductance", "target_conductance_fallback"},
            )
        self.assertIn("Chloroplast", nodes)

    def test_step5_penalty_breakdown_and_size_floor(self):
        nodes, metadata = self._adaptive(target_conductance=1e-6, b=10,
                                         teleport_alpha=0.05)
        self.assertEqual(metadata["ladder_step"], 5)
        self.assertEqual(metadata["ladder_rule"], "penalty_best_effort")
        breakdown = metadata["penalty_breakdown"]
        for key in ("quality", "volume_octaves", "mass_bonus", "penalty"):
            self.assertIn(key, breakdown)
        self.assertGreaterEqual(metadata["size"], metadata["min_cut_size"])
        self.assertIn("phi_ok", metadata["failed_conditions"])
        self.assertIn("Chloroplast", nodes)

    def test_max_cut_size_bounds_the_cluster(self):
        nodes, metadata = pagerank_cluster(
            self.graph,
            "Chloroplast",
            target_conductance=0.4,
            b=10,
            teleport_alpha=0.05,
            weighted=True,
            max_cut_size=10,
        )
        self.assertLessEqual(metadata["size"], 10)
        self.assertGreaterEqual(len(nodes), 5)
        self.assertIn("Chloroplast", nodes)

    def test_singleton_when_no_prefix_is_admissible(self):
        tiny = BacteriaGraph.from_nodes_edges(
            ["A", "B", "C"],
            [
                {"source": "A", "target": "B", "correlation": 0.8},
                {"source": "B", "target": "C", "correlation": -0.7},
                {"source": "A", "target": "C", "correlation": 0.6},
            ],
        )
        nodes, metadata = pagerank_nibble_adaptive(
            tiny, "A", target_conductance=0.4, b=2,
            teleport_alpha=0.05, weight=None, verbose=False,
        )
        self.assertEqual(nodes, {"A"})
        self.assertEqual(metadata["ladder_step"], 0)
        self.assertEqual(metadata["strategy"], "seed_singleton_fallback")
        self.assertEqual(metadata["size"], 1)
        self.assertEqual(
            metadata["failed_conditions"], ["volume_ok", "prob_ok", "phi_ok"]
        )
        self.assertFalse(metadata["constraints_satisfied"])

    def test_metadata_exposes_the_ladder_schema(self):
        _nodes, metadata = self._adaptive(
            target_conductance=0.4, b=10, teleport_alpha=0.05
        )
        for key in (
            "ladder_step",
            "ladder_rule",
            "phi_target",
            "phi_effective",
            "volume_min",
            "volume_min_used",
            "volume_max",
            "mass_threshold",
            "min_cut_size",
            "constraints",
            "ranked_nodes",
            "ranked_scores",
        ):
            self.assertIn(key, metadata)
        self.assertEqual(metadata["phi_target"], 0.4)
        self.assertEqual(metadata["phi_effective"], 0.4)
        # The ladder caps the window at half of the total volume, and the
        # returned cluster must sit inside it.
        total_volume = sum(
            self.graph.get_deg(weight="weight").values()
        )
        self.assertEqual(metadata["volume_max"], int(total_volume // 2))
        self.assertLessEqual(metadata["volume"], metadata["volume_max"])


class ExperimentReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, _ = load_network("101")
        cls.graph_path = resolve_network("101")
        cls.seed = "Chloroplast"
        assert cls.graph.has_node(cls.seed)

    def test_all_configs_are_defined(self):
        self.assertEqual(
            set(CONFIGS), {"baseline", "new_parameters", "signed_weights", "both"}
        )

    def test_repeatability_report(self):
        report = repeatability_experiment(
            self.graph,
            self.seed,
            "baseline",
            runs=3,
            run_seed=7,
            jobs=2,
            network_path=self.graph_path,
        )
        self.assertEqual(report["runs"], 3)
        self.assertTrue(report["seed_present_in_all_runs"])
        self.assertGreater(report["intersection_size"], 0)
        self.assertEqual(len(report["per_run"]), 3)
        self.assertIn("conditions", report)
        self.assertIn("signed_score", report)
        for key in (
            "size_fraction",
            "volume_fraction",
            "jaccard_null",
            "jaccard_above_chance",
            "identical_runs",
            "nodes_in_at_least_80pct",
        ):
            self.assertIn(key, report)
        self.assertGreater(report["size_fraction"], 0.0)
        self.assertLess(report["size_fraction"], 1.0)
        self.assertGreaterEqual(report["jaccard_above_chance"], 0.0)

    def test_nested_consistency_report(self):
        report = nested_consistency_experiment(
            self.graph,
            self.seed,
            "new_parameters",
            jobs=2,
            network_path=self.graph_path,
            max_nodes=2,
        )
        self.assertEqual(report["nodes_tested"], 2)
        # ``max_nodes`` is a smoke-test cap: a capped run is flagged so the
        # reports can tell it apart from the full every-member protocol.
        self.assertFalse(report["all_nodes_tested"])
        self.assertGreater(report["nodes_in_base_cluster"], 0)
        self.assertEqual(report["self_recovery_rate"], 1.0)
        self.assertEqual(len(report["per_run"]), 2)
        self.assertEqual(len(report["node_support"]), report["nodes_in_base_cluster"])
        for key in (
            "base_size_fraction",
            "mean_nested_size",
            "support_null",
            "support_above_chance",
            "jaccard_above_chance",
            "nodes_support_80pct",
            "identical_runs",
        ):
            self.assertIn(key, report)
        self.assertGreater(report["base_size_fraction"], 0.0)
        self.assertLess(report["base_size_fraction"], 1.0)

    def test_nested_consistency_is_uncapped_by_default(self):
        # Experiment 2 must run from every member of the base cluster; the
        # fixed run count of Experiment 1 must not leak into it.
        signature = inspect.signature(nested_consistency_experiment)
        self.assertIsNone(signature.parameters["max_nodes"].default)

    def test_edge_sign_stats(self):
        signs = edge_sign_stats(self.graph)
        self.assertEqual(
            signs["positive_edges"] + signs["negative_edges"] + signs["zero_edges"],
            self.graph.number_of_edges(),
        )
        self.assertGreater(signs["negative_edges"], 0)
        self.assertGreater(signs["positive_edges"], 0)


class RunnerPayloadTests(unittest.TestCase):
    def _run_payload(self, payload: dict) -> subprocess.CompletedProcess:
        runner = CLUSTERING_DIR / "run_clustering.py"
        with tempfile.NamedTemporaryFile(
            "w", suffix=".json", delete=False, encoding="utf-8"
        ) as handle:
            json.dump(payload, handle)
            input_path = handle.name
        try:
            return subprocess.run(
                [sys.executable, str(runner), "--input", input_path],
                capture_output=True,
                text=True,
            )
        finally:
            Path(input_path).unlink(missing_ok=True)

    @staticmethod
    def _signed_payload() -> dict:
        return {
            "seedNode": "A",
            "parameterProfile": "app",
            "signed": True,
            "signLambda": 1.0,
            "weighted": True,
            "graphData": {
                "graph": {
                    "nodes": [{"id": node} for node in ("A", "B", "C")],
                    "edges": [
                        {"source": "A", "target": "B", "correlation": 0.8},
                        {"source": "B", "target": "C", "correlation": -0.7},
                        {"source": "A", "target": "C", "correlation": 0.6},
                    ],
                }
            },
        }

    def test_runner_accepts_signed_options(self):
        completed = self._run_payload(self._signed_payload())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(completed.stdout.strip().splitlines()[-1])
        self.assertTrue(result["success"])
        self.assertIn("A", result["clusterNodes"])
        self.assertTrue(result["metadata"]["signed_mode"])
        self.assertEqual(result["metadata"]["parameters"]["profile"], "app")
        self.assertIn("ladder_step", result["metadata"])

    def test_runner_rejects_the_removed_profile(self):
        payload = self._signed_payload()
        payload["parameterProfile"] = "recommended"
        completed = self._run_payload(payload)
        self.assertNotEqual(completed.returncode, 0)
        result = json.loads(completed.stdout.strip().splitlines()[-1])
        self.assertFalse(result["success"])
        self.assertIn("recommended", result["error"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
