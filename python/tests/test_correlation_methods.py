"""Smoke tests for the statistical methods exposed by the Electron backend.

These test method semantics (real p-values vs StARS-selected network edges),
not only that each imported package can be called.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from python import analysis


@pytest.fixture(scope="module")
def abundance():
    rng = np.random.default_rng(12)
    n, p = 32, 5
    common = rng.normal(size=(n, 1))
    latent = 0.9 * common + 0.25 * rng.normal(size=(n, p))
    counts = np.maximum(1, np.rint(np.exp(2 + latent))).astype(float)
    counts[rng.random(counts.shape) < 0.1] = 0.0
    features = [f"taxon-{index}" for index in range(p)]
    samples = [f"sample-{index}" for index in range(n)]
    by_sample = {
        sample: {feature: float(counts[row, column]) for column, feature in enumerate(features)}
        for row, sample in enumerate(samples)
    }
    return by_sample, features, samples


def test_cclasso_integration_returns_pairwise_p_values(abundance):
    by_sample, features, samples = abundance
    correlations, p_values, diagnostics = analysis._compute_cclasso_for_matrix(
        by_sample,
        features,
        samples,
        {"cclassoBootstraps": 3, "cclassoSeed": 5},
    )

    assert diagnostics["evidence_type"] == "bootstrap_p_value"
    assert diagnostics["finite_pair_count"] == len(features) * (len(features) - 1) // 2
    values = [value for row in p_values.values() for value in row.values()]
    assert len(values) == diagnostics["finite_pair_count"] * 2
    assert all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in values)
    assert any(value not in (0.0, 1.0) for value in values)
    for left in features:
        for right, value in correlations[left].items():
            assert correlations[right][left] == value
            assert p_values[right][left] == p_values[left][right]


@pytest.mark.parametrize(
    ("method", "params", "expected_weight"),
    [
        (
            "spieceasi",
            {
                "spieceasiNLambdas": 4,
                "spieceasiNSubsamples": 4,
                "spieceasiStarsThreshold": 0.5,
                "spieceasiJobs": 1,
            },
            "selection_stability",
        ),
        (
            "spring",
            {
                "springNLambdas": 4,
                "springNSubsamples": 4,
                "springStarsThreshold": 0.8,
                "springSeed": 5,
            },
            "latent_correlation",
        ),
    ],
)
def test_network_methods_do_not_fabricate_p_values(abundance, method, params, expected_weight):
    by_sample, features, samples = abundance
    result = analysis._compute_network_for_matrix(method, by_sample, features, samples, params)
    weights, p_values, diagnostics = result

    assert diagnostics["weight_type"] == expected_weight
    assert diagnostics["evidence_type"] == "stars"
    assert diagnostics["edge_count"] > 0
    assert sum(len(row) for row in weights.values()) == diagnostics["edge_count"] * 2
    assert all(value is None for row in p_values.values() for value in row.values())

    rows = []
    seen = set()
    for left, neighbors in weights.items():
        for right, weight in neighbors.items():
            key = tuple(sorted((left, right)))
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "row": left,
                    "column": right,
                    "cor": weight,
                    "p": None,
                    "selected": True,
                    "weight_type": expected_weight,
                    "evidence_type": "stars",
                    "method": method,
                }
            )

    graph = analysis._build_graph_from_pairs(rows, p_value_threshold=None)
    assert len(graph["edges"]) == diagnostics["edge_count"]
    assert all(edge["p_value"] is None for edge in graph["edges"])
    assert graph["stats"]["p_value_threshold_applied"] is False


def test_network_edge_is_clusterable_without_a_pairwise_p_value():
    pytest.importorskip('networkx')
    from python.clustering.bacteria_graph import BacteriaGraph

    graph = BacteriaGraph.from_nodes_edges(
        ["taxon-a", "taxon-b"],
        [
            {
                "source": "taxon-a",
                "target": "taxon-b",
                "weight": 0.72,
                "p_value": None,
                "weight_type": "selection_stability",
            }
        ],
    )

    edge = graph.graph["taxon-a"]["taxon-b"]
    assert edge["p_value"] is None
    assert edge["weight"] == pytest.approx(0.72)


@pytest.mark.parametrize(
    ("method", "weight_type", "uses_p_values"),
    [
        ("spieceasi", "selection_stability", False),
        ("spring", "latent_correlation", False),
        ("cclasso", "correlation", True),
    ],
)
def test_run_summary_keeps_network_selection_distinct_from_p_values(
    tmp_path, monkeypatch, method, weight_type, uses_p_values
):
    samples = [f"s{index}" for index in range(40)]
    features = ["taxon-a", "taxon-b", "taxon-c"]
    abundance_path = tmp_path / "abundance.csv"
    metadata_path = tmp_path / "metadata.tsv"
    # taxon-b must not be an affine copy of taxon-a (index + 2 = (index + 1) + 1):
    # perfectly collinear columns drive CCLasso's ADMM to a NaN state, and the
    # port raises "missing value where TRUE/FALSE needed" instead of returning.
    abundance_path.write_text(
        "sample_id," + ",".join(features) + "\n" + "".join(
            f"{sample},{index + 1},{(index * 5) % 13 + 2},{(index * 3) % 7 + 1}\n"
            for index, sample in enumerate(samples)
        ),
        encoding="utf-8",
    )
    metadata_path.write_text(
        "sample_id\tCondition\n" + "".join(
            f"{sample}\t{'A' if index < 20 else 'B'}\n"
            for index, sample in enumerate(samples)
        ),
        encoding="utf-8",
    )

    def fake_network(
        method_arg, matrix_by_sample, features_arg, sample_ids, params, progress_range, **_kwargs
    ):
        # ``**_kwargs`` keeps this stand-in valid when the internal wrapper gains
        # provenance arguments (drop_events / scope_label).
        weights = {feature: {} for feature in features_arg}
        p_values = {feature: {} for feature in features_arg}
        for left, right in zip(features_arg, features_arg[1:]):
            weights[left][right] = weights[right][left] = 0.7
            p_values[left][right] = p_values[right][left] = None
        return weights, p_values, {
            "weight_type": weight_type,
            "evidence_type": "stars",
            "edge_count": len(features_arg) - 1,
            "selected_lambda": 0.1,
        }

    monkeypatch.setattr(analysis, "_compute_network_for_matrix", fake_network)
    result = analysis.run(
        "test-job",
        {
            "bacteria": [str(abundance_path)],
            "bacteria_metadata": [str(metadata_path)],
        },
        {
            "datasetType": "bacteria",
            "correlationMethod": method,
            "conditionColumn": "Condition",
            "conditionA": "A",
            "conditionB": "B",
            "groupingMode": "all",
            # Network methods must not use this even when it would reject every
            # ordinary pairwise p-value.
            "pValueThreshold": 0.0,
            "correlationThreshold": 0.0,
            "cclassoBootstraps": 3,
        },
    )

    taxa_report = result["summary"]["taxa_filtering"]
    assert taxa_report["meta"]["scope"] == "all_samples_in_abundance_table"
    assert taxa_report["meta"]["relative_abundance_denominator"] == "domain"
    assert isinstance(taxa_report["meta"]["method_removals"], list)

    assert result["summary"]["p_value_threshold_applied"] is uses_p_values
    assert (result["summary"]["p_value_threshold"] is not None) is uses_p_values
    for condition_graph in result["graph_pipeline"]["condition_graphs"]:
        assert condition_graph["graph"]["stats"]["weight_type"] == weight_type
        if uses_p_values:
            assert all(edge["p_value"] is not None for edge in condition_graph["graph"]["edges"])
        else:
            assert condition_graph["graph"]["edges"]
            assert all(edge["p_value"] is None for edge in condition_graph["graph"]["edges"])
        Path(condition_graph["pair_file"]["path"]).unlink(missing_ok=True)
