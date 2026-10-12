from __future__ import annotations

import itertools
from collections.abc import Callable

import numpy as np
import pytest

from spieceasi.pipeline import spiec_easi_mb
from tests.synthetic_data import (
    generate_block_signal_counts,
    generate_hub_signal_counts,
    generate_strong_signal_counts,
)
from tests.validation_config import (
    MIN_EXPANDED_SCENARIO_FAMILIES,
    MIN_EXPANDED_SCENARIOS,
    STOCHASTIC_BASELINE_TOLERANCE,
)


def _edge_set(A: np.ndarray) -> set[tuple[int, int]]:
    upper = np.triu_indices(A.shape[0], k=1)
    return {(int(i), int(j)) for i, j in zip(*upper, strict=True) if int(A[i, j]) != 0}


def _jaccard_distance(a: set[tuple[int, int]], b: set[tuple[int, int]]) -> float:
    union = a.union(b)
    if not union:
        return 0.0
    return 1.0 - (len(a.intersection(b)) / len(union))


def _degree_mae(A: np.ndarray, B: np.ndarray) -> float:
    deg_a = np.sum(A != 0, axis=1).astype(float)
    deg_b = np.sum(B != 0, axis=1).astype(float)
    return float(np.mean(np.abs(deg_a - deg_b)))


def _selected_index_from_monotone(monotone: np.ndarray, threshold: float) -> int:
    feasible = np.where(monotone <= threshold)[0]
    return int(feasible[-1]) if feasible.size > 0 else 0


def _run_r_pipeline(
    r_runtime,
    counts_df,
    *,
    nlambda: int,
    lambda_min_ratio: float,
    n_subsamples: int,
    threshold: float,
    seed: int,
) -> tuple[np.ndarray, int, int]:
    robjects = r_runtime["robjects"]
    pandas2ri = r_runtime["pandas2ri"]
    localconverter = r_runtime["localconverter"]

    robjects.r("""
        suppressMessages(library(SpiecEasi))

        run_spieceasi_baseline <- function(
            count_df, nlambda, lambda_min_ratio, rep_num, threshold, seed
        ) {
          set.seed(seed)
          X <- as.matrix(count_df)
          fit <- spiec.easi(
            X,
            method = 'mb',
            sel.criterion = 'stars',
                        sym = 'or',
            nlambda = nlambda,
            lambda.min.ratio = lambda_min_ratio,
            pulsar.params = list(rep.num = rep_num, thresh = threshold, seed = seed)
          )

          adj <- as.matrix(getRefit(fit))
          stars_summary <- as.numeric(fit$select$stars$summary)
          opt_index <- as.integer(fit$select$stars$opt.index)

          list(
            adjacency = adj,
            stars_summary = stars_summary,
            opt_index = opt_index
          )
        }
        """)
    run_fn = robjects.globalenv["run_spieceasi_baseline"]

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_counts_df = robjects.conversion.py2rpy(counts_df)

    out = run_fn(
        r_counts_df,
        nlambda,
        lambda_min_ratio,
        n_subsamples,
        threshold,
        seed,
    )
    r_adj = (np.asarray(out.rx2("adjacency"), dtype=float) != 0.0).astype(int)
    stars_summary = np.asarray(out.rx2("stars_summary"), dtype=float)
    opt_index_0b = int(np.asarray(out.rx2("opt_index"), dtype=int)[0]) - 1
    return r_adj, opt_index_0b, int(stars_summary.size)


def _expanded_scenarios() -> (
    list[tuple[str, dict[str, int | float], Callable[..., tuple[object, np.ndarray]]]]
):
    scenarios = [
        (
            "chain_small",
            dict(n_samples=160, n_taxa=16, depth=4500, rho=0.35, seed=7),
            generate_strong_signal_counts,
        ),
        (
            "chain_medium",
            dict(n_samples=190, n_taxa=20, depth=5200, rho=0.40, seed=31),
            generate_strong_signal_counts,
        ),
        (
            "block_small",
            dict(n_samples=170, n_taxa=18, depth=4800, n_blocks=3, seed=13),
            generate_block_signal_counts,
        ),
        (
            "block_medium",
            dict(n_samples=210, n_taxa=24, depth=6000, n_blocks=4, seed=23),
            generate_block_signal_counts,
        ),
        (
            "hub_small",
            dict(n_samples=170, n_taxa=18, depth=5000, n_hubs=2, seed=19),
            generate_hub_signal_counts,
        ),
        (
            "hub_medium",
            dict(n_samples=210, n_taxa=24, depth=6200, n_hubs=3, seed=29),
            generate_hub_signal_counts,
        ),
    ]

    families = {name.split("_")[0] for name, _, _ in scenarios}
    assert len(families) >= MIN_EXPANDED_SCENARIO_FAMILIES
    assert len(scenarios) >= MIN_EXPANDED_SCENARIOS
    return scenarios


@pytest.mark.r_validation
@pytest.mark.slow
def test_stochastic_variability_baseline_expanded_scenarios(r_runtime) -> None:
    nlambda = 20
    lambda_min_ratio = 1e-3
    n_subsamples = 8
    threshold = 0.05
    benchmark_seeds = [5, 11, 17, 23]

    py_py_jaccard: list[float] = []
    py_py_degree: list[float] = []
    py_py_lambda_rank: list[float] = []

    r_r_jaccard: list[float] = []
    r_r_degree: list[float] = []
    r_r_lambda_rank: list[float] = []

    py_r_jaccard: list[float] = []
    py_r_degree: list[float] = []
    py_r_lambda_rank: list[float] = []
    scenario_rows: list[str] = []

    for scenario_name, cfg, generator in _expanded_scenarios():
        counts_df, _ = generator(**cfg)

        py_by_seed: dict[int, tuple[np.ndarray, int, int]] = {}
        r_by_seed: dict[int, tuple[np.ndarray, int, int]] = {}
        scenario_py_r_jaccard: list[float] = []
        scenario_py_r_degree: list[float] = []
        scenario_py_r_lambda_rank: list[float] = []
        scenario_py_opt_rel: list[float] = []
        scenario_r_opt_rel: list[float] = []
        scenario_py_boundary0: list[int] = []
        scenario_r_boundary0: list[int] = []
        py_edge_counts: list[int] = []
        r_edge_counts: list[int] = []

        for seed in benchmark_seeds:
            py_result = spiec_easi_mb(
                counts_df,
                n_lambdas=nlambda,
                lambda_min_ratio=lambda_min_ratio,
                n_subsamples=n_subsamples,
                stars_threshold=threshold,
                sym_rule="or",
                random_state=seed,
                n_jobs=1,
            )
            py_opt_0b = _selected_index_from_monotone(
                py_result.monotone_instability,
                threshold,
            )
            py_by_seed[seed] = (
                py_result.adjacency,
                py_opt_0b,
                int(py_result.alpha_path.size),
            )

            r_by_seed[seed] = _run_r_pipeline(
                r_runtime,
                counts_df,
                nlambda=nlambda,
                lambda_min_ratio=lambda_min_ratio,
                n_subsamples=n_subsamples,
                threshold=threshold,
                seed=seed,
            )

        for s1, s2 in itertools.combinations(benchmark_seeds, 2):
            py_a, py_idx_a, py_len = py_by_seed[s1]
            py_b, py_idx_b, _ = py_by_seed[s2]
            py_py_jaccard.append(_jaccard_distance(_edge_set(py_a), _edge_set(py_b)))
            py_py_degree.append(_degree_mae(py_a, py_b))
            py_py_lambda_rank.append(abs(py_idx_a - py_idx_b) / max(py_len - 1, 1))

            r_a, r_idx_a, r_len = r_by_seed[s1]
            r_b, r_idx_b, _ = r_by_seed[s2]
            r_r_jaccard.append(_jaccard_distance(_edge_set(r_a), _edge_set(r_b)))
            r_r_degree.append(_degree_mae(r_a, r_b))
            r_r_lambda_rank.append(abs(r_idx_a - r_idx_b) / max(r_len - 1, 1))

        for seed in benchmark_seeds:
            py_adj, py_idx, py_len = py_by_seed[seed]
            r_adj, r_idx, r_len = r_by_seed[seed]
            rank_den = max(min(py_len, r_len) - 1, 1)

            # Boundary diagnostics: relative position of the StARS optimum along the
            # lambda path. 0.0 == largest lambda (empty/degenerate, R warns "optimal
            # lambda may be larger than the supplied values"); 1.0 == smallest lambda.
            scenario_py_opt_rel.append(py_idx / max(py_len - 1, 1))
            scenario_r_opt_rel.append(r_idx / max(r_len - 1, 1))
            scenario_py_boundary0.append(1 if py_idx == 0 else 0)
            scenario_r_boundary0.append(1 if r_idx == 0 else 0)

            py_edges = _edge_set(py_adj)
            r_edges = _edge_set(r_adj)
            py_edge_counts.append(len(py_edges))
            r_edge_counts.append(len(r_edges))

            jaccard = _jaccard_distance(py_edges, r_edges)
            degree = _degree_mae(py_adj, r_adj)
            lambda_rank = abs(py_idx - r_idx) / rank_den

            py_r_jaccard.append(jaccard)
            py_r_degree.append(degree)
            py_r_lambda_rank.append(lambda_rank)

            scenario_py_r_jaccard.append(jaccard)
            scenario_py_r_degree.append(degree)
            scenario_py_r_lambda_rank.append(lambda_rank)

        scenario_rows.append(
            f"{scenario_name},"
            f"{np.median(scenario_py_r_jaccard):.6f},"
            f"{np.median(scenario_py_r_degree):.6f},"
            f"{np.median(scenario_py_r_lambda_rank):.6f},"
            f"{np.median(py_edge_counts):.1f},"
            f"{np.median(r_edge_counts):.1f},"
            f"{np.median(scenario_py_opt_rel):.4f},"
            f"{np.median(scenario_r_opt_rel):.4f},"
            f"{int(np.sum(scenario_py_boundary0))},"
            f"{int(np.sum(scenario_r_boundary0))}"
        )

    baseline_jaccard = max(
        float(np.median(py_py_jaccard)), float(np.median(r_r_jaccard))
    )
    baseline_degree = max(float(np.median(py_py_degree)), float(np.median(r_r_degree)))
    baseline_lambda_rank = max(
        float(np.median(py_py_lambda_rank)),
        float(np.median(r_r_lambda_rank)),
    )

    cross_jaccard = float(np.median(py_r_jaccard))
    cross_degree = float(np.median(py_r_degree))
    cross_lambda_rank = float(np.median(py_r_lambda_rank))

    summary = (
        "metric,py_py_median,r_r_median,cross_py_r_median,baseline,margin\n"
        f"edge_jaccard_distance,{np.median(py_py_jaccard):.6f},{np.median(r_r_jaccard):.6f},"
        f"{cross_jaccard:.6f},{baseline_jaccard:.6f},"
        f"{STOCHASTIC_BASELINE_TOLERANCE.edge_jaccard_distance_margin:.6f}\n"
        f"degree_mae,{np.median(py_py_degree):.6f},{np.median(r_r_degree):.6f},"
        f"{cross_degree:.6f},{baseline_degree:.6f},"
        f"{STOCHASTIC_BASELINE_TOLERANCE.degree_mae_margin:.6f}\n"
        f"lambda_rank_gap,{np.median(py_py_lambda_rank):.6f},{np.median(r_r_lambda_rank):.6f},"
        f"{cross_lambda_rank:.6f},{baseline_lambda_rank:.6f},"
        f"{STOCHASTIC_BASELINE_TOLERANCE.lambda_rank_margin:.6f}\n"
        "scenario,py_r_jaccard_median,py_r_degree_median,py_r_lambda_rank_median,"
        "py_edges_median,r_edges_median,py_opt_rel_median,r_opt_rel_median,"
        "py_boundary0,r_boundary0\n" + "\n".join(scenario_rows)
    )

    assert cross_jaccard <= (
        baseline_jaccard + STOCHASTIC_BASELINE_TOLERANCE.edge_jaccard_distance_margin
    ), summary
    assert cross_degree <= (
        baseline_degree + STOCHASTIC_BASELINE_TOLERANCE.degree_mae_margin
    ), summary
    assert cross_lambda_rank <= (
        baseline_lambda_rank + STOCHASTIC_BASELINE_TOLERANCE.lambda_rank_margin
    ), summary
