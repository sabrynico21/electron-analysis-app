from __future__ import annotations

import numpy as np
import pytest

from spieceasi.neighborhood import mb_neighborhood_selection
from spieceasi.transforms import clr_transform
from tests.synthetic_data import generate_compositional_counts
from tests.validation_config import (
    FIXED_LAMBDA_AMBIGUOUS_EDGE_WEIGHT_MAX,
    FIXED_LAMBDA_WEIGHT_TOLERANCE,
    SUPPORT_ZERO_THRESHOLD,
)


def _r_style_lambda_path(
    X: np.ndarray,
    *,
    nlambda: int,
    lambda_min_ratio: float,
) -> np.ndarray:
    corr = np.corrcoef(X, rowvar=False)
    upper = np.triu_indices(corr.shape[0], k=1)
    lambda_max = float(np.max(np.abs(corr[upper])))
    lambda_min = lambda_max * lambda_min_ratio
    return np.exp(np.linspace(np.log(lambda_max), np.log(lambda_min), nlambda))


def _dataset_cfgs() -> list[dict[str, float | int]]:
    return [
        dict(n_samples=90, n_taxa=10, depth=2500, rho=0.25, seed=5),
        dict(n_samples=110, n_taxa=12, depth=2800, rho=0.30, seed=9),
        dict(n_samples=140, n_taxa=14, depth=3500, rho=0.35, seed=11),
        dict(n_samples=180, n_taxa=18, depth=4800, rho=0.40, seed=17),
        dict(n_samples=220, n_taxa=24, depth=6000, rho=0.45, seed=31),
    ]


def _lambda_values(X: np.ndarray) -> np.ndarray:
    alpha_path = _r_style_lambda_path(X, nlambda=25, lambda_min_ratio=1e-2)
    return alpha_path[[0, 4, 8, 12, 16, 24]]


def _register_r_helpers(robjects) -> None:
    robjects.r("""
        suppressMessages(library(SpiecEasi))

        run_sparseicov_fixed <- function(X, lambda_value, sym_rule) {
          sparse_fun <- get('sparseiCov', envir = asNamespace('SpiecEasi'))
          fit <- sparse_fun(
            X,
            method = 'mb',
            lambda = c(as.numeric(lambda_value)),
            sym = as.character(sym_rule),
            verbose = FALSE,
            cov.output = TRUE
          )

          adj <- as.matrix(fit$path[[1]])
          beta <- as.matrix(fit$beta[[1]])
          beta_sym <- as.matrix(symBeta(beta, mode = 'ave')) * (adj != 0)

          list(
            adjacency = adj,
            beta = beta,
            sym_weights = beta_sym
          )
        }

        run_spieceasi_norm <- function(count_matrix) {
          norm_fun <- get('.spiec.easi.norm', envir = asNamespace('SpiecEasi'))
          norm_fun(as.matrix(count_matrix))
        }
        """)


def _run_r_sparseicov_fixed(
    r_runtime,
    X: np.ndarray,
    *,
    alpha: float,
    sym_rule: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    robjects = r_runtime["robjects"]
    numpy2ri = r_runtime["numpy2ri"]
    localconverter = r_runtime["localconverter"]

    run_fn = robjects.globalenv["run_sparseicov_fixed"]

    with localconverter(robjects.default_converter + numpy2ri.converter):
        r_X = robjects.conversion.py2rpy(X)

    out = run_fn(r_X, float(alpha), sym_rule)
    r_adj = (np.asarray(out.rx2("adjacency"), dtype=float) != 0.0).astype(int)
    # huge stores beta with predictors on rows and responses on columns.
    # Python stores directed coefficients as [response, predictor].
    r_beta = np.asarray(out.rx2("beta"), dtype=float).T
    r_sym = np.asarray(out.rx2("sym_weights"), dtype=float)
    return r_adj, r_beta, r_sym


def _directed_presence(beta: np.ndarray, i: int, j: int) -> tuple[bool, bool]:
    return (
        bool(abs(float(beta[i, j])) > SUPPORT_ZERO_THRESHOLD),
        bool(abs(float(beta[j, i])) > SUPPORT_ZERO_THRESHOLD),
    )


def _edge_present_from_directed(forward: bool, reverse: bool, sym_rule: str) -> bool:
    if sym_rule == "or":
        return forward or reverse
    return forward and reverse


def _discordance_origin(
    py_beta: np.ndarray,
    r_beta: np.ndarray,
    i: int,
    j: int,
    *,
    sym_rule: str,
) -> str:
    py_fwd, py_rev = _directed_presence(py_beta, i, j)
    r_fwd, r_rev = _directed_presence(r_beta, i, j)

    py_edge = _edge_present_from_directed(py_fwd, py_rev, sym_rule)
    r_edge = _edge_present_from_directed(r_fwd, r_rev, sym_rule)

    if py_edge == r_edge:
        return "none"
    if (py_fwd != r_fwd) and (py_rev == r_rev):
        return "i->j"
    if (py_fwd == r_fwd) and (py_rev != r_rev):
        return "j->i"
    return "rule_or_both_directions"


@pytest.mark.r_validation
@pytest.mark.slow
def test_phase2_input_clr_and_column_order_match_r(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    numpy2ri = r_runtime["numpy2ri"]
    localconverter = r_runtime["localconverter"]

    _register_r_helpers(robjects)
    run_spieceasi_norm = robjects.globalenv["run_spieceasi_norm"]

    for cfg in _dataset_cfgs():
        counts_df, _ = generate_compositional_counts(**cfg)
        py_clr = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr

        with localconverter(robjects.default_converter + numpy2ri.converter):
            r_counts = robjects.conversion.py2rpy(counts_df.to_numpy(dtype=float))
        r_clr = np.asarray(run_spieceasi_norm(r_counts), dtype=float)

        assert py_clr.shape == r_clr.shape
        assert np.array_equal(np.arange(py_clr.shape[1]), np.arange(r_clr.shape[1]))
        assert np.allclose(py_clr, r_clr, atol=1e-12, rtol=1e-12)


@pytest.mark.r_validation
@pytest.mark.slow
def test_phase2_directed_coefficients_match_r_sparseicov_mb(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    _register_r_helpers(robjects)

    mismatches: list[str] = []

    for cfg in _dataset_cfgs():
        counts_df, _ = generate_compositional_counts(**cfg)
        X = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr
        for alpha in _lambda_values(X):
            py = mb_neighborhood_selection(
                X,
                alpha=float(alpha),
                sym_rule="or",
                standardize=True,
                max_iter=200_000,
                tol=1e-12,
                random_state=0,
                coefficient_scale="standardized",
            )
            _, r_beta, _ = _run_r_sparseicov_fixed(
                r_runtime,
                X,
                alpha=float(alpha),
                sym_rule="or",
            )

            if not np.allclose(
                py.directed_coefficients,
                r_beta,
                rtol=FIXED_LAMBDA_WEIGHT_TOLERANCE.rtol,
                atol=FIXED_LAMBDA_WEIGHT_TOLERANCE.atol,
            ):
                max_abs = float(np.max(np.abs(py.directed_coefficients - r_beta)))
                mismatches.append(
                    "directed_beta_mismatch "
                    f"n={cfg['n_samples']} p={cfg['n_taxa']} seed={cfg['seed']} "
                    f"alpha={float(alpha):.6g} max_abs={max_abs:.6g}"
                )

    assert not mismatches, "\n".join(mismatches)


@pytest.mark.r_validation
@pytest.mark.slow
@pytest.mark.parametrize("sym_rule", ["or", "and"])
def test_mb_fixed_lambda_matches_r_sparseicov(r_runtime, sym_rule: str) -> None:
    robjects = r_runtime["robjects"]
    _register_r_helpers(robjects)

    mismatch_records: list[str] = []

    for cfg in _dataset_cfgs():
        counts_df, _ = generate_compositional_counts(**cfg)
        X = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr

        for alpha in _lambda_values(X):
            py = mb_neighborhood_selection(
                X,
                alpha=float(alpha),
                sym_rule=sym_rule,
                standardize=True,
                max_iter=200_000,
                tol=1e-12,
                random_state=0,
                coefficient_scale="standardized",
            )

            r_adj, r_beta, r_weights = _run_r_sparseicov_fixed(
                r_runtime,
                X,
                alpha=float(alpha),
                sym_rule=sym_rule,
            )

            py_adj = py.adjacency
            py_beta = py.directed_coefficients
            py_weights = py.weight_matrix

            if not np.array_equal(py_adj, r_adj):
                upper = np.triu_indices(py_adj.shape[0], k=1)
                discord = np.logical_xor(py_adj[upper] != 0, r_adj[upper] != 0)
                discord_rows = upper[0][discord]
                discord_cols = upper[1][discord]

                hard_discord: list[str] = []
                for i, j in zip(discord_rows, discord_cols, strict=True):
                    py_abs = abs(float(py_weights[i, j]))
                    r_abs = abs(float(r_weights[i, j]))
                    py_ij = float(py_beta[i, j])
                    py_ji = float(py_beta[j, i])
                    r_ij = float(r_beta[i, j])
                    r_ji = float(r_beta[j, i])
                    origin = _discordance_origin(
                        py_beta,
                        r_beta,
                        int(i),
                        int(j),
                        sym_rule=sym_rule,
                    )

                    if max(py_abs, r_abs) > FIXED_LAMBDA_AMBIGUOUS_EDGE_WEIGHT_MAX:
                        hard_discord.append(
                            f"({int(i)},{int(j)}): |py|={py_abs:.6g}, |r|={r_abs:.6g}, "
                            f"dist0_py={py_abs - SUPPORT_ZERO_THRESHOLD:.6g}, "
                            f"dist0_r={r_abs - SUPPORT_ZERO_THRESHOLD:.6g}, "
                            f"py_ij={py_ij:.6g}, py_ji={py_ji:.6g}, "
                            f"r_ij={r_ij:.6g}, r_ji={r_ji:.6g}, cause={origin}"
                        )

                if hard_discord:
                    mismatch_records.append(
                        "hard_adjacency_discord "
                        f"sym={sym_rule} n={cfg['n_samples']} p={cfg['n_taxa']} "
                        f"seed={cfg['seed']} alpha={float(alpha):.6g} "
                        + "; ".join(hard_discord)
                    )

            shared = np.logical_and(py_adj != 0, r_adj != 0)
            shared_upper = np.triu_indices(shared.shape[0], k=1)
            shared_mask = shared[shared_upper]

            if np.any(shared_mask):
                py_shared = py_weights[shared_upper][shared_mask]
                r_shared = r_weights[shared_upper][shared_mask]
                if not np.allclose(
                    py_shared,
                    r_shared,
                    rtol=FIXED_LAMBDA_WEIGHT_TOLERANCE.rtol,
                    atol=FIXED_LAMBDA_WEIGHT_TOLERANCE.atol,
                ):
                    max_abs = float(np.max(np.abs(py_shared - r_shared)))
                    mismatch_records.append(
                        "weight_mismatch "
                        f"sym={sym_rule} n={cfg['n_samples']} p={cfg['n_taxa']} "
                        f"seed={cfg['seed']} alpha={float(alpha):.6g} "
                        f"shared_edges={int(np.sum(shared_mask))} max_abs={max_abs:.6g}"
                    )

    assert not mismatch_records, "\n".join(mismatch_records)
