from __future__ import annotations

import numpy as np
import pytest

from spieceasi.neighborhood import make_alpha_path
from spieceasi.pipeline import spiec_easi_mb
from spieceasi.transforms import clr_transform
from tests.synthetic_data import generate_strong_signal_counts
from tests.validation_config import (
    FIXED_LAMBDA_AMBIGUOUS_EDGE_WEIGHT_MAX,
    FIXED_LAMBDA_WEIGHT_TOLERANCE,
    MIN_SHARED_SUBSAMPLE_SCENARIOS,
    STARS_SHARED_SUBSAMPLES_TOLERANCE,
    SUPPORT_ZERO_THRESHOLD,
)


def _scenario_cfgs() -> list[dict[str, int | float]]:
    cfgs = [
        dict(n_samples=170, n_taxa=16, depth=4500, rho=0.35, seed=7),
        dict(n_samples=190, n_taxa=18, depth=5200, rho=0.40, seed=13),
        dict(n_samples=220, n_taxa=24, depth=6000, rho=0.45, seed=31),
        dict(n_samples=240, n_taxa=24, depth=6500, rho=0.42, seed=41),
        dict(n_samples=260, n_taxa=28, depth=7000, rho=0.45, seed=17),
    ]
    assert len(cfgs) >= MIN_SHARED_SUBSAMPLE_SCENARIOS
    return cfgs


def _default_subsample_size(n_samples: int) -> int:
    if n_samples > 144:
        return int(np.floor(10 * np.sqrt(n_samples)))
    return int(np.floor(0.8 * n_samples))


def _selected_index_from_monotone(monotone: np.ndarray, threshold: float) -> int:
    feasible = np.where(monotone <= threshold)[0]
    return int(feasible[-1]) if feasible.size > 0 else 0


@pytest.mark.r_validation
@pytest.mark.slow
def test_pipeline_shared_subsamples_matches_r(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    pandas2ri = r_runtime["pandas2ri"]
    numpy2ri = r_runtime["numpy2ri"]
    localconverter = r_runtime["localconverter"]

    robjects.r("""
        suppressMessages(library(SpiecEasi))
        suppressMessages(library(pulsar))

        run_pipeline_shared <- function(
            df, lambda_path, subsample_idx_1b, threshold, sym_rule
        ) {
          norm_fun <- get('.spiec.easi.norm', envir = asNamespace('SpiecEasi'))
          sparse_fun <- get('sparseiCov', envir = asNamespace('SpiecEasi'))
          stars_fun <- get('stars.stability', envir = asNamespace('pulsar'))

          X <- norm_fun(as.matrix(df))
          rep_num <- nrow(subsample_idx_1b)

          premerge <- vector('list', rep_num)
          for (r in seq_len(rep_num)) {
            idx <- as.integer(subsample_idx_1b[r, ])
            fit <- sparse_fun(
              X[idx, , drop = FALSE],
              method = 'mb',
              lambda = as.numeric(lambda_path),
              sym = as.character(sym_rule),
              verbose = FALSE,
              cov.output = TRUE
            )
            premerge[[r]] <- fit$path
          }

          premerge_reord <- lapply(seq_along(lambda_path), function(i) {
            lapply(seq_len(rep_num), function(j) premerge[[j]][[i]])
          })

          est <- stars_fun(premerge_reord, threshold, rep_num, ncol(X))
          opt_index_1b <- as.integer(est$opt.index)
          opt_lambda <- as.numeric(lambda_path[opt_index_1b])

          fit_full <- sparse_fun(
            X,
            method = 'mb',
            lambda = c(opt_lambda),
            sym = as.character(sym_rule),
            verbose = FALSE,
            cov.output = TRUE
          )

          adj <- as.matrix(fit_full$path[[1]])
          beta <- as.matrix(fit_full$beta[[1]])
          sym_weights <- as.matrix(symBeta(beta, mode = 'ave')) * (adj != 0)

          list(
            summary = as.numeric(est$summary),
            opt_index_1b = opt_index_1b,
            opt_lambda = opt_lambda,
            adjacency = adj,
            beta = beta,
            weights = sym_weights
          )
        }
        """)
    run_pipeline_shared = robjects.globalenv["run_pipeline_shared"]

    stars_threshold = 0.05
    n_subsamples = 16

    mismatch_records: list[str] = []

    for sym_rule in ("or", "and"):
        for cfg in _scenario_cfgs():
            counts_df, _ = generate_strong_signal_counts(**cfg)
            X = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr
            alpha_path = make_alpha_path(X, n_lambdas=20, lambda_min_ratio=1e-2)

            subsample_size = _default_subsample_size(X.shape[0])
            rng = np.random.default_rng(int(cfg["seed"]) + 101)
            subsamples_0b = [
                rng.choice(X.shape[0], size=subsample_size, replace=False)
                for _ in range(n_subsamples)
            ]

            py_result = spiec_easi_mb(
                counts_df,
                alpha_path=alpha_path,
                n_subsamples=n_subsamples,
                subsamples=subsamples_0b,
                stars_threshold=stars_threshold,
                sym_rule=sym_rule,
                standardize=True,
                max_iter=200_000,
                tol=1e-12,
                random_state=0,
                n_jobs=1,
            )

            subsample_idx_1b = np.vstack([idx + 1 for idx in subsamples_0b]).astype(int)
            with localconverter(robjects.default_converter + pandas2ri.converter):
                r_counts = robjects.conversion.py2rpy(counts_df)
            with localconverter(robjects.default_converter + numpy2ri.converter):
                r_lambda = robjects.conversion.py2rpy(alpha_path)
                r_subsamples = robjects.conversion.py2rpy(subsample_idx_1b)

            r_out = run_pipeline_shared(
                r_counts,
                r_lambda,
                r_subsamples,
                stars_threshold,
                sym_rule,
            )

            r_monotone = np.asarray(r_out.rx2("summary"), dtype=float)
            r_opt_0b = int(np.asarray(r_out.rx2("opt_index_1b"), dtype=int)[0]) - 1
            r_opt_lambda = float(np.asarray(r_out.rx2("opt_lambda"), dtype=float)[0])
            r_adj = (np.asarray(r_out.rx2("adjacency"), dtype=float) != 0.0).astype(int)
            r_beta = np.asarray(r_out.rx2("beta"), dtype=float).T
            r_weights = np.asarray(r_out.rx2("weights"), dtype=float)

            py_opt_0b = _selected_index_from_monotone(
                py_result.monotone_instability,
                stars_threshold,
            )

            if not np.allclose(
                py_result.monotone_instability,
                r_monotone,
                rtol=STARS_SHARED_SUBSAMPLES_TOLERANCE.rtol,
                atol=STARS_SHARED_SUBSAMPLES_TOLERANCE.atol,
            ):
                mismatch_records.append(
                    "stars_curve_mismatch "
                    f"sym={sym_rule} n={cfg['n_samples']} "
                    f"p={cfg['n_taxa']} seed={cfg['seed']}"
                )

            if py_opt_0b != r_opt_0b:
                mismatch_records.append(
                    "opt_index_mismatch "
                    f"sym={sym_rule} n={cfg['n_samples']} p={cfg['n_taxa']} "
                    f"seed={cfg['seed']} py={py_opt_0b} r={r_opt_0b}"
                )

            if not np.isclose(
                py_result.selected_lambda, r_opt_lambda, rtol=0.0, atol=1e-12
            ):
                mismatch_records.append(
                    "selected_lambda_mismatch "
                    f"sym={sym_rule} n={cfg['n_samples']} p={cfg['n_taxa']} "
                    f"seed={cfg['seed']} py={py_result.selected_lambda:.12g} "
                    f"r={r_opt_lambda:.12g}"
                )

            if not np.array_equal(py_result.adjacency, r_adj):
                upper = np.triu_indices(py_result.adjacency.shape[0], k=1)
                discord = np.logical_xor(
                    py_result.adjacency[upper] != 0,
                    r_adj[upper] != 0,
                )
                rows = upper[0][discord]
                cols = upper[1][discord]

                details: list[str] = []
                for i, j in zip(rows, cols, strict=True):
                    py_w = float(py_result.weights[i, j])
                    r_w = float(r_weights[i, j])
                    py_abs = abs(py_w)
                    r_abs = abs(r_w)
                    py_ij = float(py_result.directed_coefficients[i, j])
                    py_ji = float(py_result.directed_coefficients[j, i])
                    r_ij = float(r_beta[i, j])
                    r_ji = float(r_beta[j, i])
                    if max(py_abs, r_abs) > FIXED_LAMBDA_AMBIGUOUS_EDGE_WEIGHT_MAX:
                        details.append(
                            f"({int(i)},{int(j)}): |py|={py_abs:.6g}, |r|={r_abs:.6g}, "
                            f"dist0_py={py_abs - SUPPORT_ZERO_THRESHOLD:.6g}, "
                            f"dist0_r={r_abs - SUPPORT_ZERO_THRESHOLD:.6g}, "
                            f"py_ij={py_ij:.6g}, py_ji={py_ji:.6g}, "
                            f"r_ij={r_ij:.6g}, r_ji={r_ji:.6g}"
                        )
                if details:
                    mismatch_records.append(
                        "adjacency_mismatch "
                        f"sym={sym_rule} n={cfg['n_samples']} p={cfg['n_taxa']} "
                        f"seed={cfg['seed']} " + "; ".join(details)
                    )

            shared = np.logical_and(py_result.adjacency != 0, r_adj != 0)
            shared_upper = np.triu_indices(shared.shape[0], k=1)
            shared_mask = shared[shared_upper]
            if np.any(shared_mask):
                py_shared = py_result.weights[shared_upper][shared_mask]
                r_shared = r_weights[shared_upper][shared_mask]
                if not np.allclose(
                    py_shared,
                    r_shared,
                    rtol=FIXED_LAMBDA_WEIGHT_TOLERANCE.rtol,
                    atol=FIXED_LAMBDA_WEIGHT_TOLERANCE.atol,
                ):
                    max_abs = float(np.max(np.abs(py_shared - r_shared)))
                    mismatch_records.append(
                        "weights_mismatch "
                        f"sym={sym_rule} n={cfg['n_samples']} p={cfg['n_taxa']} "
                        f"seed={cfg['seed']} max_abs={max_abs:.6g}"
                    )

    assert not mismatch_records, "\n".join(mismatch_records)
