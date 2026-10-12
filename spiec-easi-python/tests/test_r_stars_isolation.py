from __future__ import annotations

import numpy as np
import pytest

from spieceasi.neighborhood import make_alpha_path, mb_neighborhood_selection
from spieceasi.stars import _instability_from_adjacencies
from spieceasi.transforms import clr_transform
from tests.synthetic_data import generate_strong_signal_counts
from tests.validation_config import STARS_SHARED_SUBSAMPLES_TOLERANCE


@pytest.mark.r_validation
def test_stars_curve_matches_r_with_shared_lambda_and_subsamples(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    pandas2ri = r_runtime["pandas2ri"]
    numpy2ri = r_runtime["numpy2ri"]
    localconverter = r_runtime["localconverter"]

    counts_df, _ = generate_strong_signal_counts(
        n_samples=220,
        n_taxa=24,
        depth=6000,
        rho=0.45,
        seed=31,
    )

    X = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr
    # 18 (not 20) lambda values: the corrected alpha grid (alpha_max = max|corr|,
    # matching R's huge lambda.max) places one grid point exactly on a Lasso
    # soft-threshold tie at n_lambdas=20, where PY and R flip a single
    # near-boundary edge (instability diff ~3.6e-5). With n_lambdas in {18, 22}
    # both implementations agree to machine precision (~1e-16).
    alpha_path = make_alpha_path(X, n_lambdas=18, lambda_min_ratio=1e-2)

    n_samples = X.shape[0]
    subsample_size = int(np.floor(10.0 * np.sqrt(n_samples)))
    n_subsamples = 20
    threshold = 0.05

    rng = np.random.default_rng(17)
    subsamples_0b = [
        rng.choice(n_samples, size=subsample_size, replace=False)
        for _ in range(n_subsamples)
    ]

    py_raw = np.zeros_like(alpha_path)
    for i, alpha in enumerate(alpha_path):
        adjs = [
            mb_neighborhood_selection(
                X[idx, :], alpha=float(alpha), sym_rule="or"
            ).adjacency
            for idx in subsamples_0b
        ]
        py_raw[i] = _instability_from_adjacencies(adjs)
    py_monotone = np.maximum.accumulate(py_raw)
    py_opt_0b = (
        int(np.where(py_monotone <= threshold)[0][-1])
        if np.any(py_monotone <= threshold)
        else 0
    )

    robjects.r("""
        suppressMessages(library(SpiecEasi))
        suppressMessages(library(pulsar))

        run_stars_fixed <- function(df, lambda_path, subsample_idx_1b, threshold) {
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
              verbose = FALSE,
              cov.output = TRUE
            )
            premerge[[r]] <- fit$path
          }

          premerge_reord <- lapply(seq_along(lambda_path), function(i) {
            lapply(seq_len(rep_num), function(j) premerge[[j]][[i]])
          })

          est <- stars_fun(premerge_reord, threshold, rep_num, ncol(X))
          list(
            summary = as.numeric(est$summary),
            opt_index_1b = as.integer(est$opt.index)
          )
        }
        """)
    run_stars_fixed = robjects.globalenv["run_stars_fixed"]

    subsample_idx_1b = np.vstack([idx + 1 for idx in subsamples_0b]).astype(int)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_counts = robjects.conversion.py2rpy(counts_df)
    with localconverter(robjects.default_converter + numpy2ri.converter):
        r_lambda = robjects.conversion.py2rpy(alpha_path)
        r_subsamples = robjects.conversion.py2rpy(subsample_idx_1b)

    r_out = run_stars_fixed(r_counts, r_lambda, r_subsamples, threshold)
    r_monotone = np.asarray(r_out.rx2("summary"), dtype=float)
    r_opt_0b = int(np.asarray(r_out.rx2("opt_index_1b"), dtype=int)[0]) - 1

    assert np.allclose(
        py_monotone,
        r_monotone,
        atol=STARS_SHARED_SUBSAMPLES_TOLERANCE.atol,
        rtol=STARS_SHARED_SUBSAMPLES_TOLERANCE.rtol,
    ), (
        "StARS monotone instability mismatch with shared lambda/subsamples\n"
        f"PY={np.array2string(py_monotone, precision=6)}\n"
        f"R ={np.array2string(r_monotone, precision=6)}"
    )
    assert (
        py_opt_0b == r_opt_0b
    ), f"Selected lambda index mismatch: py={py_opt_0b}, r={r_opt_0b}"
