from __future__ import annotations

import numpy as np
import pytest

from spieceasi.neighborhood import _fit_single_neighborhood
from spieceasi.transforms import clr_transform
from tests.synthetic_data import generate_compositional_counts


def _r_style_lambda_path(
    X: np.ndarray, nlambda: int, lambda_min_ratio: float
) -> np.ndarray:
    corr = np.corrcoef(X, rowvar=False)
    upper = np.triu_indices(corr.shape[0], k=1)
    lambda_max = float(np.max(np.abs(corr[upper])))
    lambda_min = lambda_max * lambda_min_ratio
    return np.exp(np.linspace(np.log(lambda_max), np.log(lambda_min), nlambda))


@pytest.mark.r_validation
def test_clr_matches_r_spieceasi_norm_transform(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    pandas2ri = r_runtime["pandas2ri"]
    localconverter = r_runtime["localconverter"]

    counts_df, _ = generate_compositional_counts(n_samples=70, n_taxa=7, seed=13)
    py_clr = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr

    robjects.r("""
        suppressMessages(library(SpiecEasi))
        run_norm <- function(df) {
          X <- as.matrix(df)
          norm_fun <- get(".spiec.easi.norm", envir = asNamespace("SpiecEasi"))
          norm_fun(X)
        }
        """)
    run_norm = robjects.globalenv["run_norm"]

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_counts = robjects.conversion.py2rpy(counts_df)
    r_clr = np.asarray(run_norm(r_counts), dtype=float)

    assert r_clr.shape == py_clr.shape
    assert np.allclose(py_clr, r_clr, atol=1e-10, rtol=1e-10)


@pytest.mark.r_validation
def test_lambda_grid_matches_r_getlampath_convention(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    numpy2ri = r_runtime["numpy2ri"]
    localconverter = r_runtime["localconverter"]

    counts_df, _ = generate_compositional_counts(n_samples=90, n_taxa=8, seed=21)
    X = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr

    nlambda = 20
    lambda_min_ratio = 1e-2
    py_path = _r_style_lambda_path(
        X, nlambda=nlambda, lambda_min_ratio=lambda_min_ratio
    )

    robjects.r("""
        suppressMessages(library(SpiecEasi))
        run_lambda_path <- function(X, nlambda, lambda_min_ratio) {
          lam_max <- getMaxCov(cor(X))
          getLamPath(lam_max, lam_max * lambda_min_ratio, nlambda, log = TRUE)
        }
        """)
    run_lambda_path = robjects.globalenv["run_lambda_path"]

    with localconverter(robjects.default_converter + numpy2ri.converter):
        r_X = robjects.conversion.py2rpy(X)
    r_path = np.asarray(run_lambda_path(r_X, nlambda, lambda_min_ratio), dtype=float)

    assert np.all(np.diff(py_path) < 0)
    assert np.allclose(py_path, r_path, atol=1e-10, rtol=1e-10)


@pytest.mark.r_validation
def test_single_neighborhood_matches_glmnet_support_pattern(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    numpy2ri = r_runtime["numpy2ri"]
    localconverter = r_runtime["localconverter"]

    counts_df, _ = generate_compositional_counts(n_samples=150, n_taxa=8, seed=8)
    X = clr_transform(counts_df, pseudocount=1.0, zero_method="additive").clr

    target_idx = 3
    alpha = 0.08

    py_row = _fit_single_neighborhood(
        X,
        target_idx,
        alpha=alpha,
        standardize=True,
        max_iter=20_000,
        tol=1e-7,
        random_state=0,
    )

    robjects.r("""
        suppressMessages(library(glmnet))
        run_glmnet_row <- function(X, target_idx_1based, lambda_value) {
          j <- as.integer(target_idx_1based)
          y <- X[, j]
          X_other <- X[, -j, drop = FALSE]
          fit <- glmnet(
            X_other,
            y,
            family = "gaussian",
            alpha = 1,
            standardize = TRUE,
            intercept = TRUE,
            lambda = c(lambda_value)
          )
          coef_vec <- as.numeric(coef(fit, s = lambda_value))[-1]
          out <- rep(0, ncol(X))
          out[-j] <- coef_vec
          out
        }
        """)
    run_glmnet_row = robjects.globalenv["run_glmnet_row"]

    with localconverter(robjects.default_converter + numpy2ri.converter):
        r_X = robjects.conversion.py2rpy(X)
    r_row = np.asarray(run_glmnet_row(r_X, target_idx + 1, alpha), dtype=float)

    py_support = np.abs(py_row) > 1e-10
    r_support = np.abs(r_row) > 1e-10
    overlap = np.sum(py_support & r_support)
    union = np.sum(py_support | r_support)
    jaccard = 1.0 if union == 0 else overlap / union

    # Solvers differ, but the selected predictor support should be similar.
    assert jaccard >= 0.7


@pytest.mark.r_validation
def test_weight_symmetrization_matches_r_symbeta_ave(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    numpy2ri = r_runtime["numpy2ri"]
    localconverter = r_runtime["localconverter"]

    rng = np.random.default_rng(10)
    beta = rng.normal(size=(6, 6))
    np.fill_diagonal(beta, 0.0)

    py_sym = 0.5 * (beta + beta.T)

    robjects.r("""
        suppressMessages(library(SpiecEasi))
        run_symbeta <- function(beta) {
          as.matrix(symBeta(beta, mode = "ave"))
        }
        """)
    run_symbeta = robjects.globalenv["run_symbeta"]

    with localconverter(robjects.default_converter + numpy2ri.converter):
        r_beta = robjects.conversion.py2rpy(beta)
    r_sym = np.asarray(run_symbeta(r_beta), dtype=float)

    assert np.allclose(py_sym, r_sym, atol=1e-12, rtol=1e-12)


@pytest.mark.r_validation
def test_optbeta_symbeta_equals_manual_average(r_runtime) -> None:
    robjects = r_runtime["robjects"]
    pandas2ri = r_runtime["pandas2ri"]
    localconverter = r_runtime["localconverter"]

    counts_df, _ = generate_compositional_counts(n_samples=180, n_taxa=12, seed=22)

    robjects.r("""
                suppressMessages(library(SpiecEasi))
                run_optbeta <- function(df) {
                    fit <- spiec.easi(
                        as.matrix(df),
                        method = "mb",
                        sel.criterion = "stars",
                        nlambda = 20,
                        lambda.min.ratio = 1e-2,
                        pulsar.params = list(rep.num = 20, thresh = 0.05, seed = 22)
                    )
                    beta_opt <- as.matrix(getOptBeta(fit))
                    sym <- as.matrix(symBeta(beta_opt, mode = "ave"))
                    list(beta = beta_opt, sym = sym)
                }
                """)
    run_optbeta = robjects.globalenv["run_optbeta"]

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_counts = robjects.conversion.py2rpy(counts_df)

    out = run_optbeta(r_counts)
    beta = np.asarray(out.rx2("beta"), dtype=float)
    r_sym = np.asarray(out.rx2("sym"), dtype=float)
    manual_sym = 0.5 * (beta + beta.T)

    assert np.allclose(manual_sym, r_sym, atol=1e-10, rtol=1e-10)
