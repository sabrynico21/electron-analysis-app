"""Core numerical routines — faithful translation of ``CCLasso-master/R/cclasso.R`` (v2.0).

Every internal function carries a ``# R: cclasso.R:<line>`` comment pointing to the
source it was translated from. Parameter names, defaults, order of operations,
stopping criteria and the RNG-injection points mirror the R code exactly.

Read ``docs/PORTING_NOTES.md`` before modifying anything here.
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = ["cclasso", "cclasso_sub", "cv_loss_cclasso", "boot_cclasso"]


# --------------------------------------------------------------------------- #
# Numeric primitives (R -> numpy mappings; see docs/PORTING_NOTES.md §2)
# --------------------------------------------------------------------------- #
def _var(x: np.ndarray) -> np.ndarray:
    """``stats::var(x)`` for a matrix: column covariance with n-1 denominator."""
    return np.cov(x, rowvar=False, ddof=1)


def _cor(x: np.ndarray) -> np.ndarray:
    """``cor(x)``: Pearson correlation, n-1 denominator."""
    return np.corrcoef(x, rowvar=False)


def _eigen_desc(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``eigen(a)`` for a symmetric matrix: eigenvalues DESCENDING (R order).

    ``np.linalg.eigh`` returns ascending eigenvalues; we reverse so that
    ``values[i]`` pairs with ``vectors[:, i]`` exactly as in R.
    """
    values, vectors = np.linalg.eigh(a)
    return values[::-1], vectors[:, ::-1]


def _lower_tri_colmajor(p: int) -> tuple[np.ndarray, np.ndarray]:
    """Indices of strictly-lower-triangle elements in R's column-major order.

    R's ``M[lower.tri(M)]`` walks columns left-to-right, rows top-to-bottom:
    ``(2,1),(3,1),...,(p,1),(3,2),...``. numpy's ``tril_indices`` is row-major,
    so we build the index pair explicitly.  See docs/DIFFERENCES.md (D7).
    """
    rows = np.concatenate([np.arange(j + 1, p) for j in range(p)])
    cols = np.concatenate([np.full(p - j - 1, j) for j in range(p)])
    return rows, cols


def _soft_threshold(a: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """``(A > l)*(A-l) + (A < -l)*(A+l)`` — strict inequalities, as in R."""
    return (a > lam) * (a - lam) + (a < -lam) * (a + lam)


# --------------------------------------------------------------------------- #
# cclasso_sub  --  R: cclasso.R:236-278
# --------------------------------------------------------------------------- #
def cclasso_sub(
    sigma: np.ndarray,
    vx: np.ndarray,
    lambda2: float,
    wd: np.ndarray,
    u_f: np.ndarray,
    u0_wd: np.ndarray,
    d0_wd: np.ndarray,
    k_max: int = 200,
    x_tol: float = 1e-4,
    relative_errors: list | None = None,
) -> tuple[np.ndarray, int]:
    """One augmented-Lagrangian (ADMM) solve for a single lambda.

    Returns ``(sigma, k)`` where ``k`` is the number of iterations performed
    (R's ``cclasso_sub`` returns only ``sigma``; ``k`` is captured for diagnostics).
    ``wd`` is accepted for signature parity but unused in v2.0 (as in the R source).
    """
    p = sigma.shape[1]
    sigma2 = sigma.copy()
    lamb = np.zeros((p, p), dtype=float)
    lambda2_mat = np.full((p, p), lambda2, dtype=float)
    np.fill_diagonal(lambda2_mat, 0.0)

    k = 0
    err = 1.0
    while True:
        # R evaluates ``while (err > x_tol && k < k_max)``.  Preserve its
        # missing-logical failure: Python would otherwise treat ``NaN > tol``
        # as False and silently return the last iterate. Test the iteration cap
        # first because R's ``NA && FALSE`` is FALSE (not an error).
        if k >= k_max:
            break
        if np.isnan(err):
            raise ValueError("missing value where TRUE/FALSE needed")
        if not err > x_tol:
            break

        # update sigma
        x_sigma = u_f.T @ ((sigma2 - vx) - lamb) @ u_f
        block = x_sigma[:-1, :-1]
        x_sigma[:-1, :-1] = u0_wd @ ((u0_wd.T @ block @ u0_wd) * d0_wd) @ u0_wd.T
        sigma_new = vx + u_f @ x_sigma @ u_f.T
        # update sigma2
        a = lamb + sigma_new
        sigma2_new = _soft_threshold(a, lambda2_mat)
        # update Lambda
        lamb = lamb + (sigma_new - sigma2_new)
        # relative error (max over all entries)
        err_sigma = np.max(np.abs(sigma_new - sigma) / (np.abs(sigma) + 1.0))
        err_sigma2 = np.max(np.abs(sigma2_new - sigma2) / (np.abs(sigma2) + 1.0))
        # R's max() propagates NA from either argument. Python's built-in max()
        # can hide a NaN depending on argument order, so propagate it explicitly.
        err = np.nan if np.isnan(err_sigma) or np.isnan(err_sigma2) else max(err_sigma, err_sigma2)
        k += 1
        sigma = sigma_new
        sigma2 = sigma2_new

    if k >= k_max:
        warnings.warn(
            f"WARNING of cclasso_sub: Maximum Iteration: {k_max} && Relative error: {err} !",
            RuntimeWarning,
            stacklevel=2,
        )

    # Diagnostic-only capture for R/Python convergence-status comparisons.
    # Do not alter the iteration loop or the returned estimate.
    if relative_errors is not None:
        relative_errors.append(float(err))

    return sigma, k


# --------------------------------------------------------------------------- #
# cv_loss_cclasso  --  R: cclasso.R:211-233
# --------------------------------------------------------------------------- #
def cv_loss_cclasso(
    lambda2: float,
    x: np.ndarray,
    k_cv: int,
    sigma: np.ndarray,
    wd: np.ndarray,
    u_f: np.ndarray,
    u0_wd: np.ndarray,
    d0_wd: np.ndarray,
    wd2: np.ndarray,
    iters: list | None = None,
    k_max: int = 200,
    relative_errors: list | None = None,
) -> tuple[float, np.ndarray]:
    """Cross-validation loss for a single lambda (R: deterministic contiguous folds).

    Folds are ``(n_b*(k-1)+1):(n_b*k)`` with ``n_b = floor(n/k_cv)``; there is no
    randomness.  ``sigma`` is warm-started across folds; the returned ``sigma`` is the
    last fold's result, exactly as in R.
    """
    n, p = x.shape
    n_b = n // k_cv
    cv_loss = 0.0
    for kk in range(1, k_cv + 1):
        itest = np.arange(n_b * (kk - 1), n_b * kk)
        vxk = _var(x[itest, :])
        vx2k = _var(np.delete(x, itest, axis=0))

        sigma, it = cclasso_sub(
            sigma,
            vx2k,
            lambda2,
            wd,
            u_f,
            u0_wd,
            d0_wd,
            k_max=k_max,
            relative_errors=relative_errors,
        )
        if iters is not None:
            iters.append(it)

        dsig = sigma - vxk
        tmp = dsig.mean(axis=1)
        dsig = dsig - tmp[:, None] - tmp[None, :] + tmp.mean()
        cv_loss += np.linalg.norm(wd2[:, None] * dsig, "fro") ** 2

    return cv_loss, sigma


# --------------------------------------------------------------------------- #
# boot_cclasso  --  R: cclasso.R:141-207
# --------------------------------------------------------------------------- #
def boot_cclasso(
    x: np.ndarray,
    sigma: np.ndarray,
    lambda2: float,
    n_boot: int,
    wd: np.ndarray,
    u_f: np.ndarray,
    u0_wd: np.ndarray,
    d0_wd: np.ndarray,
    boot_indices: np.ndarray,
    null_norm: np.ndarray,
    iters: list | None = None,
    k_max: int = 200,
    relative_errors: list | None = None,
) -> dict:
    """Bootstrap estimation + artificial-null p-values.

    ``boot_indices`` (n x n_boot, **0-based**) and ``null_norm`` (n x p, raw standard
    normal draws) replace R's ``sample()`` / ``rnorm()`` so results are reproducible
    across languages.  See docs/DIFFERENCES.md (D6).
    """
    n, p = x.shape
    rows, cols = _lower_tri_colmajor(p)

    cors_boot = np.zeros((rows.size, n_boot + 1))
    vars_boot = np.zeros((p, n_boot + 1))

    # Bootstrap procedure
    for k in range(n_boot):
        ind_samp = boot_indices[:, k]
        sigma2, it = cclasso_sub(
            sigma,
            _var(x[ind_samp, :]),
            lambda2,
            wd,
            u_f,
            u0_wd,
            d0_wd,
            k_max=k_max,
            relative_errors=relative_errors,
        )
        if iters is not None:
            iters.append(it)
        vars_boot[:, k] = np.diag(sigma2)
        inv_s = 1.0 / np.sqrt(vars_boot[:, k])
        cors_mat = inv_s[:, None] * sigma2 * inv_s[None, :]
        cors_boot[:, k] = cors_mat[rows, cols]

    inv_s = 1.0 / np.sqrt(np.diag(sigma))
    cors_mat = inv_s[:, None] * sigma * inv_s[None, :]
    cors_boot[:, n_boot] = cors_mat[rows, cols]
    vars_boot[:, n_boot] = np.diag(sigma)

    # Variance estimation via bootstrap
    vars2 = vars_boot.mean(axis=1)

    # Correlations' relationship for artificial null sample
    tol_cor = 1e-3
    sam_art0 = null_norm * np.sqrt(vars2)[None, :]
    cors_art0 = _cor(sam_art0)[rows, cols]
    sam_art = sam_art0 - np.log(np.exp(sam_art0).sum(axis=1))[:, None]
    sigma_art, it = cclasso_sub(
        sigma,
        _var(sam_art),
        lambda2,
        wd,
        u_f,
        u0_wd,
        d0_wd,
        k_max=k_max,
        relative_errors=relative_errors,
    )
    if iters is not None:
        iters.append(it)
    inv_s = 1.0 / np.sqrt(np.diag(sigma_art))
    cors_mat = inv_s[:, None] * sigma_art * inv_s[None, :]
    cors_art2 = cors_mat[rows, cols]

    cors0m2 = np.log(((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) * (1 - cors_art0)))
    tmp = np.abs(cors_art2) >= tol_cor
    # R's ``sum(tmp)`` propagates NA; its ``ifelse(NA, median(...), 0)`` then
    # propagates NA to the correction for every pair. NumPy comparisons against
    # NaN yield False, so ``tmp.any()`` alone would silently use zero bias and
    # return a partial set of finite pairs instead of R's all-nonfinite result.
    if np.isnan(cors_art2).any():
        bias02 = float("nan")
    else:
        bias02 = float(np.median(np.abs(cors0m2)[tmp])) if tmp.any() else 0.0

    # Modification of estimation for cclasso
    cors2 = np.log((1 + cors_boot) / (1 - cors_boot))
    cors2mod = (cors_boot >= tol_cor) * (cors2 + bias02) + (cors_boot <= -tol_cor) * (
        cors2 - bias02
    )
    cors2mod = 1 - (2 / (np.exp(cors2mod) + 1)).mean(axis=1)
    cors2_mat = _sym_from_lower(p, rows, cols, cors2mod)

    # P-values
    p_vals = _student_t_cdf(cors2mod * np.sqrt((n - 2) / (1 - cors2mod**2)), n - 2)
    p_vals = np.where(p_vals <= 0.5, p_vals, 1 - p_vals)
    pval_mat = _sym_from_lower(p, rows, cols, p_vals)

    return {
        "var_w": vars2,
        "cor_w": cors2_mat,
        "p_vals": pval_mat,
        "_bootstrap_variances": vars_boot,
    }


def _sym_from_lower(p: int, rows: np.ndarray, cols: np.ndarray, vals: np.ndarray):
    """Build a symmetric matrix from a column-major lower triangle + unit diagonal."""
    m = np.eye(p)
    m[rows, cols] = vals
    m[cols, rows] = vals
    return m


def _student_t_cdf(q: np.ndarray, df: float) -> np.ndarray:
    """``pt(q, df)`` — Student-t CDF (scipy)."""
    from scipy import stats

    return stats.t.cdf(q, df)


# --------------------------------------------------------------------------- #
# cclasso  --  R: cclasso.R:36-128
# --------------------------------------------------------------------------- #
def cclasso(
    x,
    counts: bool = False,
    pseudo: float = 0.5,
    k_cv: int = 3,
    lam_int=(1e-4, 1.0),
    k_max: int = 20,
    n_boot: int = 20,
    *,
    boot_indices: np.ndarray | None = None,
    null_norm: np.ndarray | None = None,
    seed: int | None = None,
) -> dict:
    """Correlation inference for compositional data through Lasso.

    Faithful translation of ``cclasso`` from the R package CCLasso v2.0.

    Parameters
    ----------
    x : array-like, shape (n, p)
        Compositional data; rows are samples, columns are variables.
    counts : bool, default False
        If True, ``x`` is treated as counts and ``(x + pseudo) / rowSums(x + pseudo)``
        is applied before the log transform.
    pseudo : float, default 0.5
        Pseudo-count, used only when ``counts=True``.
    k_cv : int, default 3
        Number of cross-validation folds (deterministic contiguous blocks).
    lam_int : (float, float), default (1e-4, 1.0)
        Tuning-parameter interval in **log10** scale (R's ``lam_int``).
    k_max : int, default 20
        Maximum golden-section iterations.
    n_boot : int, default 20
        Number of bootstrap replicates.
    boot_indices : ndarray (n, n_boot), optional
        Bootstrap row indices, **0-based**. Injected so results are reproducible.
    null_norm : ndarray (n, p), optional
        Raw standard-normal draws for the artificial null sample.
    seed : int, optional
        Used only when ``boot_indices``/``null_norm`` are not supplied.

    Returns
    -------
    dict
        ``{"var_w", "cor_w", "p_vals", "lambda", "info_cv"}`` matching the R output
        (plus, for diagnostics, nothing extra on stdout — warnings go to stderr).
    """
    from .validation import validate_input

    x = validate_input(
        x,
        counts=counts,
        pseudo=pseudo,
        k_cv=k_cv,
        lam_int=lam_int,
        n_boot=n_boot,
    )
    n, p = x.shape

    if counts:
        x = x + pseudo
        x = x / x.sum(axis=1, keepdims=True)
    x = np.log(x)
    vx2 = _var(x)

    # Diagonal weight for loss function
    rmean_vx2 = vx2.mean(axis=1)
    wd = 1.0 / (np.diag(vx2) - rmean_vx2 - rmean_vx2 + rmean_vx2.mean())
    wd2 = np.sqrt(wd)

    # Some global parameters for optimization with single lambda
    rho = 1.0
    _, u_f = _eigen_desc(np.eye(p) - 1.0 / p)
    wd_u = (u_f.T @ (wd[:, None] * u_f))[:-1, :-1]
    wd_u_eig_vals, wd_u_eig_vecs = _eigen_desc(wd_u)
    d0_wd = 1.0 / ((np.add.outer(wd_u_eig_vals, wd_u_eig_vals) / (2 * rho)) + 1.0)
    u0_wd = wd_u_eig_vecs

    # Golden section method for the selection of lambda (log10 scale)
    sigma = vx2
    lam_int2 = np.log10(np.asarray(lam_int, dtype=float))
    a1, b1 = float(lam_int2[0]), float(lam_int2[1])

    iters: list = []
    relative_errors: list = []
    lams: list = []
    fvals: list = []

    a2 = a1 + 0.382 * (b1 - a1)
    b2 = a1 + 0.618 * (b1 - a1)

    fb2_loss, fb2_sigma = cv_loss_cclasso(
        10**b2 / rho,
        x,
        k_cv,
        sigma,
        wd,
        u_f,
        u0_wd,
        d0_wd,
        wd2,
        iters,
        relative_errors=relative_errors,
    )
    lams.append(b2)
    fvals.append(fb2_loss)
    fa2_loss, fa2_sigma = cv_loss_cclasso(
        10**a2 / rho,
        x,
        k_cv,
        fb2_sigma,
        wd,
        u_f,
        u0_wd,
        d0_wd,
        wd2,
        iters,
        relative_errors=relative_errors,
    )
    lams.append(a2)
    fvals.append(fa2_loss)

    # Error tolerance for convergence
    err_lam2 = 1e-1 * max(1.0, float(lam_int2[0]), float(lam_int2[1]))
    err_fval = 1e-4

    err = b1 - a1
    k = 0
    while err > err_lam2 and k < k_max:
        fval_max = max(fa2_loss, fb2_loss)

        if fa2_loss > fb2_loss:
            a1 = a2
            a2 = b2
            fa2_loss, fa2_sigma = fb2_loss, fb2_sigma
            b2 = a1 + 0.618 * (b1 - a1)
            fb2_loss, fb2_sigma = cv_loss_cclasso(
                10**b2 / rho,
                x,
                k_cv,
                fa2_sigma,
                wd,
                u_f,
                u0_wd,
                d0_wd,
                wd2,
                iters,
                relative_errors=relative_errors,
            )
            lams.append(b2)
            fvals.append(fb2_loss)
        else:
            b1 = b2
            b2 = a2
            fb2_loss, fb2_sigma = fa2_loss, fa2_sigma
            a2 = a1 + 0.382 * (b1 - a1)
            fa2_loss, fa2_sigma = cv_loss_cclasso(
                10**a2 / rho,
                x,
                k_cv,
                fb2_sigma,
                wd,
                u_f,
                u0_wd,
                d0_wd,
                wd2,
                iters,
                relative_errors=relative_errors,
            )
            lams.append(a2)
            fvals.append(fa2_loss)

        fval_min = min(fa2_loss, fb2_loss)
        k += 1
        err = b1 - a1
        if abs(fval_max - fval_min) / (1 + fval_min) <= err_fval:
            break

    info_cv = {
        "lams": np.asarray(lams, dtype=float),
        "fvals": np.asarray(fvals, dtype=float),
        "k": k + 2,
        "lam_int": 10 ** np.array([a1, b1], dtype=float),
    }
    if a1 == lam_int2[0] or b1 == lam_int2[1]:
        warnings.warn(
            f"WARNING:\n\tOptimal lambda is near boundary! ([{10**a1},{10**b1}])\n",
            RuntimeWarning,
            stacklevel=2,
        )

    lambda_ = 10 ** ((a2 + b2) / 2)

    # Bootstrap for cclasso
    lambda2 = lambda_ / rho
    boot_indices, null_norm = _resolve_random_inputs(n, p, n_boot, boot_indices, null_norm, seed)
    info_boot = boot_cclasso(
        x,
        fb2_sigma,
        lambda2,
        n_boot,
        wd,
        u_f,
        u0_wd,
        d0_wd,
        boot_indices,
        null_norm,
        iters,
        relative_errors=relative_errors,
    )

    return {
        "var_w": info_boot["var_w"],
        "cor_w": info_boot["cor_w"],
        "p_vals": info_boot["p_vals"],
        "lambda": lambda_,
        "info_cv": info_cv,
        "_admm_iters": np.asarray(iters, dtype=int),
        "_admm_relative_errors": np.asarray(relative_errors, dtype=float),
        "_bootstrap_variances": info_boot["_bootstrap_variances"],
    }


def _resolve_random_inputs(
    n, p, n_boot, boot_indices, null_norm, seed
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    if boot_indices is None:
        # 0-based row indices in [0, n) (R's sample() is 1-based; we map at the edge).
        boot_indices = rng.integers(0, n, size=(n, n_boot))
    boot_indices = np.asarray(boot_indices, dtype=np.intp)
    if null_norm is None:
        null_norm = rng.standard_normal((n, p))
    null_norm = np.asarray(null_norm, dtype=float)
    return boot_indices, null_norm
