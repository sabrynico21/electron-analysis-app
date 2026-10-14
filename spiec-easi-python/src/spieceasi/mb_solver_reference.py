"""Reference pure-Python implementation of ``huge::mb()`` (huge_core.cpp).

This is the **reference** implementation and the automatic fallback of
:mod:`spieceasi.mb_solver`. When the compiled core (``spieceasi._mb_core``) is
importable it is used instead, and this module is then exercised by the
equivalence tests (extension vs reference). It is also what runs if the
extension cannot be built on a target platform.

It reproduces the Meinshausen-Buhlmann neighborhood-selection solver
of the ``huge`` C++ core **exactly**, including the behaviour that determines
which near-optimal vertex of the (ill-conditioned) Lasso problem is reached:

* one **warm-started** solve of the whole lambda path per node (coefficients
  carried from ``lambda[i-1]`` to ``lambda[i]``), matching ``huge_core.cpp``
  lines 679-810;
* the **sequential strong rule** ``2*lambda[i] - lambda[i-1]`` (``2*lambda[0]``
  for the first lambda);
* coordinate sweep order ``j = 0..d-1`` skipping the response ``m``;
* the same stopping criteria: inner active-set refinement with relative
  coefficient-change tolerance ``thol`` and ``MAX_ITER = 10000``, the pruning
  pass, and the full KKT certification pass that re-runs the strong loop on
  violations.

Rationale: for well-conditioned problems the Lasso solution is unique and any
correct solver agrees.  In ill-conditioned regimes (e.g. CLR-transformed data,
whose columns sum to zero and make the correlation matrix nearly singular, at
small lambda) the objective is flat and several near-optimal solutions with
**different supports** exist; the sweep order, warm start and stopping rules
then decide which one is reached.  Reproducing them is therefore required for
output parity with R.

Declared deviations from the C++ (documented in ``docs/AUDIT.md``):

* the residual inner products use ``numpy.dot`` (BLAS) instead of the scalar
  ``for`` accumulation loops in ``mb_partial_residual`` / ``mb_refine_active``;
  summation order may differ at the last-bit level.
* ``numpy`` works in double precision throughout, as does the C++ core.
"""

from __future__ import annotations

import numpy as np

# huge_core.cpp: `const double thol = 1e-4; const int MAX_ITER = 10000;`
HUGE_THOL = 1e-4
HUGE_MAX_ITER = 10_000


def soft_threshold(x: float, thr: float) -> float:
    """``threshold_l1`` from huge_core.cpp."""
    if x > thr:
        return x - thr
    if x < -thr:
        return x + thr
    return 0.0


def huge_correlation(X: np.ndarray) -> np.ndarray:
    """Centered/standardized (ddof=1) correlation matrix, diag forced to 1.

    Mirrors ``huge:::.huge_fast_cor`` (which equals ``cor(x)``).
    """
    means = np.mean(X, axis=0)
    centered = X - means
    scales = np.std(centered, axis=0, ddof=1)
    if np.any(scales <= 0):
        msg = "Each feature must have positive sample standard deviation."
        raise ValueError(msg)
    Z = centered / scales
    S = (Z.T @ Z).astype(float) / float(X.shape[0] - 1)
    np.fill_diagonal(S, 1.0)
    return S


def _residual(
    S: np.ndarray,
    m: int,
    j: int,
    idx_a: np.ndarray,
    size_a: int,
    w: np.ndarray,
) -> float:
    """``mb_partial_residual``: S[m,j] - sum_k S[j, a_k] w[a_k]."""
    if size_a:
        act = idx_a[:size_a]
        return float(S[m, j] - S[j, act].dot(w[act]))
    return float(S[m, j])


def solve_mb_path_huge_scale(
    X: np.ndarray,
    alpha_path: np.ndarray,
    *,
    max_iter: int = HUGE_MAX_ITER,
    tol: float = HUGE_THOL,
) -> np.ndarray:
    """Solve the MB nodewise Lasso path on the correlation scale, warm-started.

    Returns a ``(n_lambda, d, d)`` array where ``directed_path[i, m, k]`` is the
    coefficient of predictor ``k`` in the regression for response ``m`` at
    ``alpha_path[i]`` (the Python ``[node, predictor]`` convention; ``huge``
    stores the transpose).
    """
    X = np.asarray(X, dtype=float)
    n_samples, d = X.shape
    if d < 2:
        msg = "At least two features are required for MB neighborhood selection."
        raise ValueError(msg)
    S = huge_correlation(X)

    alpha_path = np.asarray(alpha_path, dtype=float)
    if alpha_path.ndim != 1 or alpha_path.size < 1:
        msg = "alpha_path must be a non-empty 1-D array."
        raise ValueError(msg)
    n_lambda = alpha_path.size

    max_cd_iter = min(int(max_iter), HUGE_MAX_ITER)
    if max_cd_iter < 1:
        max_cd_iter = 1
    # huge uses a fixed thol = 1e-4; the Python floor mirrors that and is a
    # documented deviation for callers requesting tol < 1e-4.
    thol = max(float(tol), HUGE_THOL)

    directed_path = np.zeros((n_lambda, d, d), dtype=float)

    for m in range(d):
        w0 = np.zeros(d, dtype=float)
        w1 = np.zeros(d, dtype=float)
        idx_a = np.empty(d, dtype=np.intp)
        size_a = 0
        idx_i = np.ones(d, dtype=np.uint8)
        idx_i[m] = 0
        grad_abs = np.abs(S[m, :]).copy()

        for i in range(n_lambda):
            lam = float(alpha_path[i])
            strong_thr = 2.0 * lam if i == 0 else 2.0 * lam - float(alpha_path[i - 1])
            strong = (idx_i == 1) & (grad_abs > strong_thr)

            iter_ext = 0
            certified = False
            while (not certified) and (iter_ext < max_cd_iter):
                gap_ext = 1
                while (gap_ext != 0) and (iter_ext < max_cd_iter):
                    size_prev = size_a

                    # --- strong-rule activation sweep (j = 0..d-1, skip m) ---
                    for j in range(d):
                        if idx_i[j] == 1 and strong[j]:
                            r = _residual(S, m, j, idx_a, size_a, w0)
                            grad_abs[j] = abs(r)
                            w_new = soft_threshold(r, lam)
                            w1[j] = w_new
                            if w_new != 0.0:
                                idx_a[size_a] = j
                                size_a += 1
                                idx_i[j] = 0
                            w0[j] = w_new
                    gap_ext = size_a - size_prev

                    # --- mb_refine_active: coordinate descent on active set --
                    gap_int = 1.0
                    iter_int = 0
                    while (gap_int > thol) and (iter_int < max_cd_iter):
                        change_sum = 0.0
                        coef_sum = 0.0
                        if size_a:
                            act = idx_a[:size_a]
                            for k in range(size_a):
                                w_idx = act[k]
                                r = (
                                    S[m, w_idx]
                                    + w0[w_idx]
                                    - S[w_idx, act].dot(w0[act])
                                )
                                w_new = soft_threshold(r, lam)
                                coef_sum += abs(w_new)
                                change_sum += abs(w_new - w0[w_idx])
                                w1[w_idx] = w_new
                                w0[w_idx] = w_new
                        gap_int = 0.0 if coef_sum == 0.0 else change_sum / coef_sum
                        iter_int += 1

                    # --- prune zeroed active coordinates ----------------------
                    junk = 0
                    for k in range(size_a):
                        w_idx = idx_a[k]
                        if w1[w_idx] == 0.0:
                            junk += 1
                            idx_i[w_idx] = 1
                            strong[w_idx] = True
                        else:
                            idx_a[k - junk] = w_idx
                    size_a -= junk
                    iter_ext += 1

                # --- KKT certification over all inactive coordinates ----------
                violations = 0
                for j in range(d):
                    if idx_i[j] == 1:
                        r = _residual(S, m, j, idx_a, size_a, w0)
                        grad_abs[j] = abs(r)
                        w_new = soft_threshold(r, lam)
                        w1[j] = w_new
                        if w_new != 0.0:
                            idx_a[size_a] = j
                            size_a += 1
                            idx_i[j] = 0
                            strong[j] = True
                            violations += 1
                        w0[j] = w_new
                certified = violations == 0

            row = np.zeros(d, dtype=float)
            if size_a:
                act = idx_a[:size_a]
                row[act] = w1[act]
            directed_path[i, m, :] = row

    return directed_path


def symmetrize_path(directed_path: np.ndarray, sym_rule: str) -> np.ndarray:
    """Symmetrize a directed coefficient path to a 0/1 adjacency path.

    Mirrors ``huge.mb``: ``path = sign(abs(beta) OP t(abs(beta)))``.
    """
    if sym_rule not in ("or", "and"):
        msg = f"Unsupported sym_rule '{sym_rule}'."
        raise ValueError(msg)
    nz = np.abs(directed_path) > 0
    if sym_rule == "or":
        adj = nz | np.transpose(nz, (0, 2, 1))
    else:
        adj = nz & np.transpose(nz, (0, 2, 1))
    adj = adj.astype(np.uint8)
    d = adj.shape[-1]
    diag = np.arange(d)
    adj[:, diag, diag] = 0
    return adj
