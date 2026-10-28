"""Faithful re-implementation of the `latentcor` Kendall-tau utilities.

R sources:
    latentcor/R/internal.R  (Kendalltau, n_x, encodeX, zratios)
"""

from __future__ import annotations

import numpy as np

__all__ = ["tie_count", "kendall_tau_b", "kendall_tau_a_lower"]


def tie_count(x: np.ndarray) -> float:
    """Return the number of tied pairs ``sum_t t(t-1)/2`` (R ``latentcor:::n_x``)."""
    _, counts = np.unique(x, return_counts=True)
    t = counts[counts > 1].astype(float)
    return float((t * (t - 1.0) / 2.0).sum())


def _tau_b_numerator(x: np.ndarray, y: np.ndarray) -> float:
    """Concordant minus discordant pairs, ignoring pairs tied in x or y.

    Implemented with a Fenwick tree over the ranks of ``y`` so ties contribute
    exactly zero, matching ``pcaPP::cor.fk`` (verified == tau-b).
    """
    n = x.size
    if n < 2:
        return 0.0
    order = np.lexsort((y, x))
    xs = x[order]
    ys = y[order]
    # compress y to ranks
    uniq, inv = np.unique(ys, return_inverse=True)
    m = uniq.size
    tree = np.zeros(m + 1, dtype=np.int64)

    def _add(i: int) -> None:
        i += 1
        while i <= m:
            tree[i] += 1
            i += i & (-i)

    def _prefix(i: int) -> int:
        # number of inserted elements with rank < i (i is 0-based rank)
        s = 0
        while i > 0:
            s += tree[i]
            i -= i & (-i)
        return s

    num = 0.0
    inserted = 0
    i = 0
    while i < n:
        j = i
        while j < n and xs[j] == xs[i]:
            j += 1
        # process the whole tie-group of x against previously inserted elements
        for k in range(i, j):
            r = int(inv[k])
            less = _prefix(r)
            leq = _prefix(r + 1)
            greater = inserted - leq
            num += less - greater
        for k in range(i, j):
            _add(int(inv[k]))
            inserted += 1
        i = j
    return num


def kendall_tau_b(x: np.ndarray, y: np.ndarray) -> float:
    """Kendall's tau-b (equals ``pcaPP::cor.fk`` for the supported cases)."""
    n = x.size
    n0 = n * (n - 1) / 2.0
    nx = tie_count(x)
    ny = tie_count(y)
    denom = np.sqrt((n0 - nx) * (n0 - ny))
    if denom == 0:
        return float("nan")
    return _tau_b_numerator(x, y) / denom


def kendall_tau_a_lower(X: np.ndarray) -> np.ndarray:
    """Kendall tau-a lower triangle of ``X`` with the latentcor tie correction.

    Returns a vector of length ``p*(p-1)/2`` in R's lower-triangle order
    (column-major: (2,1), (3,1), ..., (p,1), (3,2), ...).
    """
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    n0 = n * (n - 1) / 2.0
    n_x = np.array([tie_count(X[:, j]) for j in range(p)])
    n_x_sqd = np.sqrt(n0 - n_x)

    rows = []
    cols = []
    vals = []
    for c in range(p):
        for r in range(c + 1, p):
            kb = kendall_tau_b(X[:, r], X[:, c])
            btoa = n_x_sqd[r] * n_x_sqd[c] / n0
            rows.append(r)
            cols.append(c)
            vals.append(kb * btoa)
    return np.array(rows, dtype=int), np.array(cols, dtype=int), np.array(vals)
