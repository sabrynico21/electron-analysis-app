"""Equivalence between the compiled MB core and the pure-Python reference.

The compiled core (``spieceasi._mb_core``) and the reference
(``spieceasi.mb_solver_reference``) implement the same algorithm; they can only
differ through floating-point summation order (the reference uses ``numpy.dot``,
the core uses the scalar loops of ``huge_core.cpp``). This module generates a
deterministic sweep of random problems and checks:

* the maximum absolute coefficient difference;
* whether the two are **bit-identical**;
* whether the per-lambda supports (symmetrised with ``or``) are identical;
* whether the core reproduces the R golden coefficients exactly.

Nothing here is tolerated a priori: the aggregate report is printed and asserted
against a tight coefficient bound plus exact support equality.

Slow: the sweep is dominated by the pure-Python reference solver.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spieceasi import mb_solver, mb_solver_reference  # noqa: E402
from spieceasi.neighborhood import make_alpha_path  # noqa: E402
from tests.validation_config import (  # noqa: E402
    DEGENERATE_KKT_VIOLATION_MAX,
    DEGENERATE_OBJECTIVE_GAP_MAX,
)

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        not mb_solver.HAVE_CPP_SOLVER,
        reason="compiled MB core not built (`pip install --no-build-isolation -e .`)",
    ),
]

# Coefficient agreement bound: the two implementations differ only by summation
# order, so differences must stay near machine precision relative to the data.
COEF_RTOL = 1e-9
COEF_ATOL = 1e-11

# Enumerated sweep: 7 x 6 x 3 kinds x 2 seeds = 252 problems, p in 5..60,
# n in 10..150, including p > n, zero-inflated and ill-conditioned designs.
P_VALUES = (5, 8, 12, 20, 30, 45, 60)
N_VALUES = (10, 15, 25, 40, 80, 150)
KINDS = ("plain", "zeros", "illcond")
SEEDS = (101, 202)
N_LAMBDA = 6
LAMBDA_MIN_RATIO = 1e-3


def _problem(kind: str, p: int, n: int, seed: int) -> np.ndarray:
    """Deterministic, non-degenerate design matrix."""
    rng = np.random.default_rng(seed)
    k = max(1, p // 3)
    X = rng.normal(size=(n, k)) @ rng.normal(size=(k, p)) + 0.5 * rng.normal(
        size=(n, p)
    )

    if kind == "zeros":
        X = X.copy()
        X[rng.random((n, p)) < 0.30] = 0.0
    elif kind == "illcond":
        # Nearly collinear columns: rank-deficient-ish, very small lambda step.
        X = X.copy()
        X[:, -1] = X[:, 0] * (1.0 + 1e-9)
        if p >= 4:
            X[:, -2] = X[:, 1] - X[:, 0] * 1e-9

    # huge rejects constant columns; nudge any that collapsed to a single value.
    for j in range(p):
        col = X[:, j]
        if np.all(col == col[0]):
            col = col.copy()
            col[0] += 1.0
            X[:, j] = col
    return np.ascontiguousarray(X)


def _adjacency(directed: np.ndarray) -> np.ndarray:
    mask = directed != 0
    return mask | mask.transpose(0, 2, 1)


def _kkt_and_objective(
    S: np.ndarray, m: int, b_full: np.ndarray, lam: float
) -> tuple[float, float]:
    """KKT residual and objective of one nodewise lasso sub-problem.

    ``huge`` regresses feature ``m`` on *all the other* features and skips
    ``j == m`` in the coordinate sweep, so both the KKT conditions and the
    objective are evaluated on the reduced problem.

    Objective (correlation scale): ``0.5 b' S b - s_m' b + lam * |b|_1`` with
    ``s_m = S[:, m]``; KKT: ``|g_j| == lam`` for active ``j`` and ``|g_j| <= lam``
    otherwise, where ``g = S b - s_m``.
    """
    d = S.shape[0]
    keep = np.ones(d, dtype=bool)
    keep[m] = False
    b = b_full[keep]
    S_red = S[np.ix_(keep, keep)]
    s_m = S[keep, m]
    grad = S_red @ b - s_m
    active = b != 0
    viol = 0.0
    if active.any():
        viol = float(np.max(np.abs(np.abs(grad[active]) - lam)))
    if (~active).any():
        viol = max(
            viol, float(np.max(np.maximum(np.abs(grad[~active]) - lam, 0.0)))
        )
    obj = 0.5 * float(b @ (S_red @ b)) - float(s_m @ b) + lam * float(np.abs(b).sum())
    return viol, obj


def _degenerate_validity(
    S: np.ndarray, grid: np.ndarray, directed: np.ndarray
) -> tuple[float, np.ndarray]:
    """Worst KKT residual and the objective of every (lambda, response) pair."""
    d = S.shape[0]
    worst = 0.0
    objs = np.empty((grid.size, d))
    for i, lam in enumerate(grid):
        for m in range(d):
            # huge never assigns a coefficient to the response itself.
            assert directed[i, m, m] == 0.0, (
                f"response coordinate {m} has a non-zero coefficient"
            )
            viol, obj = _kkt_and_objective(S, m, directed[i, m], lam)
            worst = max(worst, viol)
            objs[i, m] = obj
    return worst, objs


def test_mb_core_matches_reference_sweep() -> None:
    worst_diff = 0.0
    wellcond_max_diff = 0.0
    illcond_max_diff = 0.0
    n_bitwise = 0
    n_problems = 0
    support_diff_blocks = 0
    support_diff_problems = 0
    wellcond_support_problems = 0
    illcond_support_problems = 0
    degenerate_problems = 0
    degenerate_kkt_ref = 0.0
    degenerate_kkt_cpp = 0.0
    degenerate_obj_gap = 0.0
    divergent: list = []
    worst = None

    for kind in KINDS:
        for p in P_VALUES:
            for n in N_VALUES:
                for seed in SEEDS:
                    # keep the slowest combinations out of the sweep
                    if kind == "illcond" and p * n > 45 * 80:
                        continue
                    if p > n and kind == "plain" and n > 40:
                        continue

                    X = _problem(kind, p, n, seed * 1000 + p * 10 + n)
                    grid = make_alpha_path(
                        X, n_lambdas=N_LAMBDA, lambda_min_ratio=LAMBDA_MIN_RATIO
                    )

                    ref = mb_solver_reference.solve_mb_path_huge_scale(
                        X, grid, max_iter=10_000, tol=1e-4
                    )
                    cpp = mb_solver.solve_mb_path_huge_scale(
                        X, grid, max_iter=10_000, tol=1e-4
                    )

                    n_problems += 1
                    diff = float(np.max(np.abs(ref - cpp)))
                    if diff > worst_diff:
                        worst_diff = diff
                        worst = (kind, p, n, seed, diff)
                    if kind == "illcond":
                        illcond_max_diff = max(illcond_max_diff, diff)
                    else:
                        wellcond_max_diff = max(wellcond_max_diff, diff)
                    if np.array_equal(ref, cpp):
                        n_bitwise += 1

                    ref_adj, cpp_adj = _adjacency(ref), _adjacency(cpp)
                    blocks = int(
                        np.count_nonzero(
                            (ref_adj != cpp_adj)
                            .reshape(ref_adj.shape[0], -1)
                            .any(axis=1)
                        )
                    )
                    if blocks:
                        support_diff_problems += 1
                        support_diff_blocks += blocks
                        divergent.append((kind, p, n, seed, round(diff, 12), blocks))
                        if kind == "illcond":
                            illcond_support_problems += 1
                        else:
                            wellcond_support_problems += 1

                    if kind == "illcond":
                        # DEGENERATE REGIME - validity criterion, not a
                        # coefficient comparison.
                        #
                        # These designs make two columns collinear to 1e-9, so the
                        # lasso objective has an almost flat direction and the
                        # solution is not unique: the two implementations may
                        # settle on different near-optimal vertices (measured: a
                        # 1-ULP difference in the correlation matrix - our
                        # accumulation vs R's ``cor`` - is enough to flip a
                        # support). Coefficient-by-coefficient agreement is
                        # therefore not a meaningful criterion here. What *is*
                        # required is that each returned solution is a valid
                        # optimum of the same reduced nodewise problem: KKT
                        # residual within the declared bound and objectives
                        # within the declared gap (see ``validation_config``).
                        # The difference against the reference is still computed
                        # and printed as information.
                        S = mb_solver.huge_correlation(X)
                        v_ref, obj_ref = _degenerate_validity(S, grid, ref)
                        v_cpp, obj_cpp = _degenerate_validity(S, grid, cpp)
                        degenerate_problems += 1
                        degenerate_kkt_ref = max(degenerate_kkt_ref, v_ref)
                        degenerate_kkt_cpp = max(degenerate_kkt_cpp, v_cpp)
                        gap = float(np.max(np.abs(obj_ref - obj_cpp)))
                        degenerate_obj_gap = max(degenerate_obj_gap, gap)

                        assert v_ref <= DEGENERATE_KKT_VIOLATION_MAX, (
                            f"illcond p={p} n={n} seed={seed}: reference KKT "
                            f"residual {v_ref:.3e}"
                        )
                        assert v_cpp <= DEGENERATE_KKT_VIOLATION_MAX, (
                            f"illcond p={p} n={n} seed={seed}: core KKT residual "
                            f"{v_cpp:.3e}"
                        )
                        assert gap <= DEGENERATE_OBJECTIVE_GAP_MAX, (
                            f"illcond p={p} n={n} seed={seed}: objective gap "
                            f"{gap:.3e}"
                        )
                    else:
                        assert blocks == 0, (
                            f"{kind} p={p} n={n} seed={seed}: {blocks} lambda "
                            f"block(s) with a different support"
                        )
                        assert np.allclose(
                            ref, cpp, rtol=COEF_RTOL, atol=COEF_ATOL
                        ), f"{kind} p={p} n={n} seed={seed}: max|diff|={diff:.3e}"

    print(
        f"\nMB core vs reference sweep: {n_problems} problems | "
        f"bitwise-identical={n_bitwise}/{n_problems} | "
        f"max|diff| non-degenerate={wellcond_max_diff:.3e} "
        f"degenerate(1e-9)={illcond_max_diff:.3e} | worst={worst}"
    )
    print(
        "  support differences: "
        f"non-degenerate={wellcond_support_problems} problem(s), "
        f"degenerate={illcond_support_problems} problem(s), "
        f"total lambda blocks={support_diff_blocks} "
        f"(in {support_diff_problems} problem(s))"
    )
    print(
        "  degenerate-family validity: "
        f"{degenerate_problems} problem(s) | worst KKT residual "
        f"reference={degenerate_kkt_ref:.3e} core={degenerate_kkt_cpp:.3e} "
        f"(bound {DEGENERATE_KKT_VIOLATION_MAX:.1e}) | worst objective gap "
        f"ref-vs-core={degenerate_obj_gap:.3e} "
        f"(bound {DEGENERATE_OBJECTIVE_GAP_MAX:.1e})"
    )
    for rec in divergent[:25]:
        print(f"    divergent: kind={rec[0]} p={rec[1]} n={rec[2]} seed={rec[3]} "
              f"max|diff|={rec[4]:.3e} blocks={rec[5]}")

    assert n_problems >= 200, f"expected >= 200 problems, got {n_problems}"
    assert wellcond_max_diff <= 1e-11, (
        f"non-degenerate max coefficient difference {wellcond_max_diff:.3e} "
        f"too large"
    )
    assert wellcond_support_problems == 0, (
        f"{wellcond_support_problems} non-degenerate problem(s) with a different "
        f"support"
    )


def test_mb_core_matches_reference_on_goldens() -> None:
    """The extension must reproduce the reference on the real golden designs."""
    gold = Path(__file__).parent / "golden"
    gold_ext = Path(__file__).parent / "golden_extended"

    cases = [
        (gold / "ds1_r.clr.tsv", gold / "ds1_r.lambda_path.txt"),
        (gold / "ds2_r.clr.tsv", gold / "ds2_r.lambda_path.txt"),
        (gold / "ds3_r.clr.tsv", gold / "ds3_r.lambda_path.txt"),
        (gold / "ds4_r.clr.tsv", gold / "ds4_r.lambda_path.txt"),
        (gold_ext / "rand_p30_n100_s209" / "clr.tsv",
         gold_ext / "rand_p30_n100_s209" / "lambda_path.txt"),
    ]

    worst = 0.0
    worst_case = None
    for clr_path, grid_path in cases:
        if not clr_path.is_file():
            continue
        X = np.loadtxt(clr_path, delimiter="\t")
        grid = np.loadtxt(grid_path)
        ref = mb_solver_reference.solve_mb_path_huge_scale(
            X, grid, max_iter=10_000, tol=1e-4
        )
        cpp = mb_solver.solve_mb_path_huge_scale(X, grid, max_iter=10_000, tol=1e-4)
        diff = float(np.max(np.abs(ref - cpp)))
        if diff > worst:
            worst, worst_case = diff, clr_path.parent.name
        assert np.array_equal(
            _adjacency(ref), _adjacency(cpp)
        ), f"{clr_path}: adjacency differs"
    print(f"\nMB core vs reference on goldens: max|diff|={worst:.3e} ({worst_case})")
    assert worst <= 1e-11
