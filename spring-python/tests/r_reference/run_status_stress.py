#!/usr/bin/env python3
"""Seeded R/Python status + discrete-output stress sweep for the SPRING port.

Third member of the family started by ``cclasso-python/tests/r_reference/
run_status_stress.py``; implements the predeclared 1,000-case seeded stress
suite required by ``docs/EQUIVALENCE_ACCEPTANCE_CONTRACT.md`` ("weighted toward
sparse, singular, and threshold-near inputs").

Design (predeclared; do not tune after seeing results)
------------------------------------------------------
* Master seed ``--seed`` (default 20261020) draws one per-case seed and one
  per-case draw seed from a ``numpy.random.default_rng`` stream, in case order.
* ``--cases`` cases are generated, ``case_index % 6`` selects the generator
  mode, so the six strata are equally represented:
    0. moderate Poisson-lognormal counts, ~15% zeros;
    1. sparse counts, ~65% zeros;
    2. rank-deficient design (one constant feature) + ~35% zeros;
    3. low library sizes (small Poisson mean) + ~45% zeros;
    4. p > n (p between n+1 and n+4, capped at 12 features);
    5. threshold-near design (weak/strong edge weights 0.25/0.35) + ~30% zeros.
  Dimensions are ``n in [12, 30]`` and ``p in [2, 8]`` (or the p > n stratum).
  Edge weights are capped at 0.7 / 0.6 so that the latent K parameters stay
  inside the multilinear interpolation region. A single pair that falls into
  latentcor's *bridge inversion* path costs ~20-25 s in Python (the
  deterministic 4-D normal CDF is evaluated inside Brent's method) and up to
  ~100 s per case once replicas are included, which would make a 1,000-case
  sweep take hours. That region is measured separately by
  ``compare_approx_original.R`` and the QMP goldens and is explicitly *out of
  scope* for this failure-status sweep. Therefore every candidate case (full
  data and all injected replicas) is screened with the same cutoff rule
  ``|K| > 0.9 * bound22(z1, z2)``; candidates with any cutoff pair are rejected
  and re-drawn from the same master stream (deterministically, up to
  ``MAX_REDRAWS`` attempts). The number of redraws is recorded per case and
  zero-cutoff is verified in ``summary.json``.
* Fixed SPRING settings: ``quantitative=FALSE`` (mclr path), ``nlambda=5``,
  ``lambda.min.ratio=1e-2`` with the explicit default path
  ``exp(seq(log(0.6), log(0.6*ratio), length=nlambda))``, ``thresh=0.1``,
  ``subsample.ratio=0.8``, ``rep.num=4``, ``Rmethod="approx"``,
  ``use.nearPD=TRUE``, ``nu=0.001``, ``ratio=0.9``.
  The *data-specific* lambda path is deliberately not used here: its known
  bridge-CDF lambda-max difference would be confounded with the MB/StARS
  comparison this sweep targets (it stays covered by the QMP golden cases).
* StARS subsample indices are drawn once with NumPy and injected into *both*
  languages (written 1-based for R), so the R RNG is never needed.

Classification (predeclared thresholds)
---------------------------------------
* ``status``: ``error`` (exception), ``nonfinite`` (any non-finite lambda,
  merge, summary, beta or refit entry), else ``finite``.
* ``opt_index`` disagreement is **boundary-classified** only when the
  monotonised stability curve of at least one language passes within
  ``OPT_BOUNDARY_MARGIN`` of ``thresh``.
* Refit edge-set disagreement is **boundary-classified** only when *every*
  differing edge satisfies at least one of:
  (a) its selection frequency at the selected lambda is within
      ``EDGE_BOUNDARY_FREQ`` of 0.5 (StARS tie), or
  (b) its MB coefficient is at the soft-threshold boundary in at least one
      language, ``|beta| <= COEF_BOUNDARY`` (the documented ``huge`` MB tie;
      this criterion is declared *before* the run, in response to the
      SPIEC-EASI sweep where the same tie had to be triaged afterwards).
  Anything else is an unexplained mismatch.

The runner exits nonzero if any status mismatch or unexplained discrete
mismatch exists, and keeps the complete case archive for investigation.

Usage::

    .venv/bin/python tests/r_reference/run_status_stress.py \\
        /tmp/spring-stress-1000 --cases 1000 --seed 20261020
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from spring_microbiome.latentcor.bridge import bound22
from spring_microbiome.latentcor.kendall import kendall_tau_a_lower
from spring_microbiome.mclr import mclr
from spring_microbiome.spring import default_lambdaseq, spring

MAX_REDRAWS = 20

NLAMBDA = 5
LAMBDA_MIN_RATIO = 1e-2
THRESH = 0.1
SUBSAMPLE_RATIO = 0.8
REP_NUM = 4
RMETHOD = "approx"
USE_NEARPD = True
NU = 0.001
RATIO = 0.9
RTOL = 1e-6

OPT_BOUNDARY_MARGIN = 1e-3
EDGE_BOUNDARY_FREQ = 0.05
COEF_BOUNDARY = 1e-6

STATUS_FIELDS = [
    "case",
    "status",
    "error_class",
    "opt_index",
    "n_edges",
    "edges",
    "beta_opt",
    "lambda_first",
    "stars_margin",
    "merge_opt",
    "seconds",
]


# ---------------------------------------------------------------------------
# deterministic case generation (shared design with the SPIEC-EASI sweep)
# ---------------------------------------------------------------------------


def _sparse_adj(
    rng: np.random.Generator,
    p: int,
    weak: float,
    strong: float,
    extra_frac: float = 1.0 / 3.0,
    strong_frac: float = 0.25,
) -> np.ndarray:
    """Two-tier sparse symmetric weighted adjacency."""
    w = np.zeros((p, p))
    for i in range(p - 1):
        w[i, i + 1] = weak
        w[i + 1, i] = weak
    cand = [(i, j) for i in range(p) for j in range(i + 1, p) if w[i, j] == 0.0]
    nadd = min(max(1, int(p * extra_frac)), len(cand))
    if nadd and cand:
        for k in rng.choice(len(cand), size=nadd, replace=False):
            i, j = cand[int(k)]
            w[i, j] = weak
            w[j, i] = weak
    edges = [(i, j) for i in range(p) for j in range(i + 1, p) if w[i, j] != 0.0]
    ns = min(max(1, int(len(edges) * strong_frac)), len(edges))
    if ns and edges:
        for k in rng.choice(len(edges), size=ns, replace=False):
            i, j = edges[int(k)]
            w[i, j] = strong
            w[j, i] = strong
    return w


def _gen_counts(
    rng: np.random.Generator,
    n: int,
    p: int,
    mu: float,
    weak: float,
    strong: float,
    zero_frac: float,
) -> np.ndarray:
    """Poisson-lognormal counts from a sparse, two-tier Gaussian copula."""
    w = _sparse_adj(rng, p, weak, strong)
    theta = np.diag(1.0 + w.sum(axis=1)) + w
    sigma = np.cov(np.linalg.inv(theta))
    d = np.sqrt(np.diag(sigma))
    sigma = sigma / np.outer(d, d)
    sigma = 0.5 * (sigma + sigma.T)
    eta = rng.multivariate_normal(np.zeros(p), sigma, size=n)
    lam = np.exp(np.log(mu) + eta)
    counts = rng.poisson(lam).astype(float)
    if zero_frac > 0:
        threshold = np.quantile(counts, zero_frac)
        counts[counts <= threshold] = 0.0
    return counts


def _make_case(
    case_rng: np.random.Generator, case_index: int
) -> tuple[np.ndarray, dict[str, object]]:
    mode = case_index % 6
    n = int(case_rng.integers(12, 31))
    if mode == 4:
        p = int(min(n + case_rng.integers(1, 5), 12))
    else:
        p = int(case_rng.integers(2, min(8, n + 2) + 1))

    mu, weak, strong, zero_frac = 10.0, 0.3, 0.7, 0.15
    if mode == 1:
        mu, weak, strong, zero_frac = 8.0, 0.25, 0.5, 0.65
    elif mode == 3:
        mu, weak, strong, zero_frac = 2.0, 0.25, 0.6, 0.45
    elif mode == 4:
        zero_frac = 0.30
    elif mode == 5:
        mu, weak, strong, zero_frac = 6.0, 0.25, 0.35, 0.30

    counts = _gen_counts(case_rng, n, p, mu, weak, strong, zero_frac)
    if mode == 2:
        counts[:, -1] = 1.0  # constant positive feature -> rank deficiency

    # Sample rows with zero total are not a valid SPRING input (both languages
    # reject them). Keep the sparsity pattern but guarantee one positive count
    # per row so the sweep targets the network/matrix degeneracy instead.
    for row in np.flatnonzero(counts.sum(axis=1) == 0):
        counts[row, int(case_rng.integers(0, p))] = 1.0

    params: dict[str, object] = {
        "mode": mode,
        "n": int(n),
        "p": int(p),
        "mu": float(mu),
        "weak": float(weak),
        "strong": float(strong),
        "zero_frac": float(zero_frac),
        "quantitative": False,
        "nlambda": NLAMBDA,
        "lambda_min_ratio": LAMBDA_MIN_RATIO,
        "thresh": THRESH,
        "subsample_ratio": SUBSAMPLE_RATIO,
        "rep_num": REP_NUM,
        "Rmethod": RMETHOD,
        "use_nearPD": USE_NEARPD,
        "nu": NU,
        "ratio": RATIO,
        "Rtol": RTOL,
    }
    return counts, params


# ---------------------------------------------------------------------------
# Python-side classification
# ---------------------------------------------------------------------------


def _edges_string(adjacency: np.ndarray) -> tuple[int, str]:
    idx = np.argwhere((adjacency != 0) & np.triu(np.ones_like(adjacency, bool), 1))
    edges = ";".join(f"{int(i)}:{int(j)}" for i, j in idx)
    return int(len(idx)), edges


def _beta_string(beta: np.ndarray) -> str:
    rows = []
    for i in range(beta.shape[0]):
        for j in range(beta.shape[1]):
            if beta[i, j] != 0.0:
                rows.append(f"{i}:{j}:{float(beta[i, j]):.6g}")
    return ";".join(rows)


def _merge_string(merge: np.ndarray) -> str:
    rows = []
    for i in range(merge.shape[0]):
        for j in range(i + 1, merge.shape[1]):
            rows.append(f"{i}:{j}:{float(merge[i, j]):.6g}")
    return ";".join(rows)


def _cutoff_pairs(counts: np.ndarray, indices: np.ndarray) -> int:
    """Number of pairs that would take latentcor's bridge-inversion path.

    Mirrors ``r_switch22``'s ``|K| > ratio * bound22(z1, z2)`` cutoff on the
    full data and on every injected replica (the matrices ``huge_kmb`` sees).
    """
    total = 0
    for data in [counts, *(counts[ix] for ix in indices)]:
        qdat = mclr(data)
        encoded = qdat - qdat.min(axis=0)
        zratio = (encoded == 0).mean(axis=0)
        rows, cols, k = kendall_tau_a_lower(encoded)
        ok = np.isfinite(k)
        if not ok.any():
            continue
        bound = bound22(zratio[rows][ok], zratio[cols][ok])
        total += int(np.sum(np.abs(k[ok]) > RATIO * bound))
    return total


def _python_status(counts: np.ndarray, indices: np.ndarray, case_name: str) -> dict[str, object]:
    start = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            result = spring(
                counts,
                quantitative=False,
                lambdaseq=default_lambdaseq(LAMBDA_MIN_RATIO, NLAMBDA),
                thresh=THRESH,
                subsample_ratio=SUBSAMPLE_RATIO,
                rep_num=REP_NUM,
                Rmethod=RMETHOD,
                use_nearPD=USE_NEARPD,
                nu=NU,
                ratio=RATIO,
                Rtol=RTOL,
                verbose=False,
                subsample_indices=indices,
            )
        except Exception as exc:  # status is compared, not the exception text
            return {
                "case": case_name,
                "status": "error",
                "error_class": type(exc).__name__,
                "opt_index": "",
                "n_edges": "",
                "edges": "",
                "beta_opt": "",
                "lambda_first": "",
                "stars_margin": "",
                "merge_opt": "",
                "seconds": f"{time.perf_counter() - start:.6f}",
            }
    elapsed = time.perf_counter() - start

    stars = result.output.stars
    opt = int(stars.opt_index)
    est = result.fit.est
    refit = np.asarray(result.fit.refit["stars"], dtype=float)
    beta_opt = np.asarray(est.beta[opt - 1], dtype=float)

    arrays: list[np.ndarray] = [
        np.asarray(result.lambdaseq, dtype=float),
        np.asarray(stars.summary, dtype=float),
        refit,
        beta_opt,
    ]
    arrays.extend(np.asarray(m, dtype=float) for m in stars.merge)
    arrays.extend(np.asarray(b, dtype=float) for b in est.beta)
    arrays.extend(np.asarray(p, dtype=float) for p in est.path)
    nonfinite = any(not np.all(np.isfinite(a)) for a in arrays)

    n_edges, edges = _edges_string(refit)
    stability = np.asarray(stars.summary, dtype=float)
    return {
        "case": case_name,
        "status": "nonfinite" if nonfinite else "finite",
        "error_class": "",
        "opt_index": opt,
        "n_edges": n_edges,
        "edges": "" if nonfinite else edges,
        "beta_opt": "" if nonfinite else _beta_string(beta_opt),
        "lambda_first": f"{float(result.lambdaseq[0]):.17g}",
        "stars_margin": f"{float(np.min(np.abs(stability - THRESH))):.17g}",
        "merge_opt": _merge_string(np.asarray(stars.merge[opt - 1], dtype=float)),
        "seconds": f"{elapsed:.6f}",
    }


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------


def _freq_margin(merge_string: str, edge: str) -> float:
    """Distance of an edge's selection frequency from 0.5 (large = decisive)."""
    for token in merge_string.split(";"):
        if not token:
            continue
        i, j, value = token.split(":")
        if f"{i}:{j}" == edge:
            return abs(float(value) - 0.5)
    return 0.5


def _coef_margin(beta_string: str, edge: str) -> float:
    """Smallest |coefficient| recorded for an edge (0.0 when the edge is absent)."""
    for token in str(beta_string).split(";"):
        if not token:
            continue
        i, j, value = token.split(":")
        if f"{i}:{j}" == edge or f"{j}:{i}" == edge:
            return abs(float(value))
    return 0.0


def _compare(
    python_rows: list[dict[str, object]],
    r_rows: dict[str, dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, int]]:
    unexplained: list[dict[str, object]] = []
    boundary: list[dict[str, object]] = []
    counts = {
        "status_mismatch": 0,
        "opt_index_mismatch": 0,
        "opt_index_boundary": 0,
        "edge_mismatch": 0,
        "edge_boundary": 0,
    }
    for row in python_rows:
        ref = r_rows.get(row["case"])
        if ref is None:
            unexplained.append({"case": row["case"], "reason": "missing R row"})
            counts["status_mismatch"] += 1
            continue
        if ref["status"] != row["status"]:
            counts["status_mismatch"] += 1
            unexplained.append(
                {
                    "case": row["case"],
                    "reason": "status",
                    "python": row["status"],
                    "R": ref["status"],
                    "python_error": row["error_class"],
                    "R_error": ref["error_class"],
                }
            )
            continue
        if row["status"] != "finite":
            continue

        if int(ref["opt_index"]) != int(row["opt_index"]):
            counts["opt_index_mismatch"] += 1
            margin = min(float(row["stars_margin"]), float(ref["stars_margin"]))
            entry = {
                "case": row["case"],
                "reason": "opt_index",
                "python": row["opt_index"],
                "R": ref["opt_index"],
                "stars_margin": f"{margin:.3g}",
            }
            if margin <= OPT_BOUNDARY_MARGIN:
                counts["opt_index_boundary"] += 1
                boundary.append(entry)
            else:
                unexplained.append(entry)

        py_edges = set(filter(None, str(row["edges"]).split(";")))
        r_edges = set(filter(None, ref["edges"].split(";")))
        differing = py_edges ^ r_edges
        if differing:
            counts["edge_mismatch"] += 1
            freq = min(_freq_margin(str(row["merge_opt"]), e) for e in differing)
            coef = min(
                min(_coef_margin(str(row["beta_opt"]), e), _coef_margin(ref["beta_opt"], e))
                for e in differing
            )
            entry = {
                "case": row["case"],
                "reason": "edge_set",
                "python_only": sorted(py_edges - r_edges)[:10],
                "R_only": sorted(r_edges - py_edges)[:10],
                "min_freq_margin": f"{freq:.3g}",
                "min_abs_coef": f"{coef:.3g}",
            }
            if freq <= EDGE_BOUNDARY_FREQ or coef <= COEF_BOUNDARY:
                counts["edge_boundary"] += 1
                boundary.append(entry)
            else:
                unexplained.append(entry)
    return unexplained, boundary, counts


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def _write_matrix(path: Path, values: np.ndarray, fmt: str = "%.17g") -> None:
    np.savetxt(path, np.asarray(values), delimiter=",", fmt=fmt)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--cases", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error("--cases must be positive")
    out_dir = args.out_dir
    if out_dir.exists() and any(out_dir.iterdir()):
        parser.error(f"refusing to overwrite non-empty output directory: {out_dir}")
    cases_dir = out_dir / "cases"
    cases_dir.mkdir(parents=True, exist_ok=True)

    master_rng = np.random.default_rng(args.seed)
    python_rows: list[dict[str, object]] = []
    redraw_total = 0
    residual_cutoff_cases = 0
    sweep_start = time.perf_counter()
    for i in range(args.cases):
        attempts = 0
        while True:
            case_seed = int(master_rng.integers(0, 2**32 - 1))
            draw_seed = int(master_rng.integers(0, 2**32 - 1))
            case_rng = np.random.default_rng(case_seed)
            counts, params = _make_case(case_rng, i)
            n = counts.shape[0]
            ratio = 0.8 if n <= 144 else 10.0 * np.sqrt(n) / n
            size = int(np.floor(n * ratio))
            draw_rng = np.random.default_rng(draw_seed)
            indices = np.array(
                [draw_rng.choice(n, size=size, replace=False) for _ in range(REP_NUM)],
                dtype=np.intp,
            )
            cutoff = _cutoff_pairs(counts, indices)
            attempts += 1
            if cutoff == 0 or attempts > MAX_REDRAWS:
                break

        redraw_total += attempts - 1
        residual_cutoff_cases += int(cutoff > 0)

        case_dir = cases_dir / f"case-{i:05d}"
        case_dir.mkdir()
        _write_matrix(case_dir / "counts.csv", counts)
        _write_matrix(case_dir / "subsample_indices.csv", indices + 1, fmt="%d")
        params.update(
            {
                "case_seed": case_seed,
                "draw_seed": draw_seed,
                "subsample_ratio_eff": float(ratio),
                "subsample_size": int(size),
                "redraws": int(attempts - 1),
                "bridge_cutoff_pairs": int(cutoff),
            }
        )
        (case_dir / "params.json").write_text(json.dumps(params, indent=2) + "\n", encoding="utf-8")

        python_rows.append(_python_status(counts, indices, case_dir.name))
        if (i + 1) % 100 == 0:
            print(f"[python] {i + 1}/{args.cases} cases", flush=True)

    python_seconds = [float(row["seconds"]) for row in python_rows]
    with (out_dir / "python_status.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(python_rows)

    r_status = out_dir / "r_status.csv"
    r_runner = Path(__file__).with_name("stress_status_R.R")
    r_start = time.perf_counter()
    r_run = subprocess.run(
        ["Rscript", str(r_runner), str(cases_dir), str(r_status)],
        cwd=Path(__file__).resolve().parents[2],
        capture_output=True,
        text=True,
    )
    r_elapsed = time.perf_counter() - r_start
    if r_run.returncode:
        raise RuntimeError(
            "R stress classification failed:\n"
            f"stdout: {r_run.stdout[-4000:]}\nstderr: {r_run.stderr[-4000:]}"
        )

    with r_status.open(newline="", encoding="utf-8") as stream:
        r_rows = {row["case"]: row for row in csv.DictReader(stream)}

    unexplained, boundary, counts = _compare(python_rows, r_rows)

    summary = {
        "seed": args.seed,
        "cases": args.cases,
        "design": {
            "quantitative": False,
            "nlambda": NLAMBDA,
            "lambda_min_ratio": LAMBDA_MIN_RATIO,
            "thresh": THRESH,
            "subsample_ratio": SUBSAMPLE_RATIO,
            "rep_num": REP_NUM,
            "Rmethod": RMETHOD,
            "use_nearPD": USE_NEARPD,
            "nu": NU,
            "ratio": RATIO,
            "Rtol": RTOL,
            "opt_boundary_margin": OPT_BOUNDARY_MARGIN,
            "edge_boundary_freq": EDGE_BOUNDARY_FREQ,
            "coef_boundary": COEF_BOUNDARY,
            "max_redraws": MAX_REDRAWS,
        },
        "bridge_inversion_screening": {
            "rule": "|K| > 0.9 * bound22(z1, z2) on full data and every replica",
            "rejected_candidates": redraw_total,
            "executed_cases_with_cutoff_pairs": residual_cutoff_cases,
        },
        "python_status_counts": {
            status: sum(row["status"] == status for row in python_rows)
            for status in ("finite", "nonfinite", "error")
        },
        "r_status_counts": {
            status: sum(row["status"] == status for row in r_rows.values())
            for status in ("finite", "nonfinite", "error")
        },
        "mismatch_counts": counts,
        "unexplained_count": len(unexplained),
        "unexplained_preview": unexplained[:20],
        "boundary_classified_count": len(boundary),
        "boundary_preview": boundary[:20],
        "timing": {
            "python_total_s": round(sum(python_seconds), 3),
            "python_mean_s": round(float(np.mean(python_seconds)), 6),
            "python_median_s": round(float(np.median(python_seconds)), 6),
            "python_p95_s": round(float(np.percentile(python_seconds, 95)), 6),
            "r_total_s": round(r_elapsed, 3),
            "sweep_wall_s": round(time.perf_counter() - sweep_start, 3),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    environment = {
        "python": sys.version,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "r_version": subprocess.run(
            ["Rscript", "--version"], capture_output=True, text=True
        ).stdout.splitlines()[0],
    }
    (out_dir / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")

    print(json.dumps(summary, indent=2))
    if unexplained:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
