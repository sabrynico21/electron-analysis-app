#!/usr/bin/env python3
"""Generate docs/STATISTICAL_EQUIVALENCE.md from the Phase E reference runs.

Runs the Python pipeline for every seed of every dataset in
``tests/golden_stats/`` (no subsample injection), summarises the distributions
and writes the markdown report.

Usage (from spiec-easi-python/):
    ./.venv/bin/python tools/make_statistical_report.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
sys.path.insert(0, str(ROOT))

from tests.test_statistical_equivalence import (  # noqa: E402
    CASE_IDS,
    GOLDEN_STATS,
    _python_runs,
    _read_counts,
    _summary,
)


def main() -> int:
    sections: list[str] = []
    rows: list[str] = []

    for case in CASE_IDS:
        ref = json.loads((GOLDEN_STATS / case / "r_runs.json").read_text())
        counts = _read_counts(GOLDEN_STATS / case / "counts.tsv")
        py_runs = _python_runs(counts, int(ref["n_seeds"]))
        rep = _summary(case, ref["runs"], py_runs)

        rows.append(
            "| {case} | {p} x {n} | {seeds} | {ri} | {re:.2f} / {pe:.2f} | "
            "{rc:.4f} | {mc:.4f} | {ne} |".format(
                case=case,
                p=ref["p"],
                n=ref["n"],
                seeds=rep["n_seeds"],
                ri=rep["r_opt_index_counts"],
                re=rep["r_mean_edges"],
                pe=rep["py_mean_edges"],
                rc=rep["raw_curve_max_abs_diff"],
                mc=rep["mono_curve_max_abs_diff"],
                ne=f"{rep['n_distinct_edges_r']} / {rep['n_distinct_edges_py']}",
            )
        )

        sections.append(
            "\n".join(
                [
                    f"### {case}",
                    "",
                    f"- seeds: {rep['n_seeds']} (pulsar seed 1..{rep['n_seeds']} on the R side,",
                    f"  `random_state` 1..{rep['n_seeds']} on the Python side)",
                    f"- R selected-index counts: `{rep['r_opt_index_counts']}`",
                    f"- mean edge count: R {rep['r_mean_edges']:.3f} (sd {rep['r_sd_edges']:.3f}), "
                    f"Python {rep['py_mean_edges']:.3f} (sd {rep['py_sd_edges']:.3f})",
                    f"- distinct edges seen: R {rep['n_distinct_edges_r']}, "
                    f"Python {rep['n_distinct_edges_py']}",
                    f"- per-edge selection-frequency Pearson r: {rep['edge_freq_pearson']:.4f}"
                    if rep["edge_freq_pearson"] == rep["edge_freq_pearson"]
                    else "- per-edge selection-frequency Pearson r: n/a (degenerate)",
                    f"- mean raw instability curve, max |R - Python|: "
                    f"{rep['raw_curve_max_abs_diff']:.4f} (mean {rep['raw_curve_mean_abs_diff']:.4f})",
                    f"- mean cummax instability curve, max |R - Python|: "
                    f"{rep['mono_curve_max_abs_diff']:.4f} "
                    f"(mean {rep['mono_curve_mean_abs_diff']:.4f})",
                    "",
                    "| lambda index | R raw | Python raw | R cummax | Python cummax |",
                    "|---|---|---|---|---|",
                    *[
                        f"| {i} | {r:.4f} | {p:.4f} | {rm:.4f} | {pm:.4f} |"
                        for i, (r, p, rm, pm) in enumerate(
                            zip(
                                rep["r_raw_mean_curve"],
                                rep["py_raw_mean_curve"],
                                rep["r_mono_mean_curve"],
                                rep["py_mono_mean_curve"],
                            )
                        )
                    ],
                ]
            )
        )

    doc = "\n".join(
        [
            "# Statistical equivalence with R (Phase E)",
            "",
            "R's StARS subsample indices are **not** injected here: each implementation",
            "draws its own subsamples, so the runs cannot be identical. What is compared",
            "is the *distribution* over independent seeds of:",
            "",
            "- the selected lambda index,",
            "- the StARS instability curve (raw and cummax), averaged over seeds,",
            "- the edge count and the per-edge selection frequency.",
            "",
            "Reference runs: `tests/r_reference/generate_stats_runs.R` ->",
            "`tests/golden_stats/<case>/r_runs.json`. Python side:",
            "`tests/test_statistical_equivalence.py` (marked `slow`).",
            "",
            "## Summary",
            "",
            "| case | shape | seeds | R selected index | mean edges R / Python | "
            "raw curve max|diff| | cummax curve max|diff| | distinct edges R / Python |",
            "|---|---|---|---|---|---|---|---|",
            *rows,
            "",
            "## Detail",
            *sections,
            "",
            "## Interpretation and limits",
            "",
            "- For this synthetic data family StARS selects `lambda = lambda.max` for every",
            "  seed in both implementations, so the networks are empty (or contain the single",
            "  `argmax|cor|` edge) and the selected-index statistic is degenerate. This is a",
            "  property of the data, not of the implementation.",
            "- The discriminating statistics are therefore the instability curves, which agree",
            "  up to a few hundredths in raw/cummax instability at every lambda.",
            "- This is **not** a claim of general statistical equivalence: it holds for these",
            "  datasets, these 30 seeds and this data family only.",
            "",
        ]
    )

    out = REPO / "docs" / "STATISTICAL_EQUIVALENCE.md"
    out.write_text(doc)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
