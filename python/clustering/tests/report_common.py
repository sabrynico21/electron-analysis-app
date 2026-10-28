"""Small shared helpers for the clustering result reports.

The generated reports deliberately contain only the run-specific information
(parameters, tables, findings).  The column definitions, the chance-corrected
metrics and the interpretation rules live once in ``docs/METHODOLOGY.md`` and
are linked from every report, so the same legend is never duplicated.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

#: The reports are written to ``results/``, so the shared documentation is
#: always one level up.
METHODOLOGY_LINK = "../docs/METHODOLOGY.md"


def md_table(headers: list[str], rows: Iterable[Iterable]) -> str:
    lines = [
        "| " + " | ".join(str(item) for item in headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(item) for item in row) + " |")
    return "\n".join(lines)


def pct(value: float, digits: int = 1) -> str:
    return f"{value:.{digits}%}"


def num(value: float, digits: int = 3) -> str:
    if value is None:
        return "-"
    try:
        finite = math.isfinite(value)
    except TypeError:
        return "-"
    if not finite:
        # "-inf" (ineligible) / "inf" (default ineligible) are not readable
        # scores; show an empty marker instead.
        return "-"
    return f"{value:.{digits}f}"


def size_cell(mean: float, node_count: int) -> str:
    return f"{mean:.1f} ({mean / node_count:.1%})"


def nested_size_cell(stats: dict) -> str:
    return f"{stats['mean']:.1f} ({stats['min']:.0f}-{stats['max']:.0f})"


def parameter_cell(params: dict) -> str:
    alpha = params.get("alpha", params.get("teleport_alpha"))
    phi = params.get("phi", params.get("target_conductance"))
    return f"alpha={alpha:g}, b={params['b']}, phi={phi:g}"


def latest_file(directory: Path, pattern: str) -> Path:
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No file matching {pattern!r} in {directory}")
    return files[-1]


def experiment_design(runs: int, nested_note: str) -> str:
    """The shared description of the two experiments, written once here."""

    return "\n".join([
        "",
        "## Design",
        "",
        f"- **Experiment 1 - repeatability**: {runs} runs from the same seed "
        "taxon, changing only the residual-queue order (the only source of "
        "variation: the algorithm is otherwise deterministic).",
        "- **Experiment 2 - nested consistency**: the deterministic cluster "
        "of the seed (the *base cluster*) is re-computed, and the algorithm "
        f"is then re-run from {nested_note}.",
        "- Every condition is run on the weighted graph (`|correlation|` "
        "edges) and on the unweighted graph (every edge has weight 1).",
        "",
        "Column definitions, the metrics and the interpretation rules live in "
        f"[docs/METHODOLOGY.md]({METHODOLOGY_LINK}); "
        "this report only lists the values.",
        "",
    ])


def related_files(entries: Iterable[str]) -> str:
    lines = ["## Related files", ""]
    lines.extend(f"- {entry}" for entry in entries)
    lines.append("")
    return "\n".join(lines)
