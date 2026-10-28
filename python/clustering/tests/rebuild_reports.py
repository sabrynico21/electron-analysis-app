#!/usr/bin/env python3
"""Rebuild the Markdown reports from the raw JSON, without re-running anything.

Reads the newest ``results/tuning_*.json`` and ``results/comparison_*.json``
(or the paths given with ``--tuning`` / ``--comparison``) and re-renders
``RESULTS_TUNING.md`` and ``RESULTS_COMPARISON.md`` from them. The tables,
the Design / Reading notes and the score section are reproduced exactly; the
``What ... shows`` interpretation sections are **omitted by default** so the
output is pure data (add ``--findings`` to include them).

The ``Command:`` line keeps the command that produced the data: it is read
from the ``command`` field of the JSON, never from this script.

Examples:
    python rebuild_reports.py                 # tables only, newest JSONs
    python rebuild_reports.py --findings      # include the commentary
    python rebuild_reports.py --dry-run       # print instead of writing
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

TEST_DIR = Path(__file__).resolve().parent
CLUSTERING_DIR = TEST_DIR.parent
for path in (TEST_DIR, CLUSTERING_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import comparison_experiment  # noqa: E402
import tuning_experiment  # noqa: E402
from report_common import latest_file  # noqa: E402


def _provenance(output: dict, script_name: str, json_path: Path) -> str:
    """The command that produced the data recorded in the JSON."""

    recorded = output.get("command")
    if recorded:
        return recorded
    return (
        f"python {script_name} (rebuilt from {json_path.name}; "
        "original command not recorded in this JSON)"
    )


def _render(path: Path, script_name: str, builder, include_findings: bool) -> str:
    output = json.loads(path.read_text(encoding="utf-8"))
    command = _provenance(output, script_name, path)
    return builder(output, command, include_findings=include_findings)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir", default=str(TEST_DIR / "results"),
        help="Directory holding the JSON files and the reports.",
    )
    parser.add_argument(
        "--tuning", default=None,
        help="Path to a tuning_*.json (default: newest in --results-dir).",
    )
    parser.add_argument(
        "--comparison", default=None,
        help="Path to a comparison_*.json (default: newest in --results-dir).",
    )
    parser.add_argument(
        "--findings", action="store_true",
        help="Include the 'What ... shows' interpretation sections "
             "(default: tables only, no commentary).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the reports to stdout instead of writing the files.",
    )
    parser.add_argument(
        "--only", choices=["tuning", "comparison"], default=None,
        help="Rebuild only one report.",
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    results_dir = Path(args.results_dir).expanduser()

    jobs = []
    if args.only != "comparison":
        tuning_path = (
            Path(args.tuning).expanduser() if args.tuning
            else latest_file(results_dir, "tuning_*.json")
        )
        jobs.append((
            tuning_path,
            Path(tuning_experiment.__file__).name,
            tuning_experiment.build_markdown,
            results_dir / "RESULTS_TUNING.md",
        ))
    if args.only != "tuning":
        comparison_path = (
            Path(args.comparison).expanduser() if args.comparison
            else latest_file(results_dir, "comparison_*.json")
        )
        jobs.append((
            comparison_path,
            Path(comparison_experiment.__file__).name,
            comparison_experiment.build_markdown,
            results_dir / "RESULTS_COMPARISON.md",
        ))

    for json_path, script_name, builder, target in jobs:
        markdown = _render(json_path, script_name, builder, args.findings)
        if args.dry_run:
            print(f"===== {target.name} (from {json_path.name}) =====")
            print(markdown)
            continue
        target.write_text(markdown, encoding="utf-8")
        print(f"{target} rebuilt from {json_path.name} "
              f"({'with' if args.findings else 'without'} findings)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
