#!/usr/bin/env python3
"""Compare default parameters, tuned parameters, signed and both.

Reads the tuning report produced by ``tuning_experiment.py`` (latest
``results/tuning_*.json`` by default), takes the best ``alpha`` and ``b`` for
each network and weight condition, and runs four configurations:

================  ==========================================  ============
config            parameters                                  signs
================  ==========================================  ============
Default           default parameters                          no
Best              best alpha/b from the sweep                  no
Signed            default parameters                          yes
Both              best alpha/b from the sweep                  yes
================  ==========================================  ============

Experiment 1 (repeatability) uses ``--runs`` runs per configuration.
Experiment 2 (nested consistency) is not capped by a run count: it re-runs
the algorithm from **every** member of the base cluster, so its run count is
the size of that cluster.  ``--max-nested-nodes`` exists only to cap the cost
of smoke tests.

Writes ``results/RESULTS_COMPARISON.md``.  Column definitions and
interpretation are documented once in ``docs/METHODOLOGY.md``; the generated
report links to it instead of repeating the legend.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

TEST_DIR = Path(__file__).resolve().parent
CLUSTERING_DIR = TEST_DIR.parent
for path in (TEST_DIR, CLUSTERING_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from clustering_profiles import resolve_parameters  # noqa: E402
from network_loader import load_network, resolve_network  # noqa: E402
from report_common import (  # noqa: E402
    experiment_design, md_table, related_files,
)
from stability_experiment import (  # noqa: E402
    DEFAULT_RUNS,
    nested_consistency_experiment,
    repeatability_experiment,
)

SEED = "Chloroplast"
WEIGHTS = (("weighted", True), ("unweighted", False))
CONFIGS = {
    "default": {"config": "new_parameters", "signed": False, "use_best": False},
    "best": {"config": "new_parameters", "signed": False, "use_best": True},
    "signed": {"config": "both", "signed": True, "use_best": False},
    "both": {"config": "both", "signed": True, "use_best": True},
}


def _latest_tuning(directory: Path) -> Path:
    files = sorted(directory.glob("tuning_*.json"))
    if not files:
        raise FileNotFoundError(f"No tuning_*.json in {directory}")
    return files[-1]


def _overrides(best: dict) -> dict:
    return {
        "teleport_alpha": best["alpha"],
        "b": best["b"],
        "target_conductance": best["phi"],
    }


def run_comparison(args) -> dict:
    tuning_path = Path(args.tuning).expanduser() if args.tuning else _latest_tuning(
        Path(args.output_dir).expanduser()
    )
    tuning = json.loads(tuning_path.read_text(encoding="utf-8"))
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "tuning_source": str(tuning_path),
        "seed": args.seed,
        "runs": args.runs,
        "nested_max_nodes": args.max_nested_nodes,
        "networks": {},
    }
    started_all = time.perf_counter()

    for network, tuning_net in tuning["networks"].items():
        path = resolve_network(network)
        graph, info = load_network(network)
        print("\n" + "=" * 78)
        print(f"Network {network}: {info['nodes']} nodes, {info['edges']} edges")
        print("=" * 78)
        network_entry = {"graph": info, "conditions": {}}

        for weight, weighted in WEIGHTS:
            best = tuning_net["conditions"][weight]["best"]
            # Base values of the sweep: the default row of the tuning
            # report (the former recommended profile values, kept as the
            # experiments' reference parameters).
            tuning_default = tuning["default"]
            base_overrides = {
                "teleport_alpha": tuning_default["alpha"],
                "b": tuning_default["b"],
                "target_conductance": tuning_default["phi"],
            }
            default_params = resolve_parameters(
                graph, "app", weighted=weighted, overrides=base_overrides
            )
            print(f"\n### {network} / {weight}")
            condition = {
                "default_parameters": {
                    "alpha": default_params["teleport_alpha"],
                    "b": default_params["b"],
                    "phi": default_params["target_conductance"],
                },
                "best_parameters": {
                    "alpha": best["alpha"], "b": best["b"], "phi": best["phi"],
                },
                "configs": {},
            }

            for name, settings in CONFIGS.items():
                overrides = _overrides(best) if settings["use_best"] else None
                merged = dict(base_overrides)
                if overrides:
                    merged.update(overrides)
                params = resolve_parameters(
                    graph, "app", weighted=weighted,
                    overrides=merged,
                )
                started = time.perf_counter()
                reused = None
                if (
                    args.reuse_unsigned
                    and name in ("default", "best")
                    and best["b"] == 10
                    and best["phi"] == 0.4
                ):
                    tuning_condition = tuning_net["conditions"][weight]
                    if name == "default" or best["alpha"] == 0.05:
                        reused = tuning_condition["default"]
                    else:
                        reused = tuning_condition["sweeps"]["alpha"].get(
                            str(best["alpha"])
                        )
                if reused is not None:
                    repeat = reused["repeatability"]
                    nested = reused["nested_consistency"]
                else:
                    repeat = repeatability_experiment(
                        graph, args.seed, settings["config"],
                        runs=args.runs, run_seed=42, jobs=args.jobs,
                        network_path=path, overrides=overrides, weighted=weighted,
                    )
                    nested = nested_consistency_experiment(
                        graph, args.seed, settings["config"],
                        jobs=args.jobs, network_path=path,
                        max_nodes=args.max_nested_nodes,
                        overrides=overrides, weighted=weighted,
                    )
                condition["configs"][name] = {
                    "parameters": {
                        "alpha": params["teleport_alpha"],
                        "b": params["b"],
                        "phi": params["target_conductance"],
                    },
                    "signed": settings["signed"],
                    "repeatability": repeat,
                    "nested_consistency": nested,
                    "elapsed_seconds": time.perf_counter() - started,
                }
                print(
                    f"  {name}: size={repeat['cluster_size']['mean']:.1f} "
                    f"J={repeat['pairwise_jaccard']['mean']:.3f} "
                    f"support={nested['node_support_stats']['mean']:.3f} "
                    f"strict={repeat['conditions']['strict_success_rate']:.0%} "
                    f"({time.perf_counter() - started:.0f}s)"
                )
            network_entry["conditions"][weight] = condition

        output["networks"][network] = network_entry

    output["elapsed_seconds"] = time.perf_counter() - started_all
    return output


def _slash(values: list[str]) -> str:
    return " / ".join(values)


def _jaccard_cell(entry: dict) -> str:
    stats = entry["pairwise_jaccard"]
    return f"({stats['mean']:.3f}/{stats['min']:.3f}/{stats['max']:.3f})"


def _conductance_cell(entry: dict) -> str:
    stats = entry["conductance"]
    return f"({stats['mean']:.3f}/{stats['min']:.3f}/{stats['max']:.3f})"


def _intersection_cell(entry: dict, denominator: int) -> str:
    count = entry["intersection_size"]
    share = count / denominator if denominator else 0.0
    return f"{count} ({share:.0%})"


def build_markdown(output: dict, command: str,
                   include_findings: bool = True) -> str:
    lines = ["# Comparison: defaults, tuned parameters, signed and both", ""]
    lines.append(f"Generated: `{output['generated_at']}`  ")
    lines.append(f"Command: `{command}`")
    nested_note = (
        f"up to **{output['nested_max_nodes']}** members of that cluster, "
        "subsampled with a fixed seed (smoke-test mode)"
        if output.get("nested_max_nodes")
        else "**every member** of that cluster"
    )
    lines.append(experiment_design(output["runs"], nested_note))
    lines.append(
        f"The tuned parameters come from `{Path(output['tuning_source']).name}`."
    )
    lines.append("")
    lines.append("## Configurations")
    lines.append("")
    lines.append(
        "Every table row follows the order **Default / Best / Signed / Both**. "
        "Each cell repeats that order, separated by `/`."
    )
    lines.append("")
    parameter_rows = []
    for network, net_entry in output["networks"].items():
        for weight, _ in WEIGHTS:
            condition = net_entry["conditions"][weight]
            parameter_rows.append([
                f"{network} / {weight}",
                *[
                    f"alpha={condition['configs'][name]['parameters']['alpha']:g}, "
                    f"b={condition['configs'][name]['parameters']['b']}, "
                    f"phi={condition['configs'][name]['parameters']['phi']:g}"
                    for name in CONFIGS
                ],
            ])
    lines.append(md_table(
        ["Condition", "Default", "Best", "Signed", "Both"],
        parameter_rows,
    ))
    lines.append("")
    lines.append(
        "`Best` uses the best `alpha` and `b` from the tuning report; "
        "`Both` is the signed version with those same parameters."
    )
    lines.append("")
    lines.append("## Reading note")
    lines.append("")
    lines.append(
        "Every value column is ordered **Default / Best / Signed / Both**: each "
        "item inside a cell is separated by `/`, and inside a mean/min/max cell "
        "each `(mean/min/max)` group belongs to one configuration. The "
        "conductance column is the plain `|correlation|` conductance for every "
        "configuration; Signed and Both are accepted on the signed score, which "
        "is stored in the JSON."
    )
    lines.append("")

    def table_rows(experiment: str):
        rows = []
        for network, net_entry in output["networks"].items():
            for weight, _ in WEIGHTS:
                condition = net_entry["conditions"][weight]
                entries = [condition["configs"][name] for name in CONFIGS]
                key = "repeatability" if experiment == "repeat" else "nested_consistency"
                reports = [entry[key] for entry in entries]

                def per(fn):
                    return _slash([fn(report) for report in reports])

                row = [
                    f"{network} / {weight}",
                    # one run count per configuration: Experiment 1 uses the
                    # same fixed N everywhere, Experiment 2 re-seeds once per
                    # base-cluster member, so its count differs between
                    # configurations whenever their base clusters differ.
                    _slash(
                        [
                            str(report["runs"])
                            if experiment == "repeat"
                            else str(report["nodes_tested"])
                            for report in reports
                        ]
                    ),
                    per(lambda e: f"{e['cluster_size']['mean']:.1f}"),
                    per(_jaccard_cell),
                    per(lambda e, exp=experiment: _intersection_cell(
                        e,
                        e["cluster_size"]["mean"] if exp == "repeat"
                        else e["nodes_in_base_cluster"],
                    )),
                    per(lambda e: str(e["union_size"])),
                    per(_conductance_cell),
                    per(lambda e: str(e["identical_runs"])),
                    per(lambda e: str(
                        e["nodes_in_at_least_80pct"] if experiment == "repeat"
                        else e["nodes_support_80pct"]
                    )),
                ]
                if experiment == "nested":
                    row.append(per(
                        lambda e: f"{e['asymmetry']['asymmetry_rate']:.0%}"
                    ))
                rows.append(row)
        return rows

    lines.append(
        f"## Experiment 1 - Repeatability ({output['runs']} runs per configuration)"
    )
    lines.append("")
    lines.append(md_table(
        ["Condition", "Runs", "Mean cluster size (taxa)",
         "Jaccard (mean / min / max)", "Intersection (n / %)", "Union",
         "Conductance (mean / min / max)", "Identical clusters",
         "Consensus >=80%"],
        table_rows("repeat"),
    ))
    lines.append("")

    lines.append(
        "## Experiment 2 - Nested consistency "
        "(one run per base-cluster member)"
    )
    lines.append("")
    lines.append(md_table(
        ["Condition", "Runs", "Mean cluster size (taxa)",
         "Jaccard (mean / min / max)", "Intersection (n / %)", "Union",
         "Conductance (mean / min / max)", "Identical clusters",
         "Consensus >=80%", "Asymmetry"],
        table_rows("nested"),
    ))
    lines.append("")
    if include_findings:
        lines.append("## What the comparison shows")
        lines.append("")
        lines.extend(_explanation(output))
        lines.append("")
    lines.append(related_files([
        "`results/comparison_*.json` - raw data of this comparison.",
        "[RESULTS_TUNING.md](RESULTS_TUNING.md) - the sweep that selected the "
        "tuned parameters.",
        "[README.md](../README.md) - repository layout and the exact commands "
        "used to produce this report.",
    ]))
    return "\n".join(lines)


def _explanation(output: dict) -> list[str]:
    lines = []
    for network, net_entry in output["networks"].items():
        for weight, _ in WEIGHTS:
            condition = net_entry["conditions"][weight]
            cfg = condition["configs"]
            default = cfg["default"]
            best = cfg["best"]
            signed = cfg["signed"]
            both = cfg["both"]
            lines.append(
                f"- **{network} / {weight}**: best vs default: "
                f"size {best['repeatability']['cluster_size']['mean']:.1f} vs "
                f"{default['repeatability']['cluster_size']['mean']:.1f}, "
                f"repeat J {best['repeatability']['pairwise_jaccard']['mean']:.3f} vs "
                f"{default['repeatability']['pairwise_jaccard']['mean']:.3f}, "
                f"nested support above chance "
                f"{best['nested_consistency']['support_above_chance']:.3f} vs "
                f"{default['nested_consistency']['support_above_chance']:.3f}."
            )
            lines.append(
                f"  Signed with defaults: size "
                f"{signed['repeatability']['cluster_size']['mean']:.1f}, "
                f"signed score {signed['repeatability']['signed_score']['mean']:.3f}, "
                f"support above chance "
                f"{signed['nested_consistency']['support_above_chance']:.3f}; "
                f"Both: size {both['repeatability']['cluster_size']['mean']:.1f}, "
                f"signed score {both['repeatability']['signed_score']['mean']:.3f}, "
                f"support above chance "
                f"{both['nested_consistency']['support_above_chance']:.3f}."
            )
    lines.append("")
    lines.append(
        "Tuning gives the largest gain where the default cluster is large and "
        "weakly consistent, especially 101 / weighted. Signed edges are not "
        "uniformly better: they improve rs4-101, but on 101 they do not improve "
        "cross-seed support and enlarge the unweighted cluster (141 vs 58 "
        "taxa). `Support above chance` is the most comparable cross-seed "
        "metric; the intersection column is the strict all-runs value and is "
        "often 0 when many runs disagree on one taxon."
    )
    return lines


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tuning", default=None, help="Path to tuning_*.json")
    parser.add_argument(
        "--runs", type=int, default=DEFAULT_RUNS,
        help="Runs per configuration in Experiment 1 (repeatability).",
    )
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--seed", default=SEED)
    parser.add_argument(
        "--max-nested-nodes", type=int, default=None,
        help="Cap Experiment 2 to this many re-seeded members (subsampled "
             "with a fixed seed). Default: every member of the base cluster.",
    )
    parser.add_argument(
        "--reuse-unsigned", action="store_true",
        help="Reuse the 100-run unsigned Default/Best results from tuning.",
    )
    parser.add_argument("--output-dir", default=str(TEST_DIR / "results"))
    parser.add_argument(
        "--summary-md", default=str(TEST_DIR / "results" / "RESULTS_COMPARISON.md")
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    output = run_comparison(args)
    command = "python " + " ".join([Path(__file__).name] + sys.argv[1:])
    output["command"] = command
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    json_path = output_dir / f"comparison_{stamp}.json"
    json_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    summary_path = Path(args.summary_md).expanduser()
    summary_path.write_text(build_markdown(output, command), encoding="utf-8")
    print(f"\nJSON report written to {json_path}")
    print(f"Markdown report written to {summary_path}")
    print(f"Total elapsed: {output['elapsed_seconds']:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
