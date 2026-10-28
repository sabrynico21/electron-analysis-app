#!/usr/bin/env python3
"""Default and one-at-a-time parameter sweep, weighted and unweighted.

Produces ``results/RESULTS_TUNING.md``: both stability experiments for

* weighted graph (``|correlation|``) and unweighted graph (all edges = 1);
* the default parameters;
* one parameter changed at a time (``alpha``, ``b``, ``phi``).

Experiment 1 (repeatability) uses ``--runs`` runs per value.  Experiment 2
(nested consistency) is not capped by a run count: it re-runs the algorithm
from **every** member of the base cluster, so its run count is the size of
that cluster.  ``--max-nested-nodes`` exists only to cap the cost of smoke
tests.

The "best" value is the value with the highest size-corrected stability,
provided the strict acceptance rate stays at 90% or more. Those values are
consumed by ``comparison_experiment.py``.

Column definitions and interpretation are documented once in
``docs/METHODOLOGY.md``; the generated report links to it instead of
repeating the legend.
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

from network_loader import NETWORKS, load_network, resolve_network  # noqa: E402
from report_common import (  # noqa: E402
    experiment_design, md_table, num, related_files,
)
from stability_experiment import (  # noqa: E402
    DEFAULT_RUNS,
    nested_consistency_experiment,
    repeatability_experiment,
)

CONFIG = "new_parameters"
SEED = "Chloroplast"
DEFAULT = {"alpha": 0.05, "b": 10, "phi": 0.4}
SWEEP = {
    "alpha": [0.01, 0.05, 0.10],
    "b": [8, 10, 12],
    # 0.2 is the app's bar: on network 101 nothing reaches it strictly, so the
    # row exercises the ladder's step 4 (relaxed phi) and reports the rung.
    # 0.5 was dropped: 0.4 already accepts every observed cluster, so the old
    # top value never changed anything.
    "phi": [0.20, 0.30, 0.40],
}
PARAM_KEY = {
    "alpha": "teleport_alpha",
    "b": "b",
    "phi": "target_conductance",
}
WEIGHTS = (("weighted", True), ("unweighted", False))


def _overrides(parameter: str, value: float) -> dict:
    return {PARAM_KEY[parameter]: value}


def _run_pair(graph, path, weighted: bool, overrides: Optional[dict],
              runs: int, max_nested: Optional[int], jobs: int) -> dict:
    repeat = repeatability_experiment(
        graph, SEED, CONFIG,
        runs=runs, run_seed=42, jobs=jobs,
        network_path=path, overrides=overrides, weighted=weighted,
    )
    nested = nested_consistency_experiment(
        graph, SEED, CONFIG,
        jobs=jobs, network_path=path,
        max_nodes=max_nested, overrides=overrides, weighted=weighted,
    )
    return {"repeatability": repeat, "nested_consistency": nested}


def _stability_score(report: dict) -> float:
    repeat = report["repeatability"]
    nested = report["nested_consistency"]
    if repeat["conditions"]["strict_success_rate"] < 0.9:
        return float("-inf")
    return 0.5 * (
        repeat["jaccard_above_chance"]
        + nested["support_above_chance"]
    )


def _best_values(condition: dict) -> dict:
    default_entry = condition["default"]
    best = {"alpha": DEFAULT["alpha"], "b": DEFAULT["b"], "phi": DEFAULT["phi"]}
    best_score = _stability_score(default_entry)
    for parameter in ("alpha", "b"):
        for value, entry in condition["sweeps"][parameter].items():
            score = _stability_score(entry)
            if score > best_score:
                best[parameter] = int(value) if parameter == "b" else float(value)
                best_score = score
    best["score"] = best_score
    best["default_score"] = _stability_score(default_entry)
    return best


def run_sweep(args) -> dict:
    names = list(NETWORKS) if args.network == "both" else [args.network]
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "config": CONFIG,
        "seed": SEED,
        "runs": args.runs,
        "nested_max_nodes": args.max_nested_nodes,
        "default": dict(DEFAULT),
        "sweep_values": {key: list(value) for key, value in SWEEP.items()},
        "networks": {},
    }
    total_started = time.perf_counter()

    for name in names:
        path = resolve_network(name)
        graph, info = load_network(name)
        print("\n" + "=" * 78)
        print(f"Network {name}: {info['nodes']} nodes, {info['edges']} edges")
        print("=" * 78)
        network_entry = {"graph": info, "conditions": {}}

        for label, weighted in WEIGHTS:
            print(f"\n### {name} / {label}")
            condition = {"default": None, "sweeps": {}}
            condition["default"] = _run_pair(
                graph, path, weighted, None,
                args.runs, args.max_nested_nodes, args.jobs,
            )
            print("  default: "
                  f"size={condition['default']['repeatability']['cluster_size']['mean']:.1f} "
                  f"strict={condition['default']['repeatability']['conditions']['strict_success_rate']:.0%}")

            for parameter, values in SWEEP.items():
                condition["sweeps"][parameter] = {}
                for value in values:
                    if value == DEFAULT[parameter]:
                        continue
                    started = time.perf_counter()
                    entry = _run_pair(
                        graph, path, weighted, _overrides(parameter, value),
                        args.runs, args.max_nested_nodes, args.jobs,
                    )
                    condition["sweeps"][parameter][str(value)] = entry
                    print(f"  {parameter}={value}: "
                          f"size={entry['repeatability']['cluster_size']['mean']:.1f} "
                          f"J={entry['repeatability']['pairwise_jaccard']['mean']:.3f} "
                          f"support={entry['nested_consistency']['node_support_stats']['mean']:.3f} "
                          f"({time.perf_counter() - started:.0f}s)")

            condition["best"] = _best_values(condition)
            print(f"  best: alpha={condition['best']['alpha']:g} "
                  f"b={condition['best']['b']} "
                  f"score={condition['best']['score']:.3f} "
                  f"(default {condition['best']['default_score']:.3f})")
            network_entry["conditions"][label] = condition

        output["networks"][name] = network_entry

    output["elapsed_seconds"] = time.perf_counter() - total_started
    return output


def _condition_label(network: str, weight: str) -> str:
    return f"{network} / {weight}"


def _slash(values: list[str]) -> str:
    return " / ".join(values)


def _entries(condition: dict, parameter: str) -> tuple[list, list[dict]]:
    """Values and entries of one parameter row, in sweep order.

    The row renders the values the JSON actually contains: values listed
    in ``SWEEP`` come first (sweep order, default included), and any
    extra value produced by an *older* sweep follows sorted.  This keeps
    historical reports rebuildable after the sweep ranges change
    (``rebuild_reports.py`` must not fail on a JSON that holds a value
    no longer listed, nor drop it from the table).
    """
    sweeps = condition["sweeps"][parameter]
    default_value = DEFAULT[parameter]
    values = [
        value
        for value in SWEEP[parameter]
        if value == default_value or str(value) in sweeps
    ]
    values.extend(
        sorted(
            float(key)
            for key in sweeps
            if float(key) not in SWEEP[parameter]
            and float(key) != default_value
        )
    )
    entries = [
        condition["default"] if value == default_value else sweeps[str(value)]
        for value in values
    ]
    return values, entries


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
    lines = ["# Parameter tuning: default values and one-at-a-time sweep", ""]
    lines.append(f"Generated: `{output['generated_at']}`  ")
    lines.append(f"Command: `{command}`")
    nested_note = (
        f"up to **{output['nested_max_nodes']}** members of that cluster, "
        "subsampled with a fixed seed (smoke-test mode)"
        if output["nested_max_nodes"]
        else "**every member** of that cluster"
    )
    lines.append(experiment_design(output["runs"], nested_note))
    lines.append("## Defaults and sweep ranges")
    lines.append("")
    lines.append(md_table(
        ["Parameter", "Default", "Values tested"],
        [
            [name, output["default"][name],
             " / ".join(str(value) for value in output["sweep_values"][name])]
            for name in ("alpha", "b", "phi")
        ],
    ))
    lines.append("")
    lines.append(
        "Each parameter row changes one parameter at a time; the other two stay "
        "at their default. `alpha` controls how local the random walk is, `b` "
        "sets the minimum cluster volume, and `phi` is the acceptance threshold."
    )
    lines.append("")
    lines.append("## Reading note")
    lines.append("")
    lines.append(
        "Each parameter row lists its values in the order of the `Values` "
        "column; `/` separates the values, and inside a mean/min/max cell each "
        "`(mean/min/max)` group belongs to one value. The `default` row is the "
        "reference; the `alpha`, `b` and `phi` rows repeat the default value "
        "inside their sequence."
    )
    lines.append("")

    def table_rows(experiment: str):
        rows = []
        for network, net_entry in output["networks"].items():
            for weight, _ in WEIGHTS:
                condition = net_entry["conditions"][weight]
                label = _condition_label(network, weight)
                key = "repeatability" if experiment == "repeat" else "nested_consistency"
                groups = [("default", ["-"], [condition["default"][key]])]
                for parameter in ("alpha", "b", "phi"):
                    values, entries = _entries(condition, parameter)
                    groups.append((
                        parameter, values, [entry[key] for entry in entries]
                    ))
                for parameter, values, entries in groups:
                    def per(fn):
                        return _slash([fn(entry) for entry in entries])
                    common = [
                        label, parameter,
                        _slash([str(value) for value in values]),
                        # one run count per tested value: Experiment 1 uses the
                        # same fixed N everywhere, Experiment 2 re-seeds once
                        # per base-cluster member, so its count differs
                        # whenever the base cluster differs.
                        _slash(
                            [
                                str(entry["runs"])
                                if experiment == "repeat"
                                else str(entry["nodes_tested"])
                                for entry in entries
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
                        common.append(per(
                            lambda e: f"{e['asymmetry']['asymmetry_rate']:.0%}"
                        ))
                    rows.append(common)
        return rows

    lines.append(
        f"## Experiment 1 - Repeatability ({output['runs']} runs per value)"
    )
    lines.append("")
    lines.append(md_table(
        ["Condition", "Parameter", "Values", "Runs", "Mean cluster size (taxa)",
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
        ["Condition", "Parameter", "Values", "Runs", "Mean cluster size (taxa)",
         "Jaccard (mean / min / max)", "Intersection (n / %)", "Union",
         "Conductance (mean / min / max)", "Identical clusters",
         "Consensus >=80%", "Asymmetry"],
        table_rows("nested"),
    ))
    lines.append("")

    lines.append("## Best value per condition")
    lines.append("")
    lines.append(
        "This table picks, for every network and weight mode, the `alpha` and "
        "`b` that `comparison_experiment.py` will use (the `Best` and `Both` "
        "configurations of the comparison report) and compares them with the "
        "default parameters. The comparison is a single number per "
        "experiment run, the **score**, defined below."
    )
    lines.append("")
    lines.append("### How the score is calculated")
    lines.append("")
    lines.append(
        "Each swept value produces one Experiment 1 report and one Experiment 2 "
        "report. Each report is reduced to a size-corrected value in `[0, 1]` "
        "and the two are averaged:"
    )
    lines.append("")
    lines.append(
        "**score = 0.5 x (repeatability above chance + support above chance)**"
    )
    lines.append("")
    lines.append(
        "- **repeatability above chance** (Experiment 1): the mean pairwise "
        "Jaccard `J` of the runs, corrected for cluster size. With `k` the mean "
        "cluster size and `n` the taxa of the graph, two random clusters of "
        "that size overlap with `J_null = k / (2n - k)`, so the corrected value "
        "is `(J - J_null) / (1 - J_null)`. 0 = \"no more overlap than two "
        "random clusters\", 1 = \"every run returned the same cluster\"."
    )
    lines.append(
        "- **support above chance** (Experiment 2): the fraction of base-cluster "
        "members that recover a taxon, corrected for cluster size. Random "
        "re-seeds recover it with `support_null = k / n`, so the corrected "
        "value is `(support - support_null) / (1 - support_null)`. 0 = \"no "
        "better than random re-seeds\", 1 = \"every member recovered the "
        "cluster\"."
    )
    lines.append("")
    lines.append(
        "Both inputs are corrected for size on purpose: raw overlap and raw "
        "support both grow with the cluster, so the largest cluster would win "
        "any ranking. Averaging the two views balances \"the same seed repeats "
        "its cluster\" against \"every member of that cluster finds it back\": "
        "a value scores high only when it is stable under both."
    )
    lines.append("")
    lines.append("### Eligibility")
    lines.append("")
    lines.append(
        "- Only values whose Experiment 1 **strict acceptance rate** (the share "
        "of runs that met the volume, probability and conductance conditions) "
        "is at least 90% compete; anything else is ineligible and cannot win "
        "the ranking."
    )
    lines.append(
        "- `phi` is not ranked: it changes acceptance, not the cluster shape. "
        "Only `alpha` and `b` vary in the ranking, so the winning `phi` is "
        "always its default value."
    )
    lines.append("")
    best_rows = []
    for network, net_entry in output["networks"].items():
        for weight, _ in WEIGHTS:
            best = net_entry["conditions"][weight]["best"]
            best_rows.append([
                _condition_label(network, weight),
                f"alpha={best['alpha']:g}, b={best['b']}, phi={best['phi']:g}",
                num(best["score"], 3),
                num(best["default_score"], 3),
                num(best["score"] - best["default_score"], 3),
            ])
    lines.append(md_table(
        ["Condition", "Best parameters", "Best score", "Default score", "Gap"],
        best_rows,
    ))
    lines.append("")
    lines.append("### Reading the table")
    lines.append("")
    lines.append(
        "- **Best parameters**: the winning `alpha` and `b` for this "
        "condition, with `phi` at its default."
    )
    lines.append(
        "- **Best score**: the score of that winner; these are the values "
        "consumed by the comparison report."
    )
    lines.append(
        "- **Default score**: the score of the defaults "
        "(`alpha=0.05`, `b=10`, `phi=0.4`) - the reference for the same "
        "condition."
    )
    lines.append(
        "- **Gap**: `Best score - Default score`. A positive gap means the "
        "sweep beat the defaults, a gap near 0 means the defaults sit on a "
        "plateau, and `-` means no value reached the 90% eligibility bar."
    )
    lines.append("")
    if include_findings:
        lines.append("## What the sweep shows")
        lines.append("")
        lines.extend(_explanation(output))
        lines.append("")
    lines.append(related_files([
        "`results/tuning_*.json` - raw data of this sweep.",
        "[RESULTS_COMPARISON.md](RESULTS_COMPARISON.md) - the best parameters "
        "of this sweep vs the defaults, plus the signed-edge variants.",
        "[README.md](../README.md) - repository layout and the exact commands "
        "used to produce this report.",
    ]))
    return "\n".join(lines)


def _explanation(output: dict) -> list[str]:
    lines = []
    for network, net_entry in output["networks"].items():
        for weight, _ in WEIGHTS:
            condition = net_entry["conditions"][weight]
            best = condition["best"]
            default = condition["default"]
            d_size = default["repeatability"]["cluster_size"]["mean"]
            b_size = condition["sweeps"].get("alpha", {}).get(
                str(best["alpha"]), default
            )["repeatability"]["cluster_size"]["mean"]
            improved = best["score"] > best["default_score"]
            lines.append(
                f"- **{_condition_label(network, weight)}**: best "
                f"`alpha={best['alpha']:g}`, `b={best['b']}` "
                f"(score {best['score']:.3f} vs {best['default_score']:.3f}); "
                + ("an improvement. " if improved else "no improvement over the default. ")
                + f"The default cluster has {d_size:.1f} taxa, the best-alpha "
                f"cluster has {b_size:.1f}."
            )

    lines.append("")
    lines.append(
        "`alpha` is the only parameter with a consistent effect on the "
        "cluster shape. `b` and `phi` are almost flat here: changing `b` "
        "rarely moves the selected cut, and `phi` only changes acceptance. "
        "The weighted and unweighted graphs do not prefer the same direction "
        "or the same cluster size, so the table should be read per condition "
        "rather than as one global optimum. A positive gap means the default "
        "is not the best value for that condition; a small gap means the "
        "default is on a plateau."
    )
    return lines


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", choices=["101", "rs4-101", "both"], default="both")
    parser.add_argument(
        "--runs", type=int, default=DEFAULT_RUNS,
        help="Runs per value in Experiment 1 (repeatability).",
    )
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--seed", default=SEED)
    parser.add_argument(
        "--max-nested-nodes", type=int, default=None,
        help="Cap Experiment 2 to this many re-seeded members (subsampled "
             "with a fixed seed). Default: every member of the base cluster.",
    )
    parser.add_argument("--output-dir", default=str(TEST_DIR / "results"))
    parser.add_argument(
        "--summary-md", default=str(TEST_DIR / "results" / "RESULTS_TUNING.md")
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    global SEED
    SEED = args.seed
    output = run_sweep(args)
    command = "python " + " ".join([Path(__file__).name] + sys.argv[1:])
    output["command"] = command
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    json_path = output_dir / f"tuning_{stamp}.json"
    json_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    summary_path = Path(args.summary_md).expanduser()
    summary_path.write_text(build_markdown(output, command), encoding="utf-8")
    print(f"\nJSON report written to {json_path}")
    print(f"Markdown report written to {summary_path}")
    print(f"Total elapsed: {output['elapsed_seconds']:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
