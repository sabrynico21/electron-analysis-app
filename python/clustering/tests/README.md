# PageRank-Nibble clustering experiments

Standalone, reproducible experiments for the clustering algorithm in
`python/clustering/`. They do not use the Electron app: each experiment loads
one of two exported edge CSVs, builds the same `BacteriaGraph` the app uses,
and calls the Python clustering modules directly, so the experiments
exercise the algorithm and only the algorithm.

## Repository layout

```
python/clustering/tests/
├── README.md                  # this file: what it is, layout, how to run
├── docs/
│   ├── METHODOLOGY.md         # the two experiments, columns, interpretation
│   ├── ALGORITHM.md           # parameters, signed edges, payload
│   └── SELECTION_LADDER.md    # the five-step selection ladder (methods spec)
├── data/                      # the exported edge lists the experiments use
│   ├── 101_Not_Treated_full_graph_edges.csv
│   └── RS4_101_Not_Treated_full_graph_edges.csv
├── results/                   # generated reports + raw JSON (see its README)
├── tuning_experiment.py       # Experiment driver 1: default + parameter sweep
├── comparison_experiment.py   # Experiment driver 2: defaults vs tuned vs signed
├── stability_experiment.py    # the two experiments, shared by both drivers
├── report_common.py           # shared report helpers
├── network_loader.py          # CSV -> BacteriaGraph
└── test_stability.py          # fast unit / smoke tests
```

## Networks

| alias | file (in `data/`) | nodes | edges |
|---|---|---|---|
| `101` | `101_Not_Treated_full_graph_edges.csv` | 387 | 5121 |
| `rs4-101` | `RS4_101_Not_Treated_full_graph_edges.csv` | 368 | 5223 |

## The two experiments

Every configuration (parameter set x weight mode) is evaluated with both
experiments, on the weighted graph (`|correlation|` edges) and the
unweighted graph (every edge has weight 1):

1. **Repeatability** - a fixed number of runs (`--runs`, 100 in the
   published reports) from the same seed taxon, changing only the
   residual-queue order.
2. **Nested consistency** - the deterministic cluster of the seed (the
   *base cluster*) is re-computed, and the algorithm is re-run **from every
   member of that cluster**; the run count equals the cluster size.

Definitions of every table column, the metrics and the interpretation rules
are documented once in
[docs/METHODOLOGY.md](docs/METHODOLOGY.md). The algorithm and its parameters
are documented in [docs/ALGORITHM.md](docs/ALGORITHM.md); the cluster
selection (the five-step ladder, its criteria and its defaults) is specified
in [docs/SELECTION_LADDER.md](docs/SELECTION_LADDER.md).

## Reports

Two reports are generated (into `results/`, together with the raw JSON they
are computed from):

1. [results/RESULTS_TUNING.md](results/RESULTS_TUNING.md) - the default
   parameters and a one-at-a-time sweep of `alpha`, `b` and `phi`.
2. [results/RESULTS_COMPARISON.md](results/RESULTS_COMPARISON.md) - the best
   parameters found by the sweep vs the defaults, plus the signed-edge
   variants (`Signed` and `Both`).

The raw `results/*.json` files predate the selection-ladder refactor: they
record the removed `recommended` profile as run provenance. They are
historical data, kept unchanged; the Markdown reports are re-rendered from
them with `rebuild_reports.py`, so only the static prose follows the current
code, never the numbers.

See [results/README.md](results/README.md) for the file inventory and the
exact regeneration commands.

## How to run

```bash
# 1. default + one-at-a-time sweep (both networks, both weight modes);
#    Experiment 1 uses 100 runs, Experiment 2 runs once per base-cluster member
.venv/bin/python python/clustering/tests/tuning_experiment.py --jobs 8

# 2. comparison; reads the newest results/tuning_*.json and reuses the
#    unsigned Default/Best results from it
.venv/bin/python python/clustering/tests/comparison_experiment.py \
  --reuse-unsigned --jobs 8

# quick smoke test (caps Experiment 2 to 5 members)
.venv/bin/python python/clustering/tests/tuning_experiment.py \
  --network 101 --runs 3 --max-nested-nodes 5 --jobs 2
```

Both scripts accept `--runs`, `--seed`, `--output-dir`, `--summary-md` and
`--max-nested-nodes` (smoke tests only; the default is every member of the
base cluster). `comparison_experiment.py` also accepts `--tuning
path/to/tuning_*.json`.

To re-render both reports from the JSON files already in `results/` **without
re-running any experiment** (tables only; add `--findings` for the
commentary):

```bash
.venv/bin/python python/clustering/tests/rebuild_reports.py
```

See [results/README.md](results/README.md) for the full options.

## Unit tests

```bash
.venv/bin/python -m unittest discover -s python/clustering/tests -p "test_*.py" -v
```
