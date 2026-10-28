# Results

Everything in this directory is **generated** - do not edit it by hand.
Each experiment writes a timestamped raw JSON plus the Markdown report
computed from it, so every number in a report can be recomputed from its
JSON alone.

## Reports

| File | Content |
|---|---|
| `RESULTS_TUNING.md` | Default parameters and the one-at-a-time sweep of `alpha`, `b`, `phi`, per network and weight mode. |
| `RESULTS_COMPARISON.md` | The sweep's best parameters vs the defaults, plus the signed-edge variants (`Signed`, `Both`). |

The column definitions, the metrics and the interpretation rules are
documented once in [../docs/METHODOLOGY.md](../docs/METHODOLOGY.md);
the reports link there instead of repeating the legend.

## Raw data

| Pattern | Content |
|---|---|
| `tuning_*.json` | Full output of `tuning_experiment.py`: parameters, derived thresholds, per-run rows and aggregates of both experiments for every swept value. |
| `comparison_*.json` | Full output of `comparison_experiment.py`, same structure for the four compared configurations. |

The reports always reflect the newest JSON of their kind (by timestamp in
the file name). Only the JSON behind the currently published reports is
kept; a superseded run is deleted together with the report it produced, so
no stale copy of the same information survives a protocol change.

## Protocol

* Experiment 1 (repeatability): 100 runs per value; run seeds `42..141`
  randomise the residual-queue order (the only variation).
* Experiment 2 (nested consistency): one run per member of the base cluster
  (no run cap). It needs no run seeds: the seed node itself is the variable
  and everything else is deterministic, so re-running always reproduces the
  same numbers.
* Seed taxon: `Chloroplast`; both weight modes.

## Regenerate

From the repository root:

```bash
# 1. tuning sweep -> results/tuning_<stamp>.json + results/RESULTS_TUNING.md
.venv/bin/python python/clustering/tests/tuning_experiment.py \
  --network both --runs 100 --jobs 8 \
  --output-dir python/clustering/tests/results \
  --summary-md python/clustering/tests/results/RESULTS_TUNING.md

# 2. comparison (consumes the newest tuning_*.json) ->
#    results/comparison_<stamp>.json + results/RESULTS_COMPARISON.md
.venv/bin/python python/clustering/tests/comparison_experiment.py \
  --runs 100 --jobs 8 --reuse-unsigned \
  --output-dir python/clustering/tests/results \
  --summary-md python/clustering/tests/results/RESULTS_COMPARISON.md
```

## Rebuild the tables without re-running

To re-extract both reports from the JSON files already on disk - no
experiments, a couple of seconds - run from the repository root:

```bash
# tables + Design/Reading/score sections, NO commentary (default)
.venv/bin/python python/clustering/tests/rebuild_reports.py

# include the "What ... shows" interpretation sections too
.venv/bin/python python/clustering/tests/rebuild_reports.py --findings

# preview on stdout without writing anything / rebuild only one report
.venv/bin/python python/clustering/tests/rebuild_reports.py --dry-run
.venv/bin/python python/clustering/tests/rebuild_reports.py --only tuning
```

The `Command:` line in a rebuilt report always shows the command that
**produced the data** (read from the `command` field of the JSON), never the
rebuild invocation. The full runs above record that field automatically;
older JSONs that predate it fall back to a "rebuilt from <file>" marker.
