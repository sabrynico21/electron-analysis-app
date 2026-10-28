# spring-microbiome

A faithful Python port of **SPRING** — *Semi-Parametric Rank-based approach for
INference in Graphical models* (Yoon, Gaynanova & Müller, 2019,
[*Frontiers in Genetics* 10:516](https://www.frontiersin.org/articles/10.3389/fgene.2019.00516/full)).

SPRING estimates sparse microbial association networks from (quantitative or
compositional) microbiome count data using a **truncated Gaussian copula**
latent-correlation estimate (Kendall τ based), **Meinshausen–Bühlmann**
neighborhood selection, and **StARS** stability selection.

The reference implementation is the R package SPRING
(<https://github.com/IrinaStatsLab/SPRING>; the modernised fork that uses
[`latentcor`](https://cran.r-project.org/package=latentcor) and `huge` 2.0.1).
This port aims to be **numerically equivalent** to that R code, not merely
"similar". See [`docs/EQUIVALENCE.md`](docs/EQUIVALENCE.md) for the exact
guarantees and [`docs/PORTING_NOTES.md`](docs/PORTING_NOTES.md) for the
function-by-function analysis.

> **Divergence from the SPRING paper.** The 2019 *Frontiers in Genetics* paper
> describes the latent-correlation step via `mixedCCA`. The modernised R fork —
> and therefore this port — uses `latentcor` instead. We implement the **code**
> (`latentcor`, `types="tru"`), not the paper's `mixedCCA` behaviour. This is the
> single largest paper/code divergence; it is stated up front and analysed in
> [`docs/DIFFERENCES.md`](docs/DIFFERENCES.md) §1. Results obtained with this
> package should be cited as the `latentcor`-based fork, not as the original
> `mixedCCA` method.

> **Status:** stable, scientific release `1.0.0` (2026-10-20). The equivalence
> test-suite replays R golden fixtures committed under `tests/golden/`, so
> **running the tests never requires R**. Regenerating the fixtures does not
> require an external SPRING checkout either: the R source, `QMP.rda` and the
> public-data selector are vendored under `tests/r_reference/reference/`.

## Installation

```bash
pip install spring-microbiome
```

From a checkout:

```bash
python -m venv .venv
.venv/bin/pip install -e ".[dev]"
```

Runtime dependencies are only **NumPy** and **SciPy**. The Lasso solver is a
faithful re-implementation of `huge`'s coordinate descent — **scikit-learn is
deliberately not used** (different standardisation and stopping rules).

## Usage from Python

```python
import numpy as np
from spring_microbiome import spring, mclr

# counts: n samples x p taxa (non-negative)
result = spring(counts, quantitative=False, nlambda=20, rep_num=20, seed=10010)

adj = result.adjacency  # p x p 0/1 adjacency at the optimal lambda
opt = result.output.stars.opt_index
lam = result.lambdaseq[opt - 1]
curve = result.output.stars.summary  # StARS instability curve

# The mclr transformation is exposed separately:
clr = mclr(counts)
```

### Reproducibility

StARS draws random subsamples. R's RNG and NumPy's differ, so **identical**
results across R and Python require injecting the same indices:

```python
result = spring(counts, subsample_indices=my_indices)  # list of 0-based arrays
```

Without injection the result is *statistically* equivalent, but the selected
lambda may differ. When `subsample_indices` is `None`, indices are drawn with
`numpy.random.default_rng(seed)`.

## Usage from the CLI / Electron

The CLI reads a JSON object from a file or stdin and writes a JSON object to
stdout. **Nothing but JSON goes to stdout**; errors go to stderr with a non-zero
exit code.

```bash
python -m spring_microbiome --input request.json --output results.json
echo '{"data": [[1,2,3],[4,5,6]], "nlambda": 5}' | python -m spring_microbiome
```

### Request format

```jsonc
{
  "data": [[...], ...],          // required: n x p count matrix (also accepted as "counts")
  "quantitative": false,         // false -> apply mclr first
  "nlambda": 20,
  "lambda_min_ratio": 0.01,
  "lambdaseq": null,             // null | "data-specific" | [numbers...]
  "seed": 10010,
  "thresh": 0.1,
  "subsample_ratio": 0.8,
  "rep_num": 20,
  "Rmethod": "approx",           // "approx" | "original"
  "use_nearPD": true,
  "nu": 0.001,
  "ratio": 0.9,
  "Rtol": 1e-6,
  "subsample_indices": null      // optional: list of 0-based index arrays
}
```

### Response format

```jsonc
{
  "spring_version": "1.0.0",
  "n": 50, "p": 20,
  "lambdaseq": [...],
  "opt_index": 2,
  "lambda": 0.42,
  "instability": [...],          // StARS instability curve (cummax)
  "selection_frequency": [[...]],// edge selection probability at opt_index
  "adjacency": [[...]],          // p x p 0/1 adjacency at opt_index
  "beta": [[...]],               // MB coefficients at opt_index (signed)
  "weights": [[...]],            // (beta + beta.T)/2 masked by adjacency
  "n_edges": 11,
  "params": {...}
}
```

### Calling it from Node

```js
const { spawn } = require("node:child_process");

function runSpring(payload) {
  return new Promise((resolve, reject) => {
    const py = spawn("python", ["-m", "spring_microbiome"], {
      stdio: ["pipe", "pipe", "pipe"],
    });
    let out = "";
    let err = "";
    py.stdout.on("data", (d) => (out += d));
    py.stderr.on("data", (d) => (err += d));
    py.on("close", (code) => {
      if (code !== 0) return reject(new Error(`spring failed: ${err}`));
      try { resolve(JSON.parse(out)); }
      catch (e) { reject(new Error(`bad JSON: ${err || e.message}`)); }
    });
    py.stdin.end(JSON.stringify(payload));
  });
}

// const res = await runSpring({ data: counts, nlambda: 20, rep_num: 20, seed: 10010 });
```

See [`docs/PACKAGING.md`](docs/PACKAGING.md) for bundling the CLI with Electron
(e.g. PyInstaller).

## Testing

```bash
.venv/bin/python -m pytest              # fast unit + property tests
.venv/bin/python -m pytest -m slow      # R-equivalence tests against the goldens
.venv/bin/python -m ruff check . && .venv/bin/python -m ruff format --check .
```

Golden fixtures are regenerated with R — see
[`tests/r_reference/README.md`](tests/r_reference/README.md).

## Documentation

| File | Contents |
|---|---|
| [`docs/PORTING_NOTES.md`](docs/PORTING_NOTES.md) | Full R→Python analysis, dependency map, uncertainty register |
| [`docs/EQUIVALENCE.md`](docs/EQUIVALENCE.md) | What is identical, what is identical only within tolerance, what cannot be |
| [`docs/DIFFERENCES.md`](docs/DIFFERENCES.md) | Where the code and the paper disagree (code wins) |
| [`docs/PACKAGING.md`](docs/PACKAGING.md) | Distributing with Electron |

## References

- Yoon G., Gaynanova I., Müller C.L. (2019). *Microbial networks in SPRING —
  Semi-parametric rank-based correlation and partial correlation estimation for
  quantitative microbiome data.* Frontiers in Genetics, 10:516.
- Yoon G., Müller C.L., Gaynanova I. (2021). *Fast computation of latent
  correlations.* JCGS.
- Meinshausen N., Bühlmann P. (2006). *High-dimensional graphs and variable
  selection with the lasso.* Annals of Statistics.
- Liu H., Roeder K., Wasserman L. (2010). *Stability approach to regularization
  selection (StARS).* NIPS.
- Original R implementation: <https://github.com/IrinaStatsLab/SPRING>

## License

GPL-3.0-only. SPRING is GPL-3, `latentcor` is GPL-3, `huge` is GPL-2 and
`pulsar` is GPL (≥3); a derivative work must therefore be GPL-3. See
[`LICENSE`](LICENSE).
