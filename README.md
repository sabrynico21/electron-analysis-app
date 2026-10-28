# Microbiome Co-occurrence Network Analysis

Desktop application for microbial co-occurrence network analysis from abundance tables (bacteria, fungi, or both), designed for users who want a graphical workflow.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Electron](https://img.shields.io/badge/Electron-44-47848F.svg)](https://www.electronjs.org/)

- Compute taxa associations with **Spearman**, **SparCC**, **CCLasso**,
  **SPIEC-EASI (MB)** or **SPRING**.
- Build condition-specific co-occurrence networks and compare two conditions.
- Analyse bacteria-only, fungi-only, or cross-domain bacteria-fungi relations.
- Explore the networks interactively and run seed-based community detection.
- Export tables as CSV or JSON, and reopen or compare every previous run.
- Runs entirely offline: the analysis engine is bundled inside the installer.

## Getting started

### Step 1 — Download the right file

The installers live on the **Releases** page, which
is the only place you need to look.

**https://github.com/sabrynico21/electron-analysis-app/releases/latest**

Open it and pick the file for your computer from the *Assets* list. **Nothing else has to
be installed** — the analysis engine travels inside the package.

| Your computer | Download this file |
| --- | --- |
| Windows 10/11 (64-bit) | `electron-analysis-app-<version>-win-x64.exe` |
| Mac with Apple chip (M1/M2/M3/M4) | `electron-analysis-app-<version>-mac-arm64.dmg` |
| Mac with Intel chip | `electron-analysis-app-<version>-mac-x64.dmg` |
| Linux (64-bit) | `electron-analysis-app-<version>-linux-x86_64.AppImage` |

> **macOS is not published yet.** Only the Windows and Linux builds exist at the moment,
> because a `.dmg` can only be produced on a Mac. The `mac-*.dmg` files will appear on
> the Releases page once that build is made.

Not sure which Mac you have? Open the Apple menu → *About This Mac*. If the chip line
starts with "Apple M", take the `arm64` file, otherwise take `x64`.

You do **not** need Python, `pip`, Anaconda, R, an administrator account, or an internet
connection after the download: the entire analysis engine is inside the file.

The download is large (roughly 160–200 MB) because the scientific computing libraries
are embedded. That is exactly what makes the installation effortless.

### Step 2 — Install it

#### Windows
1. Double-click the `.exe` file you downloaded.
2. Windows may show a blue **"Windows protected your PC"** box. This is normal for an
   app without a commercial signing certificate. Click **More info** → **Run anyway**.
3. Follow the wizard: choose a folder (or keep the default) and click **Install**.
4. Launch **Analysis App** from the Start menu.

#### macOS
1. Double-click the `.dmg` file.
2. Drag the **Analysis App** icon onto the **Applications** folder.
3. The first time, open it with **right-click → Open** (not a plain double click) and
   confirm **Open**. macOS shows this warning for apps that are not notarized; you only
   need to do it once.
4. On an Apple Silicon Mac, if macOS offers to install *Rosetta*, you downloaded the
   `x64` file by mistake — get the `arm64` one instead.

#### Linux
1. Make the file executable: right-click it → **Properties** → **Permissions** → tick
   *Allow executing file as program* (or run `chmod +x` on it once).
2. Double-click to run. It is a single portable file and installs nothing system-wide.
3. If it refuses to start with a sandbox error, run it from a terminal as
   `./electron-analysis-app-<version>-linux-x86_64.AppImage --no-sandbox`.

### Step 3 — First launch

The window opens straight on the analysis screen. There is no setup wizard, nothing to
configure, no account to create and no dependency to download:

- no Python to install or point at;
- no internet connection required;
- no administrator rights required.

The only optional settings live under **Settings**, and you can ignore them unless the
app explicitly asks you to point it at a FastSpar binary (see
[SparCC](#sparcc) below).

### Step 4 — Prepare your two input files

You need one **abundance table** and one **metadata** file.

**Abundance table** (CSV or TSV)

- One row per sample, one column per taxon.
- The first column holds the sample names, the header row holds the taxon names.
- The transposed layout (taxa on rows, samples on columns) is detected automatically.
- Commas and tabs both work: the delimiter is detected automatically.

```text
sample   OTU1   OTU2   OTU3
S1       10     2      12
S2       11     1      13
S3       12     3      11
```

**Metadata** (CSV or TSV)

- The first column holds the sample names, spelled exactly as in the abundance table.
- It must contain a column with the condition of each sample (for example
  `Control` / `Treated`). By default this is `Combination_treat2`, but you can pick any
  other column from the **Condition column** parameter (see below).
- At least **two different values** are required, because the analysis compares two
  conditions.

```text
sample_id   Combination_treat2
S1          Control
S2          Control
S3          Treated
```

**Minimum data requirements.** Every method needs at least **2 taxa** that vary across
samples; taxa that never vary are removed automatically. The per-condition sample
minimum depends on the method: **3** for Spearman and SparCC, **4** for SPIEC-EASI (MB),
**6** for CCLasso and SPRING (see [Correlation methods](#correlation-methods)). Below its
minimum, CCLasso or SPRING stops the run with an explicit message; the other methods log
the shortfall and produce an empty network. For SparCC, 4 or more samples per condition
are recommended.

### Step 5 — Run an analysis

1. **Start a new analysis** from the main screen.
2. **Upload the files**: drag the abundance table and the metadata file into the two
   drop areas, or use the file picker.
3. **Fill in the parameters**:
   - *Dataset type* — bacteria, fungi, or both.
   - *Correlation method* — **Spearman** is the safe default; **SparCC**,
     **CCLasso**, **SPIEC-EASI (MB)** and **SPRING** are compositional-aware
     alternatives (see [Correlation methods](#correlation-methods)).
   - *Grouping mode* — `all`, or restrict to a `phylum` / `genus` level.
   - *P-value threshold* — `0.05` is the usual choice.
   - *Correlation filter · min |r|* — `0` keeps every significant edge; raise it to
     keep only strong edges. It works on the **absolute** value, so it filters weak
     positive and weak negative edges at the same time.
   - *Condition column* — the metadata column the conditions are read from. It lists
     the columns shared by all loaded metadata files; pick the one holding the
     experimental condition (defaults to `Combination_treat2`).
   - *Condition A* / *Condition B* — the two values of the chosen condition column to
     compare.
4. **Start the analysis** and watch the progress bar. Every step is written to the log
   panel; the **Full log** button shows the complete output.
5. **Read the results**: a summary, the table of significant pairs, one network per
   condition, and — once you pick a taxon — the local community around it.
6. **Export** the table as CSV or JSON from the results header.

Each run is archived automatically under **History**, where you can reopen it, change
the parameters and re-run it, or compare two runs side by side.

#### Downloadable tables

Under the cluster graphs (both in a single analysis and in the comparison view) each
figure has a **Download CSV** button and a **Show table** toggle that reveals a preview
below the button. The available tables are:

| Table | What it contains |
| --- | --- |
| Full condition graph | every edge of the condition graph with its correlation value |
| Cluster | the edges of the cluster subgraph with the correlation value of each edge |
| Intersection | edges shared by the two clusters, with the correlation value **of each condition** (an edge with opposite signs keeps both values) |
| Complement | edges present in exactly one cluster |
| Only in … | edges exclusive to one condition |
| Nodes present in only one condition | for each node in exactly one cluster, whether it is missing only from the other *cluster* or from the other condition *graph* entirely |

#### How the comparison graphs are built

- **Intersection** — two nodes are connected when the same pair appears in both
  clusters. If the correlation has the **same sign** in both conditions the edge is a
  solid line (red = positive, blue = negative); if the sign is **opposite** the edge is
  drawn **dashed**, with a colour that encodes the direction of the flip:
  - **dark red** = positive in the first condition, negative in the second;
  - **dark blue** = negative in the first condition, positive in the second.

  A small legend under the graph explains all four cases.
- **Only in <condition>** — the nodes that belong to one cluster but not to the other,
  plus the edges between those exclusive nodes. This is the same set as the
  *"Without intersection"* switch on each cluster card, so it is downloadable right
  below the cluster table.
- **Complement** — the union of the two clusters minus the intersection (on node
  pairs): nodes shared by both clusters are dropped, and so are the edges attached to
  them.

### Correlation methods

All five methods receive the same abundance table; each one applies its own
compositional model internally.

| Method | Edge weight | Edge evidence | Min. samples per condition | Extra component |
| --- | --- | --- | --- | --- |
| Spearman | rank correlation | pairwise p-value | 3 per pair | — |
| SparCC | correlation | permutation p-value | 3 (4+ recommended) | FastSpar, bundled for Windows/Linux x64 |
| CCLasso | correlation | bootstrap p-value | 6 | — |
| SPIEC-EASI (MB) | StARS selection frequency | StARS (no p-values) | 4 | — |
| SPRING | latent correlation | StARS (no p-values) | 6 | — |

The *P-value threshold* parameter appears only for the methods that produce pairwise
p-values (Spearman, SparCC, CCLasso). SPIEC-EASI and SPRING select edges by StARS
stability, so they use the *Correlation filter · min |r|* control instead — for
SPIEC-EASI the same control is labelled *StARS edge stability · minimum*, because its
edge weight is a selection frequency rather than a correlation.

#### SparCC

SparCC needs an extra component called **FastSpar**, which is bundled for Windows and
Linux on x64 only.

- On Windows and Linux x64 it works with no extra step.
- On macOS and on ARM Linux it is not available, so the app says so in plain language and
  runs the comparison with **Spearman** instead. The substitution is recorded in the
  analysis log, so the provenance of the numbers stays clear.
- If you really need SparCC on macOS, install FastSpar yourself and point the app at it
  in **Settings → FastSpar path**.

#### SPIEC-EASI (MB)

- SPIEC-EASI (MB) is implemented directly in Python in this app (no R required); the
  port lives in `spiec-easi-python/`.
- It performs graph selection via neighborhood regression + StARS, then exports the
  selected edges as correlation-like weights.
- Practical recommendation: SPIEC-EASI MB is often used also when taxa > samples,
  but results can be less stable in very high-dimensional settings; prevalence/abundance
  filtering and taxa grouping are usually recommended.

##### Parameters used by this app

Integration parameters for SPIEC-EASI (MB), forwarded to the standalone
`spieceasi` library (`spiec-easi-python/`):

- `spieceasiNLambdas`: `20` (R default)
- `spieceasiLambdaMinRatio`: `0.001` (R default)
- `spieceasiNSubsamples`: `20` (StARS `rep.num`, R default)
- `spieceasiStarsThreshold`: `0.05` (R default)
- `spieceasiSubsampleFraction`: unset by default. When unset the StARS subsample
  ratio follows the R reference (`pulsar::.ratcheck`): `0.8` for `n <= 144`
  samples, otherwise `10*sqrt(n)/n`. Set it (e.g. `0.8`) to force a value.
- `spieceasiJobs`: worker processes for the StARS subsample loop. Unset → about
  half the logical CPU count (`os.cpu_count() // 2`). The subsamples are
  independent and aggregated in a fixed order, so results are identical for any
  worker count.
- Fixed settings: `method="mb"`, `sel_criterion="stars"`, `sym="or"`, `seed=17`.

> `spieceasiMaxIter` is no longer used. The port follows R's fixed MB solver
> iteration limit (`MAX_ITER = 10000`); the previous per-case cap was a deviation
> from R and has been removed.

Runtime behavior implemented by the app:

- Constant/all-zero taxa are removed before fitting.
- SPIEC-EASI is skipped (empty network) if fewer than 4 samples or fewer than 2
  varying taxa remain.
- For very small and high-dimensional datasets (`n_samples <= 10` and
  `n_taxa >= 200`), if the user did not explicitly override them, defaults are
  reduced to keep execution practical: `spieceasiNLambdas <= 10`,
  `spieceasiNSubsamples <= 8`.

##### Performance

The MB solver is a pure-Python port of `huge::mb` (exact coordinate descent);
there is no compiled extension. Wall-clock cost is dominated by the number of
`nlambda × rep_num × n_taxa` neighbourhood Lasso fits, and the StARS subsample
loop is parallelised across `spieceasiJobs` processes. See
`spiec-easi-python/docs/PACKAGING.md` for freezing notes.

#### CCLasso

- Python port of the R CCLasso 2.0 implementation, in `cclasso-python/`; no R and
  no extra binary are needed.
- Lasso-based correlation inference for compositional data. The app feeds it the same
  abundance table as every other method; the port applies CCLasso's own 0.5 pseudo-count
  closure internally.
- Edge weights are correlations and the p-values come from CCLasso's bootstrap +
  artificial-null procedure, so the usual *P-value threshold* applies.
- Fixed defaults: 3 cross-validation folds, 20 bootstraps, seed 17.
- Solver diagnostics (ADMM calls at the iteration cap, lambda at the search boundary)
  are reported as notes in the log, not as errors.

#### SPRING

- Python port of the modernised R SPRING fork (the `latentcor`-based one), in
  `spring-python/`; no R and no extra binary are needed.
- Truncated Gaussian-copula latent correlation (Kendall's τ) on MCLR-transformed data,
  followed by Meinshausen–Bühlmann neighbourhood selection with StARS.
- Edge weights are latent correlations; like SPIEC-EASI, SPRING selects edges by StARS
  stability and produces no pairwise p-values, so edges are filtered by
  *Correlation filter · min |r|*.
- Fixed defaults: 20 lambdas, 20 StARS subsamples, StARS threshold 0.1, seed 17.
  All-zero samples are removed before fitting.
- The 2019 paper derives the latent-correlation step from `mixedCCA`; the modern R
  fork — and therefore this port — uses `latentcor`. Cite results as the
  `latentcor`-based implementation.

### Verification against the R references

Each of the three Python ports — CCLasso, SPIEC-EASI (MB) and SPRING — is checked
step by step against its R reference (CCLasso 2.0; SpiecEasi 1.99.0 + huge 2.0.1 +
pulsar 0.3.13; the `latentcor`-based SPRING fork) on committed golden datasets (QMP,
American Gut Round 1, GlobalPatterns). R's random draws are saved and injected into
Python, so RNG differences cannot mask algorithm differences, and no comparison
tolerance was widened to make a test pass.

What is identical, what is identical within a numerical tolerance, and what cannot be
is described in `docs/NUMERICAL_EQUIVALENCE.md` and in each port's
`docs/EQUIVALENCE.md`; known deviations from R are in each port's
`docs/DIFFERENCES.md`. The results were **measured on Linux x86-64 (GCC)** only:
different compilers, libm implementations and BLAS builds can change the last bits of
floating-point results. How to regenerate and checksum the references is documented in
the [developer notes](#developer-notes).

### Where your results live

Everything stays on your computer, inside the application data folder:

| System | Folder |
| --- | --- |
| Windows | `%APPDATA%\electron-analysis-app` |
| macOS | `~/Library/Application Support/electron-analysis-app` |
| Linux | `~/.config/electron-analysis-app` |

Analysis results, the graph cache and the history archive live there. Nothing is
uploaded anywhere; copy that folder to back up your history.

### If something goes wrong

| What you see | What it means / what to do |
| --- | --- |
| A blue **"Windows protected your PC"** box | Expected for an unsigned app. *More info* → *Run anyway*. |
| macOS refuses to open the app | Right-click the app → *Open* → *Open* (first time only). |
| "SparCC unavailable" warning | Expected on macOS/ARM. The run continues with Spearman. |
| The metadata file is rejected | The selected *Condition column* must exist in every metadata file. It defaults to `Combination_treat2`. |
| The network is empty or very sparse | Relax the p-value threshold (try `0.1`), lower the correlation filter (min \|r\|), or add samples per condition. StARS methods (SPIEC-EASI, SPRING) deliberately keep sparse networks. |
| CCLasso or SPRING stops with "needs at least 6 samples" | Those methods require at least 6 samples per condition. Use Spearman/SparCC (3+) or SPIEC-EASI (4+) on smaller datasets. |
| A SPIEC-EASI or SPRING run yields no edges | StARS found no stable edge set at this sample size; add samples or relax the filter (see above). |

## Current limitations

- **macOS installers are not published yet.** Windows (x64) and Linux (x64) are
  available; a macOS build has to be produced on a Mac.
- FastSpar (SparCC) binaries ship for Windows and Linux on x64 only. On macOS and on
  ARM Linux the app detects this, warns you, and runs the comparison with Spearman
  correlation instead; the substitution is written to the analysis log.
- Multiple-testing correction is not applied automatically.
- Export formats are CSV and JSON only.
- The installers are not code-signed, so Windows SmartScreen and macOS Gatekeeper warn
  on the first launch (see [Step 2](#step-2--install-it)).

## Privacy

The app runs locally. Analysis data is not sent to cloud services.

---

## Developer notes

### Stack

- Electron 44 + React 18 (Vite 8) desktop shell
- Python 3 analysis engine, frozen with PyInstaller; NumPy, SciPy, NetworkX,
  scikit-learn

### Repository layout

```text
electron-analysis-app/
|- src/                    # Electron app: main / renderer / shared
|- python/                 # analysis engine
|  |- backend.py           # single entry point (analyze | cluster)
|  |- analysis.py          # per-method correlation pipeline
|  |- backend.spec         # PyInstaller build description
|  `- clustering/          # PageRank push + conductance-sweep cluster selection
|- cclasso-python/         # Python port of R CCLasso 2.0           (verified vs R; git submodule)
|- spiec-easi-python/      # Python port of SpiecEasi MB + pulsar   (verified vs R; git submodule)
|- spring-python/          # Python port of the SPRING fork         (verified vs R; git submodule)
|- resources/bin/<platform>/<arch>/   # bundled FastSpar binaries
|- scripts/                # build, environment-manifest, equivalence tooling
|- tests/                  # Electron unit + e2e suites
|- docs/                   # equivalence contract, dataset provenance
`- package.json
```

Each method port is a standalone source package **and a standalone git
repository**, integrated here as a git submodule (`cclasso-python`,
`spiec-easi-python`, `spring-python`): the app pins one commit of each port,
while every port keeps its own history, tags and releases. The engine imports
each port from its `src/` directory and PyInstaller freezes it into the backend
(`python/backend.spec`), so packaged releases are self-contained regardless of
the submodule layout. The upstream R reference sources are kept as read-only
provenance copies — `CCLasso-master/` and `SPRING-master/` at the repository
root (git-ignored), and the SpiecEasi/huge/pulsar sources vendored in
`spiec-easi-python/reference/` — as are the generated build directories
(`dist/`, `python-dist/`, `.pyinstaller-build/`, `.pyinstaller-venv/`,
`.pyinstaller-wine/`, `release/`).

### Working on a port

```bash
git submodule update --init            # first clone, or after git pull
cd spiec-easi-python                   # work on the port in its own repo
git checkout -b my-fix && ... && git commit
cd .. && git add spiec-easi-python     # advance the pinned commit in the app
git commit -m "spiec-easi-python: <what changed>"
```

Each port lives at `github.com/sabrynico21/<port>` and is released
independently; the submodule URLs in `.gitmodules` are relative, so they
resolve against wherever the app repository is cloned from.

### Statistical method ports and R equivalence

| Port | R reference | App method |
| --- | --- | --- |
| `cclasso-python/` | CCLasso 2.0 | CCLasso |
| `spiec-easi-python/` | SpiecEasi 1.99.0 + huge 2.0.1 + pulsar 0.3.13 | SPIEC-EASI (MB) |
| `spring-python/` | SPRING fork on `latentcor` + huge 2.0.1 | SPRING |

- Golden fixtures are generated from R by `<port>/tests/r_reference/generate_golden.R`
  and committed under `<port>/tests/golden/`; the step-by-step suites are
  `<port>/tests/test_equivalence.py` (run with `pytest`; the slow R-equivalence tests
  carry the `slow` marker).
- R's random draws are saved and injected into Python — the two RNGs are not
  bit-comparable — so only deterministic algorithm behaviour is compared.
- Tolerances were declared in advance and never widened to make a test pass.
- The comparison contract is `docs/NUMERICAL_EQUIVALENCE.md`; per-port guarantees are
  in `<port>/docs/EQUIVALENCE.md` and known deviations from R in
  `<port>/docs/DIFFERENCES.md`.
- `bash scripts/reproduce_equivalence.sh <new-output-dir>` regenerates the R
  references, runs all three suites, records R/Python/BLAS manifests and checksums the
  output. Set `PY_CCLASSO`, `PY_SPIEC` and `PY_SPRING` if the ports run in different
  Python environments; a non-empty output directory is never overwritten.
- Equivalence has been measured on Linux x86-64 (GCC) only.

### Local development

```bash
git clone --recurse-submodules https://github.com/sabrynico21/electron-analysis-app.git
cd electron-analysis-app
npm install
python3 -m venv .venv && .venv/bin/pip install -r python/requirements.txt
npm run dev
```

(If you cloned without `--recurse-submodules`, run
`git submodule update --init` before building.)

`npm run dev` runs the engine from source through the system interpreter; the three
ports are picked up from their `src/` directories automatically.

### Scripts

```bash
npm run dev            # Vite dev server + Electron with hot reload
npm start              # run Electron against the built renderer
npm run build:backend  # freeze the Python engine into python-dist/
npm run build:backend -- --target win32  # Windows engine, cross-built through wine
npm run build          # renderer + frozen backend + installer for the current OS
npm run build:app      # repackage without rebuilding the backend (faster iteration)
npm run build:win      # Windows NSIS installer (needs wine when built on Linux/macOS)
npm run build:mac      # macOS DMG, x64 + arm64 (must run on macOS)
npm run build:linux    # Linux AppImage
npm run icons          # regenerate the application icons
npm run test           # unit tests
npm run test:e2e       # end-to-end suites against a real Electron instance
npm run lint           # ESLint
```

### Tests

- `npm run test` — JS unit tests; `npm run lint` — ESLint. The Python ports have their
  own `pytest` suites inside each port directory.
- `tests/e2e/` drives a real Electron over the Chrome DevTools Protocol, including real
  drag & drop from disk — no container needed. The renderer must be built first
  (`npm run build:renderer`): the harness refuses to start without
  `dist/renderer/index.html`. `npm run test:e2e` skips the Python-backed workflow;
  `npm run test:e2e:full` runs it. The harness always uses an isolated
  `--user-data-dir`, so your own history archive is never touched.
- Headless Linux: wrap the command in `xvfb-run -a …` (Electron still needs an X
  display).

### Distribution

Packaged installs never use the user's Python: `npm run build:backend` freezes
`python/backend.py` — engine, clustering package and the three method ports — into a
self-contained executable that electron-builder ships via `extraResources`; at runtime
the main process prefers it and falls back to a system interpreter only in a
development checkout (`src/main/analysis/runtimeDetector.js`).

- PyInstaller cannot cross-compile: each build target rebuilds the engine for its
  platform. Windows installers are cross-built by running PyInstaller inside a wine
  prefix (`.pyinstaller-wine/`); macOS packages still require a macOS machine.
- `python-dist/` is recreated on every backend build, so it never mixes two platforms.
- Only `resources/bin/<platform>/<arch>` ships in each installer, so a Windows
  installer does not carry the Linux FastSpar binaries.
- Python sources are not shipped; the bundled engine carries the ports' licenses under
  `licenses/`.

## License

Released under the MIT License. See [LICENSE](LICENSE).

FastSpar components are distributed under GPL-3.0.
