# Microbiome Co-occurrence Network Analysis

Desktop application for microbial co-occurrence network analysis from abundance tables (bacteria, fungi, or both), designed for users who want a graphical workflow.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Electron](https://img.shields.io/badge/Electron-44-47848F.svg)](https://www.electronjs.org/)

- Compute taxa correlations with **Spearman**, **SparCC** or **SPIEC-EASI (MB)**.
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
[About SparCC](#about-sparcc) below).

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

**Minimum data requirements.** A condition is skipped — with a message in the log — when
it has fewer than **3 samples**, or when the two conditions share fewer than **2 taxa**.
Taxa that never vary across samples are removed automatically. For SparCC, 4 or more
samples per condition are recommended.

### Step 5 — Run an analysis

1. **Start a new analysis** from the main screen.
2. **Upload the files**: drag the abundance table and the metadata file into the two
   drop areas, or use the file picker.
3. **Fill in the parameters**:
   - *Dataset type* — bacteria, fungi, or both.
   - *Correlation method* — **Spearman** is the safe default; **SparCC** and
     **SPIEC-EASI (MB)** are compositional-aware alternatives (see
     [About SparCC](#about-sparcc)).
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

### About SparCC

SparCC needs an extra component called **FastSpar**, which is bundled for Windows and
Linux on x64 only.

- On Windows and Linux x64 it works with no extra step.
- On macOS and on ARM Linux it is not available, so the app says so in plain language and
  runs the comparison with **Spearman** instead. The substitution is recorded in the
  analysis log, so the provenance of the numbers stays clear.
- If you really need SparCC on macOS, install FastSpar yourself and point the app at it
  in **Settings → FastSpar path**.

### About SPIEC-EASI (MB)

- SPIEC-EASI (MB) is implemented directly in Python in this app (no R required).
- It performs graph selection via neighborhood regression + StARS, then exports the
  selected edges as correlation-like weights.
- Practical recommendation: SPIEC-EASI MB is often used also when taxa > samples,
  but results can be less stable in very high-dimensional settings; prevalence/abundance
  filtering and taxa grouping are usually recommended.

#### Parameters currently used in this app

Current integration defaults for SPIEC-EASI MB:

- `spieceasiNLambdas`: `20`
- `spieceasiLambdaMinRatio`: `0.001`
- `spieceasiNSubsamples`: `20`
- `spieceasiStarsThreshold`: `0.05`
- `spieceasiMaxIter`: `10000`
- `spieceasiJobs`: unset by default. When unset the worker count is the number of
  **physical** CPU cores (`spieceasi.runtime.default_n_jobs()`; `logical // 2` if
  the physical count cannot be read, `1` if nothing can be read). Set it to force
  a value.
- Parallel backend: **threads** (joblib). The compiled MB core releases the GIL,
  the subsample results are aggregated in a fixed order, and results are
  bit-identical for any worker count (verified for `n_jobs = 1…16`).
- `spieceasiSubsampleFraction`: unset by default. When unset, the StARS subsample
  fraction follows the R reference (`pulsar::.ratcheck`): `0.8` for `n <= 144`
  samples, otherwise `10*sqrt(n)/n`. Set it (e.g. `0.8`) to force a value.
- Fixed pipeline settings: `sym_rule="or"`, `random_state=17`

Runtime behavior implemented by the app:

- Constant/all-zero taxa are removed before fitting.
- The run is skipped if fewer than 4 samples or fewer than 2 varying taxa remain.
- For very small and high-dimensional datasets (`n_samples <= 10` and
  `n_taxa >= 200`), if the user did not explicitly override them, defaults are
  automatically reduced to keep execution practical:
  - `spieceasiNLambdas <= 10`
  - `spieceasiNSubsamples <= 8`
  - `spieceasiMaxIter <= 3000`

#### Performance (measured)

The MB solver is compiled (`spieceasi._mb_core`, a pybind11 translation of
`huge_core.cpp`); the pure-Python reference stays as the fallback. On the
development machine (4 physical / 8 logical cores), full pipeline, median of 3:

| case | pure Python | compiled, n_jobs=1 | compiled, n_jobs=4 | R (`spiec.easi`) |
|------|-------------|--------------------|--------------------|------------------|
| p=30, n=100 | 97.8 s | 0.52 s | 0.24 s | 2.11 s |
| p=50, n=50 | 453.4 s | 3.80 s | 1.50 s | 4.45 s |
| p=100, n=80 | 1565.8 s | 21.6 s | **8.73 s** | 14.04 s |

So with 4 workers the p=100 pipeline is faster than the R reference. The
compiled core is **not** bit-identical to the pure-Python reference (differences
at the 1e-13…1e-14 level; identical discovered structure except in a deliberately
degenerate family), it performs *the same number of coordinate-descent sweeps*
as `huge`, and the remaining serial gap to R is R's OpenMP parallelism. Details,
flags and caveats: `docs/PERFORMANCE.md`; packaging: `docs/PACKAGING.md`.

### Equivalence with R ("Equivalenza con R")

The Python implementation in `spiec-easi-python/` is a translation of the R
pipeline (SpiecEasi 1.99.0 + huge 2.0.1 + pulsar 0.3.13) and is checked against
R on generated golden datasets:

- `tests/r_reference/generate_golden.R` and `generate_golden_extended.R` produce
  17-digit R references; `tests/test_equivalence*.py` compare Python against them
  step by step (CLR, `lambda.max`, lambda grid, solver path, coefficients, StARS,
  selection, final adjacency).
- R's StARS subsample indices are injected into Python, because the two RNGs are
  not bit-comparable.
- No comparison tolerance was widened to make a test pass.

Result **on the currently tested cases** (p up to 100, n between 25 and 300,
dense and sparse data, `p > n`, `n > 144`, both the library and the application
configuration): every step matches to about `1e-15`–`1e-13`, with one
numerically indeterminate situation at `lambda = lambda.max` documented in
`docs/AUDIT.md` §2.E. This is **not** a claim of general equivalence: it holds
for the tested cases only.

Known and documented deviations (solver arithmetic, StARS RNG, the application's
reduced `max_iter` cap in the tiny/high-dimensional regime) are listed in
`docs/AUDIT.md`. Measured runtimes are in `docs/PERFORMANCE.md`.

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
| The network is empty or very sparse | Relax the p-value threshold (try `0.1`), lower the correlation filter (min \|r\|), or add samples per condition. |
| A condition is skipped | Fewer than 3 samples, or fewer than 2 taxa shared with the other condition. |

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

- Electron 44
- React 18 with Vite 8
- Python 3 analysis engine, frozen into a standalone executable with PyInstaller
- NetworkX, SciPy, NumPy

### Repository layout

```text
electron-analysis-app/
|- src/
|  |- main/
|  |- renderer/
|  `- shared/
|- python/
|  |- backend.py         # single entry point (analyze | cluster)
|  |- backend.spec       # PyInstaller build description
|  |- analysis.py
|  |- requirements.txt
|  `- clustering/
|- python-dist/          # frozen engine (generated, git-ignored)
|- resources/bin/<platform>/<arch>/   # bundled FastSpar binaries
|- scripts/build-backend.mjs
|- tests/
`- package.json
```

Generated build directories (all git-ignored): `python-dist/`, `.pyinstaller-build/`,
`.pyinstaller-venv/`, `.pyinstaller-wine/`, `release/`.

### Local development setup

```bash
git clone https://github.com/sabrynico21/electron-analysis-app.git
cd electron-analysis-app
npm install
npm run dev
```

`npm run dev` uses a system Python interpreter, so install the analysis packages in a
virtual environment first:

```bash
python3 -m venv .venv
.venv/bin/pip install -r python/requirements.txt
```

### Useful scripts

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

### How distribution works

Packaged installs do not rely on the user's Python. `npm run build:backend` freezes
`python/backend.py` (plus scipy, networkx, numpy and the clustering package) into a
self-contained `analysis-backend` executable, which electron-builder ships through
`extraResources` under `resources/python-dist`. At runtime the main process prefers
that executable and only falls back to a system interpreter in a development checkout
(`getBundledBackendPath()` in `src/main/analysis/runtimeDetector.js`).

The executable must be built for the platform you are packaging. PyInstaller has no
cross-compiler, so `build:win`, `build:mac` and `build:linux` each rebuild it for their
target before packaging — **except** Windows, which is produced from a Linux or macOS host
by installing a real Windows Python inside a wine prefix (`.pyinstaller-wine/`) and
running PyInstaller there. macOS packages still require a macOS machine.

The output directory is recreated on every backend build, so `python-dist/` always
contains exactly one platform's engine: the one being packaged.

### End-to-end tests

The suites under `tests/e2e/` drive a real Electron instance over the Chrome DevTools
Protocol: they start the app exactly as a user would and assert on what is rendered,
including real drag & drop of files from disk. Nothing about them needs a container.

```bash
npm run build:renderer   # required: the harness loads dist/renderer/index.html
npm run test:e2e         # every suite (the Python-backed one is skipped)
npm run test:e2e:full    # adds the Python-backed analysis workflow
```

Building the renderer first is not optional — `tests/e2e/support/appHarness.mjs` refuses
to start when `dist/renderer/index.html` is missing.

On a machine with a desktop session the suites run as they are. On a headless Linux box
Electron still needs an X display, so wrap them in a virtual one:

```bash
sudo pacman -S xorg-server-xvfb     # Arch; on Debian/Ubuntu: apt install xvfb
xvfb-run -a npm run test:e2e
```

`npm run test:e2e:full` is the only suite that executes the Python backend, which is why
it is skipped unless asked for explicitly. The harness always uses an isolated
`--user-data-dir`, so the tests never touch your own history archive.

To inspect the app visually on a headless machine, install `xvfb` and `x11vnc`, start them
on a virtual display and run `electron .` there — the application itself needs no container.


### Important packaging constraints

- The frozen engine is shipped via `extraResources` under `resources/python-dist`;
  the Python sources themselves are **not** included in the installer.
- FastSpar binaries are shipped per platform (`resources/bin/<platform>` →
  `bin/<platform>`), so a Windows installer does not carry the Linux binaries.
- The engine is rebuilt for the platform being packaged, so `python-dist/` never
  mixes two platforms.

## License

Released under the MIT License. See [LICENSE](LICENSE).

FastSpar components are distributed under GPL-3.0.
