# Microbiome Co-occurrence Network Analysis

Desktop application for microbial co-occurrence network analysis from abundance tables (bacteria, fungi, or both), designed for users who want a graphical workflow instead of writing code.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Electron](https://img.shields.io/badge/Electron-44-47848F.svg)](https://www.electronjs.org/)

- Compute taxa correlations with **Spearman** or **SparCC** (compositional-aware).
- Build condition-specific co-occurrence networks and compare two conditions.
- Analyse bacteria-only, fungi-only, or cross-domain bacteria-fungi relations.
- Explore the networks interactively and run seed-based community detection.
- Export tables as CSV or JSON, and reopen or compare every previous run.
- Runs entirely offline: the analysis engine is bundled inside the installer.

## Getting started

Everything in this section is meant for the person who will actually use the app. No
programming knowledge is needed: it is a download, a double click, and a form to fill in.

### Step 1 — Download the right file

The installers are **not** in the source tree: they live on the **Releases** page, which
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
- It must contain a column named `Combination_treat2` with the condition of each
  sample (for example `Control` / `Treated`).
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
   - *Correlation method* — **Spearman** is the safe default; **SparCC** is the
     compositional-aware alternative (see [About SparCC](#about-sparcc)).
   - *Grouping mode* — `all`, or restrict to a `phylum` / `genus` level.
   - *P-value threshold* — `0.05` is the usual choice.
   - *Condition A* / *Condition B* — the two `Combination_treat2` values to compare.
4. **Start the analysis** and watch the progress bar. Every step is written to the log
   panel; the **Full log** button shows the complete output.
5. **Read the results**: a summary, the table of significant pairs, one network per
   condition, and — once you pick a taxon — the local community around it.
6. **Export** the table as CSV or JSON from the results header.

Each run is archived automatically under **History**, where you can reopen it, change
the parameters and re-run it, or compare two runs side by side.

### About SparCC

SparCC needs an extra component called **FastSpar**, which is bundled for Windows and
Linux on x64 only.

- On Windows and Linux x64 it works with no extra step.
- On macOS and on ARM Linux it is not available, so the app says so in plain language and
  runs the comparison with **Spearman** instead. The substitution is recorded in the
  analysis log, so the provenance of the numbers stays clear.
- If you really need SparCC on macOS, install FastSpar yourself and point the app at it
  in **Settings → FastSpar path**.

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
| The metadata file is rejected | The column must be named exactly `Combination_treat2`. |
| The network is empty or very sparse | Relax the p-value threshold (try `0.1`), or add samples per condition. |
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

### Publishing a release

Installers are **never committed to git**: they are 150–200 MB binaries and `release/`
is git-ignored. They belong to a GitHub *Release* instead, so users get a plain download
link and the repository stays small.

```bash
npm run build:linux     # -> release/electron-analysis-app-<version>-linux-x86_64.AppImage
npm run build:win       # -> release/electron-analysis-app-<version>-win-x64.exe
```

Then create a release for the version in `package.json`, tag it (`v1.0.0`) and attach
the files listed above as assets. Either do it from the GitHub web UI
(*Releases* → *Draft a new release* → drag the files into the assets box), or with the
GitHub CLI:

```bash
gh release create v1.0.0 \
  release/electron-analysis-app-1.0.0-win-x64.exe \
  release/electron-analysis-app-1.0.0-linux-x86_64.AppImage \
  --title "1.0.0" --notes "First public release."
```

The `blockmap` and `latest*.yml` files are only needed if you later enable automatic
updates; they can be ignored for a manual release.

Each artifact must be built **on its own platform** — see the wine note below for
Windows, and note that macOS packages can only be produced on a Mac.

### Building the Windows package on Linux/macOS

`npm run build:win` needs `wine` (on Arch: `sudo pacman -S wine`). Failures here are
almost always about the wine prefix rather than the project:

- `could not load kernel32.dll, status c0000135` — the default `~/.wine` prefix was left
  half-initialised, for example by an interrupted first run. Delete it and let it rebuild
  from scratch: `rm -rf ~/.wine && wineboot -u`.
- Wine waits forever on a dialog offering to install Mono/Gecko. electron-builder builds a
  temporary installer (~570 KB, this is not a broken artifact) and runs it under wine to
  generate the uninstaller, so the prefix must be complete before packaging. `build:win`
  sets `WINEDLLOVERRIDES=mscoree,mshtml=` to keep this non-interactive; set it yourself if
  you call `electron-builder` directly.
- Do not pipe a long build through `tail`: it hides the real error until the very end.

### Important packaging constraints

- The frozen engine is shipped via `extraResources` under `resources/python-dist`;
  the Python sources themselves are **not** included in the installer.
- FastSpar binaries are shipped per platform (`resources/bin/<platform>` →
  `bin/<platform>`), so a Windows installer does not carry the Linux binaries.
- The engine is rebuilt for the platform being packaged, so `python-dist/` never
  mixes two platforms.
- Windows packages can be produced from Linux/macOS (see the section above); macOS
  packages only from macOS.
- Build output is written to `release/`.

## License

Released under the MIT License. See [LICENSE](LICENSE).

FastSpar components are distributed under GPL-3.0.
