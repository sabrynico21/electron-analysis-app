# Packaging & Electron integration

> Status: skeleton. JSON contract and PyInstaller steps are finalized in a later step.

## CLI

```bash
python -m cclasso < request.json > response.json
cclasso --input request.json --output response.json
```

- Reads a JSON request from stdin (or `--input`).
- Writes a JSON response to stdout (or `--output`). **Nothing else** is written to stdout.
- Errors/warnings go to **stderr** with a non-zero exit code on failure.

## Request schema

```json
{
  "x": [[...], [...]],
  "counts": true,
  "pseudo": 0.5,
  "k_cv": 3,
  "lam_int": [0.0001, 1.0],
  "k_max": 20,
  "n_boot": 20,
  "seed": 101
}
```

## Response schema

```json
{
  "var_w": [...],
  "cor_w": [[...]],
  "p_vals": [[...]],
  "lambda": 0.5632,
  "info_cv": {"lams": [...], "fvals": [...], "k": 8, "lam_int": [0.43, 0.72]}
}
```

## Node.js example

```js
const { spawn } = require("node:child_process");

function runCclasso(request) {
  return new Promise((resolve, reject) => {
    const py = spawn("python", ["-m", "cclasso"], { stdio: ["pipe", "pipe", "pipe"] });
    let out = "", err = "";
    py.stdout.on("data", (d) => (out += d));
    py.stderr.on("data", (d) => (err += d));
    py.on("close", (code) =>
      code === 0 ? resolve(JSON.parse(out)) : reject(new Error(err || `exit ${code}`)));
    py.stdin.write(JSON.stringify(request));
    py.stdin.end();
  });
}
```

## Bundling with PyInstaller

The CLI is a single-process, stdin→stdout program, so it bundles cleanly:

```bash
pip install pyinstaller
pyinstaller --onefile --name cclasso-backend \
  --collect-submodules scipy \
  --paths src \
  src/cclasso/__main__.py
# -> dist/cclasso-backend
```

Notes:
- `numpy` and `scipy` pull in compiled extensions; `--collect-submodules scipy`
  avoids "hidden import" failures.
- Ship `dist/cclasso-backend` in the Electron `extraResources` and spawn it the same
  way as the `python -m cclasso` example above (omit the `-m cclasso` args).
- Keep a **pure-`numpy`** fallback path in mind if you must avoid SciPy: only the
  Student-t CDF uses SciPy (`scipy.stats.t.cdf`), and it can be replaced by an
  equivalent implementation if the SciPy wheel size is a problem.
- On Windows, name the binary `cclasso-backend.exe`.
- The golden fixtures and tests are **not** shipped; only the package is packaged.

## Startup performance for Electron

- Cold start is dominated by importing numpy/scipy (~0.3–0.6 s). If the app calls
  `cclasso` repeatedly, keep the backend process alive and exchange multiple JSON
  messages, or accept the per-call process cost for simplicity.

