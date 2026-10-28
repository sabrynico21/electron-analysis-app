# Packaging for the Electron app

The CLI is the integration surface: it reads a JSON object and writes a JSON
object, so the Electron main process only needs a path to an executable and a
`child_process.spawn` call (see `README.md`).

## 1. Development / installed package

If the app ships a Python runtime, install the package and call it directly:

```bash
pip install spring-microbiome
python -m spring_microbiome --input request.json --output result.json
```

## 2. Freezing with PyInstaller

Bundle the CLI into a single executable so the Electron app does not depend on a
system Python.

```bash
pip install pyinstaller
pyinstaller \
  --name spring-microbiome \
  --onefile \
  --collect-data spring_microbiome \
  src/spring_microbiome/__main__.py
```

`--collect-data spring_microbiome` is **required**: the multilinear
interpolation table `spring_microbiome/latentcor/data/ipol_22.npz` is package
data loaded at runtime via `importlib.resources`. Without it the interpolation
table is missing and `latentcor`/`huge_kmb` will fail at import time.

The resulting executable is `dist/spring-microbiome` (Linux/macOS) or
`dist/spring-microbiome.exe` (Windows). It behaves identically to
`python -m spring_microbiome`.

### Notes

- Prefer `--onedir` over `--onefile` if startup time matters: `--onefile`
  unpacks to a temporary directory on every launch (~0.5–2 s).
- NumPy/SciPy pull in sizeable shared libraries (OpenBLAS/LAPACK); expect a
  60–120 MB bundle. Excluding SciPy is not possible (the multivariate normal CDF
  and the bridge inversion need it).
- Build on the target OS: PyInstaller does not cross-compile.

## 3. Node-side wrapper

```js
const path = require("node:path");
const { spawn } = require("node:child_process");

const CLI = path.join(process.resourcesPath, "spring-microbiome");

function runSpring(payload) {
  return new Promise((resolve, reject) => {
    const py = spawn(CLI, [], { stdio: ["pipe", "pipe", "pipe"] });
    let out = "", err = "";
    py.stdout.on("data", (d) => (out += d));
    py.stderr.on("data", (d) => (err += d));
    py.on("close", (code) => {
      if (code !== 0) return reject(new Error(err || `exit ${code}`));
      resolve(JSON.parse(out));
    });
    py.stdin.end(JSON.stringify(payload));
  });
}
```

Add the executable to `electron-builder`'s `extraResources` so it is copied into
`resources/` next to the app.

## 4. Performance notes

- The MB solver is a pure-Python/NumPy coordinate descent. It is adequate for
  the taxon counts typical of 16S genus-level tables (tens to a few hundred
  taxa); very large `p` (≳ 500) will be slow.
- Consider caching results keyed by a hash of `(data, params)` in the Electron
  main process rather than re-running the estimator on every UI change.
- Large counts arrays should be passed as JSON numbers; if serialization becomes
  a bottleneck, send a temp file path instead of inline `data`.

## 5. Licensing

The bundle contains a GPL-3 derivative (see `LICENSE`). If you distribute the
Electron app with the bundled CLI, the combined distribution must comply with
GPL-3. Until the license is finally decided with the project owner, **do not
publish** the bundle.
