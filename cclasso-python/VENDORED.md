# Vendored source notice

`cclasso-python/` is **vendored into the root `electron-analysis-app` repository**
so a clean checkout contains the complete, runnable Python port. It is no longer
an untracked directory or an embedded Git repository at this path.

## Why

The statistical release-readiness audit blocked gate 4 ("regenerate fixtures
from a clean checkout") because `cclasso-python/` was untracked in the root
repository. A clean clone therefore could not run the port or its golden
fixtures. Vendoring removes that dependency on the local working tree.

## What is pinned

| item | value |
|---|---|
| port | Python translation of CCLasso R **v2.0** |
| upstream R source | `CCLasso-master/R/cclasso.R` |
| pinned R copy | `tests/r_reference/reference/R/cclasso.R` |
| last standalone Git commit | `3ca924f9cc0f52c5952154d476656484af6e9fc5` (`main`) |
| vendoring date | 2026-10-20 |
| license | LGPL-3.0-or-later (see `LICENSE`) |

The vendored tree is the working state as of the vendoring date. The previous
standalone repository's history is preserved as a Git bundle in
`VENDOR_HISTORY.bundle` (created with `git bundle create ... --all`). It can be
restored with:

```bash
git clone cclasso-python/VENDOR_HISTORY.bundle cclasso-python-history
```

## What is not vendored

Build and cache artifacts remain ignored by the root `.gitignore`
(`.venv/`, `__pycache__/`, `*.egg-info/`, `.pytest_cache/`, `.ruff_cache/`).
The R source-of-truth package `CCLasso-master/` is still untracked; the port
itself only needs the pinned copy under `tests/r_reference/reference/`.

## Reproducing from a clean checkout

```bash
cd cclasso-python
python -m venv .venv
.venv/bin/pip install -e '.[test]'
.venv/bin/python -m pytest -q
```

Regenerating the R goldens additionally requires R plus `SpiecEasi` and
`phyloseq`; see `tests/r_reference/README.md`.
