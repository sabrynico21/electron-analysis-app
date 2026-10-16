# Licensing: upstream dependencies

This note records the licence metadata of the projects this Python code is a
derivative of. It is **informational**: it reports the facts and does not choose,
change or add any licence for this repository.

## Upstream projects

| Project | Version used as reference | Declared licence | Source |
| --- | --- | --- | --- |
| SpiecEasi | 1.99.0 | `GPL (>= 3)` | `DESCRIPTION` of the installed package and `LICENSE.md` (GNU GPL v3, 29 June 2007) |
| huge | 2.0.1 | `GPL-2` | `DESCRIPTION` of the installed package |
| pulsar | 0.3.13 | `GPL (>= 3)` | `DESCRIPTION` of the installed package |

Verified with:

```r
Rlib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", "4.6")
read.dcf(file.path(Rlib, "SpiecEasi", "DESCRIPTION"), fields = "License")
read.dcf(file.path(Rlib, "huge",     "DESCRIPTION"), fields = "License")
read.dcf(file.path(Rlib, "pulsar",   "DESCRIPTION"), fields = "License")
```

## This repository

| Artifact | Declared licence |
| --- | --- |
| `LICENSE` (repository root) | MIT License |
| `package.json` | `"license": "MIT"` |
| `spiec-easi-python/pyproject.toml` | `license = { text = "MIT" }`, classifier `MIT License` |

## What the Python code derives from

- `spiec-easi-python/src/spieceasi/mb_solver.py` re-implements the Meinshausen–
  Bühlmann graph estimation as performed by **huge** (`huge_core.cpp`, MB
  coordinate descent, warm-started lambda path).
- `spiec-easi-python/src/spieceasi/stars.py` re-implements the StARS selection
  performed by **pulsar** (`.ratcheck`, `stars.stability`, `.starsind`).
- `spiec-easi-python/src/spieceasi/pipeline.py` re-implements the MB pipeline of
  **SpiecEasi** (`spiec.easi`, `sparseiCov`).

## Compatibility note (reported, not resolved here)

The three upstream projects are distributed under the GNU GPL (GPL-2 for `huge`,
GPL-3-or-later for `SpiecEasi` and `pulsar`). The Python code in this repository
is a re-implementation / direct translation of those algorithms and the
repository currently declares the MIT licence, which is not compatible with the
GPL for a derivative work.

This document does not make any licensing decision and no licence text was added
or changed. The discrepancy is reported here so it can be resolved deliberately.
