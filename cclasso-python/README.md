# cclasso (Python)

Python port of **CCLasso** — *Correlation inference for compositional data through Lasso*
(Fang, Huang, Zhao & Deng). This is a **faithful, line-mapped translation** of the R
implementation, **version 2.0** (`CCLasso-master/R/cclasso.R`).

> Status: core port implemented and verified against R (see `docs/EQUIVALENCE.md`).
> `lams`/`lambda` are bit-identical, ADMM iteration counts exact, outputs within ~1e-13.

## What this is

- A 1:1 algorithmic port: same parameters, same defaults, same order of operations,
  same golden-section λ search, same ADMM solver, same bootstrap + artificial-null
  p-value procedure.
- Deterministic where R is deterministic; the two sources of randomness (bootstrap
  resampling and the artificial null sample) are **injectable** so results can be
  reproduced bit-for-bit against R.
- Zero extra stdout output: an error/warning is written to **stderr**, results to stdout.

## What this is *not*

- It is **not** an "improved" or re-derived CCLasso. Optimizations are allowed only when
  numerically equivalent.

## Install

```bash
pip install -e .          # from this directory
```

Runtime deps: `numpy`, `scipy`. Dev deps: `pytest`, `pytest-cov`, `ruff`.

## Use from Python

```python
import numpy as np
from cclasso import cclasso

x = np.loadtxt("counts.csv", delimiter=",")  # n samples x p variables
res = cclasso(x, counts=True, pseudo=0.5)
res["cor_w"]  # p x p correlation estimate
res["p_vals"]  # p x p p-values
res["lambda"]  # selected lambda
```

## Use from the CLI / Electron

```bash
python -m cclasso < request.json > response.json
# or
cclasso --input request.json --output response.json
```

JSON contract and a Node.js (`child_process.spawn`) example are documented in
`docs/PACKAGING.md`.

## Documentation

| file | contents |
|---|---|
| `docs/PORTING_NOTES.md` | Step-0 R→Python mapping, primitives, uncertainties |
| `docs/DIFFERENCES.md`   | deliberate, documented deviations from R |
| `docs/EQUIVALENCE.md`   | what is guaranteed identical / within-tolerance / not |
| `docs/PACKAGING.md`     | distributing with Electron (PyInstaller), JSON contract |

## Reference

- Article: Fang H., Huang C., Zhao H., Deng M. (2015) *CCLasso: correlation inference for
  compositional data through Lasso*. Bioinformatics 31(19):3172–3180.
- Original R package: CCLasso v2.0 (GNU LGPL-v3), included here under
  `tests/r_reference/reference/` for golden-fixture generation only.

## Equivalence warning

Because R's RNG (`sample`, `rnorm`) does not match NumPy's, folds/bootstrap/null draws
produce **statistically** equivalent — not identical — results unless the indices are
injected. See `docs/EQUIVALENCE.md`.

## License

LGPL-3.0-or-later (derivative work of the LGPL-v3 R package). See `LICENSE`.
