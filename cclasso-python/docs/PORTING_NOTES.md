# CCLasso R → Python — Porting Notes (Step 0)

Status: **analysis only**. No code yet. This document is the source-of-truth map
for the Python port and lists every place where equivalence is *not* yet proven.

## 0. Source of truth

- Primary source: `CCLasso-master/R/cclasso.R` (**version 2.0**, `DESCRIPTION: Version: 2.0`,
  dated 2016-01-08). This is the version being ported.
- `CCLasso-master/R/SparCC.R` is a *separate* method (SparCC). It is **not** part of the
  `cclasso` port. See "Out of scope" below.
- `CCLasso-master/v1.0/R/cclasso.R` is a **different, incompatible** algorithm
  (different lambda grid `10^seq(0,-8,-0.01)`, plain k-fold over a fixed grid,
  `nearPD` post-processing, `kmax=5000`, `x.tol=1e-6`). It is **NOT** ported.
  Differences are recorded in `docs/DIFFERENCES.md`.
- Rule applied: whenever prose/README and code disagree, **the code wins**.

### Out of scope (pending confirmation, see Questions)
`SparCC.count` / `SparCC.frac` (needs `gtools::rdirichlet`, `Matrix`), and the entire
`v1.0/` tree.

---

## 1. Public API mapping

### R (authoritative) — `R/cclasso.R:36-38`
```r
cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3,
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20)
```
| R param | default | Python param | default | notes |
|---|---|---|---|---|
| `x` | – | `x` | – | `n × p`, samples = rows, variables = columns |
| `counts` | `FALSE` | `counts` | `False` | |
| `pseudo` | `0.5` | `pseudo` | `0.5` | only used if `counts=True` |
| `k_cv` | `3` | `k_cv` | `3` | CV folds |
| `lam_int` | `c(1e-4, 1)` | `lam_int` | `(1e-4, 1)` | **log10 interval**, length 2 |
| `k_max` | `20` | `k_max` | `20` | max golden-section iterations |
| `n_boot` | `20` | `n_boot` | `20` | bootstrap replicates |

> ⚠️ The task brief guessed `lam_min_ratio`; the R source uses **`lam_int`** (a length‑2
> log10 interval `[log10(1e‑4), log10(1)]`). Per the "same names/defaults as R" rule we
> port `lam_int`.

### R output (authoritative) — `R/cclasso.R:127-128`
```r
list(var_w = info_boot$var_w, cor_w = info_boot$cor_w,
     p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv)
```
| field | shape | meaning |
|---|---|---|
| `var_w` | `p` | bootstrap-mean of `diag(sigma)` (variance estimates) |
| `cor_w` | `p×p` | modified correlation matrix, diag = 1 |
| `p_vals` | `p×p` | two-sided p-values, diag = 1 |
| `lambda` | scalar | final λ = `10^((a2+b2)/2)` |
| `info_cv` | dict | `{lams, fvals, k, lam_int}` — see §3.1 |

Python returns a dict with identical keys. `info_cv["lam_int"]` holds the final bracket.

---

## 2. R → Python primitive mapping (exact)

| R | Python | note |
|---|---|---|
| `nrow/ncol` | `x.shape` | |
| `stats::var(x)` (p×p) | `np.cov(x, rowvar=False, ddof=1)` | R uses **n−1** denominator |
| `rowMeans(M)` | `M.mean(axis=1)` | length `nrow` |
| `colMeans(M)` | `M.mean(axis=0)` | |
| `mean(v)` | `v.mean()` | |
| `cor(M)` | `np.corrcoef(M, rowvar=False)` | Pearson, n−1 |
| `diag(M)` (extract) | `np.diag(M)` | |
| `diag(v)` (build) | `np.diag(v)` | |
| `M[ind_low]`, `lower.tri` | column-major lower triangle | use **R order** (see §5.3) |
| `x[i, ]` row subset | `x[idx, :]` | 1→0 based indices |
| `x[-itest, ]` | complement rows | |
| `A > l`, `A < -l` (logical) | `A > l`, `A < -l` then cast | strict, matches `sign·max(|A|−λ,0)` |
| `t(u) %*% M %*% u` | `u.T @ M @ u` | |
| `eigen(S)` | `np.linalg.eigh(S)` **reversed to descending** | see §5.1 |
| `base::norm(M,"F")` | `np.linalg.norm(M, "fro")` | |
| `pt(q, df=n-2)` | `scipy.stats.t.cdf(q, df=n-2)` | Student-t CDF |
| `ifelse(c, a, b)` | `np.where(c, a, b)` | |
| `median(v)` | `np.median(v)` | |
| `sample(1:n, size, replace=T)` | **injected** `boot_indices` | R RNG ≠ numpy |
| `matrix(rnorm(n*p), nrow=n)` | **injected** `null_norm` (n×p, **F-order**) | R RNG ≠ numpy |
| `matrix(..., ncol=n_boot)` | fill **column-major (F-order)** | |
| `rep(v, each=k)` | `np.repeat(v, k)` | |
| `v * M` (vector×matrix) | `v[:,None] * M` | R recycles **down columns** |
| `outer(a,b,"+")` | `a[:,None] + b[None,:]` | |

---

## 3. Function-by-function

### 3.1 `cclasso` — `R/cclasso.R:36-128`

Order of operations (must be preserved exactly):

1. `n,p = shape(x)`.
2. If `counts`: `x = (x + pseudo) / rowSums(x + pseudo)`.  (`:41-44`)
   - Note: `x <- x + pseudo; x <- x / rowSums(x)`.
3. `x = log(x)`.  (`:45`)
4. `vx2 = var(x)` (p×p).  (`:46`)
5. **Diagonal weight**: (`:49-51`)
   ```r
   rmean_vx2 <- rowMeans(vx2)
   wd  <- 1/diag( vx2 - rmean_vx2 - rep(rmean_vx2, each=p) + mean(rmean_vx2) )
   wd2 <- sqrt(wd)
   ```
   → `wd = 1/(diag(vx2) - 2*rmean_vx2 + mean(rmean_vx2))`, length `p`.
6. **Basis / ADMM constants**: (`:54-60`)
   ```r
   rho = 1
   u_f     = eigen(diag(p) - 1/p)$vectors          # p×p orthogonal
   wd_u    = (t(u_f) %*% (wd * u_f))[-p, -p]        # (p-1)×(p-1), symmetric
   wd_u_eig= eigen(wd_u)
   d0_wd   = 1 / ( (rep(wd_u_eig$values, each=p-1) + wd_u_eig$values)/(2*rho) + 1 )
   u0_wd   = wd_u_eig$vectors
   ```
   These are computed **once** and reused for every `cclasso_sub` call (CV, final,
   bootstrap, artificial-null). `sigma <- vx2` (full-data covariance) is the init.
7. **Golden section on log10(λ)** — see §3.2.
8. `info_cv`, `lambda = 10^((a2+b2)/2)`, `lambda2 = lambda/rho`.
9. **Bootstrap** `boot_cclasso(sigma = fb2$sigma, ...)` — see §3.4.
10. Return dict.

### 3.2 Golden-section λ search — `R/cclasso.R:62-118`

- `lam_int2 = log10(lam_int)` = `[log10(min), log10(max)]`.
- `a1,b1 = lam_int2`; `err_lam2 = 1e-1 * max(1, lam_int2)` → with defaults `max(1,-4,0)=1`
  so `err_lam2 = 0.1`; `err_fval = 1e-4`.
- Trial points: `a2 = a1 + 0.382*(b1-a1)`, `b2 = a1 + 0.618*(b1-a1)`.
- **Evaluation order (matters, warm-start chain):**
  1. `fb2 = cv_loss(10^b2/rho, sigma=vx2)`
  2. `fa2 = cv_loss(10^a2/rho, sigma=fb2$sigma)`
- `lams = c(b2, a2)`, `fvals = c(fb2, fa2)` (append order preserved).
- `err = b1-a1`, `k = 0`.
- Loop while `err > err_lam2 and k < k_max`:
  - `fval_max = max(fa2, fb2)`.
  - If `fa2 > fb2`: `a1=a2; a2=b2; fa2=fb2; b2=a1+0.618*(b1-a1);
    fb2 = cv_loss(10^b2/rho, sigma=fa2$sigma)`; append `(b2, fb2)`.
  - Else: `b1=b2; b2=a2; fb2=fa2; a2=a1+0.382*(b1-a1);
    fa2 = cv_loss(10^a2/rho, sigma=fb2$sigma)`; append `(a2, fa2)`.
  - `fval_min = min(fa2,fb2)`; `k += 1`; `err = b1-a1`;
    if `abs(fval_max-fval_min)/(1+fval_min) <= err_fval`: break.
- `info_cv = {lams, fvals, k: k+2, lam_int: 10^[a1,b1]}`.
- Boundary warning `cat(...)` when `a1==lam_int2[1] or b1==lam_int2[2]` → **stderr** in
  Python (stdout must stay clean for the CLI). Record in DIFFERENCES.
- `lambda = 10^((a2+b2)/2)`.

### 3.3 `cv_loss_cclasso` — `R/cclasso.R:211-233`

```r
n_b <- floor(n / k_cv)
cv_loss <- 0
for k in 1:k_cv:
    itest <- (n_b*(k-1)+1):(n_b*k)          # contiguous block
    vxk  <- var(x[itest, ])                 # test covariance
    vx2k <- var(x[-itest, ])                # train covariance
    sigma <- cclasso_sub(sigma, vx=vx2k, lambda2, ...)   # warm-started across folds
    dsig <- sigma - vxk
    tmp  <- rowMeans(dsig)
    dsig <- dsig - tmp - rep(tmp, each=p) + mean(tmp)
    cv_loss <- cv_loss + norm(wd2 * dsig, "F")^2
return list(cv_loss, sigma)                 # sigma = LAST fold's result
```

- **CV folds are deterministic contiguous blocks** — there is **no RNG in CV**.
  (Trailing `n mod k_cv` rows are never in any test block but are always in training.)
  ⇒ the brief's `cv_folds` injection is unnecessary; we keep the parameter name for
  API symmetry but it does not affect R-faithful behaviour. Recorded in DIFFERENCES.
- `sigma` returned is the last fold's converged matrix and is the warm start for the
  next golden-section evaluation.

### 3.4 `cclasso_sub` (ADMM / augmented Lagrangian) — `R/cclasso.R:236-278`

```r
cclasso_sub(sigma, vx, lambda2, wd, u_f, u0_wd, d0_wd, k_max=200, x_tol=1e-4)
```
- `p = ncol(sigma)`; `sigma2 = sigma`; `LAMBDA = 0`; `lambda2 = matrix(lambda2,p,p)` with
  **`diag = 0`** (penalty only off-diagonal).
- Loop while `err > 1e-4 and k < 200`:
  1. `x_sigma = t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f`
  2. `x_sigma[-p,-p] = u0_wd %*% ((t(u0_wd) %*% x_sigma[-p,-p] %*% u0_wd) * d0_wd) %*% t(u0_wd)`
  3. `sigma_new = vx + u_f %*% x_sigma %*% t(u_f)`
  4. `A = LAMBDA + sigma_new`
     `sigma2_new = (A > lambda2)*(A - lambda2) + (A < -lambda2)*(A + lambda2)`
  5. `LAMBDA = LAMBDA + (sigma_new - sigma2_new)`
  6. `err = max( abs(sigma_new-sigma)/(abs(sigma)+1),
                abs(sigma2_new-sigma2)/(abs(sigma2)+1) )`  (elementwise→max)
  7. `k += 1; sigma = sigma_new; sigma2 = sigma2_new`
- `if k >= k_max: cat(warning)` → **stderr** in Python.
- Returns `sigma` (the last `sigma_new`).
- Note `rho` is folded into `lambda2 = lambda/rho` and does **not** appear inside
  `cclasso_sub` in v2.0 (it did in v1.0). `rho = 1`.

### 3.5 `boot_cclasso` — `R/cclasso.R:141-207`

```
n,p
cors_boot : (p(p-1)/2) × (n_boot+1)
vars_boot : p × (n_boot+1)
ind_low   : lower.tri(p×p)                     # column-major order
sam_boot  : matrix(sample(1:n, n*n_boot, replace=T), ncol=n_boot)  # RNG!
for k in 1:n_boot:
    ind_samp = sam_boot[,k]
    sigma2 = cclasso_sub(sigma, vx=var(x[ind_samp,]), lambda2, ...)
    vars_boot[,k] = diag(sigma2)
    Is = 1/sqrt(vars_boot[,k]); cors_mat = Is * sigma2 * rep(Is,each=p)
    cors_boot[,k] = cors_mat[ind_low]
# extra column = full-data sigma
Is = 1/sqrt(diag(sigma)); cors_mat = sigma * Is * rep(Is, each=p)
cors_boot[,n_boot+1] = cors_mat[ind_low]; vars_boot[,n_boot+1] = diag(sigma)

vars2 = rowMeans(vars_boot)                     # length p
tol_cor = 1e-3

# artificial null sample (RNG!)
sam_art0 = matrix(rnorm(n*p), nrow=n) * rep(sqrt(vars2), each=n)
cors_art0 = cor(sam_art0)[ind_low]
sam_art   = sam_art0 - log(rowSums(exp(sam_art0)))   # row-broadcast subtract
sigma_art = cclasso_sub(sigma, vx=var(sam_art), lambda2, ...)
Is = 1/sqrt(diag(sigma_art)); cors_mat = Is*sigma_art*rep(Is,each=p)
cors_art2 = cors_mat[ind_low]

cors0m2 = log( ((1+cors_art0)*(1-cors_art2)) / ((1+cors_art2)*(1-cors_art0)) )
tmp = abs(cors_art2) >= tol_cor
bias02 = ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0)

cors2 = log((1+cors_boot)/(1-cors_boot))
cors2mod = (cors_boot >=  tol_cor)*(cors2 + bias02) +
           (cors_boot <= -tol_cor)*(cors2 - bias02)
cors2mod = 1 - rowMeans( 2/(exp(cors2mod)+1) )      # over n_boot+1 cols
cors2_mat symmetric, diag=1, lower←cors2mod

p_vals = pt(cors2mod*sqrt((n-2)/(1-cors2mod^2)), df=n-2)
p_vals = ifelse(p_vals <= 0.5, p_vals, 1-p_vals)
pval_mat symmetric, diag=1, lower←p_vals
return {var_w: vars2, cor_w: cors2_mat, p_vals: pval_mat}
```

RNG points: `sample()` (§3.5) and `rnorm()` (artificial null). Both must be injectable.

---

## 4. Determinism / injection design

- **Only two sources of randomness**, both inside `boot_cclasso`:
  1. `sam_boot` — bootstrap row indices, `n × n_boot`, integer in `1..n`.
  2. `null_norm` — the raw `rnorm(n*p)` draws for `sam_art0`, `n × p`, **F-order**.
- Python API adds optional injections:
  - `boot_indices : (n, n_boot) int array | None`
  - `null_norm    : (n, p) float array | None`  (already the raw draws; scaled internally)
  - `seed : int | None` → `np.random.Generator(PCG64(seed))` used only if the above are `None`.
- If not injected, results are **statistically** equivalent, not identical (R RNG ≠ numpy).
- CV folds are deterministic ⇒ no injection needed there.

---

## 5. Numerical-equivalence analysis

### 5.1 `eigen` ordering & the basis-invariance result (IMPORTANT)

- R `eigen(S)` returns **decreasing** eigenvalues; `np.linalg.eigh` returns **ascending**.
  Reverse eigh output to descending so `values[i] ↔ vectors[:,i]` ordering matches R.
- `u_f = eigen(diag(p) - 1/p)$vectors`: this matrix has eigenvalue **1 with multiplicity
  p−1** and eigenvalue **0 with multiplicity 1**. The eigenvector of the zero eigenvalue is
  (up to sign) `1/√p` (the all-ones direction). Because the eigenvalue-1 eigenspace is
  degenerate, R's and numpy's eigenvectors inside it **differ**, so `u_f` cannot be made
  bit-identical. **Claim:** the algorithm only needs (a) `u_f` orthonormal with its **last
  column = ±ones/√p**, and (b) any orthonormal eigenbasis of `wd_u`.
  - In `cclasso_sub`, the special treatment `x_sigma[-p,-p] ← f(x_sigma[-p,-p])` singles out
    the **null (ones) coordinate**, and `f` (diagonalise with `u0_wd`, scale by `d0_wd`,
    rotate back) is **equivariant** under any orthogonal change of basis `Q` of the
    remaining `e^⊥`:
    `f(QᵀBQ) = Qᵀ f(B) Q` because eigenvalues are basis-invariant and
    `u0_wdᵀ B u0_wd` is invariant to the `u0_wd` basis choice.
    Hence `sigma_new = vx + u_f (·) u_fᵀ` is the same matrix regardless of which
    orthonormal basis is used **provided the last column spans `1`**.
  - ⇒ We can construct `u_f` ourselves (last column `= 1/√p`, rest any orthonormal
    complement) and still be numerically equivalent to R up to floating-point rounding.
    **Action:** implement `u_f` via reversed `eigh(diag(p)-1/p)` (structure matches R) and
    *verify numerically* against golden data. This is the single most important item to
    validate empirically.
- Signs of eigenvectors are irrelevant for the same reason (they cancel in the similarity
  transform). Only invariant combinations are compared (reconstructed matrices).
- `eigen(wd_u)`: `wd_u` is symmetric `(p−1)×(p−1)`; values are basis-invariant, and `f`
  is basis-covariant ⇒ also invariant.

### 5.2 Soft-threshold at equality
`(A > λ)*(A−λ) + (A < −λ)*(A+λ)` uses **strict** inequalities: exactly `|A| == λ` → 0.
Reproduce with strict `>`/`<`, **not** `np.sign(A)*np.maximum(|A|−λ,0)` (which keeps ±0
differently only at the boundary — numerically identical, but implement literally).

### 5.3 Lower-triangle ordering
R fills `cors_boot`/`vars_boot` rows via `M[lower.tri(M)]` in **column-major** order
`(2,1),(3,1),…,(p,1),(3,2),…`. numpy `np.tril_indices` is row-major. Use a fixed
column-major index vector when filling/extracting so intermediate golden arrays align
exactly. (`rowMeans`/`median` are order-independent, but golden comparison of the raw
`cors_boot`/`cors_art2` arrays needs matching order.)

### 5.4 ADMM stopping tolerance / max iterations
`x_tol = 1e-4`, `k_max = 200` for `cclasso_sub`; relative error uses
`|Δ|/(|old|+1)` and `max` over all entries (mixed `sigma`/`sigma2`). Iteration count is
sensitive near the boundary: a 1-iteration shift is possible from BLAS/LAPACK rounding and
will be reported as a failure at the exact step (per the test spec).

### 5.5 Golden-section `err_lam2`
`err_lam2 = 0.1 * max(1, lam_int2)` is a **scalar** (R `max` collapses the vector):
with defaults = `0.1`. Implement literally including the `max(1, ...)`.

---

## 6. Explicit uncertainties (must be resolved empirically, not assumed)

1. **`u_f` basis** (§5.1): proven invariant in theory; **must** be confirmed against golden
   R output with tight tolerance. If a mismatch appears, it will show first in
   `sigma` after the first `cclasso_sub` iteration.
2. **RNG mismatch** — folds deterministic, but bootstrap `sample()` and `rnorm()` differ:
   only injectable-index paths can be compared exactly; others only statistically.
3. **BLAS/LAPACK rounding** (~1e‑12…1e‑15) amplified by ADMM warm-start chains and the
   `fvals` comparison `fa2 > fb2` — a tie/±1‑ulp flip can change which branch golden section
   takes, changing the whole λ path. This is the most likely source of a *hard* divergence.
   Will be documented in `docs/EQUIVALENCE.md`.
4. **`pt` (Student-t CDF)** — scipy vs R's `pt`: parity expected to ~1e‑12; not proven yet.
5. **`median` with even counts** — R `median` averages the two middle values; numpy matches.
6. **`ifelse(sum(tmp), …)`** — `sum(tmp)` over a logical across `n_boot+1` columns;
   need to confirm `tmp` length (`cor` vector of `p(p−1)/2`) as implemented.
7. **`cor()` vs `np.corrcoef`** for `cors_art0` (only used through `ind_low`) — must match
   column-major extraction.
8. **Boundary `cat` warnings** — R prints to **stdout**; Python must print to **stderr**
   (CLI contract). Recorded as an intentional, documented difference.
9. **`info_cv$k = k+2`**, `lams`/`fvals` append order — must reproduce exactly for golden
   comparison.
10. **v1.0 vs v2.0** — confirm v2.0 is the target (assumed from `DESCRIPTION`).

---

## 7. Planned repository layout

```
cclasso-python/
  pyproject.toml            # deps: numpy, scipy; dev: pytest, pytest-cov, ruff
  src/cclasso/
    __init__.py
    _core.py                # cclasso, cv_loss_cclasso, cclasso_sub, boot_cclasso (+ # R: refs)
    _linalg.py              # eigen_order, lower_tri_colmajor, soft_threshold helpers
    validation.py           # input validation with clear errors
    cli.py                  # python -m cclasso  (stdin/file JSON -> stdout JSON)
    __main__.py
  tests/
    golden/                 # committed R-generated fixtures (no R needed in CI)
    r_reference/generate_golden.R
    r_reference/README.md
    test_equivalence.py
    test_properties.py
    test_cli.py
    test_validation.py
  docs/{PORTING_NOTES.md,DIFFERENCES.md,EQUIVALENCE.md,PACKAGING.md}
  README.md  CHANGELOG.md  LICENSE  .gitignore
  .github/workflows/ci.yml  # lint + pytest, py3.9–3.12, linux/win/macOS
```

## 8. Current real-data golden fixtures (`tests/r_reference/generate_golden.R`)

The committed goldens include two QMP cases plus American Gut Round 1 and
GlobalPatterns cases selected by `scripts/public_microbiome_cases.R`. Each added
case family uses three R RNG seeds. The R generator records full-precision
intermediates (`%.17g`), injected bootstrap/null draws, ADMM iteration counts and
terminal residuals, source/sample/taxon identifiers, source profiles, and
`sessionInfo()`. A separate `tests/golden/failures/` case records the
GlobalPatterns prevalence-spread R error and its matching Python error status.
Small synthetic property/recovery tests remain generated in memory; they are
not golden data files.

---

## 9. Questions for the user (blocking decisions)

1. **Target version**: port `v2.0` (root `R/cclasso.R`) — confirm (not `v1.0`).
2. **Scope**: `cclasso` only, or also `SparCC.count`/`SparCC.frac`?
3. **Repo location/name**: `electron-analysis-app/cclasso-python/` (sibling of
   `spiec-easi-python`), or a separate folder/repo? Package name `cclasso` OK?
4. **License**: CCLasso is **GNU LGPL-v3**; a faithful translation is a derivative work, so
   the Python library must be **LGPL-v3 (or later) / GPL-compatible**. `spiec-easi-python`
   in this workspace is MIT → **incompatible** if code is shared. Confirm we set the new
   library's LICENSE to LGPL-v3 and keep it independent from the MIT project. (I will not
   choose the license for you.)

---

## 10. Resolution of §6 uncertainties (after implementation)

Decisions confirmed with the user: port **v2.0**, `cclasso` only, repository at
`cclasso-python/`, license **LGPL-3.0-or-later**.

| # | uncertainty | outcome |
|---|---|---|
| 1 | `u_f` eigenvector basis | **Resolved.** The basis-invariance argument (§5.1) holds; on current QMP cases the port reproduces `lams` and `lambda` bit-identically and continuous results within the QMP maxima in `docs/EQUIVALENCE.md`. `u_f` is built from a reversed `eigh(I - 1/p)` (last column = ones direction). |
| 2 | RNG mismatch | **Resolved by injection.** `boot_indices` (0-based) and `null_norm` are injected in tests; RNG replay in the generator is asserted exact. |
| 3 | BLAS/LAPACK rounding / ADMM iteration flips | No QMP iteration-count divergence was observed on this platform. All 70 recorded QMP ADMM solves hit the cap in both languages; terminal residuals match within `3.96e-14` absolute. |
| 4 | `pt` (Student-t CDF) | **Resolved.** QMP `p_vals` match to ≤1.48e-13 absolute. |
| 5 | `median` (even counts) | **Resolved.** `np.median` matches R. |
| 6 | `ifelse(sum(tmp), …)` | **Resolved.** Implemented as `tmp.any()`; matches. |
| 7 | `cor()` for `cors_art0` | **Resolved.** `np.corrcoef` + column-major `lower.tri` matches. |
| 8 | boundary `cat` warnings | **Resolved as a difference (D3):** routed to stderr / `warnings`. |
| 9 | `info_cv$k = k+2`, append order | **Resolved.** `info_cv.k`, `lams`, `fvals` reproduce exactly. |
| 10 | v1.0 vs v2.0 | **Resolved:** v2.0 only (D4). |

Extra behaviours discovered and documented: non-convergent ADMM cases can yield
shared R/Python non-finite results, and a GlobalPatterns stress seed can raise
the same missing-logical failure in both implementations. These outcomes are
recorded, not treated as valid estimates. Private diagnostic keys `_admm_iters`,
`_admm_relative_errors`, and `_bootstrap_variances` are returned by the Python
API and stripped by the CLI.
