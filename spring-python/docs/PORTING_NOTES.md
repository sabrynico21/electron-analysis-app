# PORTING_NOTES — SPRING (R) → `spring` (Python)

**Status:** implemented and validated for release `1.0.0` (2026-10-20). The
analysis below is retained as the design record; the shipped code and the
equivalence evidence live in `src/spring_microbiome/` and `docs/EQUIVALENCE.md`.
**Author:** generated during the porting session of 16 Oct 2026.
**Purpose:** record, *before writing any code*, the exact behaviour of every R
function involved in SPRING, the external R packages it delegates to, and the
fidelity status (identical / numerically equivalent / not equivalent) of the
planned Python replacement.

> **Source of truth rule (from the task):** the R code prevails over my prior
> knowledge of SPRING and over the paper. Every time the code and the paper
> disagree, the code wins and the disagreement is listed in `DIFFERENCES.md`.

---

## 0. Sources and exact versions

The SPRING source read for this port is the **modernized fork** shipped in the
workspace, *not* the original `GraceYoon/SPRING` repository. This matters a lot:
the modern fork replaced the original `mixedCCA`-based latent correlation with
the **`latentcor`** package and rewrote the MB solver inside **`huge` 2.0.1**.

| Role | Package | Version | Origin | License |
|------|---------|---------|--------|---------|
| Port target | `SPRING` | 0.99.0 | `electron-analysis-app/SPRING-master` (URL `IrinaStatsLab/SPRING`) | GPL-3 |
| Kendall τ + latent correlation | `latentcor` | 2.0.2 | CRAN | GPL-3 |
| MB neighborhood selection / Lasso CD | `huge` | 2.0.1 | CRAN (modernized, `Gatech-Flash/huge`) | GPL-2 |
| StARS stability selection | `pulsar` | 0.3.13 | GitHub `zdk123/pulsar` @ `2bc391f` | GPL (≥3) |
| Kendall τ-a fast | `pcaPP` | 2.0.7 | CRAN (`cor.fk`) | GPL (≥2) |
| bivariate normal CDF (bridge) | `fMultivar` | 4031.84 | CRAN (`pnorm2d`) | GPL-2 |
| multivariate normal CDF (bridge) | `mnormt` | 2.1.2 | CRAN (`pmnorm`) | GPL-2/3 |
| nearPD / sparse | `Matrix` | 1.7.5 | base-recommended | GPL |
| (synthetic data only) | `mvtnorm`, `rootSolve`, `SpiecEasi` | — | CRAN / Bioconductor | GPL |

R runtime on the porting machine: **R 4.6.1 (2026-06-24)**, `x86_64-pc-linux-gnu`.

`sessionInfo()` for the golden generator is saved next to the golden data
(see §9). Dependency source trees were extracted for reading to
`/tmp/r_dep_src/{huge_src,latentcor_src,pulsar_src}`.

### 0.1 Recommended Python runtime
`numpy`, `scipy`, `scikit-learn` — but per §5 the Lasso **must not** be
`scikit-learn`; a faithful re-implementation of `huge`'s coordinate descent is
required. `scipy` is needed for `stats.norm.ppf`, `bivariate/multivariate normal
CDF`, and the optimiser used by `Rmethod="original"`.

---

## 1. Public R API to reproduce (same names, same defaults)

```r
SPRING(data, quantitative = FALSE, method = c("mb"),
       lambda.min.ratio = 1e-2, nlambda = 20,
       lambdaseq = exp(seq(log(0.6), log(0.6 * lambda.min.ratio),
                           length.out = nlambda)),
       seed = 10010, ncores = 1, thresh = 0.1, subsample.ratio = 0.8,
       rep.num = 20, Rtol = 1e-6, verbose = TRUE,
       Rmethod = "approx", use.nearPD = TRUE, nu = 0.001, ratio = 0.9, ...)
```
`SPRING.R:57` (generic, S3), `SPRING.R:98` (`SPRING.default`, the real body),
plus `SPRING.phyloseq` (`SPRING.R:78`) and `SPRING.otu_table` (`SPRING.R:88`)
which merely extract the OTU matrix via `.phy2mat` (`SPRING.R:61`).

Also exported and to be ported:
* `hugeKmb(data, lambda, type = "tru", sym = "or", verbose = TRUE,
  Rmethod = c("approx","original"), tol = 1e-6, use.nearPD = TRUE,
  nu = 0.001, ratio = 0.9)` — `helpers.R:44`
* `mclr(dat, base = exp(1), tol = 1e-16, eps = NULL, atleast = 1)` — `helpers.R:71`
* `synthData_from_ecdf(comm, mar = 2, Sigma, n, verbose = FALSE)` — `synthData.R:19`
  (see §8, open question on scope).

### 1.1 Return object (`SPRING.default`)
```
list(
  output = <pulsar object>,      # output$stars$merge, $summary, $opt.index,
                                 # $criterion, $thresh; output$est = full hugeKmb fit
  fit    = <pulsar.refit>,       # fit$est = full hugeKmb, fit$refit$stars = adjacency
  lambdaseq = <numeric vector>
)
```

---

## 2. Pipeline, in order (`SPRING.default`, `SPRING.R:98-146`)

1. `dots <- list(...)`; **any** extra argument → hard `stop("unused arguments: …")`.
2. `method <- match.arg(method)`; only `"mb"` exists.
3. `if (any(data < 0)) stop(...)` — the **only** input validation in the code.
4. `p <- ncol(data)`.
5. `quantitative=TRUE`: warn if `max(rowSums(data)) <= 1` (looks normalised);
   `qdat <- data` unchanged.
   `quantitative=FALSE`: `qdat <- mclr(data)`.
6. `gc()` (no numerical effect).
7. If `lambdaseq == "data-specific"`: compute
   `Kcor = latentcor(qdat, types="tru", method=Rmethod, tol=Rtol,
                    use.nearPD=use.nearPD, nu=nu, ratio=ratio)$R`,
   `lambda.max = max(max(Kcor-diag(p)), -min(Kcor-diag(p)))`,
   `lambda.min = lambda.min.ratio*lambda.max`,
   `lambdaseq = exp(seq(log(lambda.max), log(lambda.min), length=nlambda))`.
   Any other character → `stop`. Non-character `lambdaseq` is used as-is.
8. `fun <- hugeKmb`.
9. `out1 <- pulsar::pulsar(qdat, fun=hugeKmb,
     fargs=list(lambda=lambdaseq, Rmethod=Rmethod, tol=Rtol, verbose=verbose,
                use.nearPD=use.nearPD, nu=nu, ratio=ratio),
     rep.num=rep.num, criterion="stars", seed=seed, ncores=ncores,
     thresh=thresh, subsample.ratio=subsample.ratio)`.
10. `fit1 <- pulsar::refit(out1)`.
11. return `list(output=out1, fit=fit1, lambdaseq=lambdaseq)`.

**Note (paper vs code):** the *paper* describes the original `mixedCCA`-based
estimate. The code uses `latentcor` (`types="tru"` fixed inside `hugeKmb`).
This is a known deliberate change (see the fork's vignette) → noted in
`DIFFERENCES.md`.

---

## 3. `mclr` — modified central log-ratio (`helpers.R:71`)

```
nzero <- dat >= tol          # tol = 1e-16
LOG   <- ifelse(nzero, log(dat, base), 0)          # log base e by default
if nrow(dat) > 1:  clrdat <- ifelse(nzero, LOG - rowMeans(LOG)/rowMeans(nzero), 0)
if nrow(dat)==1:   clrdat <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0)
```
* `rowMeans(LOG)` includes the zeros (they contribute 0), and
  `rowMeans(nzero)` = (#non-zero)/p. So the subtracted constant equals
  `sum_{j: nz} log(x_j) / (#nz)` — the geometric mean of the **non-zero** part.
* `nzero` uses `>=` with `tol`, not `>`.

Branch on `eps`:
* `eps = NULL` (default): if `atleast < 0` → warn and reset `atleast <- 1`.
  `positivecst = abs(min(clrdat)) + atleast` if `min(clrdat) < 0` else `0`.
  Return `ifelse(nzero, clrdat + positivecst, 0)` — **zeros stay exactly 0**,
  all non-zero values become strictly positive.
* `eps == 0`: return `clrdat` (may contain negatives).
* `eps > 0`: return `ifelse(nzero, clrdat + eps, 0)`.
* `eps < 0`: `stop`.

**Python equivalence:** straightforward; `np.log`, row masking. Watch
`ifelse` semantics (vectorised, no short-circuit) and that `min(clrdat)` is
taken over the whole matrix (including the zeros, which are 0, so the minimum is
non-positive). Floating-point: `log` of the same double gives bit-identical
results on glibc; `rowMeans` divides by `p` (same as `sum/p`). Expected
agreement ≲ 1e-15.

**Uncertainty:** none numerically; only the "unused arguments" stop and the
`atleast<0` warning need faithful reproduction in messages.

---

## 4. `latentcor` — rank-based latent correlation (external, `types="tru"`)

`hugeKmb` (`helpers.R:44`) calls, per (sub)sample:
```r
S <- latentcor::latentcor(data, types = "tru", method = Rmethod, tol = tol,
                          use.nearPD = use.nearPD, nu = nu, ratio = ratio,
                          showplot = FALSE)$R
est <- huge::huge.mb(S, lambda, sym = "or", verbose = verbose)
```

### 4.1 `latentcor()` control flow (`latentcor.R`)
1. Validate `nu∈[0,1]`, `tol>0`, `ratio∈[0,1]`.
2. `X <- as.matrix(X)`; `p <- ncol(X)`.
3. `types` given as `"tru"` (length 1) → replicated to all `p` columns.
4. `encodeX(X, types)`: for `tru`, `X <- X - colmin` (per column). `negate=1`
   (all-ones) for `tru`, so the final `R * outer(negate,negate)` is a no-op.
5. Kendall matrix on **encoded** X (see §4.2).
6. `zratios`: proportion of zeros per column (see §4.3).
7. For each lower-triangle pair, ordered by type code so that `types_cp[1] <=
   types_cp[2]`; `comb = paste0(code1, code2)`. `"00"` (con/con) → closed form
   `sin(pi/2 * K)`; otherwise → `r_switch` (§4.4).
8. Assemble symmetric `K` (diag 1) and `R` (diag 1). Save `Rpointwise <- R`.
9. If `use.nearPD`: `R = as.matrix(Matrix::nearPD(R, corr=TRUE)$mat)`, then
   `R = (1-nu)*R + nu*I`.

Since SPRING always passes `types="tru"`, only `comb ∈ {"22"}` (tru/tru) can
occur — **all off-diagonal pairs are truncated/truncated**. This is a big
simplification: only `ipol_22`, `bound_switch("22")`, `r_sol("22")` are needed.

> ⚠️ Verify: `type` is `"tru"` for *every* column because `hugeKmb` hard-codes
> `type="tru"` and passes it through `latentcor(types=type)`. Confirmed by
> reading `helpers.R:44-49`. So `comb` is always `"22"`.

### 4.2 Kendall τ-a with tie correction (`internal.R::Kendalltau`, `n_x`)
```
n0  = n(n-1)/2
n_X = for each column: sum over tie-groups t: t(t-1)/2        # "n_x"
n_X_sqd  = sqrt(n0 - n_X)
K_b      = pcaPP::cor.fk(X)          # Kendall tau-b (numerator/denominator form)
btoa     = n_X_sqd[i]*n_X_sqd[j]/n0  # pairwise
K_a      = K_b * btoa
```
`pcaPP::cor.fk` computes τ **with the usual O(n log n) algorithm** and is *not*
exactly the textbook τ-b nor τ-a; the correction back to τ-a is done by `btoa`.
`n_x` uses `if (length(unique(x) != n))` — note the **parenthesis bug**:
`unique(x) != n` compares to a vector and `length(...)` is almost always `n`,
so the `if` is essentially always TRUE and `n_x` is computed by `rle(sort(x))`.
Faithful port: always compute tie counts via run-length of the sorted values.
(If `if` were FALSE, `n_x` would be 0 — only when *all* values are distinct,
which the `if` cannot detect reliably; but the computed value is 0 in that case
anyway, so the bug is harmless.)

**Python equivalence:** `scipy.stats.kendalltau` gives τ-b; to get τ-a we
re-derive with the same `btoa`. `pcaPP::cor.fk` uses a specific tie handling
(numerator denominator based). It should equal `(concordant - discordant) /
sqrt((n0-nx)(n0-ny))` = τ-b. **Uncertainty:** flag — must be validated against
the golden `K` matrix; if `scipy.stats.kendalltau` differs (e.g. continuity or
weighting during ties), implement `cor.fk`'s algorithm directly. Tentative
tolerance: strict (1e-12) after validation, else implement exactly.

### 4.3 `zratios` for `tru`
`colMeans(X == 0)` on the encoded matrix → proportion of zeros (after
subtracting the column min, zeros are exactly the original column minima).
`zratio1 = zratio2 = zratio` for pair (i,j).

### 4.4 Bridge inversion (`internal.R::r_sol`, `r_switch`, `bound_switch`)
For `comb="22"`:
* `bound = 1 - pmax(z1,z2)^2`.
* `r_switch("approx")`: `cutoff = abs(K) > ratio * bound`.
  * not cut → `r_ml(K/bound, z1, z2, "22")` (multilinear, §4.5) — note the
    division of K by `bound`, and the result is divided by `1e7`.
  * cut → `r_sol(K, z1, z2, "22")`: minimise `(F(r) - K)^2` over `r∈[-0.999,0.999]`
    with `optimize(tol=tol)` where F is the 4-variate normal expression using
    `mnormt::pmnorm` (two 4×4 covariance matrices `mat1`, `mat2`).
`r_sol` uses base R `optimize` (Brent, golden-section + parabolic). **Resolved:**
the port re-implements R's `Brent_fmin` verbatim from `stats/src/optimize.c`
instead of using `scipy.optimize`, because on this flat objective SciPy's Brent
locates a materially different minimiser (~1e-5). See
`spring_microbiome/latentcor/bridge.py::brent_fmin`.

**Remaining, unavoidable:** the bridge function itself needs the 4-variate
normal CDF. R's `mnormt::pmnorm` and `scipy.stats.multivariate_normal.cdf`
disagree by up to ~8e-5 on the mildly ill-conditioned 4×4 bridge matrices that
occur here. **Arbitrated** with independent 4-D Gauss–Legendre quadrature
(n=45/dim): SciPy agrees with the quadrature to ~3e-6 while `mnormt` is off by
~7.6e-5 on the worst matrix, and `mnormt`'s value does not change when
`abseps` is tightened from 1e-6 to 1e-12. So **R's value is the inaccurate one**
and we deliberately do not reproduce it; our inverted correlations are the more
accurate. Because the objective is flat, the inverted correlation still differs
from R by up to ~5e-5. Documented in `docs/EQUIVALENCE.md` §2; the discrete
pipeline outputs (`opt.index`, edge sets) are unaffected.

**Determinism gotcha:** `scipy.stats.multivariate_normal.cdf` uses randomized
quasi-Monte-Carlo and is non-deterministic by default. The port therefore passes
a *fresh* `np.random.default_rng(0)` on every call (with
`abseps=releps=1e-10, maxpts=1e6`), so each evaluation is reproducible and does
not depend on unrelated call history.

### 4.5 Multilinear approximation (`r_ml`, `ipol_*`, `ipol.c`)
* Precomputed tables `ipol_10…ipol_33` live in `latentcor/R/sysdata.rda`.
  For `"22"` the table `ipol_22` is a `41×49×49` grid (row dim = the rescaled
  `K/bound ∈ [-0.983604, 0.983604]`, other dims = `zratio ∈ [0.079656, 0.98758]`),
  stored as **integers scaled by 1e7**.
* `r_ml(K, z1, z2, "22")` → `ipol_22(rbind(K, z1, z2)) / 1e7`.
* Evaluation = **multilinear interpolation on a non-uniform Cartesian grid**
  (`C_evalmlip`, `ipol.c`), blend = `linear` (blend id 0):
  - per dimension: `gp = binsearch` = smallest index with `grid[gp] >= x`,
    clamped to `[1, n-1]` (never 0, never beyond n-1); weight
    `w = (x - grid[gp-1])/(grid[gp]-grid[gp-1])`;
  - flatten via **column-major** stride (`dim0` fastest), sum the `2^rank`
    corners weighted by products of `w`/`1-w`;
  - outside the box → **linear extrapolation** (values are *not* clipped).

**Python equivalence plan (preferred):** extract the four `ipol_22` arrays
(grid ×3 + integer values) from a locally-installed `latentcor` and ship them
as a versioned `.npz` data file, then re-implement `C_evalmlip` (≈40 lines of
numpy). This guarantees the *inputs* to interpolation are bit-identical to R, so
the only residual difference is the interpolation arithmetic itself (few ulps).
Alternative (regenerate from scratch) would depend on `mnormt`/`fMultivar`
versions → **rejected** unless the user asks otherwise.

### 4.6 nearPD + shrinkage
`Matrix::nearPD(Rpointwise, corr=TRUE)$mat` is R's Higham algorithm with
specific defaults (`corr=TRUE`, `keepDiag=FALSE`, `maxit=100`, `posd.tol=1e-8`,
`do2eigen=TRUE`, `base.matrix=FALSE`). Then ridge toward identity `nu`. 
**Python equivalence:** no bit-identical equivalent exists in numpy;
`statsmodels`/custom Higham iterations differ in details. **Plan:** port
`Matrix::nearPD` faithfully (Higham, same iteration order and tolerances) rather
than use an arbitrary PSD projection. This is a *critical* node because its output
feeds the entire MB path. **Flagged as high-risk for exactness** — will be
validated per golden case; if it is not reproducible, the discrepancy is
documented as a tolerance-level (not structural) difference.

---

## 5. MB neighborhood selection — `huge.mb` (external, huge 2.0.1)

`hugeKmb` passes the **correlation matrix** `S` (not the data) to
`huge::huge.mb(S, lambda, sym="or")`.

### 5.1 R-level wrapper (`huge::huge.mb`)
* `S` is square, symmetric, finite, diag 1 → `.huge_is_covariance_input`
  returns TRUE → `cov.input=TRUE` → `S = inp$correlation` (cov→cor conversion;
  with diag 1 it is a no-op up to rounding, clamped to `[-1,1]`).
* `lambda` supplied → validated (strictly positive, non-increasing) and used
  **verbatim** (no `lambda.max` recomputation, no `cummin`).
* Solver: `.Call("_huge_SPMBgraph", S, lambda, nlambda, d)`.
* Post-processing per λ: `beta = G_slice` (p×p sparse);
  `path = abs(beta)`; `df = colSums(path≠0)`;
  `sym="or"` → `path = sign(path + t(path))` ⇒ **unsigned {0,1} adjacency**
  (because `path≥0` before symmetrisation); `and` → `sign(path*t(path))`.
  `sparsity = sum(path)/d/(d-1)`.

### 5.2 Solver core (`huge_core.cpp::mb`, L679–L812)
Nodewise lasso objective (implicit): `min_w ½ wᵀS w − S[:,m]ᵀ w + λ|w|₁`.
Exact algorithm per response column `m` (note the response is **not** in the
active set; the diagonal `S[m,m]=1` is effectively included via the `+w0[j]`
term in the refine step):
```
thol = 1e-4 ; MAX_ITER = 10000
w0 = w1 = 0 ; active = ∅ ; inactive = all j≠m
grad_abs[j] = |S[m,j]|
for i in 0..nlambda-1:
    λ = lambda[i]
    strong_thr = (i>0) ? 2λᵢ - λᵢ₋₁ : 2λᵢ        # sequential strong rule
    strong[j] = (j inactive) and grad_abs[j] > strong_thr
    loop:
      # strong sweep (order j = 0..d-1)
      for j: if inactive and strong[j]:
          r = S[m,j] - Σ_k S[j, a_k] w0[a_k]
          grad_abs[j] = |r| ; w1[j] = soft(r, λ)
          if w1[j]!=0: activate j
          w0[j] = w1[j]
      # refine active set (cycle CD) until relative change <= thol
      mb_refine_active(...):
          repeat:
             for j in active:
                r = S[m, a_j] + w0[a_j] - Σ_k S[a_j, a_k] w0[a_k]
                w1[a_j] = soft(r, λ)
             gap = change_sum/coef_sum ; until gap <= thol or iter>=MAX_ITER
      # drop zeros from active, re-enable them (strong=1) and re-sweep
      converge when no coordinate was newly activated in the sweep
    # certification pass over all inactive coords; violators activate and the
    # strong loop re-runs; if violations==0 -> certified
    collect_sorted(...) : emit active coeffs sorted by ascending index
soft(x, λ): x>λ → x-λ ; x<-λ → x+λ ; else 0     # strict inequalities
```
* Warm start across λ (w0 carries over), `strong` recomputed each λ,
  `grad_abs` refreshed by the certification pass.
* `collect_sorted` sorts a **scratch copy** (solver sweep order preserved).
* Missing coordinates (`w1` stays 0) are simply absent from the output.

**Python equivalence:** re-implement `mb()` in pure numpy/`numexpr`-free code
following the pseudocode above *exactly* — same sweep order (0…d-1), same
`thol`, same strong-rule threshold, same strict soft-threshold, same
certification. **Do not** use `sklearn`'s Lasso. Residual differences will only
come from floating-point summation order in `r` (documented).

### 5.3 What `hugeKmb` returns (per subsample)
`list(beta, path, df, sparsity, lambda, cov.input, scr)`. `pulsar` uses only
`path` (a list of `nlambda` adjacency matrices).

---

## 6. StARS — `pulsar` (external, pulsar 0.3.13)

### 6.1 Subsample generation (`mcPulsarSelect.R::pulsar`)
```
n, p ; .lamcheck(lambda)
subsample.ratio = .ratcheck(subsample.ratio, n)
   = given value, else  0.8 if n<=144 else 10*sqrt(n)/n
set.seed(seed)                              # R RNG — NOT reproducible in numpy
ind.sample = replicate(rep.num,
    sample(1:n, floor(n*subsample.ratio), replace=FALSE), simplify=FALSE)
set.seed(NULL)
if refit: ind.sample = c(list(1:n with attr "full"), ind.sample)
estFun(ind) : do.call(hugeKmb, c(fargs, list(data[ind,]))) -> tmp$path
             (full model returns the whole object, not just path)
premerge = .try_mclapply(isamp, estFun, mc.cores=ncores, mc.preschedule=FALSE)
   -> with ncores=1 this is mclapply(..., mc.cores=1) = serial lapply.
rep.num reduced by #errors; fullmodel = premerge[[1]]; premerge = premerge[-1]
premerge.reord = .tlist(premerge, nlams, rep.num)   # [lambda][replicate]
est$stars = stars.stability(premerge.reord, thresh, rep.num, p)
```

### 6.2 `stars.stability`
```
merge[[i]] = Σ_replicates path[[i]]                       # counts
merge[[i]] = merge[[i]] / rep.num                         # selection prob.
summary[i] = 4 * Σ merge[[i]]*(1-merge[[i]]) / (p*(p-1))  # instability
summary    = cummax(summary)                              # monotone
opt.index  = .starsind(summary, thresh)
```
`.starsind(summary, thresh, offset=1)`:
`if any(summary>=thresh): max(which.max(summary>=thresh) - offset, 1)`
(= first index where cummax≥thresh, minus 1, floored at 1);
else warn and return `length(summary)`.
If `opt.index==1` pulsar emits the "Optimal lambda may be larger/smaller …"
warning (direction = "larger" if any summary≥0.1).

### 6.3 `refit.pulsar`
`est$est = obj$est` (the full-data `hugeKmb` fit computed inside `pulsar`);
`est$refit$stars = est$est$path[[opt.index]]`.

### 6.4 Where randomness lives
**Only** in `ind.sample`. Everything else is deterministic. The port will expose
`subsample_indices=None`; when provided (goldens) the exact R indices (converted
**1-based R → 0-based Python**) are reused; otherwise generated with
`numpy.random.Generator(seed)` (statistically equivalent, **not** identical).

---

## 7. `var`/`cov`/`cor`/`scale`/`sweep`/`apply`/`eigen` conventions

| R construct | Meaning | Python |
|---|---|---|
| `cov`, `var` | denominator `n-1` | `np.cov(..., ddof=1)` |
| `cor` | Pearson, `n-1` | `np.corrcoef` |
| `rowMeans(LOG)` | `sum/p` | `arr.sum(1)/p` |
| `colMeans(mask)` | `sum/n` | `arr.mean(0)` |
| `scale` | centre by mean, divide by sd (`n-1`) | manual |
| `sweep` | broadcast apply | manual |
| `apply(M, 2, f)` | per column | vectorised |
| `eigen(x, symmetric=TRUE)` | **descending** eigenvalues | `np.linalg.eigh` **ascending** → reverse |
| `lower.tri` | **column-major** order | explicit `np.tril_indices` (row-major) — must reorder to match R's `cp` indexing |
| `drop=TRUE` | 1-column matrices collapse | keep 2-D shape explicitly |
| indexing | 1-based | 0-based; careful with `[1:x]`, `which`, `seq_len` |

Critical spot: `latentcor` builds `cp = rbind(row(R)[lower.tri(R)],
col(R)[lower.tri(R)])` and iterates columns of `cp` — the **order of the lower
triangle** must match (column-major). Since the final `R` is reassembled into a
symmetric matrix, order does not affect the result, but it affects any
per-pair golden comparison, so it is documented.

---

## 8. `synthData_from_ecdf` / `qstepcdf` (only if in scope — see §10 Q3)

Uses `mvtnorm::rmvnorm` (RNG), `stats::ecdf`, `stats::pnorm`, and
`rootSolve::uniroot.all` for inverse-transform sampling, plus `SpiecEasi` graph
utilities in the examples. This is a **simulation utility**, not part of the
network-estimation path. Its output depends on R's RNG and on `uniroot.all`
numerics, so exact equivalence is impossible without RNG bridging; it can only
be "statistically equivalent". Recommendation: implement only if the user wants
it; otherwise mark "out of scope / not ported" in `EQUIVALENCE.md`.

---

## 9. Golden-test plan (matches the requested structure)

`tests/r_reference/generate_golden.R` will, for each case:
1. write `input.csv` (count matrix, full precision) + `params.json`;
2. write **intermediate** outputs at 17 significant digits:
   `mclr` output, Kendall `K`, `Rpointwise`, `R` (post nearPD+nu),
   `lambdaseq`, per-λ selection-frequency matrices `merge[[i]]`,
   `summary` (instability curve), `opt.index`, final adjacency & MB `beta`;
3. write the **exact subsample indices** used by `pulsar` (so Python reuses
   them);
4. write `sessionInfo()`.

Current committed cases use the public QMP dataset from `SPRING-master/data/QMP.rda`:
`qmp_n106_p20` (all 106 samples and 20 low-prevalence taxa) and
`qmp_n20_p30` (20 observed samples and 30 taxa, a real p>n case). Both use
`quantitative=TRUE`, `Rmethod="approx"`, and a data-specific lambda path.
The former synthetic golden files were removed; compositional MCLR, the original
bridge method, and other branches remain covered by in-memory unit/CLI tests.

`tests/test_equivalence.py` will load the goldens, inject the indices, and
compare with `np.testing.assert_allclose` and structural invariants. The QMP
p>n case uses measured data-specific lambda tolerance, exact support/adjacency,
and KKT checks rather than pointwise coefficient equality because the optimizer
is rank-deficient.

---

## 10. Non-equivalence points and open questions

### 10.1 Definitely NOT identical (documented, unavoidable)
1. **RNG.** R's `set.seed`/`sample` ≠ numpy. Identical results require injecting
   the golden indices. Without injection: statistically equivalent only, and the
   selected λ may differ. (This is the reason the indices are injected.)
2. **BLAS/LAPACK.** Different libraries → ~1e-12…1e-15 noise in `cov`, `cor`,
   `nearPD`, normal CDFs.
3. **`optimize` vs `scipy.optimize.brent`.** Same objective, different iteration
   control → agreement ≈ tol (1e-6), not bit-level.
4. **Normal CDFs.** `mnormt::pmnorm` vs `scipy.stats.multivariate_normal.cdf`
   may differ ~1e-10…1e-14 (affects only the `cutoff` region of `r_switch`, i.e.
   |K| > ratio·bound).
5. **nearPD.** Ours is a re-implementation; agreement to be measured.
6. **Lasso at λ boundaries.** A coefficient near the soft-threshold tie can be
   selected in R and not in Python (or vice-versa). Expected to be rare when the
   CD is re-implemented exactly; any occurrence must be reported, not hidden.

### 10.2 Must be identical if the port is correct
`mclr`, Kendall τ (given the same tie handling), the `comb="22"` interpolation
*inputs*, the MB active-set path (given identical CD), selection frequencies,
`summary`, `opt.index`, final adjacency — all should match up to the points in
§10.1.

### 10.3 Open questions for the user (asked before coding)
- **Q1 — PyPI name.** `spring` is already taken on PyPI (HTTP 200). `spring-microbiome`
  and `spring_microbiome` are free. Proposal: distribution `spring-microbiome`,
  import package `spring_microbiome`. OK?
- **Q2 — Interpolation tables.** Ship the `ipol_22` tables extracted from
  `latentcor 2.0.2` in the Python package (bit-faithful inputs) — confirm this is
  acceptable, vs. regenerating them in `scipy` (less faithful, R-version dependent).
- **Q3 — `synthData_from_ecdf` scope.** Port it too (needs `mvtnorm`/`rootSolve`/
  `SpiecEasi` equivalents, RNG-dependent), or declare it out of scope and rely on
  simpler synthetic generators for the recovery tests?
- **Q4 — License.** SPRING (GPL-3) + latentcor (GPL-3) + huge (GPL-2) + pulsar
  (GPL≥3) ⇒ a derivative must be **GPL-3** (or later). I will not choose the
  license: please confirm GPL-3.
- **Q5 — Original vs fork.** Confirm the target is the *modernized fork*
  (`latentcor`-based, in `SPRING-master`), and that the paper's `mixedCCA`
  description is to be recorded only as a difference, not implemented.

---

## 11. Mapping table (R → planned Python), with R source locations

| R function | R file:line | Python module (planned) | Equivalence |
|---|---|---|---|
| `SPRING.default` | `SPRING.R:98` | `spring/spring.py::spring` | orchestration, exact |
| `SPRING` (generic) | `SPRING.R:57` | `spring/spring.py::spring` | n/a (Phyloseq dropped) |
| `.phy2mat` | `SPRING.R:61` | not ported | out of scope |
| `hugeKmb` | `helpers.R:44` | `spring/huge_kmb.py::huge_kmb` | exact + latentcor |
| `mclr` | `helpers.R:71` | `spring/mclr.py::mclr` | exact (1e-15) |
| `synthData_from_ecdf` | `synthData.R:19` | `spring/synth.py` (if Q3=yes) | statistical only |
| `qstepcdf` | `synthData.R` | `spring/synth.py` | statistical only |
| `latentcor::latentcor` | latentcor pkg | `spring/latentcor/latentcor.py` | tolerance-level |
| `Kendalltau`,`n_x` | latentcor `internal.R` | `spring/latentcor/kendall.py` | verify ~1e-12 |
| `r_sol`,`r_switch`,`bound_switch` | latentcor `internal.R` | `spring/latentcor/bridge.py` | ~tol |
| `ipol_22`/`C_evalmlip` | latentcor `src/ipol.c` | `spring/latentcor/interp.py` + `data/ipol_22.npz` | few ulps |
| `Matrix::nearPD` | Matrix pkg | `spring/nearpd.py` | measure |
| `huge::huge.mb` | huge pkg | `spring/huge/mb.py` | exact CD |
| `huge::mb` (C++) | `huge_core.cpp:679` | `spring/huge/mb.py` | exact CD |
| `pulsar::pulsar` | pulsar `mcPulsarSelect.R` | `spring/pulsar.py::pulsar` | exact (indices injected) |
| `stars.stability` | pulsar `utilities.R` | `spring/pulsar.py::stars_stability` | exact |
| `.starsind`, `.ratcheck` | pulsar `utilities.R` | `spring/pulsar.py` | exact |
| `refit.pulsar` | pulsar `refit.R` | `spring/pulsar.py::refit` | exact |

---

## 12. Residual uncertainty register (explicit, per the task)

| # | Node | Why unsure | Mitigation |
|---|------|-----------|-----------|
| U1 | `pcaPP::cor.fk` vs `scipy.stats.kendalltau` | tie handling / weighting may differ | validate on golden `K`; else implement `cor.fk` |
| U2 | `Matrix::nearPD` | no bit-identical numpy equivalent | port Higham faithfully; measure |
| U3 | `optimize` vs scipy `brent` | RESOLVED: R's `Brent_fmin` ported verbatim | exact |
| U4 | `mnormt`/`fMultivar` vs scipy CDF | RESOLVED which is right: independent GL quadrature shows scipy is accurate (~3e-6) and `mnormt` is off by ~7.6e-5 on the worst bridge matrix -> we do NOT reproduce R's error | ~5e-5 in `Rpointwise` for bridge-inversion pairs only; discrete outputs unaffected |
| U5 | `huge` CD summation order | BLAS vs scalar loop | exact same loop order in numpy; residual ~1e-15 |
| U6 | `huge` covariance-input conversion | `S = correlation` branch | verify no-op for diag-1 input |
| U7 | λ boundary ties | soft-threshold equality | report any boundary flip; do not widen tol |
| U8 | `comb` set | assumed always `"22"` | assert at runtime; keep other combs unimplemented-but-erroring |
