# Differences between the R code and the paper

Per the project rule, **the R code is the source of truth**; where code and
paper disagree, the code wins and the disagreement is recorded here.

> **This divergence is declared up front in the package `README.md`** and is a
> required citation note: results from this package are the `latentcor`-based
> fork, not the original `mixedCCA` method.

## 1. Latent correlation: `latentcor` instead of `mixedCCA`

The SPRING paper (Yoon, Gaynanova & Müller, 2019) describes the latent
correlation step via the `mixedCCA` package. The R source ported here (the
modernised fork in `SPRING-master`, `Imports: latentcor`) uses
`latentcor::latentcor(..., types = "tru", method = Rmethod, ...)` instead.

**Consequence for this port:** we implement the `latentcor` behaviour, not the
paper's `mixedCCA` behaviour. This is the single largest paper/code divergence
and it is the reason `docs/PORTING_NOTES.md` analyses `latentcor` in full.

## 2. `SPRING` calls `latentcor` with `types = "tru"` unconditionally

`hugeKmb` (`helpers.R:44`) hard-codes `type = "tru"`, so inside the pipeline
**every** variable is treated as truncated/zero-inflated and the bridge function
used is always the truncated/truncated case (`comb = "22"`). The `quantitative`
flag therefore does *not* change the variable type — it only selects whether
`mclr` is applied first. The paper discusses the truncated type in the context
of zero-inflated counts, which is consistent, but readers often assume the type
is inferred per variable (`latentcor::get_types`). It is not: the pipeline never
calls `get_types`.

## 3. `path` (adjacency) is unsigned

`huge.mb` computes `path[[i]] = abs(beta)` **before** symmetrising with
`sign(path + t(path))`. The result is a **0/1 adjacency**, not a signed one,
even though `huge`'s documentation and common usage return signed paths. The
sign information lives in `beta` (the MB coefficient matrices), not in `path`.
The final `refit$stars` returned by SPRING is therefore unsigned. The README of
the fork works around this by symmetrising `beta` manually
(`(beta + t(beta))/2`) for visualisation — our `weights` field reproduces that.

## 4. Only the negative-value check exists in R

`SPRING.default` validates only `any(data < 0)`. It does not reject NaN/Inf,
all-zero columns or all-zero rows, and it does not check `p` vs `n`. This port
adds stricter validation (`spring_microbiome.validation.validate_counts`), which
is **more** restrictive than R for clearly-degenerate input. It can be disabled
with `spring(..., validate_input=False)` to reproduce R exactly.

## 5. Unused arguments

`SPRING.default` explicitly errors on any `...` argument
(`stop("unused arguments: ...")`). The Python port reproduces this naturally
through the function signature (`TypeError` on unknown keyword arguments) and,
for `huge_kmb`, through an explicit `**kwargs` check.

## 6. `Lambda.min.ratio` semantics

`lambda.min.ratio` only affects the default λ path and the `"data-specific"`
path; if the user passes an explicit `lambdaseq`, `lambda.min.ratio` and
`nlambda` are ignored by the estimator (they are, however, still used to build
the *default* sequence object). We keep that behaviour.

## 7. `ncores` does not change results

`pulsar` parallelises the estimator across subsamples. Because the estimator is
deterministic and the subsamples are independent, `ncores` cannot change the
result. The Python port accepts `ncores` for API compatibility but always runs
serially.
