# Performance notes — SPIEC-EASI MB solver

Measured on the development machine (Linux x86_64), running the Python
coordinate-descent MB solver in `spiec-easi-python`. Timings are for the full
`spiec_easi_mb` pipeline (StARS over all subsamples plus the final fit), with the
defaults `n_lambdas = 20`, `lambda_min_ratio = 1e-3`, `n_subsamples = 20`
(config `app`: `n_lambdas = 25`, `lambda_min_ratio = 1e-2`).

Everything below the "FASE 2" section is the **FASE 1 baseline**, measured before
the compiled core existed, i.e. with the pure-Python coordinate descent.

## FASE 2 — compiled MB core (`spieceasi._mb_core`)

### Environment

| | |
|---|---|
| CPU | AMD Ryzen 5 3550H, 4 cores / 8 threads |
| OS / libc | Linux 6.18.53 x86_64, glibc 2.44 |
| Python | 3.14.7, numpy 2.5.3, BLAS = scipy-openblas |
| R | 4.6.1, SpiecEasi 1.99.0, huge 2.0.1, pulsar 0.3.13 |
| Parallelism | serial on both sides (Python `n_jobs=1`, R without `pulsar` parallelism) |

### Full pipeline, serial, median of 3 runs

`benchmarks/bench_python.py --reps 3 --n-jobs 1`, same golden cases and the same
injected R subsample indices as the FASE 1 baseline:

| case | p | n | R median (s) | Python ref median (s) | Python C++ median (s) | C++ vs ref | C++ vs R |
|------|---|---|--------------|-----------------------|-----------------------|------------|----------|
| `ds2_p30_n100` | 30 | 100 | 2.108 | 97.819 | **0.572** | 171× | **0.27×** |
| `rand_p50_n50_s211` | 50 | 50 | 4.450 | 453.408 | **3.991** | 114× | **0.90×** |
| `largep_p100_n80_s215` | 100 | 80 | 14.040 | 1565.844 | **22.056** | 71× | **1.57×** |

Ratio < 1 means Python is faster than R. Versus the FASE 1 baseline the gap to R
goes from 46.6× / 100.8× / 111.8× slower to 0.27× / 0.90× / 1.57×.

*Measured 2026-10-13, Linux x86-64 (GCC), single thread on both sides, **before** the
row-major change.* Three later runs of the same cases exist and give different ratios,
which is why the repository quotes 1.57×, 1.52× and 1.54× for `p = 100`: after adopting
the row-major layout the serial pipeline is 0.514 / 3.691 / **21.348** s (ratios 0.24× /
0.83× / **1.52×**), and the `n_jobs = 1` row of the FASE 3 (a) table below is **21.619** s
(**1.54×**). They are three runs of the same case on the same machine, not a
disagreement — always quote the ratio together with the build it belongs to.

**The measured outputs are unchanged.** For all three cases the selected lambda
and the edge count are bit-for-bit identical between the reference and the C++
solver, and identical to R:

| case | R | Python ref | Python C++ |
|------|---|------------|------------|
| `ds2_p30_n100` | λ=0.2970674324372107, 0 edges | λ=0.2970674324372107, 0 edges | λ=0.2970674324372107, 0 edges |
| `rand_p50_n50_s211` | λ=0.36303491242767955, 9 edges | λ=0.36303491242767966, 9 edges | λ=0.36303491242767966, 9 edges |
| `largep_p100_n80_s215` | λ=0.33582686948865675, 22 edges | λ=0.3358268694886568, 22 edges | λ=0.3358268694886568, 22 edges |

### Test-suite wall time

| suite | before (pure Python) | after (C++) |
|-------|----------------------|-------------|
| non-slow suite (`pytest tests/ -m "not slow"`) | 662 s | **18.66 s** |

The slow suites are dominated by the R runs (rpy2) plus the Python pipeline:
after the extension the Python side is no longer the bottleneck of those runs
(e.g. `tests/test_r_variability_expanded_scenarios.py` went from being dominated
by the Python solver to 343.56 s, almost entirely R).

### Where the p = 100 time goes (measured, not inferred)

`largep_p100_n80_s215`, `n_jobs = 1`, 20 subsamples + final refit, extension and
correlation computation wrapped with timers:

| stage | time | share |
|-------|------|-------|
| compiled core `_mb_core.solve_mb_path`, 21 calls (20 subsamples + final refit) | 22.987 s | **99.7 %** |
| correlation matrix `huge_correlation`, 21 calls | 0.041 s | 0.2 % |
| everything else in Python (CLR, `symmetrize_path`, StARS aggregation, index selection, graph conversion) | 0.034 s | 0.1 % |
| total | 23.062 s | 100 % |

`cProfile` agrees: 23.739 s inside `_mb_core.solve_mb_path`, of which
`stars_select` accounts for 22.003 s (the 20 subsample paths) and the final fit
1.779 s; `symmetrize_path` 0.012 s, graph conversion 0.003 s.

So the gap to R at p = 100 (1.57×) is **entirely inside the coordinate-descent
loop**. It is not Python overhead, not adjacency conversion, not StARS
aggregation and not array copies: those together are 0.1 % of the runtime, i.e.
they cannot be optimised into anything measurable. Inside the loop the work is
identical to R's (see the loop counters below) and the difference in wall time is
R's OpenMP parallelism (see "the R comparison is not single-thread vs
single-thread" below).

### Memory layout and `-O3` (measured)

Four builds of the same source were compiled directly with `g++` (same command,
only the listed switches differ) and each was checked for **bit-identity**
against the shipped build on 7 inputs (`ds1..ds4` CLR, the `p = 100` CLR, one
`plain` and one `illcond` synthetic case), comparing the full coefficient paths:

| build | bit-identical to shipped | max \|diff\| | end-to-end p = 100 |
|-------|--------------------------|--------------|--------------------|
| `-O2`, column-major S (shipped) | — | — | 23.57 s |
| `-O2`, column-major S (rebuilt) | 7/7 | 0.0 | 23.24 s |
| `-O3`, column-major S | 7/7 | 0.0 | 23.14 s |
| `-O2`, row-major S | 7/7 | 0.0 | 22.00 s |
| `-O3`, row-major S | 7/7 | 0.0 | 21.82 s |

* **Layout.** `huge` stores S column-major and the inner loop reads a *row*
  (`S[j, a_k]` for consecutive `a_k`), i.e. stride `d` in a column-major buffer.
  Storing S row-major and adjusting the accessor accesses the *same* element with
  the *same* summation order, so the result is bit-identical, and the row becomes
  contiguous: **≈ 5.3 % faster** (23.24 s → 22.00 s).
* **`-O3`.** Worth only ≈ 0.4–0.8 % here (23.24 s → 23.14 s), but again
  bit-identical on all 7 inputs.
* **The forced C→F copy is irrelevant.** `np.asfortranarray(S)` costs 5.2 µs per
  call, i.e. 0.109 ms over the 21 calls (≈ 0.0005 % of the runtime). The
  row-major build helps through the access pattern, not through the copy.

Combined, the two switch-level changes recover ≈ 6 % of the p = 100 runtime
while remaining bit-identical.

**What was adopted.** The **row-major layout is adopted** (commit `e73ec82`):
same element access and same summation order, bit-identity verified on all 248
problems (244 sweep + 4 goldens) against the pre-change build, and the full slow
suite re-run afterwards: `50 passed, 84 deselected in 1445.48 s`. Measured on the
three benchmark cases (median of 3, serial): 0.572 → 0.514 s (p=30), 3.991 →
3.691 s (p=50), 22.056 → 21.348 s (p=100), i.e. the ratio to R becomes
0.24× / 0.83× / 1.52×.

**`-O3` is not adopted.** Rebuilt from the final source and compared on the
p = 100 pipeline (3 runs each): shipped `-O2` min 21.560 s / median 21.696 s,
`-O3` min 21.303 s / median 21.670 s — 0–1.2 %, inside run-to-run noise, while
it would widen the deviation from R's build flags. The identity gate passes
(248/248, max |diff| = 0.0), so it stays available as a one-line change if a
future build wants it.

**At very large `p` the row-major gain grows (information only, no code change).**
The `huge` source carries an explicit warning about this access pattern:

> *“NOTE: S(j, idx_a[k]) reads look stride-d, but with j advancing in the caller's
> outer loop each fixed k walks its column sequentially — size_a parallel
> prefetch streams. Rewriting to stay within column j (using symmetry) measured
> 2x SLOWER at d=2000; do not "fix" this access pattern.”*

That note is about a different rewrite (keeping the access inside column `j` by
using the symmetry of S); our change keeps the same element and the same
summation order and only changes the storage layout. Re-measured anyway, one run
each, both builds bit-identical, single-threaded:

| case | column-major | row-major | row/col | bit-identical |
|------|--------------|-----------|---------|---------------|
| p=300, n=100 | 12.697 s | **10.970 s** | 0.864 | yes |
| p=500, n=100 | 40.903 s | **30.064 s** | 0.735 | yes |

So at p = 300 and p = 500 the row-major layout is 1.16× and 1.36× faster,
i.e. the advantage *increases* with `d` on our code path — the opposite of the
effect the comment warns about. These cases are outside the verified range (the
goldens and the sweep stop at p = 100) and are reported as information only.

### The R comparison is not single-thread vs single-thread: `huge` uses OpenMP

The installed `huge` 2.0.1 is compiled **with OpenMP** (`libs/huge.so` links
`libgomp`, `src/Makevars.in` sets `OPENMP_CXXFLAGS`), and its `mb()` parallelises
the `d` nodewise regressions:

```
#pragma omp parallel
#pragma omp for schedule(dynamic)
for (int m = 0; m < d; m++) { ... }
```

Our core is deliberately single-threaded (no OpenMP over `m`), so the end-to-end
ratios versus R are not like-for-like. Isolating the solver: the same correlation
matrix S and the same lambda path, handed to huge's own C++ `mb` through
`.Call("_huge_SPMBgraph", S, lambda, nlambda, d)` (so R's `cor` is not involved),
against our core on the same input:

| case | R `huge::mb`, 1 thread | R `huge::mb`, default (8 threads) | our core, 1 thread |
|------|------------------------|-----------------------------------|--------------------|
| p=30 | 0.0274 s | **0.0176 s** (1.56×) | 0.0340 s |
| p=50 | 0.3771 s | **0.1755 s** (2.15×) | 0.3579 s |
| p=100 | 1.8621 s | **0.8732 s** (2.13×) | 1.6639 s |

Single-threaded, our core is therefore **on par with R's** (slightly faster at
p=50/p=100, 1.24× slower at p=30, single runs). The residual gap at p=100
(1.52×) is explained by R's OpenMP factor (~2.1× on the `mb` call), not by a
different amount of work — see the counters below.

### Is the work the same? Loop counters, R vs C++

Both sources were rebuilt with identical counters (they only observe the loops;
the instrumented `huge` returns results bit-identical to the official build,
verified on p=30). Counted on the same S and the same path, single-threaded:

| counter | p=30 R / ours | p=50 R / ours | p=100 R / ours |
|---------|---------------|---------------|----------------|
| active-set sweeps | 31099 / **31099** | 169646 / **169646** | 284007 / **284007** |
| coordinate updates | 862432 / **862432** | 7636682 / **7636682** | 21656514 / **21656514** |
| residual evaluations | 8457 / **8457** | 40575 / **40575** | 213966 / **213966** |
| residual terms (`sum size_a`) | 83287 / **83287** | 844344 / **844344** | 9138629 / **9138629** |
| outer (strong-rule) iterations | 692 / **692** | 1399 / **1399** | 2283 / **2283** |
| strong-rule rounds | 1031 / **1031** | 2651 / **2651** | 5387 / **5387** |
| certification passes | 692 / **692** | 1399 / **1399** | 2283 / **2283** |
| certification violations | 146 / **146** | 590 / **590** | 422 / **422** |
| lambda iterations | 600 / **600** | 1000 / **1000** | 2000 / **2000** |

Every counter matches exactly: the compiled core performs **the same work** as
`huge` on the same input. There is no "extra sweep" and therefore no first
divergence point to report in the loop structure; the remaining difference in
wall time is threads (`-O3`/layout aside).

### FASE 3 (a) — thread-parallel subsamples

The StARS subsamples are independent and the compiled core releases the GIL, so
the default joblib backend is now **`threads`** (was `processes`; the legacy
backend is still available). A thread pool cannot re-execute the interpreter,
which is what broke the frozen build with the *loky* backend.

**Determinism.** Results are **bit-identical** for `n_jobs = 1, 2, 4, 8`:
adjacency, weights, directed coefficients, instability, monotone instability,
alpha path and the selected lambda all compare equal with `np.array_equal`. The
subsample adjacency paths are stacked in submission order and summed in that
order, and no global mutable state is involved (each subsample owns its scratch).
Covered by `tests/test_n_jobs_determinism.py`; the benchmark payloads (edge count
and selected lambda) are also identical for every worker count tested (1…16).

**Speed-up** (median of 3, same cases as above; "×" is versus `n_jobs = 1`):

| n_jobs | ds2 p=30 | × | rand p=50 | × | largep p=100 | × |
|--------|----------|---|-----------|---|--------------|---|
| 1 | 0.516 s | 1.00× | 3.804 s | 1.00× | 21.619 s | 1.00× |
| 2 | 0.291 s | 1.77× | 2.202 s | 1.73× | 13.055 s | 1.66× |
| 4 | 0.239 s | 2.16× | 1.497 s | 2.54× | **8.734 s** | **2.48×** |
| 8 | 0.196 s | 2.63× | 1.250 s | 3.04× | 6.861 s | 3.15× |
| 16 | 0.204 s | 2.53× | 1.149 s | 3.31× | 6.505 s | 3.32× |

Machine: 4 physical cores / 8 logical (SMT). `n_jobs = 8` (logical count) still
beats 4; `n_jobs = 16` (twice the logical count) is *not* harmful here — the
subsample costs are heterogeneous, so a larger pool smooths the tail — but it is
over-subscription and is deliberately not a default.

**End-to-end versus R** (R runs `spiec.easi`+pulsar and its `huge` is
OpenMP-parallel internally; ours uses threads). Measured 2026-10-13, Linux x86-64,
4 physical / 8 logical cores, median of 3 runs:

| case | R | ours n_jobs=1 | ours n_jobs=4 | ours n_jobs=8 |
|------|---|---|---|---|
| p=30 | 2.108 s | 0.516 s (0.24×) | 0.239 s (0.11×) | 0.196 s (0.09×) |
| p=50 | 4.450 s | 3.804 s (0.85×) | 1.497 s (0.34×) | 1.250 s (0.28×) |
| p=100 | 14.040 s | 21.619 s (1.54×) | **8.734 s (0.62×)** | **6.861 s (0.49×)** |

So the serial-to-serial comparison was 1.54× slower than R at p=100; with threads
we are 1.6× faster at `n_jobs = 4` and 2.05× faster at `n_jobs = 8`. The
requested p=100, n=80 measurement with `n_jobs = 4` is **8.734 s**.

**Default worker count.** `spieceasi.runtime.default_n_jobs()`: physical cores
when detectable (Linux `/proc/cpuinfo`, macOS `sysctl hw.physicalcpu`, Windows
`wmic`/PowerShell CIM, or `psutil` if present); otherwise `logical // 2`
(assume SMT, stay on the safe side); otherwise 1. The desktop app uses it and
`spieceasiJobs` overrides it.

### What was measured, and what was not






* Timings above are wall-clock medians of 3 runs on an idle machine, serial.
  Individual slow-suite timings elsewhere in this repository include R and are
  not comparable to these numbers.
* `-march=native` was deliberately **not** used and LTO is explicitly **off**
  (see flags); a differently tuned build could be faster. The numbers are for
  the declared build, not for the theoretical best.

### Compilation flags — ours vs R's (declared deviations)

Read from `R CMD config` on the same machine:

| | flags |
|---|---|
| R (distro build, used for `huge`, `pulsar`, `SpiecEasi`) | `-march=x86-64 -mtune=generic -O2 -pipe -fno-plt -fexceptions -Wp,-D_FORTIFY_SOURCE=3 … -fno-omit-frame-pointer … -g -ffile-prefix-map=… -flto=auto -ffat-lto-objects` |
| `spieceasi._mb_core` | `-std=c++17 -O2 -ffp-contract=off -fno-lto` (`/O2 /fp:precise` on MSVC) |

Deliberate deviations, all conservative:

1. **No LTO.** R builds with `-flto=auto -ffat-lto-objects`; the extension is
   compiled with `-fno-lto`. LTO could reorder/merge the scalar accumulation
   loops; keeping it off makes the build reproducible from the flags alone and
   keeps the numerics of the translation traceable to the C++ source.
2. **`-ffp-contract=off`.** GCC's default is `-ffp-contract=fast`, and R's flag
   list does not override it, so R's `huge` may use FMA in the accumulation
   loops while the extension does not. This is the conservative direction: an
   FMA changes rounding, so enabling it would make the result both
   compiler-dependent and farther from the Python reference.
3. **No `-ffast-math`, no `-march=native`, no `-mtune=native`.** These reorder
   floating-point operations and/or change the ISA baseline, and would break the
   parity results. `-O2` matches R's optimisation level.
4. **Profiling deviation.** The plan asked for `line_profiler` on the `p = 100`
   subsample; `line_profiler` was run on `ds2_p30_n100` instead (the `p = 100`
   subsample was profiled with `cProfile`, which reports per-function totals and
   does not need per-line instrumentation). The `cProfile` result was the one
   that drove the decision: 63.9 M calls in a single `p = 100` subsample, with
   the coordinate-descent loop body at 76.9 s tottime, 15.96 M `ndarray.dot`,
   15.98 M `soft_threshold` and 31.7 M `abs` calls.

### Extension-vs-reference verification (sweep)

`tests/test_mb_core_equivalence.py -m slow` compares the compiled core with the
pure-Python reference on 244 problems (p = 5…60, n = 10…150; plain, with zeros,
and ill-conditioned; 6 lambdas each) plus the released goldens:

| metric | result |
|--------|--------|
| problems | 244 |
| bit-identical to the reference | **0 / 244** |
| max abs difference, non-degenerate families | 9.470e-14 |
| max abs difference, degenerate family (relative collinearity 1e-9) | 1.881e-01 (worst: `illcond p=8 n=150 seed=202`) |
| support (non-zero pattern) differences, non-degenerate families | **0 problems** |
| support differences, degenerate family | 25 problems, 40 lambda blocks |
| max abs difference vs R goldens | 1.652e-13 |

So the compiled core is **not bit-identical** to the pure-Python reference (nor
to R): the differences are at the 1e-13…1e-14 level and the discovered structure
is identical everywhere except in a deliberately degenerate family, where the
matrix is numerically singular and the lasso solution is not unique — R does not
single out either implementation there. This is reported, not hidden, and no
tolerance was widened to obtain it.

## FASE 1 baseline — end-to-end case times (`pytest tests/test_equivalence_extended.py -m slow`)

| case | p | n | config | wall time |
|------|---|---|--------|-----------|
| rand_p10_n30_s201 | 10 | 30 | r | 9.4 s |
| rand_p10_n50_s202 | 10 | 50 | r | 5.2 s |
| rand_p10_n100_s203 | 10 | 100 | r | 4.8 s |
| rand_p20_n30_s204 | 20 | 30 | r | 92.0 s |
| rand_p20_n50_s205 | 20 | 50 | r | 52.8 s |
| rand_p20_n100_s206 | 20 | 100 | r | 37.8 s |
| rand_p30_n30_s207 | 30 | 30 | r | 204.3 s |
| rand_p30_n50_s208 | 30 | 50 | r | 232.6 s |
| rand_p30_n100_s209 | 30 | 100 | r | 127.1 s |
| rand_p50_n30_s210 | 50 | 30 | r | 262.4 s |
| rand_p50_n50_s211 | 50 | 50 | r | 619.2 s |
| rand_p50_n100_s212 | 50 | 100 | r | 693.5 s |
| sparse_p30_n100_s213 | 30 | 100 | r | 118.2 s |
| sparse_p50_n60_s214 | 50 | 60 | r | 1123.1 s |
| pn_zeros_p40_n25_s219 | 40 | 25 | r | 172.7 s |
| ratcheck_p30_n200_s217 | 30 | 200 | r | 108.5 s |
| ratcheck_p50_n300_s218 | 50 | 300 | r | 464.7 s |
| app_p10_n30_s201 | 10 | 30 | app | 9.4 s |
| app_p30_n100_s209 | 30 | 100 | app | 127.3 s |
| app_ratcheck_p30_n200_s217 | 30 | 200 | app | 124.6 s |

## p = 100 (measured)

One measured subsample path (basis for the projection):

| case | dim | one subsample | projected pipeline only |
|------|-----|---------------|-------------------------|
| largep_p100_n80_s215 | 80 × 100 | 78.1 s | ≈ 27 min |
| app_largep_p100_n80_s215 | 80 × 100 | 48.6 s | ≈ 17 min |
| largep_p100_n200_s216 | 200 × 100 | 188.7 s | ≈ 66 min |

Measured end-to-end runtime of the parity test case
(`pytest tests/test_equivalence_extended.py -m slow -k p100`) on this machine:

| case | dim | status | wall time |
|------|-----|--------|-----------|
| largep_p100_n80_s215 | 80 × 100 | OK | 2425.5 s (40.4 min) |
| largep_p100_n200_s216 | 200 × 100 | OK | 4933.3 s (82.2 min) |
| app_largep_p100_n80_s215 | 80 × 100 | OK | 1367.3 s (22.8 min) |

Full run: `3 passed, 43 deselected in 8726.93 s (2:25:26)`.

**After FASE 2 — same cases, same test, compiled solver:**

| case | dim | status | wall time (C++) | speed-up |
|------|-----|--------|-----------------|----------|
| largep_p100_n80_s215 | 80 × 100 | OK | 29.80 s | 81× |
| largep_p100_n200_s216 | 200 × 100 | OK | 81.84 s | 60× |
| app_largep_p100_n80_s215 | 80 × 100 | OK | 17.51 s | 78× |

Full run: `3 passed, 67 deselected in 129.23 s (0:02:09)` — a 67.5× reduction,
with all steps still matching (`final_edges = 0`, 0 differing edges in the
per-subsample solver path).

Note: the parity test solves each case **twice** — the solver-isolation loop over
all subsamples plus the full pipeline with injected subsamples — so the end-to-end
test time is roughly twice the pipeline-only projection.

`sub0_edge_diff = 0` for all three, i.e. the per-subsample solver output already
matches R at p = 100.

## Start-up time of the frozen backend (measured)

`benchmarks/bench_startup.py`, tiny case p=5 n=20, warm median of 4–5 launches
(no other load on the machine):

| metric | FASE 1 (pure-Python solver) | + compiled core | + compiled core, lazy `scipy.stats` |
|--------|------------------------------|-----------------|--------------------------------------|
| frozen, startup only | 0.092 s | 0.108 s | **0.103 s** |
| frozen, analyze tiny | 10.740 s | 2.444 s | **1.504 s** |
| source, startup only | 0.018 s | 0.038 s | **0.037 s** |
| source, analyze tiny | 5.476 s | 1.786 s | **1.734 s** |
| frozen − source (analyze) | +5.264 s | +0.658 s | **−0.231 s** |

**Where the ~10 s went.** About **7.7 s of it was the pure-Python solver**, not
start-up: compiling the MB core brought `frozen analyze` from 10.740 s to 2.444 s.
A further **0.94 s was the eager `from scipy.stats import spearmanr`** — measured
at 0.835 s on its own — which the SPIEC-EASI path never uses; it is now imported
inside the helper that needs it, so `frozen analyze` is 1.504 s in total. There
is no longer any measurable PyInstaller penalty: the frozen run is *faster* than
the source run (−0.23 s) because pandas, which `spieceasi.transforms` imports as
an optional convenience, is correctly absent from the bundle while the source
environment has it (0.225 s of the source import time).

**Import profile of an analyze run** (`-X importtime`, source, µs cumulative):
`spieceasi` 432 · `spieceasi.pipeline` 432 · `spieceasi.transforms` 226 (all of
it pandas) · `networkx` 148 · `analysis` 121 · `numpy` 71 · `spieceasi.stars` 54
(including joblib). `scipy.stats` is no longer on the path (it was 834).

**one-file vs one-folder.** The shipped bundle is one *directory*
(`EXE(exclude_binaries=True)` + `COLLECT`). A one-file build of the identical
spec was measured on the same machine and tiny case:

| layout | startup only | analyze tiny |
|--------|--------------|--------------|
| one-folder (shipped) | **0.094 s** | **1.499 s** |
| one-file (58.5 MB single binary) | 0.783 s | 2.200 s |

One-file pays **+0.69 s per launch** for unpacking itself into a temporary
directory, i.e. +46 % on this job. One-folder stays.

**What is still avoidable.** Only the residual import cost of `numpy` (71 ms),
`networkx` (148 ms) and the `spieceasi` package (432 ms, of which 225 ms is
pandas in the source environment only). Removing them would mean dropping
functionality, so no further start-up work is planned. Overall: of the original
≈10 s frozen launch-and-run, **≈9.2 s is gone** and ~1.5 s remains, of which
0.10 s is the PyInstaller bootloader itself.

## Known limitations / future work

State after FASE 2: the coordinate-descent solver is no longer the dominant cost
of the Python pipeline. The remaining gap to R at `p = 100` (1.57×) is inside the
solver loop; the rest of the pipeline (CLR, correlations, lambda path, StARS
subsampling, selection, refit) is still pure Python and is now a visible share of
the total for small problems.

Known, still-open items:

1. **Parallelism over subsamples — done (FASE 3 a).** The default backend is now
   `threads`, results are bit-identical for `n_jobs = 1…16`, and the p=100 case
   goes from 21.6 s serial to 8.7 s at four workers; see "FASE 3 (a)" above.
   A future option (not exposed) is to default to the *logical* core count on SMT
   machines: measured 1.27× more at p=100 (6.86 s at 8 workers vs 8.73 s at 4).
2. **OpenMP over neighbourhoods.** Not implemented — see the dedicated
   "Future work: OpenMP" section below for the proposal, its bit-identity gates
   and the per-platform packaging work it needs.
3. **`collect_sorted` / sequential strong rule.** `huge` collects and sorts the
   active-set candidates; the translation follows the simpler variable-order
   scan. Same fixed point, different order — it is a declared deviation, not a
   bug, but it is why results are not bit-identical to R. The loop counters show
   the *amount* of work is identical, so this does not cost time on our side.
4. **Start-up time of the packaged app.** See the startup section below for the
   measured breakdown and what was removed.

Earlier hypothesis that is now obsolete: replacing BLAS with scalar loops was
considered to remove the `lambda = lambda.max` indeterminacy. It is not needed —
at `lambda = lambda.max` the result is indeterminate in R as well (see
`docs/AUDIT.md` §2.E), and Python returns the empty model there by choice.

Any further optimisation must keep the R parity results in
`tests/test_equivalence.py`, `tests/test_equivalence_extended.py` and
`tests/test_app_integration.py` unchanged.

## Future work: OpenMP over neighbourhoods (proposal — not implemented)

**Why.** `huge` parallelises the `d` nodewise regressions with
`#pragma omp parallel` + `#pragma omp for schedule(dynamic)` and gets 1.56× (p=30)
to 2.13× (p=100) on the `mb` call itself; the measured loop counters are
*identical* to ours, so this is the only remaining source of the serial gap. The
`d` problems are independent — no reduction crosses `m` — so the output should be
bit-identical; the only shared state is `hit_max_iter` (a single bool) and the
per-`m` scratch, which must become thread-local instead of per-`m`.

**Bit-identity gates (explicit, not assumed).** All three must run 1 vs 2/4/8
OpenMP threads and compare with `np.array_equal` (and against the pure-Python
reference):

1. the 244-problem sweep plus the ds1–ds4 goldens through the core;
2. the end-to-end golden pipeline (ds1–ds4, configs `r` and `app`) and the p=100
   case;
3. the whole slow suite, run twice (`OMP=1` and `OMP=8`).

The thread count is set programmatically (`spieceasi.openmp.set_threads(n)` /
`openmp_threads()`, backed by `omp_set_num_threads` so it does not depend on
environment variables), and a new slow test module drives the gates.

**Packaging per platform.**

| platform | OpenMP runtime | proposal |
|----------|----------------|----------|
| Linux (manylinux/musllinux) | `libgomp` (GCC) | link it **statically** (`-Wl,-Bstatic -lgomp -Wl,-Bdynamic`) so the wheel does not depend on the system copy; otherwise let `auditwheel` bundle it |
| macOS Intel + arm64 | `libomp` (Apple clang has no `-fopenmp`) | explicit `-I/-L .../libomp -lomp` against the keg-only Homebrew install; arm64 wheel must be built on an arm64 runner (no cross-compiling); `delocate-wheel` copies `libomp.dylib` into the wheel |
| Windows | MSVC `/openmp` → `vcomp140.dll` | bundle it with `delvewheel` (or ship it next to the app) |

**Mandatory serial fallback.** `setup.py` compiles a probe with the OpenMP flag;
if it fails, or if `SPIECEASI_NO_OPENMP=1`, the **same source** is built without
OpenMP (`#ifdef _OPENMP` around the pragmas) and `openmp_threads()` returns 1. A
CI job builds that variant and runs the extension-vs-reference tests, so the
serial path is always exercised. For PyInstaller, the static-link option means
nothing extra to collect; a dynamic runtime would need a hook
(`collect_dynamic_libs("spieceasi")`) and is caught by the existing frozen test.

**Anti-oversubscription rule.** `n_jobs × openmp_threads ≤ physical cores`, with
the product clamped (and logged) if a user forces inconsistent values. Proposed
default: the app exposes a single worker budget (`spieceasiJobs`, default =
physical cores) and derives `n_jobs = P, threads = 1`; presets `(P,1)`, `(P/2,2)`,
`(P/4,4)`, `(1,P)`, plus `spieceasiOpenmpThreads` for advanced use. The matrix to
measure is **4×1, 2×2, 1×4** (all = 4 ≤ P) and **2×4, 4×4** (over-subscribed,
expected to degrade), on p=100 and p=50, with `OMP=1, n_jobs=8` as today's
reference point. R's own configuration is effectively **1×8** (pulsar serial,
`huge` with 8 OpenMP threads).

**Possible future option (not exposed).** On SMT machines the *logical* core
count is a better worker budget for the thread-parallel subsamples alone:
measured 1.27× more at p=100 (6.861 s with 8 workers vs 8.734 s with 4) on this
4-physical/8-logical CPU. It is deliberately not a default because it
over-subscribes the FP units and because with OpenMP enabled the product rule
above starts from the physical count.

