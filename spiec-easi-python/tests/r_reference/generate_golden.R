#!/usr/bin/env Rscript
# =============================================================================
# generate_golden.R — Phase 2 golden-file generator for SPIEC-EASI (MB).
#
# Runs the REFERENCE R stack (SpiecEasi 1.99.0 + huge 2.0.1 + pulsar 0.3.13)
# on deterministic synthetic count matrices and writes every intermediate step
# needed to audit the Python translation (spiec-easi-python/src/spieceasi).
#
# Usage:
#   Rscript tests/r_reference/generate_golden.R [output_dir]
#   (default output_dir = tests/golden relative to the current directory)
#
# Precision: all floating-point values are written with %.17g (round-trips an
# IEEE-754 double exactly). Subsample indices are written 0-based (Python
# convention); R uses 1-based internally.
#
# Datasets: p=10/n=50 with zeros; p=30/n=100 with zeros; p>n (p=30/n=15);
#           p=10/n=50 WITHOUT zeros.
# Configs:  "r"   = pure R defaults (nlambda=20, lambda.min.ratio=1e-3, rep.num=20)
#           "app" = application params (nlambda=25, lambda.min.ratio=1e-2, rep.num=20)
# =============================================================================

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(pulsar)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1] else "tests/golden"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
seldir <- file.path(outdir, "selection")
dir.create(seldir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Precision-safe writers
## ---------------------------------------------------------------------------
write_num_matrix <- function(m, path) {
  m <- as.matrix(m)
  txt <- matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m))
  write.table(txt, path, sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}
write_num_vec <- function(v, path) {
  writeLines(sprintf("%.17g", as.numeric(v)), path)
}
write_int_matrix <- function(m, path) {
  write.table(m, path, sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}
write_params <- function(kv, path) {
  writeLines(paste0(names(kv), "=", vapply(kv, as.character, character(1))), path)
}

## ---------------------------------------------------------------------------
## Deterministic synthetic data
## ---------------------------------------------------------------------------
make_counts <- function(n, p, depth, seed, want_zeros) {
  set.seed(seed)
  k <- max(1L, p %/% 4L)
  W <- matrix(abs(rnorm(p * k, sd = 1)), p, k)
  Z <- matrix(rnorm(n * k), n, k)
  logmu <- scale(Z %*% t(W))
  mu <- exp(as.numeric(logmu) * 1.0)
  mu <- matrix(pmax(mu, 1e-6), n, p)
  counts <- matrix(0L, n, p)
  for (i in seq_len(n)) {
    counts[i, ] <- as.integer(rmultinom(1L, depth, mu[i, ] / sum(mu[i, ])))
  }
  if (want_zeros) {
    mask <- matrix(runif(n * p) < 0.15, n, p)
    counts[mask] <- 0L
  } else {
    counts[counts < 1L] <- 1L
  }
  # Guarantee every column has variation (huge rejects constant columns).
  for (j in seq_len(p)) {
    if (length(unique(counts[, j])) == 1L) {
      counts[1L, j] <- counts[1L, j] + 1L
    }
  }
  counts
}

datasets <- list(
  ds1 = list(p = 10L, n = 50L,  want_zeros = TRUE,  depth = 3000L,   seed = 101L),
  ds2 = list(p = 30L, n = 100L, want_zeros = TRUE,  depth = 6000L,   seed = 102L),
  ds3 = list(p = 30L, n = 15L,  want_zeros = TRUE,  depth = 800L,    seed = 103L),
  ds4 = list(p = 10L, n = 50L,  want_zeros = FALSE, depth = 200000L, seed = 104L)
)

configs <- list(
  r   = list(nlambda = 20L, ratio = 1e-3, rep.num = 20L),
  app = list(nlambda = 25L, ratio = 1e-2, rep.num = 20L)
)

threshold <- 0.05
stars_seed <- 17L

## Reproduce the EXACT subsample indices pulsar draws for a given seed.
## Mirrors pulsar::pulsar:
##   subsample.ratio <- .ratcheck(subsample.ratio, n)
##   set.seed(seed); replicate(rep.num, sample(1:n, floor(n*ratio), replace=FALSE))
pulsar_subsamples <- function(n, rep.num, seed) {
  ratio <- pulsar:::.ratcheck(NULL, n)
  set.seed(seed)
  ind <- replicate(
    rep.num,
    sample(c(1L:n), floor(n * ratio), replace = FALSE),
    simplify = FALSE
  )
  do.call(rbind, ind)
}

## Cross-check the reproduced indices against those actually consumed by pulsar
## via an instrumented `fun` (best-effort: may be skipped if mclapply forks).
capture_pulsar_subsamples <- function(X, fargs, rep.num, seed) {
  rec <- new.env()
  rec$idx <- list()
  cap_fun <- function(x, ...) {
    rec$idx[[length(rec$idx) + 1L]] <- as.integer(rownames(x))
    SpiecEasi::sparseiCov(x, ...)
  }
  ok <- tryCatch({
    suppressWarnings(suppressMessages(
      pulsar::pulsar(
        data = X, fun = cap_fun, fargs = fargs,
        criterion = "stars", thresh = threshold,
        rep.num = rep.num, seed = seed, ncores = 1
      )
    ))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  n <- nrow(X)
  idx <- Filter(function(v) length(v) != n, rec$idx)
  if (length(idx) != rep.num) return(NULL)
  do.call(rbind, idx)  # 1-based
}

## ---------------------------------------------------------------------------
## Main loop
## ---------------------------------------------------------------------------
manifest <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  X <- make_counts(ds$n, ds$p, ds$depth, ds$seed, ds$want_zeros)
  rownames(X) <- as.character(seq_len(nrow(X)))
  colnames(X) <- paste0("t", seq_len(ncol(X)))
  n <- nrow(X); p <- ncol(X)

  X_clr <- SpiecEasi:::.spiec.easi.norm(X)   # t(clr(X + 1, 1))

  for (cfg_name in names(configs)) {
    cfg <- configs[[cfg_name]]
    case <- sprintf("%s_%s", ds_name, cfg_name)
    message("== generating ", case, " (n=", n, ", p=", p, ") ==")

    ## --- lambda.max and grid (spiec.easi.default convention) --------------
    lambda_max <- pulsar::getMaxCov(cor(X_clr))
    lambda_path <- pulsar::getLamPath(
      lambda_max, lambda_max * cfg$ratio, cfg$nlambda, log = TRUE
    )

    ## --- fixed lambda: isolate the lasso (no StARS) -----------------------
    fixed_index <- max(1L, cfg$nlambda %/% 2L)      # 1-based
    fixed_lambda <- lambda_path[[fixed_index]]
    iso <- suppressWarnings(suppressMessages(
      SpiecEasi::sparseiCov(X_clr, method = "mb", lambda = fixed_lambda,
                            verbose = FALSE)
    ))
    # huge stores beta as [predictor, node] (see SPMBgraph.cpp: column = node m,
    # row-within-block = predictor j). We transpose to [node, predictor] so the
    # golden matches the Python `directed_coefficients` convention
    # (directed[node, predictor]) and the "coefficients per node" reading.
    coef_fixed <- t(as.matrix(iso$beta[[1]]))       # d x d, [node, predictor]
    adj_fixed  <- as.matrix(iso$path[[1]])          # symmetrized ("or")

    ## --- full pipeline via spiec.easi ------------------------------------
    fit <- suppressWarnings(suppressMessages(
      spiec.easi(
        X, method = "mb", sel.criterion = "stars",
        nlambda = cfg$nlambda, lambda.min.ratio = cfg$ratio,
        pulsar.params = list(rep.num = cfg$rep.num, thresh = threshold,
                             seed = stars_seed),
        verbose = FALSE
      )
    ))

    stars <- fit$select$stars
    grid <- as.numeric(fit$lambda)
    opt_index <- as.integer(stars$opt.index)        # 1-based
    opt_lambda <- grid[[opt_index]]
    monotone <- as.numeric(stars$summary)           # already cummax'd by pulsar
    merge_list <- stars$merge                       # length nlambda, d x d P

    ## raw instability per lambda, reconstructed from the averaged P matrices
    raw_instability <- vapply(merge_list, function(Pm) {
      Pm <- as.matrix(Pm)
      4 * sum(Pm * (1 - Pm)) / (p * (p - 1))
    }, numeric(1))

    final_adj <- as.matrix(fit$refit$stars)
    final_adj[is.na(final_adj)] <- 0

    ## --- subsample indices (exclude the full model) ----------------------
    ind_rng <- pulsar_subsamples(n, cfg$rep.num, stars_seed)
    ind_cap <- capture_pulsar_subsamples(
      X, fargs = list(method = "mb", lambda = grid),
      rep.num = cfg$rep.num, seed = stars_seed
    )
    ind <- ind_rng
    if (!is.null(ind_cap)) {
      if (!identical(unname(ind_cap), unname(ind_rng))) {
        warning(sprintf("%s: captured pulsar subsamples differ from RNG reproduction", case))
        ind <- ind_cap
      }
    }

    ## --- write everything -------------------------------------------------
    write_int_matrix(X, file.path(outdir, paste0(case, ".counts.tsv")))
    write_num_matrix(X_clr, file.path(outdir, paste0(case, ".clr.tsv")))
    writeLines(sprintf("%.17g", lambda_max), file.path(outdir, paste0(case, ".lambda_max.txt")))
    write_num_vec(grid, file.path(outdir, paste0(case, ".lambda_path.txt")))
    writeLines(sprintf("%.17g", fixed_lambda), file.path(outdir, paste0(case, ".fixed_lambda.txt")))
    write_int_matrix(ind - 1L, file.path(outdir, paste0(case, ".subsample_indices.tsv")))
    write_num_matrix(coef_fixed, file.path(outdir, paste0(case, ".coef_fixed.tsv")))
    write_num_matrix(adj_fixed, file.path(outdir, paste0(case, ".adj_fixed.tsv")))
    write_num_vec(raw_instability, file.path(outdir, paste0(case, ".stars_raw.txt")))
    write_num_vec(monotone, file.path(outdir, paste0(case, ".stars_monotone.txt")))
    writeLines(as.character(opt_index), file.path(outdir, paste0(case, ".opt_index.txt")))
    writeLines(sprintf("%.17g", opt_lambda), file.path(outdir, paste0(case, ".opt_lambda.txt")))
    write_int_matrix(final_adj, file.path(outdir, paste0(case, ".final_adj.tsv")))
    write_params(list(
      dataset = ds_name, n = n, p = p, want_zeros = ds$want_zeros,
      depth = ds$depth, data_seed = ds$seed,
      nlambda = cfg$nlambda, lambda_min_ratio = cfg$ratio,
      rep.num = cfg$rep.num, threshold = threshold, stars_seed = stars_seed,
      fixed_index = fixed_index
    ), file.path(outdir, paste0(case, ".params.txt")))

    manifest[[case]] <- c(ds_name = ds_name, cfg = cfg_name, n = n, p = p,
                          nlambda = cfg$nlambda, ratio = cfg$ratio,
                          opt_index = opt_index)
  }
}

## ---------------------------------------------------------------------------
## Selection-rule boundary cases (pulsar:::.starsind) — synthetic vectors.
## Each file: line 1 = thresh, line 2 = R opt index, rest = monotone summary.
## ---------------------------------------------------------------------------
starsind <- pulsar:::.starsind
sel_cases <- list(
  normal_cross   = list(s = c(0.01, 0.02, 0.03, 0.06, 0.10), thresh = 0.05),
  equality_mid   = list(s = c(0.01, 0.02, 0.05, 0.07, 0.10), thresh = 0.05),
  equality_first = list(s = c(0.05, 0.06, 0.08, 0.10),        thresh = 0.05),
  all_below      = list(s = c(0.01, 0.02, 0.03, 0.04),        thresh = 0.05),
  all_above      = list(s = c(0.06, 0.07, 0.08, 0.09),        thresh = 0.05),
  exactly_zero   = list(s = c(0.00, 0.00, 0.05, 0.06),        thresh = 0.05)
)
for (nm in names(sel_cases)) {
  sc <- sel_cases[[nm]]
  idx <- suppressWarnings(starsind(sc$s, sc$thresh, 1))
  writeLines(c(
    sprintf("%.17g", sc$thresh),
    as.character(as.integer(idx)),
    sprintf("%.17g", sc$s)
  ), file.path(seldir, paste0(nm, ".txt")))
}

## ---------------------------------------------------------------------------
## Solver-isolation goldens (ds3, ds4): one warm-started `huge` solve per pulsar
## subsample over the FULL lambda path (exactly what pulsar's `fun` does).
## Saved so the Python warm-started path solver can be compared BEFORE StARS,
## isolating "solver divergence" from "StARS divergence".
## ---------------------------------------------------------------------------
for (case in c("ds3_r", "ds3_app", "ds4_r", "ds4_app")) {
  X <- as.matrix(read.table(file.path(outdir, paste0(case, ".counts.tsv")), sep = "\t"))
  X_clr <- SpiecEasi:::.spiec.easi.norm(X)
  grid <- as.numeric(readLines(file.path(outdir, paste0(case, ".lambda_path.txt"))))
  ind0 <- as.matrix(read.table(file.path(outdir, paste0(case, ".subsample_indices.tsv")), sep = "\t"))
  nsub <- nrow(ind0); L <- length(grid); d <- ncol(X_clr)

  rows <- vector("list", nsub * L)
  k <- 0
  for (s in seq_len(nsub)) {
    Xs <- X_clr[ind0[s, ] + 1L, , drop = FALSE]
    res <- suppressWarnings(suppressMessages(
      SpiecEasi::sparseiCov(Xs, method = "mb", lambda = grid, verbose = FALSE)))
    for (li in seq_len(L)) {
      k <- k + 1L
      m <- as.matrix(res$path[[li]])       # symmetric 0/1
      rows[[k]] <- as.numeric(t(m))        # row-major flatten
    }
  }
  write.table(t(sapply(rows, as.numeric)),
              file.path(outdir, paste0(case, ".solver_adj_path.tsv")),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

  # coefficients over the full path for the first subsample (node x predictor)
  Xs <- X_clr[ind0[1, ] + 1L, , drop = FALSE]
  res <- suppressWarnings(suppressMessages(
    SpiecEasi::sparseiCov(Xs, method = "mb", lambda = grid, verbose = FALSE)))
  crow <- vector("list", L)
  for (li in seq_len(L)) {
    m <- t(as.matrix(res$beta[[li]]))      # node x predictor
    crow[[li]] <- as.numeric(t(m))         # row-major flatten
  }
  write.table(t(sapply(crow, as.numeric)),
              file.path(outdir, paste0(case, ".solver_coef_path.tsv")),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

  write_params(list(d = d, n_lambda = L, n_sub = nsub),
               file.path(outdir, paste0(case, ".solver_dims.txt")))
  message("  solver-path golden written for ", case)
}

## ---------------------------------------------------------------------------
## Environment provenance
## ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
writeLines(c(
  paste0("SpiecEasi=", as.character(packageVersion("SpiecEasi"))),
  paste0("huge=", as.character(packageVersion("huge"))),
  paste0("pulsar=", as.character(packageVersion("pulsar"))),
  paste0("glmnet=", as.character(packageVersion("glmnet"))),
  paste0("R=", as.character(getRversion()))
), file.path(outdir, "versions.txt"))

message("Golden files written to: ", normalizePath(outdir))
