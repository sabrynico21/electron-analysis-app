#!/usr/bin/env Rscript
# =============================================================================
# generate_r_golden.R - R reference generator for the real-data verification
# battery (docs/VERIFICA_DATI_REALI.md).
#
# For one (case, seed) it produces the SAME intermediate artifacts as
# tests/r_reference/generate_golden_extended.R (so the Python comparator can
# reuse the same 9 steps), under:
#     tests/real_data/golden/<case>/seed_<seed>/
#
# The case table below is the authoritative definition of the battery cases.
#
# Usage:
#   Rscript tests/real_data/generate_r_golden.R --list
#   Rscript tests/real_data/generate_r_golden.R --case R1 --seed 1001
#
# Config (default everywhere): SpiecEasi 1.99.0 / huge 2.0.1 / pulsar 0.3.13,
#   nlambda=20, lambda.min.ratio=1e-3, rep.num=20, thresh=0.05, sym rule "or".
# =============================================================================

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(pulsar)
})

THIS <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
HERE <- dirname(normalizePath(THIS))
REPO <- normalizePath(file.path(HERE, "..", "..", ".."))
DATA_DIR <- file.path(HERE, "data")
genera_path <- file.path(DATA_DIR, "genere_table.txt")
meta_path <- file.path(DATA_DIR, "sample_metadatab.txt")

NSEEDS_10 <- 10L
NSEEDS_5 <- 5L

## --- case table (authoritative) --------------------------------------------
# kind: real | limit | synthetic | pathological
mk_real <- function(id, group, thr, nseeds = NSEEDS_10)
  list(id = id, kind = "real", group = group, thr = thr, nseeds = nseeds)
mk_limit <- function(id, mode, thr = NA_real_, nsamp = NA_integer_, ngenus = NA_integer_,
                     nseeds = NSEEDS_10)
  list(id = id, kind = "limit", mode = mode, thr = thr, nsamp = nsamp,
       ngenus = ngenus, nseeds = nseeds)
mk_syn <- function(id, p, n, signal, depth = 5000L, zero_frac = 0, nseeds = NSEEDS_10)
  list(id = id, kind = "synthetic", p = p, n = n, signal = signal,
       depth = depth, zero_frac = zero_frac, nseeds = nseeds)
mk_path <- function(id, variant, nseeds = NSEEDS_10)
  list(id = id, kind = "pathological", variant = variant, nseeds = nseeds)

CASES <- list(
  mk_real("R1", "all", 0.80),
  mk_real("R2", "all", 0.90),
  mk_real("R3", "all", 0.70),
  mk_real("R4", "Treated", 0.80),
  mk_real("R5", "Not_Treated", 0.80),
  mk_real("R6", "all", 0.20, NSEEDS_5),
  mk_limit("L1", "random_samples", thr = 0.20, nsamp = 24L),
  mk_limit("L2", "random_genera", ngenus = 100L),
  mk_real("L3", "all", 0.05, NSEEDS_5),
  mk_limit("L4", "top_abundance", ngenus = 100L),
  mk_syn("S1", 50L, 100L, 0.15),
  mk_syn("S2", 50L, 100L, 0.30),
  mk_syn("S3", 50L, 100L, 0.50),
  mk_syn("S4", 60L, 200L, 0.30),
  mk_path("P1", 1L),
  mk_path("P2", 2L),
  mk_path("P3", 3L)
)
names(CASES) <- vapply(CASES, function(c) c$id, character(1))

args <- commandArgs(trailingOnly = TRUE)
if ("--list" %in% args) {
  for (c in CASES) cat(sprintf("%s %d\n", c$id, c$nseeds))
  quit(status = 0)
}
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}
case_id <- getarg("--case")
seed <- as.integer(getarg("--seed", "1001"))
if (is.null(case_id) || !case_id %in% names(CASES)) stop("unknown or missing --case")
cs <- CASES[[case_id]]

outdir <- file.path(HERE, "golden", case_id, sprintf("seed_%d", seed))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## --- writers ----------------------------------------------------------------
write_num_matrix <- function(m, path) {
  m <- as.matrix(m)
  txt <- matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m))
  write.table(txt, path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
write_num_vec <- function(v, path) writeLines(sprintf("%.17g", as.numeric(v)), path)
write_int_matrix <- function(m, path)
  write.table(m, path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

## --- data builders ----------------------------------------------------------
load_genera <- function() as.matrix(read.delim(genera_path, row.names = 1, check.names = FALSE))
load_group_samples <- function(group) {
  meta <- read.delim(meta_path, check.names = FALSE)
  ids <- meta[[1]][meta[["Condition"]] == group]
  ids
}

build_real <- function(cs, seed) {
  counts <- load_genera()  # genera x samples
  samples <- if (cs$group == "all") colnames(counts) else {
    ids <- load_group_samples(cs$group); intersect(colnames(counts), ids)
  }
  sub <- counts[, samples, drop = FALSE]
  prev <- rowMeans(sub > 0)
  keep <- prev >= cs$thr
  t(sub[keep, , drop = FALSE])  # samples x genera
}

build_limit <- function(cs, seed) {
  if (cs$mode == "random_samples") {
    counts <- load_genera()
    set.seed(seed)
    idx <- sort(sample(seq_len(ncol(counts)), cs$nsamp, replace = FALSE))
    sub <- counts[, idx, drop = FALSE]
    prev <- rowMeans(sub > 0)
    keep <- prev >= cs$thr
    return(t(sub[keep, , drop = FALSE]))
  }
  counts <- load_genera()
  if (cs$mode == "random_genera") {
    set.seed(seed)
    g <- sort(sample(seq_len(nrow(counts)), cs$ngenus, replace = FALSE))
    return(t(counts[g, , drop = FALSE]))
  }
  if (cs$mode == "top_abundance") {
    g <- order(rowSums(counts), decreasing = TRUE)[seq_len(cs$ngenus)]
    return(t(counts[g, , drop = FALSE]))
  }
  stop("unknown limit mode")
}

build_synthetic <- function(cs, seed) {
  set.seed(seed)
  p <- cs$p; n <- cs$n; s <- cs$signal
  Theta <- diag(p)
  for (i in seq_len(p - 1L)) { Theta[i, i + 1L] <- -s; Theta[i + 1L, i] <- -s }
  Sigma <- solve(Theta)
  ev <- eigen(Sigma, symmetric = TRUE)
  if (min(ev$values) <= 0) { Sigma <- Sigma + (abs(min(ev$values)) + 0.01) * diag(p) }
  L <- t(chol(Sigma))
  Z <- matrix(rnorm(n * p), n, p) %*% t(L)
  mu <- exp(Z)
  mu <- mu / rowSums(mu)
  counts <- matrix(0L, n, p)
  for (i in seq_len(n)) counts[i, ] <- as.integer(rmultinom(1L, cs$depth, mu[i, ]))
  if (cs$zero_frac > 0) {
    mask <- matrix(runif(n * p) < cs$zero_frac, n, p)
    counts[mask] <- 0L
  }
  for (j in seq_len(p)) if (length(unique(counts[, j])) == 1L) counts[1L, j] <- counts[1L, j] + 1L
  for (i in seq_len(n)) if (sum(counts[i, ]) == 0L) counts[i, 1L] <- counts[i, 1L] + 1L
  counts
}

build_pathological <- function(cs, seed) {
  set.seed(seed)
  p <- 40L; n <- 60L; depth <- 3000L
  mu <- matrix(exp(rnorm(n * p, sd = 0.7)), n, p)
  mu <- mu / rowSums(mu)
  counts <- matrix(0L, n, p)
  for (i in seq_len(n)) counts[i, ] <- as.integer(rmultinom(1L, depth, mu[i, ]))
  if (cs$variant == 1L) {
    counts[, p] <- counts[, 1L] * 2L                     # proportional duplicate
  } else if (cs$variant == 2L) {
    counts[, p] <- 0L; counts[1L, p] <- depth            # present in one sample only
  } else if (cs$variant == 3L) {
    counts[1L, ] <- 0L; counts[1L, 1L] <- 3L             # one sample with very few reads
  }
  for (j in seq_len(p)) if (length(unique(counts[, j])) == 1L) counts[1L, j] <- counts[1L, j] + 1L
  for (i in seq_len(n)) if (sum(counts[i, ]) == 0L) counts[i, 1L] <- counts[i, 1L] + 1L
  counts
}

## --- main -------------------------------------------------------------------
build_Y <- function(cs, seed) switch(
  cs$kind,
  real = , limit = if (cs$kind == "real") build_real(cs, seed) else build_limit(cs, seed),
  synthetic = build_synthetic(cs, seed),
  pathological = build_pathological(cs, seed)
)

errfile <- file.path(outdir, "error.txt")
if (file.exists(errfile)) unlink(errfile)

res <- tryCatch({
  Y <- build_Y(cs, seed)
  Y <- matrix(as.numeric(Y), nrow = nrow(Y))
  storage.mode(Y) <- "integer"
  p <- ncol(Y); n <- nrow(Y)
  # Guard: huge rejects constant columns.
  const <- which(apply(Y, 2, function(col) length(unique(col)) == 1L))
  if (length(const)) {
    writeLines(sprintf("constant_column:%s", paste(const, collapse = ",")), errfile)
    quit(status = 0)
  }
  X_clr <- SpiecEasi:::.spiec.easi.norm(Y)
  lambda_max <- pulsar::getMaxCov(cor(X_clr))
  grid <- pulsar::getLamPath(lambda_max, lambda_max * 1e-3, 20L, log = TRUE)
  ratio <- pulsar:::.ratcheck(NULL, n)
  set.seed(seed)
  ind <- do.call(rbind, replicate(20L, sample(seq_len(n), floor(n * ratio), replace = FALSE),
                                  simplify = FALSE))
  fit <- suppressWarnings(suppressMessages(spiec.easi(
    Y, method = "mb", sel.criterion = "stars",
    nlambda = 20L, lambda.min.ratio = 1e-3,
    pulsar.params = list(rep.num = 20L, thresh = 0.05, seed = seed), verbose = FALSE)))
  stars <- fit$select$stars
  grid_fit <- as.numeric(fit$lambda)
  opt_index <- as.integer(stars$opt.index)
  monotone <- as.numeric(stars$summary)
  raw <- vapply(stars$merge, function(P) {
    P <- as.matrix(P); 4 * sum(P * (1 - P)) / (p * (p - 1))
  }, numeric(1))
  final_adj <- as.matrix(fit$refit$stars); final_adj[is.na(final_adj)] <- 0

  L <- length(grid_fit); nsub <- nrow(ind)
  rows <- vector("list", nsub * L); k <- 0L
  for (s in seq_len(nsub)) {
    Xs <- X_clr[ind[s, ], , drop = FALSE]
    r <- suppressWarnings(suppressMessages(
      SpiecEasi::sparseiCov(Xs, method = "mb", lambda = grid_fit, verbose = FALSE)))
    for (li in seq_len(L)) { k <- k + 1L; rows[[k]] <- as.numeric(t(as.matrix(r$path[[li]]))) }
  }
  r1 <- suppressWarnings(suppressMessages(
    SpiecEasi::sparseiCov(X_clr[ind[1, ], , drop = FALSE], method = "mb",
                          lambda = grid_fit, verbose = FALSE)))
  crow <- vector("list", L)
  for (li in seq_len(L)) { m <- t(as.matrix(r1$beta[[li]])); crow[[li]] <- as.numeric(t(m)) }

  write_int_matrix(Y, file.path(outdir, "counts.tsv"))
  write_num_matrix(X_clr, file.path(outdir, "clr.tsv"))
  writeLines(sprintf("%.17g", lambda_max), file.path(outdir, "lambda_max.txt"))
  write_num_vec(grid_fit, file.path(outdir, "lambda_path.txt"))
  write_int_matrix(ind - 1L, file.path(outdir, "subsample_indices.tsv"))
  write_num_vec(raw, file.path(outdir, "stars_raw.txt"))
  write_num_vec(monotone, file.path(outdir, "stars_monotone.txt"))
  writeLines(as.character(opt_index), file.path(outdir, "opt_index.txt"))
  writeLines(sprintf("%.17g", grid_fit[[opt_index]]), file.path(outdir, "opt_lambda.txt"))
  write_int_matrix(final_adj, file.path(outdir, "final_adj.tsv"))
  write.table(t(sapply(rows, as.numeric)), file.path(outdir, "solver_adj_path.tsv"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(t(sapply(crow, as.numeric)), file.path(outdir, "solver_coef_path.tsv"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  writeLines(c(paste0("d=", p), paste0("n_lambda=", L), paste0("n_sub=", nsub),
               paste0("p=", p), paste0("n=", n)),
             file.path(outdir, "solver_dims.txt"))
  writeLines("ok", file.path(outdir, "done.txt"))
  cat(sprintf("RGEN case=%s seed=%d p=%d n=%d opt_index=%d ok\n", case_id, seed, p, n, opt_index))
}, error = function(e) {
  writeLines(sprintf("r_error:%s", conditionMessage(e)), errfile)
  cat(sprintf("RGEN case=%s seed=%d ERROR %s\n", case_id, seed, conditionMessage(e)))
})
