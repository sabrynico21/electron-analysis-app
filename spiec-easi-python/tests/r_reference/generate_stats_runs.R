#!/usr/bin/env Rscript
# =============================================================================
# generate_stats_runs.R — Phase E statistical-equivalence reference runs.
#
# NO subsample injection here: R and Python each draw their own StARS
# subsamples, so this measures *distributional* agreement (selected lambda,
# edge counts, per-edge selection frequency) rather than bit-level parity.
#
# For each dataset and each of N seeds: run the full spiec.easi MB+StARS
# pipeline with pulsar.params$seed = s and record the selected index, the
# selected lambda, the edge count and the edge list.
#
# Output: tests/golden_stats/<id>/{counts.tsv,r_runs.json}
#
# Usage: Rscript tests/r_reference/generate_stats_runs.R [outdir] [n_seeds]
# =============================================================================

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(pulsar)
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1] else "tests/golden_stats"
n_seeds <- if (length(args) >= 2L) as.integer(args[2]) else 30L
id_filter <- if (length(args) >= 3L) args[3] else ".*"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_counts <- function(n, p, depth, seed, zero_frac) {
  set.seed(seed)
  k <- max(1L, p %/% 4L)
  W <- matrix(abs(rnorm(p * k, sd = 1)), p, k)
  Z <- matrix(rnorm(n * k), n, k)
  mu <- exp(as.numeric(scale(Z %*% t(W))))
  mu <- matrix(pmax(mu, 1e-6), n, p)
  counts <- matrix(0L, n, p)
  for (i in seq_len(n)) {
    counts[i, ] <- as.integer(rmultinom(1L, depth, mu[i, ] / sum(mu[i, ])))
  }
  if (zero_frac > 0) {
    mask <- matrix(runif(n * p) < zero_frac, n, p)
    counts[mask] <- 0L
  }
  for (j in seq_len(p)) {
    if (length(unique(counts[, j])) == 1L) counts[1L, j] <- counts[1L, j] + 1L
  }
  for (i in seq_len(n)) {
    if (sum(counts[i, ]) == 0L) counts[i, 1L] <- counts[i, 1L] + 1L
  }
  counts
}

datasets <- list(
  list(id = "stats_p10_n60_s401",  p = 10L, n = 60L,  zeros = 0.15, depth = 5000L, seed = 401L),
  list(id = "stats_p10_n100_s402", p = 10L, n = 100L, zeros = 0.15, depth = 5000L, seed = 402L)
)

nlambda <- 20L
ratio <- 1e-3
rep_num <- 20L
threshold <- 0.05

for (ds in datasets) {
  if (!grepl(id_filter, ds$id)) next
  cdir <- file.path(outdir, ds$id)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  X <- make_counts(ds$n, ds$p, ds$depth, ds$seed, ds$zeros)
  rownames(X) <- sprintf("s%03d", seq_len(nrow(X)))
  colnames(X) <- sprintf("f%02d", seq_len(ncol(X)))
  write.table(X, file.path(cdir, "counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)

  runs <- vector("list", n_seeds)
  for (s in seq_len(n_seeds)) {
    fit <- suppressWarnings(suppressMessages(spiec.easi(
      X, method = "mb", sel.criterion = "stars",
      nlambda = nlambda, lambda.min.ratio = ratio,
      pulsar.params = list(rep.num = rep_num, thresh = threshold, seed = s),
      verbose = FALSE)))
    grid <- as.numeric(fit$lambda)
    opt <- as.integer(fit$select$stars$opt.index)
    monotone <- as.numeric(fit$select$stars$summary)
    raw <- vapply(fit$select$stars$merge, function(P) {
      P <- as.matrix(P); 4 * sum(P * (1 - P)) / (ds$p * (ds$p - 1))
    }, numeric(1))
    A <- as.matrix(fit$refit$stars); A[is.na(A)] <- 0
    pairs <- which(upper.tri(A) & A != 0, arr.ind = TRUE)
    edges <- if (nrow(pairs) == 0) list() else
      lapply(seq_len(nrow(pairs)), function(k)
        c(colnames(X)[pairs[k, 1]], colnames(X)[pairs[k, 2]]))
    runs[[s]] <- list(
      seed = s, opt_index = opt, opt_lambda = grid[[opt]],
      n_edges = length(edges), edges = edges,
      raw_instability = raw, monotone_instability = monotone
    )
    message("  ", ds$id, " seed=", s, " opt=", opt, " edges=", length(edges))
  }

  # per-edge selection frequency over seeds
  all_keys <- character(0)
  for (r in runs) for (e in r$edges) all_keys <- c(all_keys, paste(e, collapse = "|"))
  freq <- list()
  for (key in sort(unique(all_keys))) {
    freq[[key]] <- sum(vapply(runs, function(r)
      any(vapply(r$edges, function(e) paste(e, collapse = "|") == key, logical(1))),
      logical(1))) / n_seeds
  }

  out <- list(
    id = ds$id, p = as.integer(ds$p), n = as.integer(ds$n), seed = as.integer(ds$seed),
    nlambda = nlambda, lambda_min_ratio = ratio, rep_num = rep_num,
    threshold = threshold, n_seeds = n_seeds,
    runs = runs, edge_frequency = freq
  )
  writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = 17, pretty = TRUE),
             file.path(cdir, "r_runs.json"))
  message("wrote ", ds$id)
}
