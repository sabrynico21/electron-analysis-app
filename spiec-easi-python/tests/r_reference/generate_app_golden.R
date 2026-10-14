#!/usr/bin/env Rscript
# =============================================================================
# generate_app_golden.R — Phase D app-integration golden generator.
#
# Produces reference outputs for the *application* entry point
#   python/analysis.py :: _compute_spieceasi_for_matrix
# using the same reference stack as the other generators
#   (SpiecEasi 1.99.0 / huge 2.0.1 / pulsar 0.3.13, 17-digit output).
#
# Two small cases:
#   app_plain_p12_n40_s301  — 12 varying features
#   app_const_p13_n40_s301  — same data + 1 constant feature (the app drops it,
#                             R would error, so the golden is computed on the
#                             reduced matrix the app keeps)
#
# Output: tests/golden_app/<case_id>/*  (+ manifest.json)
#
# Usage: Rscript tests/r_reference/generate_app_golden.R [outdir]
# =============================================================================

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(pulsar)
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1] else "tests/golden_app"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_num_matrix <- function(m, path) {
  m <- as.matrix(m)
  txt <- matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m))
  write.table(txt, path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
write_num_vec <- function(v, path) writeLines(sprintf("%.17g", as.numeric(v)), path)
write_int_matrix <- function(m, path) {
  write.table(m, path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

## ---------------------------------------------------------------------------
## Deterministic synthetic counts (same construction as generate_golden_extended.R)
## ---------------------------------------------------------------------------
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

## ---------------------------------------------------------------------------
## Case parameters (match the application defaults)
## ---------------------------------------------------------------------------
n_lambda   <- 20L
ratio      <- 1e-3
rep_num    <- 20L
threshold  <- 0.05
stars_seed <- 4242L
n          <- 40L
p          <- 12L
depth      <- 5000L
seed       <- 301L
zero_frac  <- 0.15

X <- make_counts(n, p, depth, seed, zero_frac)
rownames(X) <- sprintf("s%02d", seq_len(n))
colnames(X) <- sprintf("f%02d", seq_len(p))
X_clr <- SpiecEasi:::.spiec.easi.norm(X)

## constant feature column appended for the const case
const_feature <- "fconst"
Xfull <- cbind(X, fconst = rep(100L, n))
Xfull_clr <- cbind(X_clr, fconst = rep(0, n))  # not used by R, informational

lambda_max <- pulsar::getMaxCov(cor(X_clr))
grid <- pulsar::getLamPath(lambda_max, lambda_max * ratio, n_lambda, log = TRUE)

ratio_sub <- pulsar:::.ratcheck(NULL, nn <- nrow(X))
set.seed(stars_seed)
ind <- do.call(rbind, replicate(rep_num, sample(c(1L:nn), floor(nn * ratio_sub), replace = FALSE),
                                simplify = FALSE))

fit <- suppressWarnings(suppressMessages(spiec.easi(
  X, method = "mb", sel.criterion = "stars",
  nlambda = n_lambda, lambda.min.ratio = ratio,
  pulsar.params = list(rep.num = rep_num, thresh = threshold, seed = stars_seed),
  verbose = FALSE)))
stars <- fit$select$stars
grid_fit <- as.numeric(fit$lambda)
opt_index <- as.integer(stars$opt.index)
monotone <- as.numeric(stars$summary)
raw <- vapply(stars$merge, function(P) {
  P <- as.matrix(P); 4 * sum(P * (1 - P)) / (p * (p - 1))
}, numeric(1))
final_adj <- as.matrix(fit$refit$stars); final_adj[is.na(final_adj)] <- 0
dimnames(final_adj) <- list(colnames(X), colnames(X))

emit <- function(case_id, counts_full, kept, note) {
  cdir <- file.path(outdir, case_id)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  write.table(counts_full, file.path(cdir, "counts_full.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  writeLines(kept, file.path(cdir, "expected_kept_features.txt"))
  write_num_matrix(X_clr, file.path(cdir, "clr.tsv"))
  writeLines(sprintf("%.17g", lambda_max), file.path(cdir, "lambda_max.txt"))
  write_num_vec(grid_fit, file.path(cdir, "lambda_path.txt"))
  write_int_matrix(ind - 1L, file.path(cdir, "subsample_indices.tsv"))
  write_num_vec(raw, file.path(cdir, "stars_raw.txt"))
  write_num_vec(monotone, file.path(cdir, "stars_monotone.txt"))
  writeLines(as.character(opt_index), file.path(cdir, "opt_index.txt"))
  writeLines(sprintf("%.17g", grid_fit[[opt_index]]), file.path(cdir, "opt_lambda.txt"))
  write.table(final_adj, file.path(cdir, "final_adj.tsv"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  writeLines(c(paste0("d=", p), paste0("n_lambda=", length(grid_fit)),
               paste0("n_sub=", nrow(ind)), paste0("p=", p), paste0("n=", n),
               paste0("nlambda=", n_lambda), paste0("lambda_min_ratio=", format(ratio)),
               paste0("rep_num=", rep_num), paste0("threshold=", threshold)),
             file.path(cdir, "meta.txt"))
  list(id = case_id, p = p, n = n, seed = seed, nlambda = n_lambda,
       lambda_min_ratio = ratio, rep_num = rep_num, threshold = threshold,
       n_kept_features = length(kept), n_sub = nrow(ind), note = note)
}

manifest <- list()
manifest[["app_plain_p12_n40_s301"]] <- emit(
  "app_plain_p12_n40_s301", X, colnames(X), "no constant feature")
manifest[["app_const_p13_n40_s301"]] <- emit(
  "app_const_p13_n40_s301", Xfull, colnames(X),
  "one constant feature (fconst) appended; app must drop it and match the plain golden")

json <- jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE)
writeLines(json, file.path(outdir, "manifest.json"))
message("Wrote ", length(manifest), " app case(s) to ", outdir)
