#!/usr/bin/env Rscript
# =============================================================================
# validate_solver_helper.R - recompute R's correlation matrix S and R's MB
# coefficients for selected subsamples of a case (used by the L2/L3 validity
# check; does NOT touch the golden files).
#
# Usage:
#   Rscript validate_solver_helper.R --case L2 --seed 1001 --subs 0,1 --outdir DIR
#
# Writes, per subsample s:
#   s_<s>.tsv      R's correlation matrix S = cor(X_clr[ind[s],])
#   beta_<s>.tsv   R's huge MB coefficients over the path (n_lambda, d, d) row-major
# =============================================================================

suppressPackageStartupMessages({ library(SpiecEasi); library(pulsar) })

THIS <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
HERE <- dirname(normalizePath(THIS))
DATA_DIR <- file.path(HERE, "data")
genera_path <- file.path(DATA_DIR, "genere_table.txt")

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[[i + 1L]]
}
case_id <- getarg("--case")
seed <- as.integer(getarg("--seed"))
subs <- as.integer(strsplit(getarg("--subs", "0"), ",")[[1]])
outdir <- getarg("--outdir")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_num_matrix <- function(m, path) {
  m <- as.matrix(m)
  write.table(matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m)),
              path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

counts <- as.matrix(read.delim(genera_path, row.names = 1, check.names = FALSE))
if (case_id == "L2") {
  set.seed(seed)
  g <- sort(sample(seq_len(nrow(counts)), 100L, replace = FALSE))
  Y <- t(counts[g, , drop = FALSE])
} else if (case_id == "L3") {
  prev <- rowMeans(counts > 0); Y <- t(counts[prev >= 0.05, , drop = FALSE])
} else stop("unsupported case")

p <- ncol(Y); n <- nrow(Y)
X_clr <- SpiecEasi:::.spiec.easi.norm(Y)
lambda_max <- pulsar::getMaxCov(cor(X_clr))
grid <- pulsar::getLamPath(lambda_max, lambda_max * 1e-3, 20L, log = TRUE)
ratio <- pulsar:::.ratcheck(NULL, n)
set.seed(seed)
ind <- do.call(rbind, replicate(20L, sample(seq_len(n), floor(n * ratio), replace = FALSE),
                                simplify = FALSE))

for (s in subs) {
  Xs <- X_clr[ind[s + 1L, ], , drop = FALSE]
  S <- cor(Xs)
  write_num_matrix(S, file.path(outdir, sprintf("s_%d.tsv", s)))
  res <- suppressWarnings(suppressMessages(
    SpiecEasi::sparseiCov(Xs, method = "mb", lambda = grid, verbose = FALSE)))
  L <- length(grid); crow <- vector("list", L)
  for (li in seq_len(L)) { m <- t(as.matrix(res$beta[[li]])); crow[[li]] <- as.numeric(t(m)) }
  write.table(t(sapply(crow, as.numeric)), file.path(outdir, sprintf("beta_%d.tsv", s)),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
cat(sprintf("VHELPER case=%s seed=%d subs=%s p=%d n=%d ok\n", case_id, seed,
            paste(subs, collapse = ","), p, n))
