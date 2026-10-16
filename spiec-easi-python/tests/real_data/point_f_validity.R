#!/usr/bin/env Rscript
# =============================================================================
# point_f_validity.R - p=1135: R's S, the real lambda path, and R's huge MB
# coefficients at lambda indices 1..5 (0-based), for the KKT/objective check.
#
# Writes results/point_f_validity/: S.tsv, lambda_path.txt, beta_idx1_5.tsv
# (5 x d x d, row-major, node x predictor). No timing comparison.
# =============================================================================

suppressPackageStartupMessages({ library(SpiecEasi); library(pulsar); library(huge) })

THIS <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
HERE <- dirname(normalizePath(THIS))
OUT <- file.path(HERE, "results", "point_f_validity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(SpiecEasi))
counts <- as.matrix(read.delim(file.path(HERE, "data", "genere_table.txt"),
                               row.names = 1, check.names = FALSE))
Y <- t(counts)
X_clr <- SpiecEasi:::.spiec.easi.norm(Y)
lambda_max <- pulsar::getMaxCov(cor(X_clr))
grid <- pulsar::getLamPath(lambda_max, lambda_max * 1e-3, 20L, log = TRUE)
ratio <- pulsar:::.ratcheck(NULL, nrow(Y))
set.seed(1001L)
idx <- sort(sample(seq_len(nrow(Y)), floor(nrow(Y) * ratio), replace = FALSE))
Xs <- X_clr[idx, , drop = FALSE]
S <- cor(Xs)
d <- ncol(Xs)

write_num_matrix <- function(m, path) {
  m <- as.matrix(m)
  write.table(matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m)),
              path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
write_num_matrix(S, file.path(OUT, "S.tsv"))
writeLines(sprintf("%.17g", grid), file.path(OUT, "lambda_path.txt"))

res <- suppressWarnings(suppressMessages(huge(Xs, method = "mb", lambda = grid, verbose = FALSE)))
idx01 <- 1:5  # 0-based indices 1..5  -> 1-based 2..6
rows <- vector("list", length(idx01))
for (k in seq_along(idx01)) {
  m <- t(as.matrix(res$beta[[idx01[k] + 1L]]))
  rows[[k]] <- as.numeric(t(m))
}
write.table(t(sapply(rows, as.numeric)), file.path(OUT, "beta_idx1_5.tsv"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
cat(sprintf("POINT_F_VALIDITY p=%d d=%d idx=1..5 written\n", d, d))
