#!/usr/bin/env Rscript
# =============================================================================
# point_f.R - Point F of docs/VERIFICA_DATI_REALI.md (p=1135).
#
# One subsample of the FULL dataset (all 1135 genera, all 96 samples), the REAL
# lambda path (lambda.max + pulsar::getLamPath), and R's own huge::huge(mb) run
# both with OMP_NUM_THREADS=1 and with the default thread count. Writes S, the
# lambda path and the per-lambda adjacency (0/1, "or"-symmetrised) for the Python
# side to reproduce. NO timing comparison.
# =============================================================================

suppressPackageStartupMessages({ library(SpiecEasi); library(pulsar); library(huge) })

THIS <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
HERE <- dirname(normalizePath(THIS))
REPO <- normalizePath(file.path(HERE, "..", "..", ".."))
OUT <- file.path(HERE, "results", "point_f")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

seed <- 1001L
counts <- as.matrix(read.delim(file.path(HERE, "data", "genere_table.txt"),
                               row.names = 1, check.names = FALSE))
Y <- t(counts)                        # samples x genera
X_clr <- SpiecEasi:::.spiec.easi.norm(Y)
lambda_max <- pulsar::getMaxCov(cor(X_clr))
grid <- pulsar::getLamPath(lambda_max, lambda_max * 1e-3, 20L, log = TRUE)
ratio <- pulsar:::.ratcheck(NULL, nrow(Y))
set.seed(seed)
idx <- sort(sample(seq_len(nrow(Y)), floor(nrow(Y) * ratio), replace = FALSE))
Xs <- X_clr[idx, , drop = FALSE]
S <- cor(Xs)
p <- ncol(Xs)

write_num_matrix <- function(m, path) {
  m <- as.matrix(m)
  write.table(matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m)),
              path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
writeLines(sprintf("%.17g", lambda_max), file.path(OUT, "lambda_max.txt"))
writeLines(sprintf("%.17g", grid), file.path(OUT, "lambda_path.txt"))
write_num_matrix(S, file.path(OUT, "S.tsv"))
write.table(idx, file.path(OUT, "subsample_indices.tsv"), sep = "\t",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

adj_of <- function(res) {
  L <- length(res$path)
  rows <- vector("list", L)
  for (li in seq_len(L)) {
    A <- as.matrix(res$path[[li]] != 0)
    A <- A | t(A)                     # "or" symmetrisation
    rows[[li]] <- as.numeric(t(A))
  }
  matrix(unlist(rows), nrow = L, byrow = TRUE)
}

run_mb <- function(threads) {
  if (is.na(threads)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = threads)
  suppressWarnings(suppressMessages(huge(Xs, method = "mb", lambda = grid, verbose = FALSE)))
}

res1 <- run_mb("1")
res8 <- run_mb(NA)
adj1 <- adj_of(res1); adj8 <- adj_of(res8)
thread_same <- identical(adj1, adj8)

write.table(adj1, file.path(OUT, "spmb_adj_r.tsv"), sep = "\t",
            quote = FALSE, col.names = FALSE, row.names = FALSE)
edges <- rowSums(adj1) / 2
writeLines(sprintf("%d", as.integer(edges)), file.path(OUT, "spmb_edges_r.txt"))
cat(sprintf("POINT_F p=%d n_subsample=%d lambda_max=%.10g grid[1]=%.10g grid[20]=%.10g\n",
            p, nrow(Xs), lambda_max, grid[1], grid[20]))
cat(sprintf("        threads_1_vs_default_identical=%s  edges_per_lambda=%s\n",
            thread_same, paste(as.integer(edges), collapse = ",")))
