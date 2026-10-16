#!/usr/bin/env Rscript
# =============================================================================
# calibrate_r_time.R - Part A micro-calibration of R wall-clock time per case.
#
# Runs the FULL default protocol on a real-data sub-case (a prevalence-threshold
# genus subset of the real 1135x96 matrix) and prints the elapsed seconds, so
# that the per-case R time of the battery can be estimated.
#
# Protocol: spiec.easi(method="mb"), nlambda=20, lambda.min.ratio=1e-3,
#           pulsar.params=list(rep.num=20, thresh=0.05, seed=<seed>).
#
# Usage:
#   Rscript spy.../calibrate_r_time.R <threshold> [group] [seed]
#     threshold : prevalence threshold in [0,1] on the group's samples
#     group     : all | Treated | Not_Treated (default all)
#     seed      : StARS seed (default 17)
# =============================================================================

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(pulsar)
})

args <- commandArgs(trailingOnly = TRUE)
threshold <- if (length(args) >= 1L) as.numeric(args[1]) else 0.95
group <- if (length(args) >= 2L) args[2] else "all"
stars_seed <- if (length(args) >= 3L) as.integer(args[3]) else 17L

repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "..", "..", ".."))
data_dir <- file.path(repo_root, "spiec-easi-python", "tests", "real_data", "data")
genera_path <- file.path(data_dir, "genere_table.txt")
meta_path <- file.path(data_dir, "sample_metadatab.txt")

counts <- as.matrix(read.delim(genera_path, row.names = 1, check.names = FALSE))
meta <- read.delim(meta_path, check.names = FALSE)
cond_col <- "Condition"
id_col <- names(meta)[1]

if (group == "all") {
  samples <- colnames(counts)
} else {
  ids <- meta[[id_col]][meta[[cond_col]] == group]
  samples <- intersect(colnames(counts), ids)
}
sub <- counts[, samples, drop = FALSE]

present <- (sub > 0)
prev <- rowMeans(present)
keep <- prev >= threshold
Y <- t(sub[keep, , drop = FALSE])  # samples x genera
p <- ncol(Y); n <- nrow(Y)

X_clr <- SpiecEasi:::.spiec.easi.norm(Y)

t0 <- Sys.time()
fit <- suppressWarnings(suppressMessages(spiec.easi(
  Y, method = "mb", sel.criterion = "stars",
  nlambda = 20L, lambda.min.ratio = 1e-3,
  pulsar.params = list(rep.num = 20L, thresh = 0.05, seed = stars_seed),
  verbose = FALSE)))
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

opt_index <- as.integer(fit$select$stars$opt.index)
cat(sprintf("CALIB threshold=%.2f group=%s seed=%d n=%d p=%d zero_frac=%.4f opt_index=%d elapsed_s=%.2f\n",
            threshold, group, stars_seed, n, p, mean(sub[keep, ] == 0), opt_index, dt))
