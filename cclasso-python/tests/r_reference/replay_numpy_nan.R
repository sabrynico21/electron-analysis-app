#!/usr/bin/env Rscript
# Replay the exact NumPy input, bootstrap indices, and null draws through the
# pinned CCLasso v2.0 R source. This answers whether the known non-finite result
# is shared by R or is specific to the Python port.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript replay_numpy_nan.R <diagnostic_dir> <output_dir>")
}
in_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("tests", "r_reference", "reference", "R", "cclasso.R"))
# Keep this harness base-R-only like the pinned reference implementation. These
# values match make_numpy_nan_replay.py; params.json is retained as provenance.
params <- list(counts = FALSE, pseudo = 0.5, k_cv = 3L, n_boot = 6L)
x <- as.matrix(read.csv(file.path(in_dir, "input.csv"), header = FALSE))
boot <- as.matrix(read.csv(file.path(in_dir, "boot_indices.csv"), header = FALSE))
null_norm <- as.matrix(read.csv(file.path(in_dir, "null_norm.csv"), header = FALSE))
stopifnot(identical(dim(boot), c(nrow(x), params$n_boot)))
stopifnot(identical(dim(null_norm), dim(x)))

# Supply the same arrays at the only two stochastic calls in boot_cclasso().
# The returned vectors are consumed by R's matrix() in column-major order.
base_sample <- base::sample
base_rnorm <- stats::rnorm
sample <- function(x, size, replace = FALSE, prob = NULL) {
  if (size == length(boot) && isTRUE(replace)) return(as.vector(boot))
  base_sample(x, size = size, replace = replace, prob = prob)
}
rnorm <- function(n, mean = 0, sd = 1) {
  if (n == length(null_norm) && mean == 0 && sd == 1) return(as.vector(null_norm))
  base_rnorm(n, mean = mean, sd = sd)
}

trace("boot_cclasso", exit = quote({
  REPLAY$vars_boot <- vars_boot
}), print = FALSE)
REPLAY <- new.env(parent = emptyenv())
result <- cclasso(
  x = x,
  counts = params$counts,
  pseudo = params$pseudo,
  k_cv = params$k_cv,
  n_boot = params$n_boot
)
untrace("boot_cclasso")
sample <- base_sample
rnorm <- base_rnorm

write.csv(REPLAY$vars_boot, file.path(out_dir, "r_bootstrap_variances.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(result$var_w, file.path(out_dir, "r_var_w.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(result$cor_w, file.path(out_dir, "r_cor_w.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(result$p_vals, file.path(out_dir, "r_p_vals.csv"),
          row.names = FALSE, quote = FALSE)
status <- data.frame(
  output = c("var_w", "cor_w", "p_vals"),
  finite = c(all(is.finite(result$var_w)), all(is.finite(result$cor_w)),
             all(is.finite(result$p_vals))),
  nonfinite_count = c(sum(!is.finite(result$var_w)),
                       sum(!is.finite(result$cor_w)),
                       sum(!is.finite(result$p_vals))),
  na_count = c(sum(is.na(result$var_w)), sum(is.na(result$cor_w)),
                sum(is.na(result$p_vals))),
  nan_count = c(sum(is.nan(result$var_w)), sum(is.nan(result$cor_w)),
                sum(is.nan(result$p_vals)))
)
write.csv(status, file.path(out_dir, "r_status.csv"), row.names = FALSE, quote = FALSE)
writeLines(c(
  paste("negative bootstrap variance count:", sum(REPLAY$vars_boot < 0)),
  paste("non-finite var_w:", sum(!is.finite(result$var_w))),
  paste("non-finite cor_w:", sum(!is.finite(result$cor_w))),
  paste("non-finite p_vals:", sum(!is.finite(result$p_vals))),
  paste("R version:", R.version.string)
), file.path(out_dir, "r_diagnostic.txt"))
print(status)
cat(readLines(file.path(out_dir, "r_diagnostic.txt")), sep = "\n")
