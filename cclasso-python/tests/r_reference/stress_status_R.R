#!/usr/bin/env Rscript
# Batch R-side convergence/failure classification for the CCLasso stress runner.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript stress_status_R.R <cases_dir> <status.csv>")
}
cases_dir <- args[[1]]
status_path <- args[[2]]
source(file.path("tests", "r_reference", "reference", "R", "cclasso.R"))

STRESS <- new.env(parent = emptyenv())
base_sample <- base::sample
base_rnorm <- stats::rnorm
sample <- function(x, size, replace = FALSE, prob = NULL) {
  if (exists("boot", STRESS, inherits = FALSE) &&
      size == length(STRESS$boot) && isTRUE(replace)) {
    return(as.vector(STRESS$boot))
  }
  base_sample(x, size = size, replace = replace, prob = prob)
}
rnorm <- function(n, mean = 0, sd = 1) {
  if (exists("null_norm", STRESS, inherits = FALSE) &&
      n == length(STRESS$null_norm) && mean == 0 && sd == 1) {
    return(as.vector(STRESS$null_norm))
  }
  base_rnorm(n, mean = mean, sd = sd)
}

case_dirs <- sort(list.dirs(cases_dir, recursive = FALSE, full.names = TRUE))
rows <- vector("list", length(case_dirs))
for (i in seq_along(case_dirs)) {
  dir <- case_dirs[[i]]
  params <- jsonlite::fromJSON(file.path(dir, "params.json"))
  x <- as.matrix(read.csv(file.path(dir, "input.csv"), header = FALSE))
  STRESS$boot <- as.matrix(read.csv(file.path(dir, "boot_indices.csv"), header = FALSE))
  STRESS$null_norm <- as.matrix(read.csv(file.path(dir, "null_norm.csv"), header = FALSE))
  result <- tryCatch(
    cclasso(
      x,
      counts = params$counts,
      pseudo = params$pseudo,
      k_cv = params$k_cv,
      lam_int = as.numeric(params$lam_int),
      k_max = params$k_max,
      n_boot = params$n_boot
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    rows[[i]] <- data.frame(
      case = basename(dir), status = "error", error_class = class(result)[[1]],
      nonfinite_var = NA_integer_, nonfinite_cor = NA_integer_,
      nonfinite_p = NA_integer_, stringsAsFactors = FALSE
    )
  } else {
    finite_counts <- c(
      sum(!is.finite(result$var_w)),
      sum(!is.finite(result$cor_w)),
      sum(!is.finite(result$p_vals))
    )
    rows[[i]] <- data.frame(
      case = basename(dir),
      status = if (any(finite_counts > 0)) "nonfinite" else "finite",
      error_class = "",
      nonfinite_var = finite_counts[[1]],
      nonfinite_cor = finite_counts[[2]],
      nonfinite_p = finite_counts[[3]],
      stringsAsFactors = FALSE
    )
  }
}
utils::write.csv(do.call(rbind, rows), status_path, row.names = FALSE, na = "")
