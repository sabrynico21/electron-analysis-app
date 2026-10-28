#!/usr/bin/env Rscript
# Batch R-side status + discrete-output classification for the SPRING stress
# runner (tests/r_reference/run_status_stress.py).
#
# For every archived case it recomputes, from the injected counts and the
# injected StARS subsample indices (no R RNG is used):
#   mclr -> explicit lambda path -> latentcor + huge.mb (full data)
#        -> per-replica latentcor + huge.mb -> StARS merge/summary/opt index
#        -> refit adjacency and MB coefficients at the selected lambda
# and writes status (finite/nonfinite/error), opt_index, refit edge set, the
# non-zero coefficients at the selected lambda, and per-edge merge frequencies.
#
# Usage: Rscript stress_status_R.R <cases_dir> <status.csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript stress_status_R.R <cases_dir> <status.csv>")
}
cases_dir <- args[[1]]
status_path <- args[[2]]

suppressWarnings(suppressMessages({
  library(jsonlite)
  library(latentcor)
  library(huge)
}))

find_spring <- function() {
  cands <- c(
    file.path("..", "SPRING-master", "R"),
    file.path("..", "..", "SPRING-master", "R"),
    Sys.getenv("SPRING_R_HOME")
  )
  for (c in cands) {
    if (nzchar(c) && file.exists(file.path(c, "helpers.R"))) return(normalizePath(c))
  }
  stop("Cannot find SPRING-master/R. Set SPRING_R_HOME or run from spring-python/.")
}
spring_r_dir <- find_spring()
source(file.path(spring_r_dir, "helpers.R"))

empty_row <- function(case) {
  data.frame(
    case = case, status = "error", error_class = "unknown",
    opt_index = NA_integer_, n_edges = NA_integer_, edges = "", beta_opt = "",
    lambda_first = NA_real_, stars_margin = NA_real_, merge_opt = "",
    seconds = NA_real_, stringsAsFactors = FALSE
  )
}

classify_case <- function(dir) {
  case <- basename(dir)
  start <- proc.time()[["elapsed"]]
  params <- fromJSON(file.path(dir, "params.json"))
  counts <- as.matrix(read.csv(file.path(dir, "counts.csv"), header = FALSE))
  idx <- as.matrix(read.csv(file.path(dir, "subsample_indices.csv"), header = FALSE))
  storage.mode(counts) <- "double"
  storage.mode(idx) <- "integer"

  p <- ncol(counts)
  reps <- nrow(idx)

  qdat <- mclr(counts)
  lams <- exp(seq(log(0.6), log(0.6 * params$lambda_min_ratio),
                  length.out = as.integer(params$nlambda)))

  fargs <- list(lambda = lams, verbose = FALSE,
                Rmethod = params$Rmethod, tol = params$Rtol,
                use.nearPD = isTRUE(params$use_nearPD),
                nu = params$nu, ratio = params$ratio)

  full <- do.call(hugeKmb, c(fargs, list(data = qdat)))

  merge <- lapply(seq_along(lams), function(i)
    Matrix::Matrix(0, p, p, sparse = TRUE))
  for (k in seq_len(reps)) {
    qsub <- qdat[idx[k, ], , drop = FALSE]
    hk <- do.call(hugeKmb, c(fargs, list(data = qsub)))
    for (i in seq_along(lams)) merge[[i]] <- merge[[i]] + hk$path[[i]]
  }

  summary_raw <- vapply(seq_along(lams), function(i) {
    m <- as.matrix(merge[[i]]) / reps
    4 * sum(m * (1 - m)) / (p * (p - 1))
  }, numeric(1))
  mono <- cummax(summary_raw)
  thresh <- params$thresh
  opt_index <- if (any(mono >= thresh)) {
    max(which.max(mono >= thresh)[1] - 1, 1)
  } else {
    length(mono)
  }

  refit <- as.matrix(full$path[[opt_index]])
  ei <- which(refit != 0 & upper.tri(refit), arr.ind = TRUE)
  edges <- if (nrow(ei)) {
    paste(sprintf("%d:%d", ei[, 1] - 1L, ei[, 2] - 1L), collapse = ";")
  } else {
    ""
  }

  beta <- as.matrix(full$beta[[opt_index]])
  bi <- which(beta != 0, arr.ind = TRUE)
  beta_str <- if (nrow(bi)) {
    paste(sprintf("%d:%d:%.6g", bi[, 1] - 1L, bi[, 2] - 1L, beta[bi]),
          collapse = ";")
  } else {
    ""
  }

  mo <- as.matrix(merge[[opt_index]]) / reps
  mi <- which(upper.tri(mo), arr.ind = TRUE)
  merge_str <- if (nrow(mi)) {
    paste(sprintf("%d:%d:%.6g", mi[, 1] - 1L, mi[, 2] - 1L, mo[mi]),
          collapse = ";")
  } else {
    ""
  }

  values <- c(lams, mono, refit, beta)
  for (i in seq_along(lams)) values <- c(values, as.matrix(merge[[i]]))
  values <- c(values, as.matrix(full$path[[opt_index]]))
  nonfinite <- any(!is.finite(values))

  data.frame(
    case = case,
    status = if (nonfinite) "nonfinite" else "finite",
    error_class = "",
    opt_index = as.integer(opt_index),
    n_edges = as.integer(nrow(ei)),
    edges = if (nonfinite) "" else edges,
    beta_opt = if (nonfinite) "" else beta_str,
    lambda_first = as.numeric(lams[[1]]),
    stars_margin = as.numeric(min(abs(mono - thresh))),
    merge_opt = merge_str,
    seconds = as.numeric(proc.time()[["elapsed"]] - start),
    stringsAsFactors = FALSE
  )
}

case_dirs <- sort(list.dirs(cases_dir, recursive = FALSE, full.names = TRUE))
rows <- vector("list", length(case_dirs))
for (i in seq_along(case_dirs)) {
  dir <- case_dirs[[i]]
  rows[[i]] <- tryCatch(
    suppressWarnings(classify_case(dir)),
    error = function(e) {
      row <- empty_row(basename(dir))
      row$error_class <- class(e)[[1]]
      row
    }
  )
}
utils::write.csv(do.call(rbind, rows), status_path, row.names = FALSE, na = "")
