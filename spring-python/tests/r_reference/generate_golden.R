#!/usr/bin/env Rscript
# =============================================================================
# generate_golden.R — deterministic golden fixtures for the SPRING Python port.
# =============================================================================
# For every QMP and public-dataset case this script writes, at FULL DOUBLE
# PRECISION (17 significant
# digits, round-trip safe), enough information for the Python test suite to
# replay the exact computation *without* R:
#
#   input.csv              raw count matrix (n x p), samples in rows
#   params.json            SPRING parameters + subsample indices (0-based too)
#   qdat.csv               matrix actually fed to pulsar (raw or mclr output)
#   lambdaseq.csv          lambda path used
#   full_K.csv             full-data Kendall tau matrix (latentcor)
#   full_Rpointwise.csv    full-data pointwise latent correlation (pre nearPD)
#   full_R.csv             full-data latent correlation after nearPD + nu
#   merge_<i>.csv          selection-frequency matrix at lambda i (i = 1..L)
#   summary.csv            StARS instability curve (after cummax)
#   opt_index.txt          selected lambda index
#   beta_full_<i>.csv      MB coefficient matrix at lambda i (full data)
#   adjacency_full_<i>.csv MB adjacency at lambda i (full data, sym = "or")
#   subsample_indices.csv  1-based R indices used by pulsar (rep.num rows)
#   subsample_R_<k>.csv    latentcor R for subsample k (k = 1..rep.num)
#   subsample_paths_long.csv  long format: k,lambda,row,col  (MB adjacency nonzeros)
#   subsample_beta_long.csv   long format: k,lambda,row,col,value (MB coeff nonzeros)
#
# The indices are generated with the *same* R expressions pulsar uses, so the
# Python side can inject them and remove the RNG difference entirely.
#
# Usage:
#   Rscript tests/r_reference/generate_golden.R [out_dir]
#
# Self-contained: the SPRING source, QMP data and the public-data selector are
# vendored under tests/r_reference/reference/. Only CRAN/Bioconductor packages
# (latentcor, huge, pulsar, and for the public cases SpiecEasi/phyloseq) are
# external. SPRING_R_HOME can still override the vendored source.
# =============================================================================

suppressMessages({
  library(latentcor)
  library(huge)
  library(pulsar)
})

# --- vendored-reference locations --------------------------------------------
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
ref_candidates <- c(
  file.path("tests", "r_reference", "reference"),
  file.path("reference"),
  file.path(script_dir, "reference")
)
ref_root <- ref_candidates[file.exists(ref_candidates)][1]
if (is.na(ref_root)) ref_root <- file.path("tests", "r_reference", "reference")

public_helper <- file.path(ref_root, "public_microbiome_cases.R")
if (!file.exists(public_helper)) stop("vendored public-data selector not found: ", public_helper)
source(public_helper)
public_cases <- load_public_microbiome_cases()

options(digits = 17, warn = 1)
set.seed(20261016)

args <- commandArgs(trailingOnly = TRUE)
out_root <- if (length(args) >= 1) args[[1]] else "tests/golden"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# --- load the SPRING source under test (vendored copy first) ------------------
find_spring <- function() {
  cands <- c(
    file.path(ref_root, "R"),
    file.path("..", "SPRING-master", "R"),
    file.path("..", "..", "SPRING-master", "R"),
    Sys.getenv("SPRING_R_HOME")
  )
  for (c in cands) {
    if (nzchar(c) && file.exists(file.path(c, "SPRING.R"))) return(normalizePath(c))
  }
  stop("Cannot find SPRING/R. Set SPRING_R_HOME or keep the vendored reference.")
}
spring_r_dir <- find_spring()
for (f in c("helpers.R", "SPRING.R", "synthData.R")) {
  source(file.path(spring_r_dir, f))
}

# --- full-precision writers --------------------------------------------------
write_mat <- function(x, path) {
  x <- as.matrix(x)
  s <- sprintf("%.17g", x)
  dim(s) <- dim(x)
  write.table(s, path, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
}
write_vec <- function(x, path) write_mat(matrix(x, ncol = 1), path)
write_long <- function(df, path) {
  if (nrow(df) == 0) {
    write.table(df, path, sep = ",", row.names = FALSE, col.names = TRUE,
                quote = FALSE)
    return(invisible())
  }
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], function(v) sprintf("%.17g", v))
  write.table(df, path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = FALSE)
}

# --- data generators (deterministic; result is written to CSV) ---------------
# Build a p x p precision matrix with a random network and return a count
# matrix with controllable zero-inflation. Counts are always non-negative.
make_counts <- function(n, p, n_edges, zero_frac, seed,
                        strength = 0.7, scale = 1.3, mu = 2.5) {
  set.seed(seed)
  A <- matrix(0, p, p)
  if (p > 1) {
    for (e in seq_len(n_edges)) {
      i <- sample.int(p, 1); j <- sample.int(p, 1)
      while (j == i) j <- sample.int(p, 1)
      A[i, j] <- A[j, i] <- 1
    }
  }
  Omega <- A * strength
  diag(Omega) <- 1
  eig <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) <= 0) diag(Omega) <- diag(Omega) + (abs(min(eig)) + 0.3)
  Sigma <- solve(Omega)
  d <- sqrt(diag(Sigma))
  Sigma <- Sigma / outer(d, d)
  Z <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
  eta <- mu + scale * Z
  counts <- matrix(0L, n, p)
  for (j in seq_len(p)) {
    lam <- exp(eta[, j])
    counts[, j] <- as.integer(round(lam))
  }
  # zero-inflation: knock out the smallest values per column
  if (zero_frac > 0) {
    for (j in seq_len(p)) {
      k <- floor(zero_frac * n)
      if (k > 0) {
        idx <- order(counts[, j])[seq_len(k)]
        counts[idx, j] <- 0L
      }
    }
  }
  counts
}

# --- one case -----------------------------------------------------------------
run_case <- function(name, counts, params) {
  cat("== case:", name, "  dim:", paste(dim(counts), collapse = "x"), "\n")
  case_dir <- file.path(out_root, name)
  dir.create(case_dir, showWarnings = FALSE, recursive = TRUE)

  counts <- as.matrix(counts)
  feature_names <- colnames(counts)
  sample_ids <- rownames(counts)
  if (is.null(feature_names)) feature_names <- paste0("t", seq_len(ncol(counts)))
  if (is.null(sample_ids)) sample_ids <- sprintf("sample-%06d", seq_len(nrow(counts)))
  storage.mode(counts) <- "double"
  colnames(counts) <- paste0("t", seq_len(ncol(counts)))
  write_mat(counts, file.path(case_dir, "input.csv"))
  writeLines(feature_names, file.path(case_dir, "taxa.txt"))
  writeLines(sample_ids, file.path(case_dir, "sample_ids.txt"))

  n <- nrow(counts); p <- ncol(counts)

  # --- parameters with SPRING defaults ---------------------------------------
  p_lambda_min_ratio <- params$lambda.min.ratio %||% 1e-2
  p_nlambda          <- params$nlambda %||% 20
  p_seed             <- params$seed %||% 10010
  p_thresh           <- params$thresh %||% 0.1
  p_subsample_ratio  <- params$subsample.ratio %||% 0.8
  p_rep_num          <- params$rep.num %||% 20
  p_Rtol             <- params$Rtol %||% 1e-6
  p_Rmethod          <- params$Rmethod %||% "approx"
  p_use_nearPD       <- params$use.nearPD %||% TRUE
  p_nu               <- params$nu %||% 0.001
  p_ratio            <- params$ratio %||% 0.9
  p_quantitative     <- params$quantitative %||% FALSE
  p_lambdaseq_spec   <- params$lambdaseq %||% "explicit"
  p_lambdaseq        <- if (is.character(p_lambdaseq_spec)) {
    exp(seq(log(0.6), log(0.6 * p_lambda_min_ratio), length.out = p_nlambda))
  } else {
    p_lambdaseq_spec
  }

  # data-specific lambda sequence, if requested
  if (is.character(p_lambdaseq_spec) && p_lambdaseq_spec == "data-specific") {
    qdat_ds <- if (p_quantitative) counts else mclr(counts)
    Kcor <- suppressMessages(
      latentcor::latentcor(qdat_ds, types = "tru", method = p_Rmethod,
                           tol = p_Rtol, use.nearPD = p_use_nearPD,
                           nu = p_nu, ratio = p_ratio, showplot = FALSE)$R)
    lambda_max <- max(max(Kcor - diag(p)), -min(Kcor - diag(p)))
    lambda_min <- p_lambda_min_ratio * lambda_max
    p_lambdaseq <- exp(seq(log(lambda_max), log(lambda_min), length = p_nlambda))
  }

  # --- qdat (mclr) and full-data latentcor intermediates ---------------------
  qdat <- if (p_quantitative) counts else mclr(counts)
  write_mat(qdat, file.path(case_dir, "qdat.csv"))

  lc <- suppressMessages(
    latentcor::latentcor(qdat, types = "tru", method = p_Rmethod, tol = p_Rtol,
                         use.nearPD = p_use_nearPD, nu = p_nu, ratio = p_ratio,
                         showplot = FALSE))
  write_mat(lc$K, file.path(case_dir, "full_K.csv"))
  write_mat(lc$Rpointwise, file.path(case_dir, "full_Rpointwise.csv"))
  write_mat(lc$R, file.path(case_dir, "full_R.csv"))
  write_vec(p_lambdaseq, file.path(case_dir, "lambdaseq.csv"))

  # --- subsample indices: identical to what pulsar will draw -----------------
  ratio <- pulsar:::.ratcheck(p_subsample_ratio, n)
  set.seed(p_seed)
  ind.sample <- replicate(p_rep_num, sample(c(1L:n), floor(n * ratio),
                                            replace = FALSE), simplify = FALSE)
  set.seed(NULL)
  write_mat(t(vapply(ind.sample, as.numeric, numeric(floor(n * ratio)))),
            file.path(case_dir, "subsample_indices.csv"))

  # --- per-subsample latentcor R and MB paths --------------------------------
  fargs <- list(lambda = p_lambdaseq, Rmethod = p_Rmethod, tol = p_Rtol,
                verbose = FALSE, use.nearPD = p_use_nearPD, nu = p_nu,
                ratio = p_ratio)
  path_long <- list(); beta_long <- list()
  for (k in seq_len(p_rep_num)) {
    rows <- ind.sample[[k]]
    qsub <- qdat[rows, , drop = FALSE]
    Rs <- suppressMessages(
      latentcor::latentcor(qsub, types = "tru", method = p_Rmethod, tol = p_Rtol,
                           use.nearPD = p_use_nearPD, nu = p_nu,
                           ratio = p_ratio, showplot = FALSE)$R)
    write_mat(Rs, file.path(case_dir, sprintf("subsample_R_%02d.csv", k)))
    fit <- do.call(hugeKmb, c(fargs, list(data = qsub)))
    for (i in seq_along(p_lambdaseq)) {
      adj <- as.matrix(fit$path[[i]])
      idx <- which(adj != 0, arr.ind = TRUE)
      if (nrow(idx) > 0) {
        path_long[[length(path_long) + 1]] <- data.frame(
          k = k, lambda = i, row = idx[, 1] - 1L, col = idx[, 2] - 1L)
      }
      bet <- as.matrix(fit$beta[[i]])
      bidx <- which(bet != 0, arr.ind = TRUE)
      if (nrow(bidx) > 0) {
        beta_long[[length(beta_long) + 1]] <- data.frame(
          k = k, lambda = i, row = bidx[, 1] - 1L, col = bidx[, 2] - 1L,
          value = bet[bidx])
      }
    }
  }
  write_long(if (length(path_long)) do.call(rbind, path_long) else
               data.frame(k = integer(), lambda = integer(), row = integer(),
                          col = integer()),
             file.path(case_dir, "subsample_paths_long.csv"))
  write_long(if (length(beta_long)) do.call(rbind, beta_long) else
               data.frame(k = integer(), lambda = integer(), row = integer(),
                          col = integer(), value = numeric()),
             file.path(case_dir, "subsample_beta_long.csv"))

  # --- the real SPRING run ---------------------------------------------------
  fit_list <- suppressWarnings(SPRING(
    counts, quantitative = p_quantitative, lambda.min.ratio = p_lambda_min_ratio,
    nlambda = p_nlambda, lambdaseq = if (identical(p_lambdaseq_spec, "data-specific"))
      "data-specific" else p_lambdaseq,
    seed = p_seed, ncores = 1, thresh = p_thresh,
    subsample.ratio = p_subsample_ratio, rep.num = p_rep_num, Rtol = p_Rtol,
    verbose = FALSE, Rmethod = p_Rmethod, use.nearPD = p_use_nearPD,
    nu = p_nu, ratio = p_ratio))

  stars <- fit_list$output$stars
  for (i in seq_along(stars$merge)) {
    write_mat(stars$merge[[i]], file.path(case_dir, sprintf("merge_%02d.csv", i)))
  }
  write_vec(stars$summary, file.path(case_dir, "summary.csv"))
  writeLines(sprintf("%d", stars$opt.index), file.path(case_dir, "opt_index.txt"))

  full_est <- fit_list$fit$est
  for (i in seq_along(full_est$beta)) {
    write_mat(as.matrix(full_est$beta[[i]]),
              file.path(case_dir, sprintf("beta_full_%02d.csv", i)))
    write_mat(as.matrix(full_est$path[[i]]),
              file.path(case_dir, sprintf("adjacency_full_%02d.csv", i)))
  }
  write_mat(as.matrix(fit_list$fit$refit$stars),
            file.path(case_dir, "adjacency_refit.csv"))

  # --- params.json -----------------------------------------------------------
  params_out <- list(
    case = name, dataset = params$dataset %||% "QMP",
    n = n, p = p,
    quantitative = p_quantitative,
    lambda.min.ratio = p_lambda_min_ratio,
    nlambda = p_nlambda,
    lambdaseq = as.numeric(p_lambdaseq),
    lambdaseq_spec = if (is.character(p_lambdaseq_spec) &&
                         identical(p_lambdaseq_spec, "data-specific")) {
      "data-specific"
    } else if (is.numeric(p_lambdaseq_spec)) {
      "injected_explicit"
    } else {
      "explicit"
    },
    seed = p_seed, thresh = p_thresh,
    subsample.ratio = p_subsample_ratio,
    rep.num = p_rep_num, Rtol = p_Rtol,
    Rmethod = p_Rmethod, use.nearPD = p_use_nearPD,
    nu = p_nu, ratio = p_ratio,
    opt.index = stars$opt.index,
    source = params$source %||% "",
    source_n = params$source_n %||% n,
    source_p = params$source_p %||% p,
    available_n = params$available_n %||% n,
    available_p = params$available_p %||% p,
    zero_fraction = params$zero_fraction %||% mean(counts == 0),
    library_size_min = params$library_size_min %||% min(rowSums(counts)),
    library_size_median = params$library_size_median %||% median(rowSums(counts)),
    library_size_max = params$library_size_max %||% max(rowSums(counts))
  )
  json <- paste0(
    "{\n",
    paste(vapply(names(params_out), function(k) {
      v <- params_out[[k]]
      val <- if (is.character(v)) paste0("\"", v, "\"")
             else if (is.logical(v)) tolower(as.character(v))
             else if (length(v) == 1) sprintf("%.17g", v)
             else paste0("[", paste(sprintf("%.17g", v), collapse = ", "), "]")
      sprintf("  \"%s\": %s", k, val)
    }, character(1)), collapse = ",\n"),
    "\n}\n")
  cat(json, file = file.path(case_dir, "params.json"))

  cat("   opt.index =", stars$opt.index,
      " edges =", sum(as.matrix(fit_list$fit$refit$stars)) / 2, "\n")
  cat(json)
  invisible(NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- minimal JSON writer for the dedicated latentcor case --------------------
write_params_json <- function(x, path) {
  json <- paste0(
    "{\n",
    paste(vapply(names(x), function(k) {
      v <- x[[k]]
      val <- if (is.character(v)) paste0("\"", v, "\"")
             else if (is.logical(v)) tolower(as.character(v))
             else sprintf("%.17g", v)
      sprintf("  \"%s\": %s", k, val)
    }, character(1)), collapse = ",\n"),
    "\n}\n")
  cat(json, file = path)
}

# Rmethod="original" inverts every pair with Brent over the 4-D normal CDF
# (~20 s per pair), so a full StARS golden would be unusably slow. This
# dedicated case compares the *latentcor stage* with R directly instead.
run_latentcor_case <- function(name, counts, params) {
  cat("== latentcor case:", name, "\n")
  case_dir <- file.path(out_root, name)
  dir.create(case_dir, showWarnings = FALSE, recursive = TRUE)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  colnames(counts) <- paste0("t", seq_len(ncol(counts)))
  rownames(counts) <- sprintf("sample-%06d", seq_len(nrow(counts)))
  write_mat(counts, file.path(case_dir, "input.csv"))
  p_quant <- params$quantitative %||% FALSE
  p_Rmethod <- params$Rmethod %||% "original"
  p_Rtol <- params$Rtol %||% 1e-6
  p_near <- params$use.nearPD %||% TRUE
  p_nu <- params$nu %||% 0.001
  p_ratio <- params$ratio %||% 0.9
  qdat <- if (p_quant) counts else mclr(counts)
  write_mat(qdat, file.path(case_dir, "qdat.csv"))
  lc <- suppressMessages(latentcor::latentcor(
    qdat, types = "tru", method = p_Rmethod, tol = p_Rtol,
    use.nearPD = p_near, nu = p_nu, ratio = p_ratio, showplot = FALSE))
  write_mat(lc$K, file.path(case_dir, "full_K.csv"))
  write_mat(lc$Rpointwise, file.path(case_dir, "full_Rpointwise.csv"))
  write_mat(lc$R, file.path(case_dir, "full_R.csv"))
  write_params_json(list(
    case = name, kind = "latentcor", dataset = params$dataset %||% "synthetic",
    n = nrow(qdat), p = ncol(qdat), quantitative = p_quant,
    Rmethod = p_Rmethod, Rtol = p_Rtol, use.nearPD = p_near,
    nu = p_nu, ratio = p_ratio
  ), file.path(case_dir, "params.json"))
  cat("   done\n")
}

# =============================================================================
# QMP real-data cases
# =============================================================================
qmp_rda <- Sys.getenv("QMP_RDA", file.path(ref_root, "data", "QMP.rda"))
if (!file.exists(qmp_rda)) {
  qmp_rda <- "../SPRING-master/data/QMP.rda"
}
if (!file.exists(qmp_rda)) stop("QMP data not found: ", qmp_rda)
load(qmp_rda)
qmp <- as.matrix(QMP)
prevalence_order <- order(colSums(qmp > 0), colnames(qmp))
qmp_full <- qmp[, prevalence_order[seq_len(20)], drop = FALSE]

small_taxa <- prevalence_order[
  unique(round(seq(1, length(prevalence_order), length.out = 30)))
]
qmp_wide <- qmp[, small_taxa, drop = FALSE]
covered <- rep(FALSE, ncol(qmp_wide))
sample_rows <- integer()
while (!all(covered)) {
  candidates <- setdiff(seq_len(nrow(qmp_wide)), sample_rows)
  gains <- vapply(candidates, function(i)
    sum(qmp_wide[i, ] > 0 & !covered), integer(1))
  next_row <- candidates[which.max(gains)]
  sample_rows <- c(sample_rows, next_row)
  covered <- covered | (qmp_wide[next_row, ] > 0)
}
if (length(sample_rows) > 20L) stop("Could not construct a 20-sample QMP p>n case")
sample_rows <- c(sample_rows,
                 head(setdiff(seq_len(nrow(qmp_wide)), sample_rows),
                      20L - length(sample_rows)))
qmp_pgt_n <- qmp_wide[sample_rows, , drop = FALSE]

# QMP is a quantitative microbiome profile in the SPRING package examples.
# The first case retains all 106 samples; the second is a real p > n subset.
cases <- list(
  list(name = "qmp_n106_p20", counts = qmp_full,
       params = list(nlambda = 3, rep.num = 2, seed = 1101,
                     quantitative = TRUE, thresh = 0.8,
                     lambdaseq = "data-specific")),
  list(name = "qmp_n20_p30", counts = qmp_pgt_n,
       params = list(nlambda = 3, rep.num = 2, seed = 1102,
                     quantitative = TRUE, thresh = 0.8,
                     lambdaseq = "data-specific")),
  # --- SPRING DEFAULTS: nlambda = 20, rep.num = 20, thresh = 0.1 ------------
  # The reduced StARS settings used by the other fixtures never exercised the
  # package defaults. This case runs the real default configuration with the
  # *default* lambda path (lambdaseq left unspecified -> "explicit" in the
  # generator, which Python mirrors with lambdaseq=None).
  list(name = "qmp_defaults_n106_p8", counts = qmp_full[, seq_len(8), drop = FALSE],
       params = list(nlambda = 20L, rep.num = 20L, seed = 2101,
                     quantitative = TRUE, thresh = 0.1))
)

for (public_case in public_cases) {
  public_qdat <- mclr(public_case$input)
  public_R <- suppressMessages(latentcor::latentcor(
    public_qdat, types = "tru", method = "approx", tol = 1e-6,
    use.nearPD = TRUE, nu = 0.001, ratio = 0.9, showplot = FALSE
  )$R)
  public_lambda_max <- max(
    max(public_R - diag(ncol(public_R))),
    -min(public_R - diag(ncol(public_R)))
  )
  # Inject the R-derived numeric path for the historical public-data cases.
  # The separate `*_dataseq` cases below exercise Python's own
  # latentcor -> lambda-max -> grid chain on compositional (mclr) data.
  public_lambdaseq <- exp(seq(
    log(public_lambda_max), log(public_lambda_max * 1e-2), length.out = 3L
  ))

  meta <- list(
    dataset = public_case$dataset,
    source = public_case$source,
    taxa_strategy = public_case$taxa_strategy,
    source_n = public_case$raw_n,
    source_p = public_case$raw_p,
    available_n = public_case$available_n,
    available_p = public_case$available_p,
    zero_fraction = public_case$zero_fraction,
    library_size_min = public_case$library_size_min,
    library_size_median = public_case$library_size_median,
    library_size_max = public_case$library_size_max
  )

  for (seed in c(3101L, 3102L, 3103L)) {
    cases[[length(cases) + 1L]] <- list(
      name = sprintf("%s_seed%d", public_case$name, seed),
      counts = public_case$input,
      params = c(meta, list(
        nlambda = 3L,
        rep.num = 2L,
        seed = seed,
        quantitative = FALSE,
        thresh = 0.8,
        lambdaseq = public_lambdaseq
      ))
    )
  }

  # --- data-specific lambda path on compositional public data ---------------
  if (public_case$name %in% c("amgut_r1_n80_p20", "globalpatterns_n26_p15")) {
    cases[[length(cases) + 1L]] <- list(
      name = sprintf("%s_dataseq", public_case$name),
      counts = public_case$input,
      params = c(meta, list(
        nlambda = 5L,
        rep.num = 5L,
        seed = 3201L,
        quantitative = FALSE,
        thresh = 0.1,
        lambdaseq = "data-specific"
      ))
    )
  }

  # --- SPRING defaults on compositional public data -------------------------
  # American Gut has no interpolation/bridge-boundary pairs in any replica, so
  # the 20x20 default StARS run stays cheap (GlobalPatterns is covered by the
  # `_dataseq` case below).
  if (identical(public_case$name, "amgut_r1_n80_p20")) {
    cases[[length(cases) + 1L]] <- list(
      name = sprintf("amgut_r1_n80_p8_defaults"),
      counts = public_case$input[, seq_len(8), drop = FALSE],
      params = c(meta, list(
        nlambda = 20L,
        rep.num = 20L,
        seed = 2103L,
        quantitative = FALSE,
        thresh = 0.1
      ))
    )
  }
}

session_path <- file.path(out_root, "sessionInfo.txt")
sink(session_path)
cat("== SPRING golden generator ==\n")
cat("source of truth:", spring_r_dir, "\n\n")
print(sessionInfo())
sink()

for (cs in cases) {
  run_case(cs$name, cs$counts, cs$params)
}

run_latentcor_case(
  "original_latentcor_n60_p3",
  make_counts(n = 60, p = 3, n_edges = 2, zero_frac = 0.05, seed = 424242),
  list(dataset = "synthetic", Rmethod = "original")
)

cat("\nDone. Goldens in:", normalizePath(out_root), "\n")
