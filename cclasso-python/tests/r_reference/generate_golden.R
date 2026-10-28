#!/usr/bin/env Rscript
################################################################################
# File: generate_golden.R
# Aim : Generate deterministic golden fixtures for the Python cclasso port.
#
# The script sources the PINNED reference R source in
#   tests/r_reference/reference/R/cclasso.R   (CCLasso v2.0)
# and writes, for each case, the input, the parameters, the injected random
# indices, all significant intermediate quantities and the final outputs, with
# full double precision (%.17g).
#
# The generated fixtures are committed under tests/golden/ so the Python test
# suite never needs R.
#
# Usage (from the cclasso-python/ directory):
#     Rscript tests/r_reference/generate_golden.R [out_dir]
# Default out_dir: tests/golden
#
# NOTE: The bootstrap indices and the artificial-null draws are *replayed* with
# the same seed that is used for the real run (see .replay_indices). The R RNG
# stream is consumed in exactly the order sample() -> rnorm() inside
# boot_cclasso, and nothing before it consumes RNG, so the replay is exact; this
# is asserted for every case.
################################################################################

args <- commandArgs(trailingOnly = TRUE)
out_root <- if (length(args) >= 1) args[1] else file.path("tests", "golden")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

ref_r <- file.path("tests", "r_reference", "reference", "R", "cclasso.R")
if (!file.exists(ref_r)) {
  stop("Reference R source not found at: ", ref_r,
       "\nRun this script from the cclasso-python/ directory.")
}
source(ref_r)

## ------------------------------------------------------------------ helpers --
write_mat <- function(m, path) {
  m <- as.matrix(m)
  s <- matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m))
  write.table(s, path, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
}

write_vec <- function(v, path) {
  writeLines(sprintf("%.17g", as.numeric(v)), path)
}

.json_val <- function(x) {
  if (is.null(x)) return("null")
  if (is.character(x)) {
    return(paste0("\"", gsub("\"", "\\\\\"", x), "\""))
  }
  if (is.logical(x)) {
    if (length(x) == 1) return(if (x) "true" else "false")
    return(paste0("[", paste(ifelse(x, "true", "false"), collapse = ","), "]"))
  }
  # numeric
  if (length(x) == 1) return(trimws(sprintf("%.17g", x)))
  return(paste0("[", paste(trimws(sprintf("%.17g", as.numeric(x))), collapse = ","), "]"))
}

write_json <- function(obj, path) {
  stopifnot(is.list(obj), !is.null(names(obj)))
  items <- vapply(names(obj), function(nm)
    paste0("\"", nm, "\":", .json_val(obj[[nm]])), character(1))
  writeLines(paste0("{", paste(items, collapse = ","), "}"), path)
}

## --------------------------------------------------------------- data gen --
# Logistic-normal fractions, optionally with a sparse known precision structure.
make_case_data <- function(n, p, seed, structure = FALSE, counts = FALSE,
                           low_totals = FALSE) {
  set.seed(seed)
  Theta <- diag(p)
  if (structure && p >= 6) {
    edges <- rbind(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6))
    edges <- edges[edges[, 1] <= p & edges[, 2] <= p, , drop = FALSE]
    for (e in seq_len(nrow(edges))) {
      i <- edges[e, 1]; j <- edges[e, 2]
      Theta[i, j] <- 0.45; Theta[j, i] <- 0.45
    }
    diag(Theta) <- diag(Theta) + 1.0
  }
  Sigma <- solve(Theta)
  E <- chol(Sigma)
  z <- matrix(rnorm(n * p), nrow = n) %*% E
  frac <- exp(z) / rowSums(exp(z))
  if (!counts) return(frac)
  tot <- if (low_totals) sample(15:60, n, replace = TRUE) else
    sample(2000:5000, n, replace = TRUE)
  x <- t(sapply(seq_len(n), function(i) rmultinom(1, size = tot[i], prob = frac[i, ])))
  storage.mode(x) <- "double"
  x
}

## ------------------------------------------------ trace / replay machinery --
GOLDEN <- new.env(parent = emptyenv())

# Replay the RNG draws consumed inside boot_cclasso, in call order:
#   sample(1:n, n*n_boot, replace = TRUE)  then  rnorm(n*p)
# Nothing before boot_cclasso consumes RNG in cclasso().
replay_indices <- function(seed, n, p, n_boot) {
  set.seed(seed)
  boot <- matrix(sample(1:n, size = n * n_boot, replace = TRUE), ncol = n_boot)
  null_norm <- matrix(rnorm(n * p), nrow = n)
  list(boot = boot, null_norm = null_norm)
}

run_traced <- function(seed, x, counts, pseudo, k_cv, lam_int, k_max, n_boot) {
  GOLDEN$boot <- NULL
  GOLDEN$sam_art0 <- NULL
  GOLDEN$iters <- NULL
  GOLDEN$relative_errors <- NULL

  trace("boot_cclasso", exit = quote({
    GOLDEN$boot <- sam_boot
    GOLDEN$sam_art0 <- sam_art0
    GOLDEN$vars2 <- vars2
  }), print = FALSE)
  trace("cclasso_sub", exit = quote({
    GOLDEN$iters <- c(GOLDEN$iters, k)
    GOLDEN$relative_errors <- c(GOLDEN$relative_errors, err)
  }), print = FALSE)
  on.exit({
    untrace("boot_cclasso"); untrace("cclasso_sub")
  }, add = TRUE)

  set.seed(seed)
  res <- cclasso(x = x, counts = counts, pseudo = pseudo, k_cv = k_cv,
                 lam_int = lam_int, k_max = k_max, n_boot = n_boot)
  res
}

## --------------------------------------------------------------- one case --
run_case <- function(id, n, p, seed, structure, counts, low_totals,
                     pseudo, k_cv, lam_int, k_max, n_boot, input = NULL,
                     dataset = "QMP", provenance = list()) {
  x <- if (is.null(input)) {
    make_case_data(n = n, p = p, seed = seed, structure = structure,
                   counts = counts, low_totals = low_totals)
  } else {
    as.matrix(input)
  }
  n <- nrow(x)
  p <- ncol(x)
  cat(sprintf("[case %s] n=%d p=%d counts=%s\n", id, n, p, counts))

  reps <- replay_indices(seed, n, p, n_boot)
  res <- tryCatch(
    run_traced(seed, x, counts, pseudo, k_cv, lam_int, k_max, n_boot),
    error = function(e) e
  )

  case_params <- c(
    list(case = id, dataset = dataset, counts = counts,
         pseudo = pseudo, k_cv = k_cv,
         lam_int = lam_int, k_max = k_max, n_boot = n_boot,
         seed = seed, n = n, p = p, structure_adj = structure),
    provenance
  )

  if (inherits(res, "error")) {
    # Preserve R-side failures as first-class equivalence fixtures so one failed
    # seed never aborts generation of the remaining public-data cases.
    failure_dir <- file.path(out_root, "failures", id)
    dir.create(failure_dir, showWarnings = FALSE, recursive = TRUE)
    write_mat(x, file.path(failure_dir, "input.csv"))
    if (!is.null(colnames(x))) writeLines(colnames(x), file.path(failure_dir, "taxa.txt"))
    if (!is.null(rownames(x))) writeLines(rownames(x), file.path(failure_dir, "sample_ids.txt"))
    write_mat(reps$boot, file.path(failure_dir, "boot_indices.csv"))
    write_mat(reps$null_norm, file.path(failure_dir, "null_norm.csv"))
    write_vec(GOLDEN$iters, file.path(failure_dir, "completed_admm_iters.csv"))
    write_vec(GOLDEN$relative_errors,
              file.path(failure_dir, "completed_admm_relative_errors.csv"))
    write_json(c(case_params, list(expected_status = "error",
                                   error_class = class(res)[1],
                                   error_message = conditionMessage(res))),
               file.path(failure_dir, "params.json"))
    writeLines(c("R_status=error", paste0("R_error=", conditionMessage(res))),
               file.path(failure_dir, "status.txt"))
    warning(sprintf("case %s: R reference failed (%s); saved under %s",
                    id, conditionMessage(res), failure_dir))
    return(invisible(list(status = "error", case = id)))
  }

  # verify the replay assumptions (RNG consumed only as sample() then rnorm())
  if (!identical(dim(GOLDEN$boot), dim(reps$boot)) ||
        !all(GOLDEN$boot == reps$boot)) {
    stop("Replay mismatch for bootstrap indices in case ", id)
  }
  recon_null <- reps$null_norm * rep(sqrt(GOLDEN$vars2), each = n)
  if (!isTRUE(all.equal(recon_null, GOLDEN$sam_art0, tolerance = 1e-15))) {
    stop("Replay mismatch for artificial-null draws in case ", id)
  }

  dir <- file.path(out_root, id)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  write_mat(x, file.path(dir, "input.csv"))
  if (!is.null(colnames(x))) writeLines(colnames(x), file.path(dir, "taxa.txt"))
  if (!is.null(rownames(x))) writeLines(rownames(x), file.path(dir, "sample_ids.txt"))

  write_json(case_params, file.path(dir, "params.json"))

  write_mat(reps$boot, file.path(dir, "boot_indices.csv"))
  write_mat(reps$null_norm, file.path(dir, "null_norm.csv"))

  # Intermediates: transform, log-covariance, weights
  xt <- x
  if (counts) { xt <- xt + pseudo; xt <- xt / rowSums(xt) }
  xt <- log(xt)
  write_mat(xt, file.path(dir, "x_log.csv"))
  vx2 <- stats::var(xt)
  write_mat(vx2, file.path(dir, "vx2.csv"))
  rmean_vx2 <- rowMeans(vx2)
  wd <- 1 / diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2))
  write_vec(wd, file.path(dir, "wd.csv"))

  # Cross-validation information
  write_vec(res$info_cv$lams, file.path(dir, "lams.csv"))
  write_vec(res$info_cv$fvals, file.path(dir, "fvals.csv"))
  write_json(list(k = res$info_cv$k), file.path(dir, "info_cv.json"))
  write_json(list(lambda = res$lambda, lam_lo = res$info_cv$lam_int[1],
                  lam_hi = res$info_cv$lam_int[2]),
             file.path(dir, "lambda.json"))

  # Final outputs
  write_vec(res$var_w, file.path(dir, "var_w.csv"))
  write_mat(res$cor_w, file.path(dir, "cor_w.csv"))
  write_mat(res$p_vals, file.path(dir, "p_vals.csv"))

  # ADMM iteration counts, in call order
  write_vec(GOLDEN$iters, file.path(dir, "admm_iters.csv"))
  write_vec(GOLDEN$relative_errors, file.path(dir, "admm_relative_errors.csv"))

  invisible(NULL)
}

## ------------------------------------------------------------------- cases --
# QMP is vendored under tests/r_reference/data/ so fixture regeneration works
# from a clean checkout. QMP_RDA can still override it, and the historical
# ../SPRING-master location is used as a fallback for developer trees.
qmp_rda <- Sys.getenv("QMP_RDA", file.path("tests", "r_reference", "data", "QMP.rda"))
if (!file.exists(qmp_rda)) {
  legacy_qmp <- file.path("..", "SPRING-master", "data", "QMP.rda")
  if (file.exists(legacy_qmp)) {
    qmp_rda <- legacy_qmp
  } else {
    stop("QMP data not found at ", qmp_rda, " or ", legacy_qmp)
  }
}
load(qmp_rda)
qmp <- as.matrix(QMP)
prevalence_order <- order(colSums(qmp > 0), colnames(qmp))

# Full-sample real-data case: 20 least-prevalent taxa, retaining all samples.
qmp_full <- qmp[, prevalence_order[seq_len(20)], drop = FALSE]

# Real p > n case: 30 taxa distributed across prevalence ranks and a deterministic
# 20-sample subset containing at least one observation for every selected taxon.
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

cases <- list(
  list(id = "qmp_n106_p20", n = nrow(qmp_full), p = ncol(qmp_full),
       seed = 1101, structure = FALSE, counts = TRUE, low_totals = FALSE,
       pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20,
       n_boot = 5, input = qmp_full,
       provenance = list(
         source = "SPRING::QMP pinned copy at tests/r_reference/data/QMP.rda",
         source_package = "SPRING",
         source_package_version = "0.99.0",
         source_license = "GPL-3",
         source_citation = paste0(
           "Vandeputte et al. (2017) Nature 551:507-511; ",
           "Yoon, Gaynanova and Muller (2019) Front Genet 10:516"
         ),
         retrieval_date = "2026-10-20",
         redistribution_cleared = FALSE
       )),
  list(id = "qmp_n20_p30", n = nrow(qmp_pgt_n), p = ncol(qmp_pgt_n),
       seed = 1102, structure = FALSE, counts = TRUE, low_totals = FALSE,
       pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20,
       n_boot = 3, input = qmp_pgt_n,
       provenance = list(
         source = "SPRING::QMP pinned copy at tests/r_reference/data/QMP.rda",
         source_package = "SPRING",
         source_package_version = "0.99.0",
         source_license = "GPL-3",
         source_citation = paste0(
           "Vandeputte et al. (2017) Nature 551:507-511; ",
           "Yoon, Gaynanova and Muller (2019) Front Genet 10:516"
         ),
         retrieval_date = "2026-10-20",
         redistribution_cleared = FALSE
       ))
)

for (cs in cases) {
  do.call(run_case, cs)
}

## -------------------------------------------------- public independent cases --
public_helper <- file.path("..", "scripts", "public_microbiome_cases.R")
if (!file.exists(public_helper)) stop("shared public-data selector not found: ", public_helper)
source(public_helper)
public_cases <- load_public_microbiome_cases()
public_seeds <- c(3101L, 3102L, 3103L)
for (public_case in public_cases) {
  for (seed in public_seeds) {
    id <- sprintf("%s_seed%d", public_case$name, seed)
    run_case(
      id = id,
      n = public_case$n,
      p = public_case$p,
      seed = seed,
      structure = FALSE,
      counts = TRUE,
      low_totals = FALSE,
      pseudo = 0.5,
      k_cv = 3,
      lam_int = c(1e-4, 1),
      k_max = 20,
      n_boot = 3,
      input = public_case$input,
      dataset = public_case$dataset,
      provenance = list(
        source = public_case$source,
        source_package = public_case$source_package,
        source_package_version = public_case$source_package_version,
        source_license = public_case$source_license,
        source_citation = public_case$source_citation,
        retrieval_date = public_case$retrieval_date,
        redistribution_cleared = public_case$redistribution_cleared,
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
    )
  }
}

# Additional sparse-feature stress case: retained separately because the R
# reference may error in its bootstrap ADMM condition for this seed. The runner
# records an explicit failure artifact instead of aborting all golden generation.
global_sparse <- load_public_microbiome_cases("spread")[[4]]
run_case(
  id = "globalpatterns_n26_p40_spread_seed3102",
  n = global_sparse$n,
  p = global_sparse$p,
  seed = 3102L,
  structure = FALSE,
  counts = TRUE,
  low_totals = FALSE,
  pseudo = 0.5,
  k_cv = 3,
  lam_int = c(1e-4, 1),
  k_max = 20,
  n_boot = 3,
  input = global_sparse$input,
  dataset = global_sparse$dataset,
  provenance = list(
    source = global_sparse$source,
    source_package = global_sparse$source_package,
    source_package_version = global_sparse$source_package_version,
    source_license = global_sparse$source_license,
    source_citation = global_sparse$source_citation,
    retrieval_date = global_sparse$retrieval_date,
    redistribution_cleared = global_sparse$redistribution_cleared,
    source_n = global_sparse$raw_n,
    source_p = global_sparse$raw_p,
    available_n = global_sparse$available_n,
    available_p = global_sparse$available_p,
    zero_fraction = global_sparse$zero_fraction,
    library_size_min = global_sparse$library_size_min,
    library_size_median = global_sparse$library_size_median,
    library_size_max = global_sparse$library_size_max,
    taxa_strategy = global_sparse$taxa_strategy
  )
)

## ---------------------------------------------- synthetic coverage cases --
# Close R-golden coverage gaps that the real-data families leave open:
#   * the counts = FALSE branch (no pseudo-count / row-normalisation),
#   * non-default k_cv, n_boot and lam_int values.
# Inputs come from the deterministic logistic-normal generator above, so these
# cases carry no third-party data-provenance constraints.
synthetic_cases <- list(
  # counts = FALSE baseline; n %/% k_cv leaves a remainder (37 %/% 3 = 12).
  list(id = "synth_frac_n37_p8_k3_b4_seed4101", n = 37, p = 8, seed = 4101,
       structure = TRUE, counts = FALSE, low_totals = FALSE, pseudo = 0.5,
       k_cv = 3, lam_int = c(1e-4, 1), k_max = 20, n_boot = 4),
  # counts = FALSE with non-default k_cv.
  list(id = "synth_frac_n45_p8_k5_b3_seed4106", n = 45, p = 8, seed = 4106,
       structure = TRUE, counts = FALSE, low_totals = FALSE, pseudo = 0.5,
       k_cv = 5, lam_int = c(1e-4, 1), k_max = 20, n_boot = 3),
  # counts = FALSE with a non-default lam_int interval.
  list(id = "synth_frac_n34_p8_k3_b5_lam_seed4103", n = 34, p = 8, seed = 4103,
       structure = TRUE, counts = FALSE, low_totals = FALSE, pseudo = 0.5,
       k_cv = 3, lam_int = c(1e-3, 1e-1), k_max = 20, n_boot = 5),
  # counts = TRUE with non-default k_cv, n_boot and lam_int together.
  list(id = "synth_counts_n36_p8_k4_b7_lam_seed4104", n = 36, p = 8, seed = 4104,
       structure = TRUE, counts = TRUE, low_totals = FALSE, pseudo = 0.5,
       k_cv = 4, lam_int = c(1e-5, 1), k_max = 20, n_boot = 7)
)
for (cs in synthetic_cases) {
  do.call(run_case, c(cs, list(
    dataset = "Synthetic logistic-normal fractions",
    provenance = list(
      source = "deterministic generator in tests/r_reference/generate_golden.R",
      source_license = "n/a (synthetic)",
      redistribution_cleared = TRUE
    )
  )))
}

## ------------------------------------------------------------- session info --
si <- file.path(out_root, "sessionInfo.txt")
sink(si)
cat("Generated by tests/r_reference/generate_golden.R\n")
cat("R version:", R.version.string, "\n\n")
print(sessionInfo())
sink()

cat("Done. Fixtures written to:", normalizePath(out_root), "\n")
