#!/usr/bin/env Rscript
# Quantify the R latentcor approximation-vs-inversion sensitivity on existing
# QMP fixtures. No SPRING/MB or StARS stages are run here.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript compare_approx_original.R <golden_dir> <out_dir>")
}
golden_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(latentcor))
params_json <- requireNamespace("jsonlite", quietly = TRUE)
if (!params_json) stop("This diagnostic uses jsonlite to read fixture parameters")

cases <- c(
  "qmp_n106_p20", "qmp_n20_p30",
  "amgut_r1_n80_p20_seed3101", "amgut_r1_n20_p30_seed3101",
  "globalpatterns_n26_p15_seed3101", "globalpatterns_n26_p40_seed3101"
)
for (case in cases) {
  case_dir <- file.path(golden_dir, case)
  params <- jsonlite::fromJSON(file.path(case_dir, "params.json"))
  qdat <- as.matrix(read.csv(file.path(case_dir, "qdat.csv"), header = FALSE))
  fit <- function(method) latentcor::latentcor(
    qdat,
    types = "tru",
    method = method,
    tol = params$Rtol,
    use.nearPD = params$`use.nearPD`,
    nu = params$nu,
    ratio = params$ratio,
    showplot = FALSE
  )
  approx <- fit("approx")
  original <- fit("original")
  p <- ncol(qdat)
  indices <- which(lower.tri(matrix(FALSE, p, p)), arr.ind = TRUE)
  zratio <- colMeans(sweep(qdat, 2L, apply(qdat, 2L, min), `==`))
  bound <- 1 - pmax(zratio[indices[, 1]], zratio[indices[, 2]])^2
  inversion <- abs(approx$K[indices]) > params$ratio * bound
  rows <- data.frame(
    i = indices[, 1],
    j = indices[, 2],
    approx_uses_inversion = inversion,
    K = approx$K[indices],
    bound = bound,
    r_approx = approx$Rpointwise[indices],
    r_original = original$Rpointwise[indices]
  )
  rows$abs_method_difference <- abs(rows$r_approx - rows$r_original)
  write.csv(rows, file.path(out_dir, paste0(case, "_bridge_sensitivity.csv")),
            row.names = FALSE, quote = FALSE)
  cat(case, "pairs", nrow(rows),
      "approx inversion region", sum(inversion),
      "max raw Rmethod difference", max(rows$abs_method_difference), "\n")
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
