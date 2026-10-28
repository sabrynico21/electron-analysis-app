#!/usr/bin/env Rscript
# Write reference-runtime versions and numerical-library metadata for a run.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript write_r_environment.R <output.json>")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")

packages <- c(
  "SpiecEasi", "phyloseq", "latentcor", "huge", "pulsar", "mnormt",
  "Matrix", "jsonlite", "pcaPP", "fMultivar", "Rcpp", "RcppArmadillo"
)
package_versions <- lapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NULL
  }
})
names(package_versions) <- packages

manifest <- list(
  r_version = R.version.string,
  platform = R.version$platform,
  os = as.list(Sys.info()[c("sysname", "release", "version", "machine")]),
  numerical_libraries = as.list(extSoftVersion()),
  rng_kind = as.list(RNGkind()),
  packages = package_versions,
  session_info = capture.output(sessionInfo())
)
out <- args[[1]]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(manifest, out, pretty = TRUE, auto_unbox = TRUE, null = "null")
