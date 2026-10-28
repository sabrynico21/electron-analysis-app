#!/usr/bin/env Rscript
# Extract the precomputed multilinear-interpolation tables from the installed
# `latentcor` package and dump them in a language-neutral form so that the
# Python port can ship bit-identical interpolation *inputs*.
#
# SPRING always calls latentcor with types = "tru" for every column, therefore
# every off-diagonal pair is truncated/truncated -> the only table that is ever
# used is `ipol_22`. We still export the metadata of the others for reference.
#
# Usage:
#   Rscript scripts/extract_ipol_tables.R [out_dir]
#
# Outputs (little-endian, R column-major flat order preserved):
#   <out_dir>/ipol_22_grid_K.f64   double[41]
#   <out_dir>/ipol_22_grid_z1.f64  double[49]
#   <out_dir>/ipol_22_grid_z2.f64  double[49]
#   <out_dir>/ipol_22_values.i32   int32[41*49*49]
#   <out_dir>/ipol_meta.json       dims, scale, source versions

suppressMessages(library(latentcor))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "src/spring_microbiome/latentcor/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

get_parts <- function(name) {
  o <- get(name, envir = asNamespace("latentcor"))
  f <- get("fun", envir = environment(o))
  e <- environment(f)
  list(grid = get("grid", envir = e), val = get("val", envir = e))
}

write_f64 <- function(x, path) {
  con <- file(path, "wb")
  writeBin(as.double(x), con, size = 8L, endian = "little")
  close(con)
}

write_i32 <- function(x, path) {
  con <- file(path, "wb")
  writeBin(as.integer(x), con, size = 4L, endian = "little")
  close(con)
}

p22 <- get_parts("ipol_22")
grid <- p22$grid
val <- p22$val
dims <- as.integer(sapply(grid, length))
stopifnot(identical(dims, c(41L, 49L, 49L)))
stopifnot(length(val) == prod(dims))
# latentcor builds these tables with `int = TRUE` (values * 1e7 rounded), but
# depending on the package version the array may be stored as double-typed
# whole numbers rather than integer. Accept either, require integral values.
stopifnot(all(val == round(val), na.rm = TRUE))

write_f64(grid[[1]], file.path(out_dir, "ipol_22_grid_K.f64"))
write_f64(grid[[2]], file.path(out_dir, "ipol_22_grid_z1.f64"))
write_f64(grid[[3]], file.path(out_dir, "ipol_22_grid_z2.f64"))
write_i32(val, file.path(out_dir, "ipol_22_values.i32"))

# Metadata for the other tables (not exported, documented for completeness).
other <- c("ipol_10", "ipol_11", "ipol_20", "ipol_21",
           "ipol_30", "ipol_31", "ipol_32", "ipol_33")
other_dims <- lapply(other, function(nm) {
  p <- get_parts(nm)
  as.integer(sapply(p$grid, length))
})
names(other_dims) <- other

meta <- sprintf(
  '{\n  "table": "ipol_22",\n  "dims": [%d, %d, %d],\n  "scale": 10000000,\n  "layout": "column-major",\n  "latentcor_version": "%s",\n  "r_version": "%s",\n  "other_table_dims": {%s}\n}\n',
  dims[1], dims[2], dims[3],
  as.character(packageVersion("latentcor")),
  paste(R.version$major, R.version$minor, sep = "."),
  paste(sprintf('"%s": [%s]', names(other_dims),
                vapply(other_dims, paste, "", collapse = ", ")), collapse = ", ")
)
cat(meta, file = file.path(out_dir, "ipol_meta.json"))

message("Wrote ipol_22 tables to ", out_dir)
message("latentcor ", as.character(packageVersion("latentcor")),
        " | R ", R.version.string)
