#!/usr/bin/env Rscript
# FASE 1 baseline — timing of the *complete* MB+StARS pipeline (R, serial).
#
# Runs SpiecEasi::spiec.easi with the same settings as the Python benchmark
# (nlambda=20, lambda.min.ratio=1e-3, rep.num=20, thresh=0.05) and the same
# pulsar seed used by the golden generator, so the subsample indices coincide.
# pulsar is used serially (ncores = 1, no batchtools registry).
#
# Usage: Rscript benchmarks/bench_r.R [reps] [id_filter] [out.json]

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(pulsar)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1L) as.integer(args[1]) else 3L
id_filter <- if (length(args) >= 2L) args[2] else ".*"
out_json <- if (length(args) >= 3L) args[3] else NA_character_

repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
gold <- file.path(repo, "spiec-easi-python", "tests", "golden")
gold_ext <- file.path(repo, "spiec-easi-python", "tests", "golden_extended")

cases <- list(
  list(id = "ds2_p30_n100", counts = file.path(gold, "ds2_r.counts.tsv"),
       subs = file.path(gold, "ds2_r.subsample_indices.tsv"),
       nlambda = 20L, ratio = 1e-3, rep = 20L, seed = 17L),
  list(id = "rand_p50_n50_s211", counts = file.path(gold_ext, "rand_p50_n50_s211", "counts.tsv"),
       subs = file.path(gold_ext, "rand_p50_n50_s211", "subsample_indices.tsv"),
       nlambda = 20L, ratio = 1e-3, rep = 20L, seed = 17L),
  list(id = "largep_p100_n80_s215", counts = file.path(gold_ext, "largep_p100_n80_s215", "counts.tsv"),
       subs = file.path(gold_ext, "largep_p100_n80_s215", "subsample_indices.tsv"),
       nlambda = 20L, ratio = 1e-3, rep = 20L, seed = 17L)
)

cat("== environment ==\n")
cat("  platform:", R.version$platform, "\n")
cat("  R:", R.version.string, "\n")
cat("  SpiecEasi:", as.character(packageVersion("SpiecEasi")),
    " huge:", as.character(packageVersion("huge")),
    " pulsar:", as.character(packageVersion("pulsar")), "\n")
cat("  cores:", parallel::detectCores(), "\n")
cat(sprintf("== R pipeline (serial), reps=%d ==\n", reps))

results <- list()
for (cs in cases) {
  if (!grepl(id_filter, cs$id)) next
  X <- as.matrix(read.table(cs$counts))
  n <- nrow(X)

  # verify the pulsar subsample indices equal the golden ones (identical work)
  ratio <- pulsar:::.ratcheck(NULL, n)
  set.seed(cs$seed)
  ind <- do.call(rbind, replicate(cs$rep, sample(seq_len(n), floor(n * ratio), replace = FALSE),
                                  simplify = FALSE))
  golden_ind <- as.matrix(read.table(cs$subs))
  same_indices <- identical(unname(ind), unname(golden_ind) + 1L)

  times <- numeric(reps)
  edges <- NA_integer_
  opt_lambda <- NA_real_
  for (r in seq_len(reps)) {
    t0 <- Sys.time()
    fit <- suppressWarnings(suppressMessages(spiec.easi(
      X, method = "mb", sel.criterion = "stars",
      nlambda = cs$nlambda, lambda.min.ratio = cs$ratio,
      pulsar.params = list(rep.num = cs$rep, thresh = 0.05, seed = cs$seed),
      verbose = FALSE)))
    times[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    A <- as.matrix(fit$refit$stars); A[is.na(A)] <- 0
    edges <- as.integer(sum(A) / 2)
    grid <- as.numeric(fit$lambda)
    opt_lambda <- grid[[as.integer(fit$select$stars$opt.index)]]
    cat(sprintf("  %s rep %d/%d: %.1f s\n", cs$id, r, reps, times[r]))
  }
  results[[cs$id]] <- list(
    min = min(times), median = median(times), times = times,
    edge_count = edges, selected_lambda = opt_lambda,
    n_samples = n, n_taxa = ncol(X), indices_match_golden = same_indices
  )
  cat(sprintf("  => %s: min=%.1fs median=%.1fs edges=%d lambda=%.6g indices_match=%s\n",
              cs$id, min(times), median(times), edges, opt_lambda, same_indices))
}

if (!is.na(out_json)) {
  writeLines(jsonlite::toJSON(results, auto_unbox = TRUE, digits = 17, pretty = TRUE), out_json)
  cat("wrote ", out_json, "\n", sep = "")
}
