#!/usr/bin/env Rscript
# =============================================================================
# generate_golden_extended.R — Phase B extended golden generator.
#
# Same reference stack and conventions as tests/r_reference/generate_golden.R:
#   SpiecEasi 1.99.0 / huge 2.0.1 / pulsar 0.3.13, 17-digit output, pulsar
#   subsample indices exported WITHOUT the full model, intermediate steps saved.
#
# Output: tests/golden_extended/<case_id>/* with a top-level manifest.json.
#
# Usage:
#   Rscript tests/r_reference/generate_golden_extended.R [outdir] [id_regex]
#   (id_regex filters case ids, e.g. "^rand_p10"; default = all)
# =============================================================================

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(pulsar)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1] else "tests/golden_extended"
id_filter <- if (length(args) >= 2L) args[2] else ".*"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Writers
## ---------------------------------------------------------------------------
write_num_matrix <- function(m, path) {
  m <- as.matrix(m)
  txt <- matrix(sprintf("%.17g", as.numeric(m)), nrow = nrow(m), ncol = ncol(m))
  write.table(txt, path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
write_num_vec <- function(v, path) writeLines(sprintf("%.17g", as.numeric(v)), path)
write_int_matrix <- function(m, path) {
  write.table(m, path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

## ---------------------------------------------------------------------------
## Deterministic synthetic data
## ---------------------------------------------------------------------------
make_counts <- function(n, p, depth, seed, zero_frac) {
  set.seed(seed)
  k <- max(1L, p %/% 4L)
  W <- matrix(abs(rnorm(p * k, sd = 1)), p, k)
  Z <- matrix(rnorm(n * k), n, k)
  mu <- exp(as.numeric(scale(Z %*% t(W))))
  mu <- matrix(pmax(mu, 1e-6), n, p)
  counts <- matrix(0L, n, p)
  for (i in seq_len(n)) {
    counts[i, ] <- as.integer(rmultinom(1L, depth, mu[i, ] / sum(mu[i, ])))
  }
  if (zero_frac > 0) {
    mask <- matrix(runif(n * p) < zero_frac, n, p)
    counts[mask] <- 0L
  }
  # Guarantee: no constant column (huge rejects it) and no all-zero row.
  for (j in seq_len(p)) {
    if (length(unique(counts[, j])) == 1L) counts[1L, j] <- counts[1L, j] + 1L
  }
  for (i in seq_len(n)) {
    if (sum(counts[i, ]) == 0L) counts[i, 1L] <- counts[i, 1L] + 1L
  }
  counts
}

## ---------------------------------------------------------------------------
## Case definitions (pure-R config unless config=="app")
## ---------------------------------------------------------------------------
configs <- list(
  r   = list(nlambda = 20L, ratio = 1e-3, rep.num = 20L),
  app = list(nlambda = 25L, ratio = 1e-2, rep.num = 20L)
)
threshold <- 0.05
stars_seed <- 17L

mk <- function(id, p, n, zeros, depth, seed, cfg = "r", desc) {
  list(id = id, p = as.integer(p), n = as.integer(n), target_zero = zeros,
       depth = as.integer(depth), seed = as.integer(seed), cfg = cfg, desc = desc)
}

cases <- list(
  # --- 12 random datasets: p in {10,20,30,50}, n in {30,50,100} ---
  mk("rand_p10_n30_s201",  10,  30, 0.15, 5000,  201, "r", "random p=10 n=30"),
  mk("rand_p10_n50_s202",  10,  50, 0.15, 5000,  202, "r", "random p=10 n=50"),
  mk("rand_p10_n100_s203", 10, 100, 0.15, 5000,  203, "r", "random p=10 n=100"),
  mk("rand_p20_n30_s204",  20,  30, 0.15, 5000,  204, "r", "random p=20 n=30"),
  mk("rand_p20_n50_s205",  20,  50, 0.15, 5000,  205, "r", "random p=20 n=50"),
  mk("rand_p20_n100_s206", 20, 100, 0.15, 5000,  206, "r", "random p=20 n=100"),
  mk("rand_p30_n30_s207",  30,  30, 0.15, 5000,  207, "r", "random p=30 n=30"),
  mk("rand_p30_n50_s208",  30,  50, 0.15, 5000,  208, "r", "random p=30 n=50"),
  mk("rand_p30_n100_s209", 30, 100, 0.15, 5000,  209, "r", "random p=30 n=100"),
  mk("rand_p50_n30_s210",  50,  30, 0.15, 8000,  210, "r", "random p=50 n=30"),
  mk("rand_p50_n50_s211",  50,  50, 0.15, 8000,  211, "r", "random p=50 n=50"),
  mk("rand_p50_n100_s212", 50, 100, 0.15, 8000,  212, "r", "random p=50 n=100"),
  # --- very sparse: zero fraction > 70% ---
  mk("sparse_p30_n100_s213", 30, 100, 0.75, 50000, 213, "r", "sparse 75% zeros p=30 n=100"),
  mk("sparse_p50_n60_s214",  50,  60, 0.75, 50000, 214, "r", "sparse 75% zeros p=50 n=60"),
  # --- large p ---
  mk("largep_p100_n80_s215",  100,  80, 0.10, 20000, 215, "r", "p>n p=100 n=80"),
  mk("largep_p100_n200_s216", 100, 200, 0.10, 20000, 216, "r", "large p=100 n=200"),
  # --- n > 144 (exercises .ratcheck: 10*sqrt(n)/n) ---
  mk("ratcheck_p30_n200_s217", 30, 200, 0.10, 20000, 217, "r", "n>144 p=30 n=200"),
  mk("ratcheck_p50_n300_s218", 50, 300, 0.10, 20000, 218, "r", "n>144 p=50 n=300"),
  # --- p>n with zeros ---
  mk("pn_zeros_p40_n25_s219", 40, 25, 0.30, 2000, 219, "r", "p>n with zeros p=40 n=25"),
  # --- app config subset ---
  mk("app_p10_n30_s201",   10,  30, 0.15, 5000,  201, "app", "random p=10 n=30 (app cfg)"),
  mk("app_p30_n100_s209",  30, 100, 0.15, 5000,  209, "app", "random p=30 n=100 (app cfg)"),
  mk("app_largep_p100_n80_s215", 100, 80, 0.10, 20000, 215, "app", "p>n p=100 n=80 (app cfg)"),
  mk("app_ratcheck_p30_n200_s217", 30, 200, 0.10, 20000, 217, "app", "n>144 p=30 n=200 (app cfg)")
)

sel <- cases[vapply(cases, function(c) grepl(id_filter, c$id), logical(1))]
message("Selected ", length(sel), " case(s).")

manifest <- list()

for (cs in sel) {
  t0 <- Sys.time()
  cfg <- configs[[cs$cfg]]
  message("== ", cs$id, " (p=", cs$p, ", n=", cs$n, ", zeros~", cs$target_zero,
          ", cfg=", cs$cfg, ") ==")
  cdir <- file.path(outdir, cs$id)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)

  X <- make_counts(cs$n, cs$p, cs$depth, cs$seed, cs$target_zero)
  rownames(X) <- as.character(seq_len(nrow(X)))
  actual_zero <- mean(X == 0L)
  X_clr <- SpiecEasi:::.spiec.easi.norm(X)

  lambda_max <- pulsar::getMaxCov(cor(X_clr))
  grid <- pulsar::getLamPath(lambda_max, lambda_max * cfg$ratio, cfg$nlambda, log = TRUE)

  # pulsar subsample indices (exclude full model)
  ratio <- pulsar:::.ratcheck(NULL, n <- nrow(X))
  set.seed(stars_seed)
  ind <- do.call(rbind, replicate(cfg$rep.num,
    sample(c(1L:n), floor(n * ratio), replace = FALSE), simplify = FALSE))

  # full pipeline
  fit <- suppressWarnings(suppressMessages(spiec.easi(
    X, method = "mb", sel.criterion = "stars",
    nlambda = cfg$nlambda, lambda.min.ratio = cfg$ratio,
    pulsar.params = list(rep.num = cfg$rep.num, thresh = threshold, seed = stars_seed),
    verbose = FALSE)))
  stars <- fit$select$stars
  grid_fit <- as.numeric(fit$lambda)
  opt_index <- as.integer(stars$opt.index)
  monotone <- as.numeric(stars$summary)
  raw <- vapply(stars$merge, function(P) {
    P <- as.matrix(P); 4 * sum(P * (1 - P)) / (cs$p * (cs$p - 1))
  }, numeric(1))
  final_adj <- as.matrix(fit$refit$stars); final_adj[is.na(final_adj)] <- 0

  # solver isolation: one warm-started huge solve per subsample over the full path
  d <- ncol(X_clr); L <- length(grid_fit); nsub <- nrow(ind)
  rows <- vector("list", nsub * L); k <- 0L
  for (s in seq_len(nsub)) {
    Xs <- X_clr[ind[s, ], , drop = FALSE]
    res <- suppressWarnings(suppressMessages(
      SpiecEasi::sparseiCov(Xs, method = "mb", lambda = grid_fit, verbose = FALSE)))
    for (li in seq_len(L)) {
      k <- k + 1L
      rows[[k]] <- as.numeric(t(as.matrix(res$path[[li]])))  # row-major 0/1
    }
  }
  # coefficient path for subsample 1 (node x predictor)
  res1 <- suppressWarnings(suppressMessages(
    SpiecEasi::sparseiCov(X_clr[ind[1, ], , drop = FALSE], method = "mb",
                          lambda = grid_fit, verbose = FALSE)))
  crow <- vector("list", L)
  for (li in seq_len(L)) {
    m <- t(as.matrix(res1$beta[[li]]))
    crow[[li]] <- as.numeric(t(m))
  }

  ## write outputs
  write_int_matrix(X, file.path(cdir, "counts.tsv"))
  write_num_matrix(X_clr, file.path(cdir, "clr.tsv"))
  writeLines(sprintf("%.17g", lambda_max), file.path(cdir, "lambda_max.txt"))
  write_num_vec(grid_fit, file.path(cdir, "lambda_path.txt"))
  write_int_matrix(ind - 1L, file.path(cdir, "subsample_indices.tsv"))
  write_num_vec(raw, file.path(cdir, "stars_raw.txt"))
  write_num_vec(monotone, file.path(cdir, "stars_monotone.txt"))
  writeLines(as.character(opt_index), file.path(cdir, "opt_index.txt"))
  writeLines(sprintf("%.17g", grid_fit[[opt_index]]), file.path(cdir, "opt_lambda.txt"))
  write_int_matrix(final_adj, file.path(cdir, "final_adj.tsv"))
  write.table(t(sapply(rows, as.numeric)), file.path(cdir, "solver_adj_path.tsv"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(t(sapply(crow, as.numeric)), file.path(cdir, "solver_coef_path.tsv"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  writeLines(c(paste0("d=", d), paste0("n_lambda=", L), paste0("n_sub=", nsub),
               paste0("p=", cs$p), paste0("n=", cs$n)),
             file.path(cdir, "solver_dims.txt"))

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  manifest[[cs$id]] <- list(
    id = cs$id, p = cs$p, n = cs$n, seed = cs$seed, config = cs$cfg,
    nlambda = cfg$nlambda, lambda_min_ratio = cfg$ratio, rep_num = cfg$rep.num,
    target_zero_fraction = cs$target_zero, actual_zero_fraction = round(actual_zero, 4),
    description = cs$desc, gen_seconds = round(dt, 1)
  )
  message(sprintf("   done in %.1fs (actual zeros %.3f)", dt, actual_zero))
}

## ---------------------------------------------------------------------------
## manifest.json (written manually to avoid a jsonlite dependency)
## ---------------------------------------------------------------------------
json_escape <- function(s) gsub('"', '\\"', s, fixed = TRUE)
entries <- lapply(manifest, function(m) {
  sprintf(paste0('    {"id": "%s", "p": %d, "n": %d, "seed": %d, "config": "%s", ',
                 '"nlambda": %d, "lambda_min_ratio": %s, "rep_num": %d, ',
                 '"target_zero_fraction": %s, "actual_zero_fraction": %s, ',
                 '"gen_seconds": %s, "description": "%s"}'),
          json_escape(m$id), m$p, m$n, m$seed, m$config, m$nlambda,
          format(m$lambda_min_ratio), m$rep_num, format(m$target_zero_fraction),
          format(m$actual_zero_fraction), format(m$gen_seconds),
          json_escape(m$description))
})
writeLines(c("[", paste(unlist(entries), collapse = ",\n"), "]"),
           file.path(outdir, "manifest.json"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
writeLines(c(
  paste0("SpiecEasi=", as.character(packageVersion("SpiecEasi"))),
  paste0("huge=", as.character(packageVersion("huge"))),
  paste0("pulsar=", as.character(packageVersion("pulsar"))),
  paste0("glmnet=", as.character(packageVersion("glmnet"))),
  paste0("R=", as.character(getRversion()))
), file.path(outdir, "versions.txt"))

message("Extended golden written to: ", normalizePath(outdir))
