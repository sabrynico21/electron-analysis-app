# Deterministic source-data preprocessing shared by the three R golden
# generators.  Matrices returned here have samples in rows and taxa in columns.
# This script records no personal/sample metadata beyond anonymized row labels.

.evenly_spaced_indices <- function(size, total) {
  if (size < 1L || size > total) stop("requested subset size is out of range")
  as.integer(round(seq(1, total, length.out = size)))
}

.normalize_abundance_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(rownames(x))) rownames(x) <- sprintf("sample-%06d", seq_len(nrow(x)))
  if (is.null(colnames(x))) colnames(x) <- sprintf("taxon-%06d", seq_len(ncol(x)))
  raw_dimensions <- dim(x)
  x <- x[order(rownames(x)), order(colnames(x)), drop = FALSE]
  x <- x[rowSums(x) > 0, colSums(x) > 0, drop = FALSE]
  if (any(!is.finite(x)) || any(x < 0)) stop("abundances must be finite and non-negative")
  attr(x, "raw_dimensions") <- raw_dimensions
  x
}

.public_case <- function(name, dataset, source, x, n, p,
                         taxa_strategy = c("spread", "prevalent"),
                         source_package = NA_character_,
                         source_package_version = NA_character_,
                         source_license = NA_character_,
                         source_citation = NA_character_,
                         retrieval_date = "2026-10-20",
                         redistribution_cleared = FALSE) {
  taxa_strategy <- match.arg(taxa_strategy)
  sample_idx <- .evenly_spaced_indices(n, nrow(x))
  samples <- x[sample_idx, , drop = FALSE]
  # Select a reproducible spread of taxa across observed prevalence, excluding
  # constant features in the selected sample set. Ties are broken by taxon ID.
  varying <- which(apply(samples, 2, function(v) length(unique(v)) > 1L))
  if (length(varying) < p) {
    stop(sprintf("%s has only %d varying taxa; requested %d", name, length(varying), p))
  }
  varying_names <- colnames(samples)[varying]
  prevalence <- colSums(samples[, varying, drop = FALSE] > 0)
  abundance <- colSums(samples[, varying, drop = FALSE])
  if (taxa_strategy == "prevalent") {
    prevalence_order <- order(-prevalence, -abundance, varying_names)
    chosen <- varying[prevalence_order[seq_len(p)]]
  } else {
    prevalence_order <- order(prevalence, abundance, varying_names)
    # Preserve a deterministic spread across observed prevalence ranks.
    chosen <- varying[prevalence_order[.evenly_spaced_indices(p, length(varying))]]
  }
  counts <- samples[, chosen, drop = FALSE]
  list(
    name = name,
    dataset = dataset,
    source = source,
    source_package = source_package,
    source_package_version = source_package_version,
    source_license = source_license,
    source_citation = source_citation,
    retrieval_date = retrieval_date,
    redistribution_cleared = redistribution_cleared,
    taxa_strategy = taxa_strategy,
    input = counts,
    raw_n = unname(attr(x, "raw_dimensions")[1]),
    raw_p = unname(attr(x, "raw_dimensions")[2]),
    available_n = nrow(x),
    available_p = ncol(x),
    n = nrow(counts),
    p = ncol(counts),
    zero_fraction = mean(counts == 0),
    library_size_min = min(rowSums(counts)),
    library_size_median = stats::median(rowSums(counts)),
    library_size_max = max(rowSums(counts))
  )
}

load_public_microbiome_cases <- function(globalpatterns_taxa_strategy = "prevalent") {
  if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
    stop("SpiecEasi R data package is required for the American Gut cases")
  }
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq is required for the GlobalPatterns cases")
  }

  # Public American Gut Project count matrix, shipped with SpiecEasi. Load from
  # the installed package namespace before the SpiecEasi reference source is
  # loaded by that package's fixture generator.
  agp_env <- new.env(parent = emptyenv())
  utils::data("amgut1.filt", package = "SpiecEasi", envir = agp_env)
  if (!exists("amgut1.filt", envir = agp_env, inherits = FALSE)) {
    stop("SpiecEasi package does not expose amgut1.filt")
  }
  agp <- .normalize_abundance_matrix(get("amgut1.filt", envir = agp_env))

  # Public 16S OTU counts from Caporaso et al. (2011), distributed with
  # phyloseq::GlobalPatterns. Drop metadata/phylogeny and use the OTU table only.
  global_env <- new.env(parent = emptyenv())
  utils::data("GlobalPatterns", package = "phyloseq", envir = global_env)
  if (!exists("GlobalPatterns", envir = global_env, inherits = FALSE)) {
    stop("phyloseq package does not expose GlobalPatterns")
  }
  gp <- get("GlobalPatterns", envir = global_env)
  gp_otu <- methods::as(phyloseq::otu_table(gp), "matrix")
  if (phyloseq::taxa_are_rows(gp)) gp_otu <- t(gp_otu)
  global_patterns <- .normalize_abundance_matrix(gp_otu)

  spiec_version <- as.character(utils::packageVersion("SpiecEasi"))
  phyloseq_version <- as.character(utils::packageVersion("phyloseq"))
  agp_citation <- paste0(
    "SpiecEasi::amgut1.filt (American Gut Project Round 1); ",
    "documented source http://humanfoodproject.com/americangut/"
  )
  gp_citation <- paste0(
    "Caporaso et al. (2011) PNAS 108(Suppl 1):4516-4522, ",
    "PMCID PMC3063599; OTU table supplied for phyloseq::GlobalPatterns"
  )

  list(
    .public_case(
      "amgut_r1_n80_p20", "American Gut Project Round 1",
      "SpiecEasi::amgut1.filt", agp, n = 80L, p = 20L,
      source_package = "SpiecEasi",
      source_package_version = spiec_version,
      source_license = "GPL (>= 3)",
      source_citation = agp_citation,
      redistribution_cleared = FALSE
    ),
    .public_case(
      "amgut_r1_n20_p30", "American Gut Project Round 1",
      "SpiecEasi::amgut1.filt", agp, n = 20L, p = 30L,
      source_package = "SpiecEasi",
      source_package_version = spiec_version,
      source_license = "GPL (>= 3)",
      source_citation = agp_citation,
      redistribution_cleared = FALSE
    ),
    .public_case(
      "globalpatterns_n26_p15", "GlobalPatterns / Caporaso et al. 2011",
      "phyloseq::GlobalPatterns", global_patterns, n = nrow(global_patterns),
      p = 15L, taxa_strategy = globalpatterns_taxa_strategy,
      source_package = "phyloseq",
      source_package_version = phyloseq_version,
      source_license = "AGPL-3",
      source_citation = gp_citation,
      redistribution_cleared = FALSE
    ),
    .public_case(
      "globalpatterns_n26_p40", "GlobalPatterns / Caporaso et al. 2011",
      "phyloseq::GlobalPatterns", global_patterns, n = nrow(global_patterns),
      p = 40L, taxa_strategy = globalpatterns_taxa_strategy,
      source_package = "phyloseq",
      source_package_version = phyloseq_version,
      source_license = "AGPL-3",
      source_citation = gp_citation,
      redistribution_cleared = FALSE
    )
  )
}
