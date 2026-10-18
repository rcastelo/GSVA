# Benchmarks the scoring stage (gsvaColScores) of the currently installed
# version of the GSVA package.
#
# Usage: Rscript benchmarks/benchmark.R <label> [--quick]
#   --quick: fast smoke test; benchmarks only the T1 (Hallmark) tier on the
#   first 1000 cells with 1-2 iterations. Timings are not representative.

suppressPackageStartupMessages({
  library(bench)
  library(GSEABase)
  library(GSVA)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(SingleCellExperiment)
})

cache_dir <- "benchmarks/cache"

main <- function() {
  parsed <- parse_args()
  label <- parsed$label
  quick <- parsed$quick

  cat("GSVA scoring benchmark\n")
  cat(sprintf("  state: %s\n", label))
  cat(sprintf("  GSVA version: %s | R version: %s\n",
              as.character(packageVersion("GSVA")),
              as.character(getRversion())))
  cat(sprintf("  mode: %s\n", if (quick) "quick" else "full"))

  se_path <- file.path(cache_dir, "se_full.rds")
  if (!file.exists(se_path)) {
    stop("Missing benchmarks/cache/se_full.rds. Run `Rscript benchmarks/setup.R` first.")
  }

  se_full <- readRDS(se_path)
  if (quick) {
    se_full <- se_full[, seq_len(min(1000, ncol(se_full)))]
    cat("  quick mode: subsampled to the first 1000 cells\n")
  }
  cat(sprintf("\nLoaded se_full: %d genes x %d cells\n", nrow(se_full), ncol(se_full)))

  ranks <- compute_native_ranks(se_full)
  rm(se_full)

  all_tiers <- benchmark_tiers(quick)
  cat(sprintf("\n%d tiers to benchmark\n", length(all_tiers)))

  outfile <- file.path("benchmarks", paste0(label, "_results.csv"))
  results <- load_prev_results(outfile)

  for (tier in all_tiers) {
    if (tier$id %in% results$tier) {
      cat(sprintf("\n── %s: %s ── (already done, skipping)\n", tier$id, tier$name))
      next
    }
    cat(sprintf("\n── %s: %s ──\n", tier$id, tier$name))
    results <- rbind(results, run_tier(tier, ranks))
    write.csv(results, file = outfile, row.names = FALSE)
    gc(verbose = FALSE)
  }

  print_summary(label, results)
  write.csv(results, file = outfile, row.names = FALSE)
  cat(sprintf("\nResults written to %s\n", outfile))
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  values <- args[!args %in% "--quick"]
  if (length(values) < 1) {
    stop("usage: Rscript benchmark.R <label> [--quick]")
  }
  list(label = values[1], quick = any(args == "--quick"))
}

compute_native_ranks <- function(se_full) {
  cat("\nComputing version-native ranks (untimed) ...\n")
  h <- msigdbr(species = "human", collection = "H")
  gsets_h <- split(h$gene_symbol, h$gs_name)
  par <- gsvaParam(se_full, gsets_h, verbose = FALSE)
  norm <- gsvaRowNorm(par, verbose = FALSE)
  rm(par)
  gc(verbose = FALSE)

  ranks <- gsvaColRanks(norm, verbose = FALSE)
  rm(norm)
  gc(verbose = FALSE)

  cat(sprintf("  ranks: %d genes x %d cells\n", nrow(ranks), ncol(ranks)))
  ranks
}

benchmark_tiers <- function(quick = FALSE) {
  if (quick) {
    return(list(
      list(id = "T1", name = "Hallmark", code = "H", max_gsets = NULL,
           min_iter = 1, max_iter = 2, extra_codes = character())
    ))
  }
  list(
    list(id = "T1", name = "Hallmark",                 code = "H",  max_gsets = NULL, min_iter = 3, max_iter = 5, extra_codes = character()),
    list(id = "T2", name = "C7 (top 1000)",            code = "C7", max_gsets = 1000, min_iter = 3, max_iter = 5, extra_codes = character()),
    list(id = "T3", name = "C2 (top 1000)",            code = "C2", max_gsets = 1000, min_iter = 3, max_iter = 5, extra_codes = character()),
    list(id = "T4", name = "C5:BP (top 2000)",         code = "C5", max_gsets = 2000, min_iter = 2, max_iter = 5, extra_codes = character()),
    list(id = "T5", name = "Hallmark + C5 (top 2000)", code = "C5", max_gsets = 2000, min_iter = 2, max_iter = 5, extra_codes = "H")
  )
}

load_prev_results <- function(outfile) {
  if (!file.exists(outfile)) {
    return(new_results_df())
  }

  prev <- read.csv(outfile)
  need <- c("tier", "name", "n_gsets", "n_cells", "median_ms",
            "n_itr", "ms_per_gset", "us_per_cell_gset", "throughput_kops_s")
  if (nrow(prev) > 0 && all(need %in% names(prev))) {
    cat(sprintf("  Resuming: skipping %d already-done tier(s): %s\n",
                nrow(prev), paste(prev$tier, collapse = ", ")))
    return(prev)
  }

  new_results_df()
}

new_results_df <- function() {
  data.frame(
    tier              = character(),
    name              = character(),
    n_gsets           = integer(),
    n_cells           = integer(),
    median_ms         = numeric(),
    n_itr             = integer(),
    ms_per_gset       = numeric(),
    us_per_cell_gset  = numeric(),
    throughput_kops_s = numeric()
  )
}

run_tier <- function(tier, ranks) {
  n_cells <- ncol(ranks)

  mapped_gsets <- build_mapped_gene_sets(
    tier$code, tier$max_gsets,
    row_names = rownames(ranks),
    extra_codes = tier$extra_codes
  )
  n_gsets <- length(mapped_gsets)
  if (n_gsets == 0) {
    stop(sprintf("  %s has 0 gene sets after mapping", tier$name))
  }
  cat(sprintf("  Benchmarking %d gene sets x %d cells ...\n", n_gsets, n_cells))

  # Swap gene sets
  ranks_tier <- ranks
  mdata <- metadata(ranks_tier)
  mdata$gsvaParam$geneSets <- mapped_gsets
  mdata$gsvaParam$annotation <- NullIdentifier()
  metadata(ranks_tier) <- mdata

  # Warmup (untimed)
  invisible(gsvaColScores(ranks_tier, verbose = FALSE))

  # Timed benchmark
  bnc <- tryCatch({
    bench::mark(
      cpu = gsvaColScores(ranks_tier, verbose = FALSE),
      min_iterations = tier$min_iter,
      max_iterations = tier$max_iter,
      memory         = FALSE,
      check          = NULL,
      filter_gc      = FALSE
    )
  }, error = function(e) {
    stop(sprintf("  ERROR: %s", conditionMessage(e)), call. = FALSE)
  })

  r <- bnc[1, ]
  median_ms <- round(as.numeric(r$median * 1e3), 1)
  n_its     <- r$n_itr
  ms_per_gset      <- round(median_ms / n_gsets, 3)
  us_per_cell_gset <- round((median_ms * 1e3) / (n_gsets * n_cells), 3)
  throughput       <- round((n_gsets * n_cells) / median_ms, 1)

  cat(sprintf("  median: %.1f ms | ms/gset: %.2f | us/(cell*gset): %.2f | throughput: %.1f kops/s | %d iterations\n",
              median_ms, ms_per_gset, us_per_cell_gset, throughput, n_its))

  data.frame(
    tier              = tier$id,
    name              = tier$name,
    n_gsets           = n_gsets,
    n_cells           = n_cells,
    median_ms         = median_ms,
    n_itr             = n_its,
    ms_per_gset       = ms_per_gset,
    us_per_cell_gset  = us_per_cell_gset,
    throughput_kops_s = throughput
  )
}

build_mapped_gene_sets <- function(coll_code, max_gsets = NULL,
                                   min_genes = 5, row_names,
                                   extra_codes = character()) {
  cat(sprintf("  Loading %s (species=human)...", coll_code))
  msig_df <- msigdbr(species = "human", collection = coll_code)

  if (length(extra_codes) > 0) {
    for (extra in extra_codes) {
      cat(sprintf("+%s", extra))
      extra_df <- msigdbr(species = "human", collection = extra)
      extra_df <- extra_df[!extra_df$gs_name %in% msig_df$gs_name, ]
      msig_df <- rbind(msig_df, extra_df)
    }
  }

  if (!is.null(max_gsets)) {
    unique_names <- unique(msig_df$gs_name)
    if (length(unique_names) > max_gsets) {
      unique_names <- head(unique_names, max_gsets)
      msig_df <- msig_df[msig_df$gs_name %in% unique_names, ]
    }
  }

  gsets_list <- split(msig_df$gene_symbol, msig_df$gs_name)
  n_raw <- length(gsets_list)
  all_symbols <- unique(unlist(lapply(gsets_list, identity)))
  lookup <- suppressMessages(mapIds(org.Hs.eg.db, keys = all_symbols, column = "ENSEMBL",
                                    keytype = "SYMBOL", multiVals = "first"))

  # ENSEMBL -> row index
  ens_to_idx <- setNames(seq_along(row_names), row_names)
  gsets_idx <- lapply(gsets_list, function(symbols) {
    ens <- lookup[symbols]
    idx <- ens_to_idx[unique(ens[!is.na(ens)])]
    idx[!is.na(idx)]
  })
  names(gsets_idx) <- names(gsets_list)

  # Drop overly small gene sets
  gsets_idx <- Filter(function(gs) length(gs) >= min_genes, gsets_idx)
  cat(sprintf(" %d raw -> %d kept (>= %d genes)\n",
              n_raw, length(gsets_idx), min_genes))
  gset_sizes <- lengths(gsets_idx)
  cat(sprintf("  Sizes: min=%d, median=%g, max=%d\n",
              min(gset_sizes), median(gset_sizes), max(gset_sizes)))

  gsets_idx
}

print_summary <- function(label, results) {
  cat("\n")
  cat("================================================================\n")
  cat(sprintf("  GSVA scoring benchmark | state: %s | GSVA %s\n",
              label, as.character(packageVersion("GSVA"))))
  cat(sprintf("  %d tiers completed\n", nrow(results)))
  cat("  Timing covers the scoring stage (gsvaColScores) only.\n")
  cat("================================================================\n")

  for (i in seq_len(nrow(results))) {
    row <- results[i, ]
    tname <- row$name
    if (nchar(tname) > 40) tname <- paste0(substr(tname, 1, 37), "...")
    cat(sprintf("  %-4s  %-40s  %4d gsets  %5d cells\n",
                row$tier, tname, row$n_gsets, row$n_cells))
    cat(sprintf("    median=%8.1fms  ms/gset=%7.2f  us/(cell*gset)=%7.2f  throughput=%8.1f kops/s  %d itr\n",
                row$median_ms, row$ms_per_gset,
                row$us_per_cell_gset, row$throughput_kops_s, row$n_itr))
  }
  cat("================================================================\n")
}

main()
