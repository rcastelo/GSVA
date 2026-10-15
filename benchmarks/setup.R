# One-time data preparation for the benchmark
#
# Produces benchmarks/cache/se_full.rds: the preprocessed pbmc68k
# SingleCellExperiment.
#
# The dataset is downloaded on the first run using TENxPBMCData, and
# is later cached.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(TENxPBMCData)
  library(scuttle)
  library(GSVA)
  library(GSEABase)
})

cat("Data preparation for the benchmark\n")
cat("R version:", as.character(getRversion()), "\n")

cat("\n[1/3] Loading PBMC68K...\n")
sce <- TENxPBMCData(dataset = "pbmc68k")

cat("\n[2/3] Preprocessing...\n")
is_mito <- grepl("^MT-", rowData(sce)$Symbol_Tenx)
sce <- quickPerCellQC(sce, subsets = list(Mito = is_mito),
                      sub.fields = "subsets_Mito_percent")

cntxgene <- rowSums(assays(sce)$counts) + 1
sce <- sce[cntxgene >= 100, ]
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

# Load in-memory from HDF5
assays(sce) <- list(logcounts = as.matrix(assay(sce, "logcounts")))
gsvaAnnotation(sce) <- ENSEMBLIdentifier("org.Hs.eg.db")

n_cells <- ncol(sce)
n_genes <- nrow(sce)
cat(sprintf("  Expression matrix: %d genes x %d cells\n", n_genes, n_cells))

cat("[3/3] Saving...\n")

cache_dir <- "benchmarks/cache"
unlink(cache_dir, recursive = TRUE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, file.path(cache_dir, "se_full.rds"))
cat(sprintf("  %s/se_full.rds (%d genes x %d cells)\n",
            cache_dir, nrow(sce), ncol(sce)))

cat("\nDone. Run `Rscript benchmarks/benchmark.R <label>` to benchmark.\n")
