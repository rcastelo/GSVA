#!/usr/bin/env Rscript
## writes, from the dense HDF5 dataset written by prepare_calibration_data.R
## with the option --dense, the genes whose fraction of nonzero values in the
## first 20,000 cells is below <maxdensity>, for the first <ncells> cells, in
## another dense HDF5 dataset with the same chunks,
## TENxBrainDataMouseE18_processed_dense_sparsegenes.sce. its fraction of
## nonzero values, about 3% with genes below 6%, is much smaller than the
## 17.6% of the whole data, so that reading it takes more memory per nonzero
## value, through the dense form of its chunks, see 'densereadfactor' in
## README.md, which validate_maxmem_sweep.sh, or validate_maxmem_sweep_slurm.sh
## on a computer cluster, validates with <ondisk> set to "no", e.g.:
##   Rscript prepare_sparse_dense_data.R 0.06 300000
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 2)
    stop("usage: Rscript prepare_sparse_dense_data.R [maxdensity] [ncells]")
maxdensity <- if (length(args) > 0) as.numeric(args[1]) else 0.06
ncells <- if (length(args) > 1) as.integer(args[2]) else 300000L

suppressPackageStartupMessages({
    library(HDF5Array)
    library(SparseArray)
    library(SummarizedExperiment)
})

indir <- file.path("data", "TENxBrainDataMouseE18_processed_dense.sce")
outdir <- file.path("data", "TENxBrainDataMouseE18_processed_dense_sparsegenes.sce")
sce <- loadHDF5SummarizedExperiment(indir)
X <- assay(sce, "logcounts", withDimnames=FALSE)
ncells <- min(ncells, ncol(sce))

## fraction of nonzero values of each gene in the first 20,000 cells
x <- as(X[, seq_len(min(20000L, ncol(X)))], "SVT_SparseMatrix")
density <- rowSums(x != 0) / ncol(x)
genes <- which(density < maxdensity)
cat(sprintf("%d genes with fewer than %.1f%% of nonzero values, %.1f%% on average\n",
            length(genes), 100 * maxdensity, 100 * mean(density[genes])))

saveHDF5SummarizedExperiment(sce[genes, seq_len(ncells)], dir=outdir,
                             replace=TRUE, chunkdim=chunkdim(X),
                             verbose=TRUE)
cat("written", outdir, "\n")
