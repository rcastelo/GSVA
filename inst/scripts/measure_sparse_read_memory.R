#!/usr/bin/env Rscript
## reads into main memory, as sparse, one block of the first <ncols> columns
## of the processed data stored on disk in the layout <layout>, the CSC layout
## (TENxMatrix) or the dense HDF5 dataset, as GSVA's loader of sparse on-disk
## data does with each of its blocks, see .load_sparse_by_blocks(), printing
## the number of nonzero values of the block, the memory of the result and of
## its dense form, and the peak memory allocated by R while reading it, together with the times at
## which the reading starts and ends, from which
## measure_sparse_read_memory.sh takes the peak memory of the process. it
## calibrates 'sparsereadbytes' and 'densereadfactor', the memory of a block
## while it is read, see README.md. the dense layout needs the option --dense
## of prepare_calibration_data.R. run it through
## measure_sparse_read_memory.sh, from the calibration directory, see
## README.md:
##   Rscript measure_sparse_read_memory.R <csc|dense> <ncols>
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2)
    stop("usage: Rscript measure_sparse_read_memory.R <csc|dense> <ncols>")
layout <- match.arg(args[1], c("csc", "dense"))
ncols <- as.integer(args[2])

suppressPackageStartupMessages({
    library(GSVA)
    library(HDF5Array)
    library(SparseArray)
})
cat("GSVA", as.character(packageVersion("GSVA")), "from", find.package("GSVA"),
    "HDF5Array", as.character(packageVersion("HDF5Array")),
    "SparseArray", as.character(packageVersion("SparseArray")), "\n")

if (layout == "csc") {
    X <- TENxMatrix("data/TENxBrainDataMouseE18_processed_csc.sce/logcounts.h5",
                    group="matrix")
} else
    X <- HDF5Array("data/TENxBrainDataMouseE18_processed_dense.sce/assays.h5",
                   "assay001", as.sparse=TRUE)
if (ncols > ncol(X))
    stop(sprintf("the data has only %d columns", ncol(X)))
## the first block of a grid of blocks of 'ncols' columns, read as GSVA's
## loader reads each of its blocks
grid <- DelayedArray::colAutoGrid(X, ncol=ncols)
readblock <- function()
    as(DelayedArray::read_block(X, grid[[1L]], as.sparse=TRUE),
       "SVT_SparseMatrix")
if (exists(".read_sparse_block", envir=asNamespace("GSVA")))
    readblock <- function() GSVA:::.read_sparse_block(1L, X, grid)

## the peak memory allocated by R, in bytes, since the last reset
rmaxused <- function() {
    g <- gc()
    sum(g[, max(which(colnames(g) == "(Mb)"))]) * 1024^2
}
invisible(gc(reset=TRUE))
rstart <- rmaxused()

cat("PHASE start", format(as.numeric(Sys.time()), nsmall=3), "\n")
t <- system.time(res <- readblock())
cat("PHASE end", format(as.numeric(Sys.time()), nsmall=3), "\n")

cat("NNZ", format(nzcount(res), scientific=FALSE), "\n")
cat("RESULT", format(as.numeric(object.size(res)), scientific=FALSE), "\n")
## the size of the dense form of the block, which the reading of the dense
## layout goes through, chunk by chunk
cat("DENSE", format(as.numeric(nrow(res)) * ncol(res) *
                    (if (type(res) == "integer") 4 else 8),
                    scientific=FALSE), "\n")
cat("RPEAK", format(rmaxused() - rstart, scientific=FALSE), "\n")
cat("ELAPSED", t[["elapsed"]], "\n")
