#!/usr/bin/env Rscript
## prepares the data on which the memory model of GSVA is calibrated, see
## README.md: the 1.3 million brain cells of the TENxBrainData package,
## fetched through ExperimentHub, filtered for low-quality cells and lowly
## expressed genes, and normalized into log-counts with the scrapper package,
## together with the M8 collection of mouse cell type signatures of MSigDB.
## run it from the calibration directory, where it writes, in 'data':
##   - TENxBrainDataMouseE18_processed_csc.sce, the log-counts in the
##     compressed sparse column (CSC) layout of 10x Genomics, which stores only
##     the nonzero values and is read by blocks of columns several times faster
##     than a dense dataset, and the rest of the object, without assays, in an
##     RDS file, loaded with:
##       sce <- readRDS(file.path(dir, "se.rds"))
##       assay(sce, "logcounts", withDimnames=FALSE) <-
##           TENxMatrix(file.path(dir, "logcounts.h5"), group="matrix")
##   - with the option --dense, TENxBrainDataMouseE18_processed_dense.sce, the
##     same log-counts in a dense HDF5 dataset written with
##     saveHDF5SummarizedExperiment(), which takes about two hours and 26 GB
##     of disk space, loaded with loadHDF5SummarizedExperiment()
##   - m8.all.v2025.1.Mm.symbols.gmt, the gene sets, downloaded from MSigDB
## usage:
##   Rscript prepare_calibration_data.R [--dense]
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1 || (length(args) == 1 && args[1] != "--dense"))
    stop("usage: Rscript prepare_calibration_data.R [--dense]")
savedense <- length(args) == 1

suppressPackageStartupMessages({
    library(ExperimentHub)
    library(HDF5Array)
    library(SingleCellExperiment)
    library(SparseArray)
    library(scrapper)
    library(TENxBrainData)
    library(GSEABase)
    library(GSVA)
})
dir.create("data", showWarnings=FALSE)

## the gene sets
gmtfile <- file.path("data", "m8.all.v2025.1.Mm.symbols.gmt")
if (!file.exists(gmtfile))
    download.file(paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/",
                         "release/2025.1.Mm/m8.all.v2025.1.Mm.symbols.gmt"),
                  gmtfile)

hub <- ExperimentHub()

## fetch row and column data
rdat <- query(hub, "TENxBrainData/1M_neurons_filtered_gene_bc_matrices_h5_rowData.rds")[[1]]
cdat <- query(hub, "TENxBrainData/1M_neurons_filtered_gene_bc_matrices_h5_colData.rds")[[1]]

## 1.3 Million Brain Cell dataset (sparse version)
h5sparsefname <- hub[["EH1039"]]
h5array <- TENxMatrix(h5sparsefname, group="mm10")
Sarray <- SparseArray(h5array) ## load it into memory as a SparseArray object for speed

stopifnot(identical(rownames(Sarray), rdat$Ensembl)) ## QC
stopifnot(identical(colnames(Sarray), rownames(cdat))) ## QC

sce <- SingleCellExperiment(assays=list(counts=Sarray), rowData=rdat, colData=cdat)

## discard low-quality cells
is_mito <- grepl("^mt-", rowData(sce)$Symbol)
stopifnot(sum(is_mito) == 13) ## QC
sce <- quickRnaQc.se(sce, subsets=list(mito=is_mito))
sce <- sce[, sce$keep]

## discard genes expressed in less than 1% of the cells
cellsxgene <- rowSums(counts(sce) > 0)
sce <- sce[cellsxgene >= floor(ncol(sce) * 0.01), ]

## normalize data
sce <- normalizeRnaCounts.se(sce, size.factors=sce$sum)

## add gene identifier metadata
gsvaAnnotation(sce) <- ENSEMBLIdentifier("org.Mm.eg.db")

## remove raw counts to save space
assays(sce)$counts <- NULL

## save the normalized values, which normalizeRnaCounts.se() gives as a
## delayed operation on the counts, in the CSC layout, and the rest of the
## object in an RDS file. the object stores no path to the HDF5 file, so the
## directory can be moved or copied
dir <- file.path("data", "TENxBrainDataMouseE18_processed_csc.sce")
dir.create(dir, showWarnings=FALSE)
system.time(writeTENxMatrix(logcounts(sce), file.path(dir, "logcounts.h5"),
                            group="matrix"))
assays(sce) <- list()
saveRDS(sce, file.path(dir, "se.rds"))

## the same normalized values in a dense dataset, which is read by blocks of
## columns more slowly, unless the blocks span whole chunks of the dataset
if (savedense) {
    assay(sce, "logcounts", withDimnames=FALSE) <-
        TENxMatrix(file.path(dir, "logcounts.h5"), group="matrix")
    saveHDF5SummarizedExperiment(sce,
        dir=file.path("data", "TENxBrainDataMouseE18_processed_dense.sce"),
        replace=TRUE)
}
