
#'
#' @useDynLib GSVA, .registration=TRUE
#'
#' @import methods
#'
#' @importMethodsFrom Biobase featureNames phenoData experimentData exprs annotation
#' @importMethodsFrom S4Vectors metadata "metadata<-"
#' @importMethodsFrom IRanges match
#' @importMethodsFrom SummarizedExperiment assays assayNames colData
#' @importMethodsFrom GSEABase geneIds incidence
#' @importMethodsFrom BiocParallel bpiterate "bpworkers<-" bplapply "bpprogressbar<-"
#'
#' @importFrom graphics plot
#' @importFrom stats ecdf na.omit rnorm rpois sd
#' @importFrom utils installed.packages setTxtProgressBar txtProgressBar
#' read.csv write.csv capture.output
#' @importFrom S4Vectors SimpleList
#' @importFrom IRanges CharacterList
#' @importFrom  GSEABase AnnoOrEntrezIdentifier mapIdentifiers getGmt
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom parallel splitIndices
#' @importFrom BiocParallel SerialParam MulticoreParam multicoreWorkers bpnworkers
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom sparseMatrixStats colRanks
#' @importFrom DelayedArray rowAutoGrid colAutoGrid defaultAutoGrid
#' currentBlockId read_block gridReduce write_block close t colSums
#' @importFrom HDF5Array HDF5RealizationSink writeHDF5Array
#' @importFrom DelayedMatrixStats rowSds
#' @importFrom BiocSingular runRandomSVD
NULL
