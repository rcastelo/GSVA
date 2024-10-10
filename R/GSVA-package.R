
#'
#' @useDynLib GSVA, .registration=TRUE
#'
#' @import methods
#'
#' @importMethodsFrom Biobase featureNames phenoData experimentData exprs
#' annotation "annotation<-"
#' @importMethodsFrom S4Vectors metadata "metadata<-"
#' @importMethodsFrom IRanges match
#' @importMethodsFrom SummarizedExperiment assays assayNames colData rowData
#' @importMethodsFrom SpatialExperiment imgData spatialCoords
#' @importMethodsFrom GSEABase geneIds incidence GeneSetCollection GeneSet
#' @importMethodsFrom BiocParallel bpiterate "bpworkers<-" bplapply "bpprogressbar<-"
#'
#' @importFrom graphics plot
#' @importFrom stats ecdf na.omit rnorm rpois sd
#' @importFrom utils installed.packages setTxtProgressBar txtProgressBar head tail
#' read.csv write.csv capture.output
#' @importFrom Matrix nnzero
#' @importFrom Biobase selectSome
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom IRanges CharacterList
#' @importFrom GSEABase AnnoOrEntrezIdentifier mapIdentifiers NullIdentifier
#' NullCollection geneIdType SymbolIdentifier ENSEMBLIdentifier EntrezIdentifier
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom parallel splitIndices
#' @importFrom BiocParallel SerialParam MulticoreParam multicoreWorkers bpnworkers
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom sparseMatrixStats colRanks rowRanges
#' @importFrom DelayedArray rowAutoGrid colAutoGrid defaultAutoGrid
#' currentBlockId read_block gridReduce write_block close t colSums
#' @importFrom HDF5Array HDF5RealizationSink writeHDF5Array
#' @importFrom BiocSingular runRandomSVD
NULL
