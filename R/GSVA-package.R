
#'
#' @useDynLib GSVA, .registration=TRUE
#'
#' @import methods
#'
#' @importMethodsFrom Biobase featureNames phenoData experimentData exprs
#' @importMethodsFrom S4Vectors metadata "metadata<-"
#' @importMethodsFrom SummarizedExperiment assays assayNames colData rowData
#' @importMethodsFrom SpatialExperiment imgData spatialCoords
#' @importMethodsFrom GSEABase geneIds incidence GeneSetCollection GeneSet
#' @importMethodsFrom BiocParallel bpiterate "bpworkers<-" bplapply "bpprogressbar<-"
#'
#' @importFrom stats ecdf na.omit rnorm rpois sd
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom BiocParallel SerialParam MulticoreParam multicoreWorkers bpnworkers
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom DelayedArray rowAutoGrid colAutoGrid defaultAutoGrid
#' currentBlockId read_block gridReduce write_block close t colSums
#' @importFrom HDF5Array HDF5RealizationSink writeHDF5Array
#' @importFrom cli cli_abort cli_alert_info cli_alert_warning cli_alert_success
NULL
