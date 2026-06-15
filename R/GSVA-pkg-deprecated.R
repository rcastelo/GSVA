## GSVA-pkg-deprecated.R
#' @title Deprecated functions in package `GSVA`.
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned.
#' @name GSVA-pkg-deprecated
#' @keywords internal
NULL

#' @description The `gsvaRanks()` method is deprecated. Please use `gsvaRowNorm()`
#' and `gsvaColRanks()` instead.
#'
#' @aliases gsvaRanks,gsvaParam-method
#' @name gsvaRanks
#' @rdname GSVA-pkg-deprecated
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @exportMethod gsvaRanks
setMethod("gsvaRanks", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {
              
              .Deprecated(new="gsvaRowNorm() and gsvaColRanks()",
                          package="GSVA",
                          msg=paste("The 'gsvaRanks()' method is deprecated.",
                                    "Please use 'gsvaRowNorm()' and",
                                    "'gsvaColRanks()'."))

              if (verbose && gsva_global$show_start_and_end_messages) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
              }

              .check_bpparam(BPPARAM)

              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              maxmem <- .check_maxmem(param, maxmem=maxmem, verbose=verbose)
              ondisk <- .check_ondisk(param, maxmem=maxmem, first=NA, last=NA,
                                      whdim=2, verbose=verbose)

              dataMatrix <- .check_sparse_load_input_expr(dataMatrix, "GSVA",
                                                          first=NA, last=NA,
                                                          whdim=2, ondisk,
                                                          verbose)

              filtDataMatrix <- dataMatrix
              BPPARAM <- .check_open_parallelism(filtDataMatrix, BPPARAM,
                                                 minparrows=100, minparcols=100,
                                                 verbose)

              if (.get_filterRows(param))
                  filtDataMatrix <- .filterGenes(dataMatrix, anyNA(param),
                                                 removeConstant=TRUE,
                                                 removeNzConstant=TRUE,
                                                 errorOnTooFewRows=TRUE,
                                                 verbose=verbose,
                                                 BPPARAM=BPPARAM, maxmem=maxmem)
              else if (verbose) {
                  msg <- "Skipping filtering of constant rows (filterRows=FALSE)"
                  cli_alert_warning(msg)
              }
              
              if (verbose)
                  cli_alert_info(sprintf("Calculating GSVA ranks"))

              kcdfminssize <- .get_kcdfNoneMinSampleSize(param)
              gsvarnorm <- .compute_row_norm(expr=filtDataMatrix,
                                             kcdf=.get_kcdf(param),
                                             kcdf.min.ssize=kcdfminssize,
                                             sparse=.get_sparse(param),
                                             any_na=anyNA(param),
                                             na_use=.get_NAuse(param),
                                             verbose=verbose,
                                             BPPARAM=BPPARAM,
                                             maxmem=maxmem)

              gsvarnks <- .compute_gsva_ranks(Z=gsvarnorm,
                                              verbose=verbose,
                                              BPPARAM=BPPARAM,
                                              maxmem=maxmem)

              rownames(gsvarnks) <- rownames(filtDataMatrix)
              colnames(gsvarnks) <- colnames(filtDataMatrix)

              rnkscontainer <- wrapData(get_exprData(param), gsvarnks, param,
                                        "gsvaranks", first=NA, last=NA, whdim=2,
                                        dropAssays=FALSE)
              rval <- new("gsvaRanksParam",
                          exprData=rnkscontainer, geneSets=get_geneSets(param),
                          assay="gsvaranks", annotation=get_annotation(param),
                          minSize=get_minSize(param), maxSize=get_maxSize(param),
                          kcdf=.get_kcdf(param),
                          kcdfNoneMinSampleSize=.get_kcdfNoneMinSampleSize(param),
                          tau=.get_tau(param), maxDiff=.get_maxDiff(param),
                          absRanking=.get_absRanking(param),
                          sparse=.get_sparse(param), checkNA=.get_checkNA(param),
                          didCheckNA=.get_didCheckNA(param), anyNA=anyNA(param),
                          use=.get_NAuse(param), filterRows=.get_filterRows(param),
                          nzcount=nzcount(param), ondisk=.get_ondisk(param))

              if (verbose && gsva_global$show_start_and_end_messages)
                  cli_alert_success("Calculations finished")

              return(rval)
          })


#' @description The `gsvaScores()` method is deprecated. Please use
#' `gsvaColScores()` instead.
#'
#' @param param A parameter object of the [`gsvaRanksParam-class`] class.
#'
#' @aliases gsvaScores,gsvaRanksParam-method
#' @name gsvaScores
#' @rdname GSVA-pkg-deprecated
#'
#' @importFrom S4Arrays is_sparse
#' @importFrom cli cli_alert_info cli_alert_success
#' @exportMethod gsvaScores
setMethod("gsvaScores", signature(param="gsvaRanksParam"),
          function(param, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {

              .Deprecated(new="gsvaColScores()",
                          package="GSVA",
                          msg=paste("The 'gsvaScores()' method is deprecated.",
                                    "Please use 'gsvaColScores()'."))

              if (verbose && gsva_global$show_start_and_end_messages) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
              }

              .check_bpparam(BPPARAM)

              ## assuming rows in the rank data have been already filtered
              exprData <- get_exprData(param)
              filtDataMatrix <- unwrapData(exprData, get_assay(param))

              filtMappedGeneSets <- .filterAndMapGeneSets(param=param,
                                           filteredDataMatrix=filtDataMatrix,
                                           verbose=verbose)

              sparse <- .get_sparse(param)
              if (sparse && !is_sparse(filtDataMatrix))
                  sparse <- FALSE

              if (verbose) {
                if (sparse)
                    cli_alert_info("GSVA sparse algorithm")
                  else
                    cli_alert_info("GSVA dense (classical) algorithm")
              }

              maxmem <- .check_maxmem(param, maxmem=maxmem, verbose=verbose)
              ondisk <- .check_ondisk(param, maxmem=maxmem, first=NA, last=NA,
                                      whdim=2, verbose=verbose)

              filtDataMatrix <- .check_sparse_load_input_expr(filtDataMatrix,
                                                              "GSVA", first=NA,
                                                              last=NA, whdim=2,
                                                              ondisk, verbose)

              BPPARAM <- .check_open_parallelism(filtDataMatrix, BPPARAM,
                                                 minparrows=100, minparcols=100,
                                                 verbose)

              ondisk <- .check_es_memory_requirements(filtDataMatrix,
                                                      filtMappedGeneSets,
                                                      ondisk, maxmem)
              if (verbose) {
                  n <- length(filtMappedGeneSets)
                  cli_alert_info("Calculating GSVA scores for {n} gene sets")
              }

              gsva_es <- .processMatrixCols(filtDataMatrix,
                                            FUN=.compute_gsva_scores,
                                            geneSetsIdx=filtMappedGeneSets,
                                            tau=.get_tau(param),
                                            maxDiff=.get_maxDiff(param),
                                            absRanking=.get_absRanking(param),
                                            sparse=sparse, any_na=anyNA(param),
                                            na_use=.get_NAuse(param),
                                            minSize=get_minSize(param),
                                            ondisk=ondisk, verbose=verbose,
                                            minparrows=100, minparcols=100,
                                            BPPARAM=BPPARAM,
                                            maxmem=ceiling(maxmem/100)) ## use
                                            ## of memory increases here about
                                            ## 10-fold over block size memory

              rownames(gsva_es) <- names(filtMappedGeneSets)
              colnames(gsva_es) <- colnames(filtDataMatrix)

              gs <- .geneSetsIndices2Names(indices=filtMappedGeneSets,
                                           names=rownames(filtDataMatrix))
              rval <- wrapData(get_exprData(param), gsva_es, param, "es",
                               first=NA, last=NA, whdim=2, dropAssays=FALSE, gs)

              if (verbose && gsva_global$show_start_and_end_messages)
                  cli_alert_success("Calculations finished")

              return(rval)
          })

#' @description The `saveHDF5GSVAranks()` function is deprecated. Please use
#' `saveHDF5GSVA()` instead.
#'
#' @param rankExprData A column-rank expression data set obtained with
#' [`gsvaColRanks`]. Must be one of the classes supported by
#' [`GsvaExprData-class`]. For a list of these classes, see its help page
#' using `help(GsvaExprData)`.
#'
#' @return For `saveHDF5GSVAranks()`, the path to the directory where the data
#' has been saved is returned invisibly. For `loadHDF5GSVAranks()`, an object
#' is returned containing the corresponding loaded GSVA row-normalized or rank
#' expression values, and their corresponding metadata. If the saved GSVA
#' output was originally stored in a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment] object
#' or one of its derived classes, then the returned object will be a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment].
#' Otherwise, the returned object will be a
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix] object.
#'
#' @name saveHDF5GSVAranks
#' @rdname GSVA-pkg-deprecated
#'
#' @export
saveHDF5GSVAranks <- function(rankExprData, dir, ...) {
    .Deprecated(new="saveHDF5GSVA", old="saveHDF5GSVAranks",
                msg=paste("The 'saveHDF5GSVAranks()' function is deprecated.",
                          "Please use 'saveHDF5GSVA()' instead."))

    saveHDF5GSVA(rankExprData, dir, assay="gsvaranks", ...)
}

#' @description The `loadHDF5GSVAranks()` function is deprecated. Please use
#' `loadHDF5GSVA()` instead.
#'
#' @param dir The path to the directory where to save or load the GSVA rank
#' values.
#'
#' @param ... Additional arguments to be passed to the underlying HDF5
#' saving/loading functions
#' [`saveHDF5SummarizedExperiment`][HDF5Array::saveHDF5SummarizedExperiment]
#' and [`loadHDF5SummarizedExperiment`][HDF5Array::loadHDF5SummarizedExperiment],
#' respectively.
#'
#' @name loadHDF5GSVAranks
#' @rdname GSVA-pkg-deprecated
#'
#' @export
loadHDF5GSVAranks <- function(dir, ...) {

    .Deprecated(new="loadHDF5GSVA", old="loadHDF5GSVAranks",
                msg=paste("The 'loadHDF5GSVAranks()' function is deprecated.",
                          "Please use 'loadHDF5GSVA()' instead."))

    rankscontainer <- loadHDF5GSVA(dir, assay="gsvaranks", ...)

    return(rankscontainer)
}

