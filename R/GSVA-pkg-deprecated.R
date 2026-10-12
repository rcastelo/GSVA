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

              showmsg <- verbose && gsva_global$show_start_and_end_messages
              if (showmsg) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
              }
              ## the steps do not show their own start and end messages
              oldmsg <- gsva_global$show_start_and_end_messages
              gsva_global$show_start_and_end_messages <- FALSE
              on.exit(gsva_global$show_start_and_end_messages <- oldmsg,
                      add=TRUE)

              .check_bpparam(BPPARAM)

              ## the ranks are calculated as gsva() does, with the row
              ## normalization and the column ranks steps, which decide where
              ## to process their data and check their memory requirements
              gsvarnorm <- gsvaRowNorm(param=param, verbose=verbose,
                                       dropExistingAssays=TRUE,
                                       errorOnTooFewRows=TRUE,
                                       BPPARAM=BPPARAM, maxmem=maxmem)

              ## the input data in main memory remains allocated while the
              ## column ranks are calculated, which their blocks cannot take
              oldheld <- gsva_global$heldmem
              gsva_global$heldmem <- oldheld + .input_held_mem(param)
              on.exit(gsva_global$heldmem <- oldheld, add=TRUE)

              gsvaranks <- gsvaColRanks(rowNormExprData=gsvarnorm,
                                        verbose=verbose,
                                        dropExistingAssays=TRUE,
                                        BPPARAM=BPPARAM, maxmem=maxmem)
              rm(gsvarnorm)
              gsvarnks <- unwrapData(gsvaranks, "gsvaranks")

              rnkscontainer <- wrapData(get_exprData(param), gsvarnks, param,
                                        "gsvaranks", first=NA, last=NA, rem=NA,
                                        whdim=2, dropAssays=FALSE)
              rval <- new("gsvaRanksParam",
                          exprData=rnkscontainer, geneSets=get_geneSets(param),
                          assay="gsvaranks", annotation=get_annotation(param),
                          minSize=get_minSize(param), maxSize=get_maxSize(param),
                          rowNorm=.get_rowNorm(param), kcdf=.get_kcdf(param),
                          kcdfNoneMinSampleSize=.get_kcdfNoneMinSampleSize(param),
                          tau=.get_tau(param), maxDiff=.get_maxDiff(param),
                          absRanking=.get_absRanking(param),
                          sparse=.get_sparse(param), checkNA=.get_checkNA(param),
                          didCheckNA=.get_didCheckNA(param), anyNA=anyNA(param),
                          use=.get_NAuse(param), filterRows=.get_filterRows(param),
                          nzcount=nzcount(param), ondisk=.get_ondisk(param))

              if (showmsg)
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
#' @exportMethod gsvaScores
setMethod("gsvaScores", signature(param="gsvaRanksParam"),
          function(param, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {

              .Deprecated(new="gsvaColScores()",
                          package="GSVA",
                          msg=paste("The 'gsvaScores()' method is deprecated.",
                                    "Please use 'gsvaColScores()'."))

              ## the scores are calculated by gsvaColScores() on the ranks,
              ## with the gene sets and the other parameters of 'param'
              gsvaColScores(rankExprData=.ranks_param_data(param),
                            verbose=verbose, BPPARAM=BPPARAM, maxmem=maxmem)
          })

## column ranks of the parameter object 'param' given by gsvaRanks(), with its
## metadata set to the parameters of 'param', which may have changed after
## gsvaRanks() gave it, e.g., its gene sets, keeping the assay of the input
## data from which the ranks were calculated
#' @importFrom S4Vectors metadata metadata<-
.ranks_param_data <- function(param) {
    x <- get_exprData(param)
    p <- .gsvaParam_as_list(param)
    if (is(x, "SummarizedExperiment")) {
        p$assay <- metadata(x)$gsvaParam$assay
        metadata(x)$gsvaParam <- p
    } else {
        p$assay <- attr(x, "gsvaParam")$assay
        attr(x, "gsvaParam") <- p
    }

    x
}

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
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @name saveHDF5GSVAranks
#' @rdname GSVA-pkg-deprecated
#'
#' @export
saveHDF5GSVAranks <- function(rankExprData, dir, ...) {
    .Deprecated(new="saveHDF5GSVA", old="saveHDF5GSVAranks",
                msg=paste("The 'saveHDF5GSVAranks()' function is deprecated.",
                          "Please use 'saveHDF5GSVA()' instead."))

    ## saves the ranks in the directory format of the HDF5Array package, as
    ## before saveHDF5GSVA() saved GSVA output to a single HDF5 file
    se <- .gsva_output_to_se(rankExprData, "gsvaranks")$se
    saveHDF5SummarizedExperiment(se, dir, ...)

    invisible(dir)
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
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment assayNames
#' @name loadHDF5GSVAranks
#' @rdname GSVA-pkg-deprecated
#'
#' @export
loadHDF5GSVAranks <- function(dir, ...) {

    .Deprecated(new="loadHDF5GSVA", old="loadHDF5GSVAranks",
                msg=paste("The 'loadHDF5GSVAranks()' function is deprecated.",
                          "Please use 'loadHDF5GSVA()' instead."))

    ## loads ranks saved in the directory format of the HDF5Array package, by
    ## saveHDF5GSVAranks() or by saveHDF5GSVA() in previous versions of GSVA
    rankscontainer <- loadHDF5SummarizedExperiment(dir, ...)
    assay <- .check_assay_ranks_rnorm(assayNames(rankscontainer), "gsvaranks")
    rankscontainer <- .se_to_gsva_output(rankscontainer, assay)

    return(rankscontainer)
}

