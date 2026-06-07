#' @title Save/load GSVA output to disk using HDF5 format
#'
#' @description The functions `saveHDF5GSVA()` and `loadHDF5GSVA()` allow one
#' to save and load the output from GSVA to/from disk. The `saveHDF5GSVA()`
#' function takes the output of [`gsvaRowNorm`] or [`gsvaColRanks`] as input,
#' and saves the output from these methods with the relevant metadata to a
#' specified directory. The `loadHDF5GSVA()` function reads the saved data from
#' the specified directory and returns an object with the corresponding
#' GSVA row-normalized or rank expression values, and their corresponding
#' metadata.
#'
#' @param gsvaExprData An obtained either with [`gsvaRowNorm`] or
#' [`gsvaColRanks`]. Must be one of the classes supported by
#' [`GsvaExprData-class`].  For a list of these classes, see its help page
#' using `help(GsvaExprData)`.
#'
#' @param dir The path to the directory where to save or load the GSVA output
#' data.
#'
#' @param assay A single character string specifying the assay that contains
#' the GSVA output to be saved or loaded. By default, `assay="auto"`, which in
#' the case of saving and a `gsvaExprData` object that is a
#' `SummarizedExperiment` or one of its derivatives, it will look for an assay
#' named `gsvaranks`, and if not found, it will look for an assay named
#' `gsvarnorm'. If `gsvaExprData` is not a `SummarizedExperiment` or one of its
#' derivatives, then the assay to be saved will be determined by the `assay`
#' attribute of the `gsvaExprData` object.
#'
#' @param ... Additional arguments to be passed to the underlying HDF5
#' saving/loading functions
#' [`saveHDF5SummarizedExperiment`][HDF5Array::saveHDF5SummarizedExperiment]
#' and [`loadHDF5SummarizedExperiment`][HDF5Array::loadHDF5SummarizedExperiment],
#' respectively.
#'
#' @return For `saveHDF5GSVArnorm()` and `saveHDF5GSVAranks()`, the path to the
#' directory where the data has been saved is returned invisibly. For
#' `loadHDF5GSVArnorm()` and `loadHDF5GSVAranks()`, an object is returned
#' containing the corresponding loaded GSVA row-normalized or rank expression
#' values, and their corresponding metadata. If the saved GSVA output was
#' originally stored in a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment] object
#' or one of its derived classes, then the returned object will be a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment].
#' Otherwise, the returned object will be a
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix] object.
#'
#' @examples
#'
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#' nGrp1 <- 15 ## number of samples in group 1
#' nGrp2 <- n - nGrp1 ## number of samples in group 2
#'
#' ## consider three disjoint gene sets
#' geneSets <- list(gset1=paste0("g", 1:3),
#'                  gset2=paste0("g", 4:6),
#'                  gset3=paste0("g", 7:10))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate row-normalized expression values
#' gsvarownorm <- gsvaRowNorm(gsvapar)
#'
#' ## calculate GSVA column ranks
#' gsvacolranks <- gsvaColRanks(gsvarownorm)
#'
#' ## calculate GSVA scores
#' es <- gsvaColScores(gsvacolranks)
#'
#' ## save the GSVA ranks to disk
#' dir <- tempfile()
#' saveHDF5GSVAranks(gsvacolranks, dir)
#'
#' ## load the GSVA ranks from disk
#' loaded_gsvacolranks <- loadHDF5GSVAranks(dir)
#'
#' ## check that the loaded ranks provide the
#' ## same scores as the original ranks
#' loaded_es <- gsvaColScores(loaded_gsvacolranks)
#' identical(es, loaded_es)
#'
#' @importFrom cli cli_abort
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @importFrom S4Vectors metadata "metadata<-"
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames "assay<-"
#'
#' @rdname gsva-serialization
#'
#' @export
saveHDF5GSVA <- function(gsvaExprData, dir, assay="auto", ...) {
    if (!is(gsvaExprData, "GsvaExprData")) {
        msg <- paste("The input object in 'gsvaExprData' must a subclass of",
                     "'GsvaExprData'. See 'help(GsvaExprData)' for details.")
        cli_abort(msg)
    }

    se <- gsvaExprData
    if (!is(se, "SummarizedExperiment")) {
        if (!is.null(attributes(gsvaExprData)$assay))
            assay <- attributes(gsvaExprData)$assay
        if (!assay %in% c("gsvarnorm", "gsvaranks"))
            cli_abort(c("x"=paste("The object in 'gsvaExprData' does not have",
                                  "GSVA output.")))
        param <- .pull_param(gsvaExprData, "dummy")
        first <- last <- NA_real_
        whdim <- NA_integer_
        annot <- NULL
        if (!is.null(attributes(gsvaExprData)$annotation) &&
            is(attributes(gsvaExprData)$annotation, "GeneIdentifierType")) {
            annot <- attributes(gsvaExprData)$annotation
            attributes(gsvaExprData)$annotation <- NULL
        }
        if (!is.null(attributes(gsvaExprData)$restrict)) {
            first <- attributes(gsvaExprData)$restrict$first
            last <- attributes(gsvaExprData)$restrict$last
            whdim <- attributes(gsvaExprData)$restrict$whdim
            attributes(gsvaExprData)$restrict <- NULL
        }
        se <- SummarizedExperiment(assays=list(dummy=gsvaExprData))
        if (!is.null(annot))
            gsvaAnnotation(se) <- annot
        se <- wrapData(se, gsvaExprData, param, assay, first=first,
                       last=last, whdim=whdim, dropAssays=TRUE)

    } else {
        ## 'SummarizedExperiment' object, remove all assays except 'gsvarnorm'
        an <- assayNames(se)
        assay <- .check_assay_ranks_rnorm(an, assay)

        if (is.null(metadata(se)$gsvaParam)) {
            msg <- paste("Cannot find the GSVA parameters in the metadata of the",
                         "input object given in the 'gsvaExprData' parameter.")
            cli_abort(c("x"=msg))
        }

        for (a in an) ## remove all assays except the selected one
            if (a != assay)
                assay(se, a) <- NULL
    }

  saveHDF5SummarizedExperiment(se, dir, ...)

  invisible(dir)
}


#' @importFrom cli cli_abort
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames
#'
#' @rdname gsva-serialization
#'
#' @export
loadHDF5GSVA <- function(dir, assay="auto", ...) {

    gsvacontainer <- loadHDF5SummarizedExperiment(dir, ...)

    assay <- .check_assay_ranks_rnorm(assayNames(gsvacontainer), assay)

    if (is.null(metadata(gsvacontainer)$gsvaParam)) {
        msg <- paste("Cannot find the GSVA parameters in the metadata of the",
                     "loaded GSVA container object.")
        cli_abort(c("x"=msg))
    }

    if (is.null(metadata(gsvacontainer)$gsvaParam$originalClassWasSE)) {
        msg <- paste("Cannot find the GSVA parameters in the metadata of the",
                     "loaded GSVA container object.")
        cli_abort(c("x"=msg))
    }

    if (!metadata(gsvacontainer)$gsvaParam$originalClassWasSE) {
        gsvapar <- metadata(gsvacontainer)$gsvaParam
        restrict <- metadata(gsvacontainer)$restrict
        annotation <- metadata(gsvacontainer)$annotation
        gsvacontainer <- unwrapData(gsvacontainer, assay)
        attr(gsvacontainer, "gsvaParam") <- gsvapar
        attr(gsvacontainer, "assay") <- assay
        if (!is.null(annotation))
            attr(gsvacontainer, "geneIdType") <- annotation
        if (!is.null(restrict))
            attr(gsvacontainer, "restrict") <- restrict
    }

    return(gsvacontainer)
}

#' @importFrom SummarizedExperiment assayNames
.check_assay_ranks_rnorm <- function(an, assay) {
    if (assay == "auto") {
        if ("gsvaranks" %in% an)
            assay <- "gsvaranks"
        else if ("gsvarnorm" %in% an)
            assay <- "gsvarnorm"
        else
            cli_abort(c("x"="Cannot find a GSVA assay in the object."))
    } else
        if (!assay %in% an)
            cli_abort(c("x"="Cannot find assay {assay} in the object."))

    assay
}
