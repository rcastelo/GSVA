#' @title Save/load GSVA rank values to disk using HDF5 format
#'
#' @description The functions `saveHDF5GSVAranks` and `loadHDF5GSVAranks` can
#' be used to save and load the GSVA rank values to/from disk, respectively.
#' The `saveHDF5GSVAranks()` function takes the output of [`gsvaColRanks`] as
#' input, and saves the rank values along with the relevant metadata to a
#' specified directory. The `loadHDF5GSVAranks()` function reads the saved data
#' from the specified directory and returns an object with the GSVA rank values
#' and their corresponding metadata.
#'
#' @param rankExprData A column-rank expression data set obtained with
#' [`gsvaColRanks`]. Must be one of the classes supported by
#' [`GsvaExprData-class`]. For a list of these classes, see its help page
#' using `help(GsvaExprData)`.
#'
#' @param dir The path to the directory where to save or load the GSVA ranks
#' data.
#'
#' @param ... Additional arguments to be passed to the underlying HDF5
#' saving/loading functions
#' [`saveHDF5SummarizedExperiment`][HDF5Array::saveHDF5SummarizedExperiment]
#' and [`loadHDF5SummarizedExperiment`][HDF5Array::loadHDF5SummarizedExperiment],
#' respectively.
#'
#' @return For `saveHDF5GSVAranks`, the path to the directory where the data
#' has been saved is returned invisibly. For `loadHDF5GSVAranks`, a
#' an object is returned containing the loaded GSVA rank values and their
#' corresponding metadata. If the saved ranks were originally stored in a
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
#' @rdname gsvaRanks_serialization
#'
#' @export
saveHDF5GSVAranks <- function(rankExprData, dir, ...) {
    if (!is(rankExprData, "GsvaExprData")) {
        msg <- paste("The input object in 'rankExprData' must a subclass of",
                     "'GsvaExprData'. See 'help(GsvaExprData)' for details.")
        cli_abort(msg)
    }

    se <- rankExprData
    if (!is(se, "SummarizedExperiment")) {
        param <- .pull_param(rankExprData, "gsvaranks")
        first <- last <- NA_real_
        whdim <- NA_integer_
        annot <- NULL
        if (!is.null(attributes(rankExprData)$annotation) &&
            is(attributes(rankExprData)$annotation, "GeneIdentifierType")) {
            annot <- attributes(rankExprData)$annotation
            attributes(rankExprData)$annotation <- NULL
        }
        if (!is.null(attributes(rankExprData)$restrict)) {
            first <- attributes(rankExprData)$restrict$first
            last <- attributes(rankExprData)$restrict$last
            whdim <- attributes(rankExprData)$restrict$whdim
            attributes(rankExprData)$restrict <- NULL
        }
        se <- SummarizedExperiment(assays=list(dummy=rankExprData))
        if (!is.null(annot))
            gsvaAnnotation(se) <- annot
        se <- wrapData(se, rankExprData, param, "gsvaranks", first=first,
                       last=last, whdim=whdim, dropAssays=TRUE)

    } else { ## 'SummarizedExperiment' object, remove all assays except 'gsvaranks'
        an <- assayNames(se)
        if (!"gsvaranks" %in% an) {
            msg <- paste("Cannot find the ranks in the input object given in the",
                         "'rankExprData' parameter.")
            cli_abort(c("x"=msg))
        }
        if (is.null(metadata(se)$gsvaParam)) {
            msg <- paste("Cannot find the GSVA parameters in the metadata of the",
                         "input object given in the 'rankExprData' parameter.")
            cli_abort(c("x"=msg))
        }

        for (a in an) ## remove all assays except the one with the ranks
            if (a != "gsvaranks")
                assay(se, a) <- NULL
    }

  saveHDF5SummarizedExperiment(se, dir, ...)

  invisible(dir)
}

#' @importFrom cli cli_abort
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames "assay<-"
#'
#' @rdname gsvaRanks_serialization
#'
#' @export
loadHDF5GSVAranks <- function(dir, ...) {

    rankscontainer <- loadHDF5SummarizedExperiment(dir, ...)

    an <- assayNames(rankscontainer)
    if (!"gsvaranks" %in% an)
            cli_abort("Cannot find the ranks in the loaded object.")

    if (is.null(metadata(rankscontainer)$gsvaParam)) {
        msg <- paste("Cannot find the GSVA parameters in the metadata of the",
                     "loaded GSVA ranks object.")
        cli_abort(c("x"=msg))
    }

    if (is.null(metadata(rankscontainer)$gsvaParam$originalClassWasSE)) {
        msg <- paste("Cannot find the GSVA parameters in the metadata of the",
                     "loaded GSVA ranks object.")
        cli_abort(c("x"=msg))
    }

    if (!metadata(rankscontainer)$gsvaParam$originalClassWasSE) {
        gsvapar <- metadata(rankscontainer)$gsvaParam
        restrict <- metadata(rankscontainer)$restrict
        annotation <- metadata(rankscontainer)$annotation
        rankscontainer <- unwrapData(rankscontainer, "gsvaranks")
        attr(rankscontainer, "gsvaParam") <- gsvapar
        attr(rankscontainer, "assay") <- "gsvaranks"
        if (!is.null(annotation))
            attr(rankscontainer, "geneIdType") <- annotation
        if (!is.null(restrict))
            attr(rankscontainer, "restrict") <- restrict
    }

    return(rankscontainer)
}
