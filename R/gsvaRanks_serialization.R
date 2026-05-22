#' @title Save/load GSVA rank values to disk using HDF5 format
#'
#' @description The functions `saveHDF5GSVAranks` and `loadHDF5GSVAranks` can
#' be used to save and load the GSVA rank values to/from disk, respectively.
#' The `saveHDF5GSVAranks` function takes a `gsvaRanksParam` object and saves
#' the rank values along with the relevant metadata to a specified directory.
#' The `loadHDF5GSVAranks` function reads the saved data from the specified
#' directory and reconstructs the `gsvaRanksParam` object with the rank values
#' and their corresponding metadata.
#'
#' @param x A [`gsvaRanksParam-class`] object to save to disk.
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
#' `gsvaRanksParam` object is returned containing the loaded GSVA rank values
#' and their corresponding metadata.
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
#' ## calculate GSVA ranks
#' gsvarankspar <- gsvaRanks(gsvapar)
#'
#' ## calculate GSVA scores
#' es <- gsvaScores(gsvarankspar)
#'
#' ## save the GSVA ranks to disk
#' dir <- tempfile()
#' saveHDF5GSVAranks(gsvarankspar, dir)
#'
#' ## load the GSVA ranks from disk
#' loaded_gsvarankspar <- loadHDF5GSVAranks(dir)
#'
#' ## check that the loaded ranks provide the
#' ## same scores as the original ranks
#' loaded_es <- gsvaScores(loaded_gsvarankspar)
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
saveHDF5GSVAranks <- function(x, dir, ...) {
    if (!inherits(x, "gsvaRanksParam"))
      cli_abort("The input object in 'x' must be of class 'gsvaRanksParam'")

    edata <- get_exprData(x)
    if (is(edata, "SummarizedExperiment")) {
        an <- assayNames(edata)
        if (!"gsvaranks" %in% an)
            cli_abort("Cannot find the ranks in the input object 'x'")

        for (a in an)
            if (a != "gsvaranks")
                assay(edata, a) <- NULL
    } else {
        annot <- NULL
        if (!is.null(attributes(edata)$annotation) &&
            is(attributes(edata)$annotation, "GeneIdentifierType")) {
          annot <- attributes(edata)$annotation
          attributes(edata)$annotation <- NULL
        }
        edata <- SummarizedExperiment(assays=list(gsvaranks=edata))
        if (!is.null(annot))
          gsvaAnnotation(edata) <- annot
    }

    metadata(edata) <- c(metadata(edata),
                         list(gsvaRanksParam=.gsvaParam_as_list(x)))

  saveHDF5SummarizedExperiment(edata, dir, ...)

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

    x <- loadHDF5SummarizedExperiment(dir, ...)
    rnksmdata <- metadata(x)$gsvaRanksParam
    if (is.null(rnksmdata)) {
        msg <- "The given directory does not contain valid GSVA ranks data"
        cli_abort(c("x"=msg))
    }
    md <- metadata(x)
    md$gsvaRanksParam <- NULL
    metadata(x) <- md
    if (is.null(rnksmdata$originalClassWasSE))
        cli_abort("Metadata is missing the original class information")

    rnkscontainer <- x
    if (!rnksmdata$originalClassWasSE) {
        if (!"gsvaranks" %in% assayNames(x)) {
            msg <- "The given directory does not contain valid GSVA ranks data"
            cli_abort(c("x"=msg))
        }
        rnkscontainer <- assay(x, "gsvaranks")
        if (!is.null(gsvaAnnotation(x)))
            gsvaAnnotation(rnkscontainer) <- gsvaAnnotation(x)
    }

    new("gsvaRanksParam",
        exprData=rnkscontainer,
        geneSets=rnksmdata$geneSets,
        assay=rnksmdata$assay,
        annotation=rnksmdata$annotation,
        minSize=rnksmdata$minSize,
        maxSize=rnksmdata$maxSize,
        kcdf=rnksmdata$kcdf,
        kcdfNoneMinSampleSize=rnksmdata$kcdfNoneMinSampleSize,
        tau=rnksmdata$tau,
        maxDiff=rnksmdata$maxDiff,
        absRanking=rnksmdata$absRanking,
        sparse=rnksmdata$sparse,
        checkNA=rnksmdata$checkNA,
        didCheckNA=rnksmdata$didCheckNA,
        anyNA=rnksmdata$anyNA,
        use=rnksmdata$use,
        filterRows=rnksmdata$filterRows,
        nzcount=rnksmdata$nzcount,
        ondisk=rnksmdata$ondisk)
}
