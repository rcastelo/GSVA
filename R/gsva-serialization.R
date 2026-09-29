#' @title Save/load GSVA output to disk using HDF5 or Apache Parquet format
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
#' The functions `saveParquetGSVA()` and `loadParquetGSVA()` do the same using
#' a single file in Apache Parquet format, which is organized to efficiently
#' read blocks of columns, as [`gsvaColRanks`] and [`gsvaColScores`] do. These
#' two functions require the package
#' [arrow](https://cran.r-project.org/package=arrow).
#'
#' @details `saveParquetGSVA()` stores sparse input, such as a
#' [`dgCMatrix`][Matrix::dgCMatrix-class] or an
#' [`SVT_SparseMatrix`][SparseArray::SVT_SparseMatrix-class] object, keeping
#' only its non-zero values, and stores GSVA ranks as integer values. The
#' object returned by `loadParquetGSVA()` holds its values in a
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix] object that reads them from
#' the file only when needed, by blocks of columns, and [`gsvaColScores`]
#' processes ranks stored in this way from disk. `loadParquetGSVA()` can also
#' read files stored in Amazon S3 or Google Cloud Storage, given as `s3://` or
#' `gs://` URIs, as long as the installed arrow package supports them; see
#' [`arrow_with_s3`][arrow::arrow_with_s3].
#'
#' @param gsvaExprData An object obtained with [`gsvaRowNorm`] or
#' [`gsvaColRanks`]. Must be one of the classes supported by
#' [`GsvaExprData-class`].  For a list of these classes, see its help page
#' using `help(GsvaExprData)`.
#'
#' @param dir The path to the directory where to save or load the GSVA output
#' data.
#'
#' @param file The path to the file where to save the GSVA output data in
#' Apache Parquet format or, for loading it, either that path or an `s3://` or
#' `gs://` URI.
#'
#' @param colsPerRowGroup Either `"auto"` (default), or the number of matrix
#' columns stored in each row group of the Apache Parquet file, which is the
#' smallest unit of data that can be read from it. With `"auto"`, dense data is
#' stored in row groups of at most 2^20 values, which is also the maximum
#' number of values allowed, while sparse data is stored in row groups of 500
#' columns.
#'
#' @param replace Logical vector of length 1. When `TRUE`, an existing file in
#' `file` is replaced. By default, `replace=FALSE`.
#'
#' @param verbose Logical vector of length 1. When `TRUE`, a progress bar is
#' shown while writing the file. By default, `verbose=FALSE`.
#'
#' @param assay A single character string specifying the assay that contains
#' the GSVA output to be saved or loaded. By default, `assay="auto"`, which in
#' the case of saving and a `gsvaExprData` object that is a
#' `SummarizedExperiment` or one of its derivatives, it will look for an assay
#' named `gsvaranks`, and if not found, it will look for an assay named
#' `gsvarnorm`. If `gsvaExprData` is not a `SummarizedExperiment` or one of its
#' derivatives, then the assay to be saved will be determined by the `assay`
#' attribute of the `gsvaExprData` object.
#'
#' @param ... Only for `saveHDF5GSVA()` and `loadHDF5GSVA()`, additional
#' arguments to be passed to the underlying HDF5 saving/loading functions
#' [`saveHDF5SummarizedExperiment`][HDF5Array::saveHDF5SummarizedExperiment]
#' and [`loadHDF5SummarizedExperiment`][HDF5Array::loadHDF5SummarizedExperiment],
#' respectively.
#'
#' @return For `saveHDF5GSVA()` the path to the directory where the data has
#' been saved is returned invisibly, and for `saveParquetGSVA()` the path to
#' the file. For `loadHDF5GSVA()` and `loadParquetGSVA()`, an object is returned
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
#' ## save the GSVA row-normalized expression values to disk
#' rnormdir <- tempfile()
#' saveHDF5GSVA(gsvarownorm, rnormdir)
#'
#' ## save the GSVA rank values to disk
#' ranksdir <- tempfile()
#' saveHDF5GSVA(gsvacolranks, ranksdir)
#'
#' ## load the GSVA row-normalized values from disk               
#' loaded_gsvarownorm <- loadHDF5GSVA(rnormdir)
#'
#' ## check that the loaded row-normalized values provide the
#' ## same ranks as the ones calculated from the original values
#' gsvacolranks_from_loaded_gsvarownorm <- gsvaColRanks(loaded_gsvarownorm)
#'
#' identical(gsvacolranks, gsvacolranks_from_loaded_gsvarownorm)
#'
#' ## load the GSVA ranks from disk
#' loaded_gsvacolranks <- loadHDF5GSVA(ranksdir)
#'
#' ## check that the loaded ranks provide the
#' ## same scores as the original ranks
#' loaded_es <- gsvaColScores(loaded_gsvacolranks)
#' identical(es, loaded_es)
#'
#' ## the same using Apache Parquet files
#' if (requireNamespace("arrow", quietly=TRUE)) {
#'     ranksfile <- tempfile(fileext=".parquet")
#'     saveParquetGSVA(gsvacolranks, ranksfile)
#'     loaded_gsvacolranks <- loadParquetGSVA(ranksfile)
#'     loaded_es <- gsvaColScores(loaded_gsvacolranks)
#'     all.equal(es, loaded_es, check.attributes=FALSE)
#' }
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
    se <- .gsva_output_to_se(gsvaExprData, assay)$se

    saveHDF5SummarizedExperiment(se, dir, ...)

    invisible(dir)
}


## put the GSVA output in 'gsvaExprData' into a 'SummarizedExperiment' object
## with a single assay and the GSVA metadata, as required for serialization;
## returns a list with that object in 'se' and the name of the assay in 'assay'
#' @importFrom cli cli_abort
#' @importFrom S4Vectors metadata "metadata<-"
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames "assay<-"
.gsva_output_to_se <- function(gsvaExprData, assay) {
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
        param <- .pull_param(gsvaExprData)
        first <- last <- NA_real_
        rem <- 0
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
            if (!is.null(attributes(gsvaExprData)$restrict$rem))
                rem <- attributes(gsvaExprData)$restrict$rem
            attributes(gsvaExprData)$restrict <- NULL
        }
        se <- SummarizedExperiment(assays=list(dummy=gsvaExprData))
        if (!is.null(annot))
            gsvaAnnotation(se) <- annot
        se <- wrapData(se, gsvaExprData, param, assay, first=first,
                       last=last, rem=rem, whdim=whdim, dropAssays=TRUE)

    } else {
        ## 'SummarizedExperiment' object, remove all assays except the selected GSVA assay
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

    list(se=se, assay=assay)
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

    .se_to_gsva_output(gsvacontainer, assay)
}


## turn a 'SummarizedExperiment' object with GSVA output, read from disk, back
## into the class of the object that was saved: if it was not originally a
## 'SummarizedExperiment', return the matrix in the assay 'assay' with the
## GSVA metadata stored as attributes
#' @importFrom cli cli_abort
#' @importFrom S4Vectors metadata
.se_to_gsva_output <- function(gsvacontainer, assay) {
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
        else if ("es" %in% an)
            assay <- "es"
        else
            cli_abort(c("x"="Cannot find a GSVA assay in the object."))
    } else
        if (!assay %in% an)
            cli_abort(c("x"="Cannot find assay {assay} in the object."))

    assay
}


#' @importFrom cli cli_abort
#' @importFrom BiocGenerics type
#' @importFrom SummarizedExperiment assay "assay<-"
#'
#' @rdname gsva-serialization
#'
#' @export
saveParquetGSVA <- function(gsvaExprData, file, assay="auto",
                            colsPerRowGroup="auto", replace=FALSE,
                            verbose=FALSE) {
    if (!is.character(file) || length(file) != 1L || is.na(file))
        cli_abort(c("x"="'file' must be a single character string."))
    if (!is.logical(replace) || length(replace) != 1L || is.na(replace))
        cli_abort(c("x"="'replace' must be either TRUE or FALSE."))

    .require_arrow()

    if (.is_uri(file))
        cli_abort(c("x"=paste("GSVA output in Apache Parquet format can only",
                              "be saved to local files.")))
    if (dir.exists(file))
        cli_abort(c("x"="{.file {file}} is a directory."))
    if (file.exists(file) && !replace)
        cli_abort(c("x"=paste("The file {.file {file}} already exists; use",
                              "'replace=TRUE' to replace it.")))
    if (!dir.exists(dirname(file)))
        cli_abort(c("x"=paste("The directory {.file {dirname(file)}} does not",
                              "exist.")))

    cnt <- .gsva_output_to_se(gsvaExprData, assay)
    X <- assay(cnt$se, cnt$assay, withDimnames=FALSE)

    ## GSVA ranks are integer values, even when they are stored as doubles,
    ## such as in a 'dgCMatrix' object
    type <- if (cnt$assay == "gsvaranks") "integer" else type(X)
    if (!type %in% c("integer", "double"))
        cli_abort(c("x"=paste("Values of type {.val {type}} cannot be stored",
                              "in Apache Parquet format.")))

    ## the rest of the container, including dimnames and GSVA metadata
    shell <- cnt$se
    assay(shell, cnt$assay) <- NULL
    md <- list(gsva_assay=cnt$assay, gsva_shell=.encode_r_object(shell))

    .write_gsva_parquet(X, file, type, colsPerRowGroup, md, verbose)

    invisible(file)
}


#' @importFrom cli cli_abort
#' @importFrom BiocGenerics path
#' @importFrom DelayedArray DelayedArray
#' @importFrom SummarizedExperiment assayNames "assay<-"
#'
#' @rdname gsva-serialization
#'
#' @export
loadParquetGSVA <- function(file, assay="auto") {
    seed <- GsvaParquetSeed(file)
    md <- .parquet_reader(path(seed))$GetSchema()$metadata
    if (is.null(md$gsva_assay) || is.null(md$gsva_shell))
        cli_abort(c("x"=paste("The file {.file {file}} does not contain GSVA",
                              "output saved with 'saveParquetGSVA()'.")))

    gsvacontainer <- .decode_r_object(md$gsva_shell)
    X <- DelayedArray(seed)
    dimnames(X) <- dimnames(gsvacontainer)
    assay(gsvacontainer, md$gsva_assay) <- X

    assay <- .check_assay_ranks_rnorm(assayNames(gsvacontainer), assay)

    .se_to_gsva_output(gsvacontainer, assay)
}


## load GSVA output given as a path in the arguments 'rowNormExprData' of
## gsvaColRanks() or 'rankExprData' of gsvaColScores(), named in 'argname',
## which can be either a directory with GSVA output saved with saveHDF5GSVA(),
## or a file or an 's3://' or 'gs://' URI with GSVA output saved with
## saveParquetGSVA(); 'assay' is the name of the GSVA assay to load
#' @importFrom cli cli_abort cli_alert_info
.load_gsva_path <- function(path, assay, argname, verbose) {
    if (length(path) != 1L || is.na(path))
        cli_abort(c("x"="'{argname}' must be a single character string."))

    if (.is_uri(path))
        parquet <- TRUE
    else if (dir.exists(path))
        parquet <- FALSE
    else if (file.exists(path)) {
        if (!.is_parquet_file(path))
            cli_abort(c("x"=paste("{.file {path}} is neither a directory with",
                                  "GSVA output saved with 'saveHDF5GSVA()',",
                                  "nor a file with GSVA output saved with",
                                  "'saveParquetGSVA()'.")))
        parquet <- TRUE
    } else
        cli_abort(c("x"="{path} cannot be found in the filesystem"))

    if (verbose) {
        if (.is_uri(path))
            cli_alert_info("Loading {path}")
        else
            cli_alert_info("Loading {basename(path)} from disk")
    }

    if (parquet)
        loadParquetGSVA(path, assay=assay)
    else
        loadHDF5GSVA(path, assay=assay)
}
