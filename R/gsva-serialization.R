#' @title Save/load GSVA output to disk using HDF5 or Apache Parquet format
#'
#' @description The functions `saveHDF5GSVA()` and `loadHDF5GSVA()` allow one
#' to save and load the output from GSVA to/from disk. The `saveHDF5GSVA()`
#' function takes the output of [`gsvaRowNorm`], [`gsvaColRanks`] or
#' [`gsvaColScores`] as input, and saves the output from these methods with the
#' relevant metadata to a single file in HDF5 format. The `loadHDF5GSVA()`
#' function reads the saved data from that file and returns an object with
#' the corresponding GSVA row-normalized or rank expression values, or GSVA
#' scores, and their corresponding metadata.
#'
#' The functions `saveParquetGSVA()` and `loadParquetGSVA()` do the same using
#' a single file in Apache Parquet format, which is organized to efficiently
#' read blocks of columns, as [`gsvaColRanks`] and [`gsvaColScores`] do. These
#' two functions require having installed additional packages by the user. More
#' concretely, loading Parquet files from the local filesystem or lazily over
#' Amazon S3 or Google Cloud Storage buckets, and saving them, require the CRAN
#' package [arrow](https://cran.r-project.org/package=arrow), while loading
#' Parquet files lazily over HTTP(S) URLs requires the CRAN packages
#' [duckdb](https://cran.r-project.org/package=duckdb) and
#' [DBI](https://cran.r-project.org/package=DBI).
#'
#' @details `saveHDF5GSVA()` stores the GSVA output values in a single HDF5
#' file, keeping only the non-zero values of sparse data, such as sparse GSVA
#' ranks, in the compressed sparse column layout used by 10x Genomics, which
#' takes less space and is read faster by blocks of columns, as
#' [`gsvaColRanks`] and [`gsvaColScores`] do, than storing all values. Dense
#' data is stored in a dense HDF5 dataset whose chunks span all its rows. The
#' rest of the object storing the GSVA output, such as its row and column
#' names and data, and its GSVA metadata, is stored in the same file. The file
#' does not store its own location, so it can be moved or copied. The object
#' returned by `loadHDF5GSVA()` holds its values in a
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix] object that reads them from
#' the file only when needed. GSVA output loaded with `loadHDF5GSVA()` can be
#' saved in other formats of the HDF5Array package, e.g., with
#' [`saveHDF5SummarizedExperiment`][HDF5Array::saveHDF5SummarizedExperiment].
#'
#' `saveParquetGSVA()` stores sparse input, such as a
#' [`dgCMatrix`][Matrix::dgCMatrix-class] or an
#' [`SVT_SparseMatrix`][SparseArray::SVT_SparseMatrix-class] object, keeping
#' only its non-zero values, and stores GSVA ranks as integer values. The
#' object returned by `loadParquetGSVA()` holds its values in a
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix] object that reads them from
#' the file only when needed, by blocks of columns, and [`gsvaColScores`]
#' processes ranks stored in this way from disk. `loadParquetGSVA()` can also
#' read files stored in Amazon S3 or Google Cloud Storage, given as `s3://` or
#' `gs://` URIs, as long as the installed arrow package supports them; see
#' [`arrow_with_s3`][arrow::arrow_with_s3]. Files in public Google Cloud Storage
#' buckets should be given with URIs of the form
#' `gs://anonymous@<bucket>/<path>`, which request anonymous access, since
#' otherwise arrow first looks for Google Cloud credentials, and opening the
#' file fails if they are not available. Files in private buckets require
#' setting up such credentials, e.g., with
#' `gcloud auth application-default login`. Requests to Google Cloud Storage
#' are retried for at most 15 seconds, unless the URI sets another limit with
#' the query parameter `retry_limit_seconds`, e.g.,
#' `gs://<bucket>/<path>?retry_limit_seconds=60`. Files in private Amazon S3
#' buckets require AWS credentials, which arrow reads from the environment
#' variables `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` (and
#' `AWS_SESSION_TOKEN` for temporary credentials), or from the file
#' `~/.aws/credentials`, whose profile can be selected with the environment
#' variable `AWS_PROFILE`. Although arrow also accepts credentials within the
#' URI, e.g., `s3://<key>:<secret>@<bucket>/<path>`, this is not recommended,
#' because they may end up stored in scripts and in the R history; GSVA hides
#' them in its messages.
#'
#' `loadParquetGSVA()` can also lazily read files available at public `http://`
#' or `https://` URLs using the API of the package
#' [duckdb](https://cran.r-project.org/package=duckdb), which only downloads
#' the parts of the file that are needed, as with `s3://` and `gs://` URIs.
#' This requires that the HTTP(S) server supports range requests, which most do.
#' The first time, DuckDB downloads its extension `httpfs` from its own servers,
#' which GSVA installs in the directory
#' `file.path(tools::R_user_dir("GSVA", which="cache"), "duckdb")`. If DuckDB
#' cannot verify the certificate of an HTTPS server, e.g., because the server
#' does not send its intermediate certificates, a file with the certificates
#' of the certification authorities can be given with
#' `options(GSVA.ca_cert_file="<file>")`. This file replaces the certificates of
#' the certification authorities that DuckDB trusts, so it should only contain
#' certificates from trustworthy sources, including the root certificates
#' needed by every HTTPS server accessed in the same R session.
#'
#' Loading GSVA output in HDF5 or Apache Parquet format restores R objects
#' stored in the file, as [`readRDS`][base::readRDS] does, and therefore, as
#' with files read with `readRDS()`, files should only be loaded from trusted
#' sources.
#'
#' @param gsvaExprData An object obtained with [`gsvaRowNorm`],
#' [`gsvaColRanks`] or [`gsvaColScores`]. Must be one of the classes supported by
#' [`GsvaExprData-class`].  For a list of these classes, see its help page
#' using `help(GsvaExprData)`.
#'
#' @param file The path to the file where to save the GSVA output data in HDF5
#' or Apache Parquet format or, for loading it, that path or, in Apache Parquet
#' format, also an `s3://` or `gs://` URI, or an `http://` or `https://` URL.
#'
#' @param colsPerRowGroup Either `"auto"` (default), or the number of matrix
#' columns stored in each row group of the Apache Parquet file, which is the
#' smallest unit of data that can be read from it. With `"auto"`, dense data is
#' stored in row groups of at most 2^20 values, which is also the maximum
#' number of values allowed, while sparse data is stored in row groups of 500
#' columns. Values larger than the number of columns of the matrix store all
#' its columns in a single row group.
#'
#' @param replace Logical vector of length 1. When `TRUE`, an existing file in
#' `file` is replaced. By default, `replace=FALSE`.
#'
#' @param verbose Logical vector of length 1. In `saveHDF5GSVA()` and
#' `saveParquetGSVA()`, when `TRUE`, a progress bar is shown while writing the
#' file, and by default `verbose=FALSE`. In `loadHDF5GSVA()` and
#' `loadParquetGSVA()`, when `TRUE` (default), messages inform about the steps
#' of loading the file.
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
#' @return For `saveHDF5GSVA()` and `saveParquetGSVA()`, the path to the file
#' where the data has been saved is returned invisibly. For `loadHDF5GSVA()` and `loadParquetGSVA()`, an object is returned
#' containing the corresponding loaded GSVA row-normalized or rank expression
#' values, and their corresponding metadata. If the saved GSVA output was
#' originally stored in a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment] object
#' or one of its derived classes, then the returned object will be a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment].
#' Otherwise, the returned object will be a
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix] object. This is also the case
#' for GSVA output originally stored in an
#' [`ExpressionSet`][Biobase::ExpressionSet] object, of which only the values
#' and the GSVA metadata are saved.
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
#' rnormfile <- tempfile(fileext=".h5")
#' saveHDF5GSVA(gsvarownorm, rnormfile)
#'
#' ## save the GSVA rank values to disk
#' ranksfile <- tempfile(fileext=".h5")
#' saveHDF5GSVA(gsvacolranks, ranksfile)
#'
#' ## load the GSVA row-normalized values from disk
#' loaded_gsvarownorm <- loadHDF5GSVA(rnormfile)
#'
#' ## check that the loaded row-normalized values provide the
#' ## same ranks as the ones calculated from the original values
#' gsvacolranks_from_loaded_gsvarownorm <- gsvaColRanks(loaded_gsvarownorm)
#'
#' identical(gsvacolranks, gsvacolranks_from_loaded_gsvarownorm)
#'
#' ## load the GSVA ranks from disk
#' loaded_gsvacolranks <- loadHDF5GSVA(ranksfile)
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
#' @importFrom BiocGenerics type
#' @importFrom S4Vectors metadata "metadata<-"
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#' @importFrom SummarizedExperiment "assay<-"
#'
#' @rdname gsva-serialization
#'
#' @export
saveHDF5GSVA <- function(gsvaExprData, file, assay="auto", replace=FALSE,
                         verbose=FALSE) {
    if (!.isCharLength1(file))
        cli_abort(c("x"="'file' must be a single character string."))
    if (!is.logical(replace) || length(replace) != 1L || is.na(replace))
        cli_abort(c("x"="'replace' must be either TRUE or FALSE."))
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
                              "in the HDF5 file.")))

    ## the rest of the container, including dimnames and GSVA metadata
    shell <- cnt$se
    assay(shell, cnt$assay) <- NULL

    ## the file is written with a temporary name and renamed when complete,
    ## so that an interrupted saving does not leave an incomplete file
    tmpname <- .unique_tmpname(file)
    on.exit(unlink(tmpname), add=TRUE)
    .write_gsva_h5(X, tmpname, cnt$assay, type == "integer", shell, verbose)
    if (file.exists(file))
        unlink(file)
    if (!file.rename(tmpname, file))
        cli_abort(c("x"="Cannot write the file {.file {file}}."))

    invisible(file)
}

## write the GSVA output 'X', of the assay named 'assayname', integer or not
## ('int'), with the rest of its container 'shell', into the new HDF5 file
## 'fname', see the details of saveHDF5GSVA(): the serialized 'shell' in the
## dataset "/gsva/shell", the name of the assay in "/gsva/assayname", and
## its values in "/gsva/assay", in CSC layout when they are sparse, or
## otherwise in a dense dataset whose chunks span its rows. its values are
## written by blocks of columns, read from 'X' in main memory or on disk
#' @importFrom rhdf5 h5createFile h5createGroup h5createDataset h5write
#' @importFrom S4Arrays is_sparse
#' @importFrom DelayedArray getAutoBlockLength
#' @importFrom HDF5Array getHDF5DumpChunkLength getHDF5DumpCompressionLevel
#' @importFrom IRanges IRanges
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
.write_gsva_h5 <- function(X, fname, assayname, int, shell, verbose) {
    h5createFile(fname)
    h5createGroup(fname, "gsva")
    r <- serialize(shell, connection=NULL)
    h5createDataset(fname, "gsva/shell", length(r), storage.mode="raw",
                    chunk=max(1L, min(length(r), 2^20)),
                    level=getHDF5DumpCompressionLevel())
    h5write(r, fname, "gsva/shell")
    rm(r)
    sparse <- is_sparse(X)
    h5write(assayname, fname, "gsva/assayname")
    h5write(if (sparse) "csc" else "dense", fname, "gsva/layout")

    nr <- nrow(X)
    nc <- ncol(X)
    ## blocks of columns of the default block size, which, for dense values,
    ## are a multiple of the width of the chunks, so that each block fills
    ## whole chunks
    chunkcols <- max(1, floor(getHDF5DumpChunkLength() / max(1, nr)))
    width <- max(1, floor(getAutoBlockLength(if (int) "integer" else "double") /
                          max(1, nr)))
    if (!sparse) {
        width <- max(chunkcols, chunkcols * floor(width / chunkcols))
        .h5_part_create(fname, c(nr, nc), int, c(nr, chunkcols),
                        name="gsva/assay")
    } else
        st <- .h5_csc_create(fname, int, group="gsva/assay")

    idpb <- NULL
    if (verbose)
        idpb <- cli_progress_bar("Writing the HDF5 file", total=nc)
    for (first in seq(1, max(1, nc), by=width)) {
        if (nc == 0)
            break
        last <- min(nc, first + width - 1)
        block <- if (is(X, "DelayedArray"))
                     .read_block_range(X, 2L, IRanges(first, last))
                 else
                     X[, first:last, drop=FALSE]
        if (sparse)
            st <- .h5_csc_append(st, block)
        else {
            block <- as.matrix(block)
            attributes(block) <- list(dim=dim(block))
            if (int)
                storage.mode(block) <- "integer"
            h5write(block, fname, "gsva/assay", start=c(1, first),
                    count=dim(block))
        }
        if (verbose)
            cli_progress_update(id=idpb, inc=last - first + 1)
    }
    if (sparse)
        .h5_csc_close(st, c(nr, nc))
    if (verbose)
        cli_progress_done(id=idpb)

    invisible(fname)
}



## put the GSVA output in 'gsvaExprData' into a 'SummarizedExperiment' object
## with a single assay and the GSVA metadata, as required for serialization;
## returns a list with that object in 'se' and the name of the assay in 'assay'.
## GSVA scores not stored in a 'SummarizedExperiment' have their gene sets in
## the attribute 'geneSets', which are stored in the row data column 'gs', as
## done for GSVA scores stored in a 'SummarizedExperiment'. Of an
## 'ExpressionSet' object, only its values and GSVA metadata are stored.
#' @importFrom cli cli_abort
#' @importFrom S4Vectors metadata "metadata<-"
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames "assay<-"
#' @importMethodsFrom Biobase exprs
.gsva_output_to_se <- function(gsvaExprData, assay) {
    if (!is(gsvaExprData, "GsvaExprData")) {
        msg <- paste("The input object in 'gsvaExprData' must a subclass of",
                     "'GsvaExprData'. See 'help(GsvaExprData)' for details.")
        cli_abort(msg)
    }

    if (is(gsvaExprData, "ExpressionSet")) {
        gsvaattr <- c("gsvaParam", "assay", "restrict", "geneSets")
        gsvaattr <- attributes(gsvaExprData)[intersect(gsvaattr,
                                                names(attributes(gsvaExprData)))]
        gsvaExprData <- exprs(gsvaExprData)
        for (a in names(gsvaattr))
            attr(gsvaExprData, a) <- gsvaattr[[a]]
    }

    se <- gsvaExprData
    if (!is(se, "SummarizedExperiment")) {
        if (!is.null(attributes(gsvaExprData)$assay))
            assay <- attributes(gsvaExprData)$assay
        if (!assay %in% c("gsvarnorm", "gsvaranks", "es"))
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
        geneSets <- attributes(gsvaExprData)$geneSets
        attributes(gsvaExprData)$geneSets <- NULL
        se <- SummarizedExperiment(assays=list(dummy=gsvaExprData))
        if (!is.null(annot))
            gsvaAnnotation(se) <- annot
        ## 'se' holds only the rows or columns of the chunk given by 'first'
        ## and 'last', so wrapData() should not subset it by them, and the
        ## 'restrict' metadata is added afterwards, as wrapData() would do
        if (assay == "es" && !is.null(geneSets))
            se <- wrapData(se, gsvaExprData, param, assay, first=NA_real_,
                           last=NA_real_, rem=rem, whdim=whdim,
                           dropAssays=TRUE, geneSets=geneSets)
        else
            se <- wrapData(se, gsvaExprData, param, assay, first=NA_real_,
                           last=NA_real_, rem=rem, whdim=whdim,
                           dropAssays=TRUE)
        if (!is.na(first) || !is.na(last))
            metadata(se)$restrict <- list(first=first, last=last, rem=rem,
                                          whdim=whdim)
        ## keep the version of GSVA that produced the object, instead of the
        ## one that is saving it, which wrapData() records
        if (!is.null(attributes(gsvaExprData)$gsvaVersion))
            metadata(se)$gsvaVersion <- attributes(gsvaExprData)$gsvaVersion

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


#' @importFrom cli cli_abort cli_alert_info
#' @importFrom rhdf5 H5Fis_hdf5 h5ls h5read
#' @importFrom HDF5Array HDF5Array H5SparseMatrix
#' @importFrom DelayedArray DelayedArray
#' @importFrom SummarizedExperiment assayNames "assay<-"
#'
#' @rdname gsva-serialization
#'
#' @export
loadHDF5GSVA <- function(file, assay="auto", verbose=TRUE) {
    if (!.isCharLength1(file))
        cli_abort(c("x"="'file' must be a single character string."))
    if (dir.exists(file))
        cli_abort(c("x"=paste("{.file {file}} is a directory, not a file",
                              "with GSVA output saved with 'saveHDF5GSVA()'.")))
    if (!file.exists(file))
        cli_abort(c("x"="The file {.file {file}} does not exist."))
    if (!.is_gsva_h5_file(file))
        cli_abort(c("x"=paste("The file {.file {file}} does not contain GSVA",
                              "output saved with 'saveHDF5GSVA()'.")))

    file <- normalizePath(file)
    if (verbose)
        cli_alert_info("Reading the metadata of {.file {file}}")
    gsvacontainer <- .restore_gsva_shell(unserialize(as.raw(h5read(file,
                                                       "gsva/shell"))), file)
    assayname <- as.character(h5read(file, "gsva/assayname"))
    layout <- as.character(h5read(file, "gsva/layout"))

    X <- if (layout == "csc") H5SparseMatrix(file, "gsva/assay") else
         HDF5Array(file, "gsva/assay")
    X <- DelayedArray(X)
    if (verbose) {
        cls <- class(gsvacontainer)[1]
        cli_alert_info(paste("Restored a {cls} object with {nrow(X)} rows",
                             "and {ncol(X)} columns"))
    }
    dimnames(X) <- dimnames(gsvacontainer)
    assay(gsvacontainer, assayname) <- X

    assay <- .check_assay_ranks_rnorm(assayNames(gsvacontainer), assay)

    .se_to_gsva_output(gsvacontainer, assay)
}

## the container of GSVA output 'x', without its assay, restored from the file
## 'file' by loadHDF5GSVA() or loadParquetGSVA(), updated to the current
## definition of its class, which may have changed since it was saved, as
## loadHDF5SummarizedExperiment() of the HDF5Array package does. its validity
## is not checked, because it lacks its assay, which is added afterwards
#' @importFrom BiocGenerics updateObject
#' @importFrom cli cli_abort
.restore_gsva_shell <- function(x, file) {
    x <- updateObject(x, check=FALSE)
    if (!is(x, "SummarizedExperiment"))
        cli_abort(c("x"=paste("The file {.file {file}} does not contain a",
                              "SummarizedExperiment object or derivative",
                              "with GSVA output.")))
    x
}

## whether 'file' is an HDF5 file with GSVA output saved with saveHDF5GSVA()
#' @importFrom rhdf5 H5Fis_hdf5 h5ls
.is_gsva_h5_file <- function(file) {
    if (!isTRUE(suppressWarnings(H5Fis_hdf5(file))))
        return(FALSE)
    content <- tryCatch(h5ls(file, recursive=2L), error=function(e) NULL)
    !is.null(content) &&
        all(c("shell", "assayname", "layout", "assay") %in%
            content$name[content$group == "/gsva"])
}


## turn a 'SummarizedExperiment' object with GSVA output, read from disk, back
## into the class of the object that was saved: if it was not originally a
## 'SummarizedExperiment', return the matrix in the assay 'assay' with the
## GSVA metadata stored as attributes, including the gene sets of GSVA scores
#' @importFrom cli cli_abort
#' @importFrom S4Vectors metadata
#' @importMethodsFrom SummarizedExperiment rowData
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
        geneSets <- NULL
        if (assay == "es" && !is.null(rowData(gsvacontainer)$gs)) {
            geneSets <- as.list(rowData(gsvacontainer)$gs)
            names(geneSets) <- rownames(gsvacontainer)
        }
        gsvaversion <- metadata(gsvacontainer)$gsvaVersion
        ranksnrow <- metadata(gsvacontainer)$ranksNrow
        gsvacontainer <- unwrapData(gsvacontainer, assay)
        attr(gsvacontainer, "gsvaParam") <- gsvapar
        attr(gsvacontainer, "assay") <- assay
        if (!is.null(annotation))
            attr(gsvacontainer, "geneIdType") <- annotation
        if (!is.null(restrict))
            attr(gsvacontainer, "restrict") <- restrict
        if (!is.null(geneSets))
            attr(gsvacontainer, "geneSets") <- geneSets
        if (!is.null(gsvaversion))
            attr(gsvacontainer, "gsvaVersion") <- gsvaversion
        if (!is.null(ranksnrow))
            attr(gsvacontainer, "ranksNrow") <- ranksnrow
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


#' @importFrom cli cli_abort cli_alert_info
#' @importFrom BiocGenerics path
#' @importFrom DelayedArray DelayedArray
#' @importFrom SummarizedExperiment assayNames "assay<-"
#'
#' @rdname gsva-serialization
#'
#' @export
loadParquetGSVA <- function(file, assay="auto", verbose=TRUE) {
    if (!.isCharLength1(file))
        cli_abort(c("x"="'file' must be a single character string."))

    if (verbose) {
        dpath <- .display_path(file)
        if (!.is_uri(file))
            cli_alert_info("Reading the metadata of {.file {dpath}}")
        else if (!.remote_parquet_info_cached(file))
            cli_alert_info("Downloading the metadata of {.val {dpath}}")
    }

    seed <- GsvaParquetSeed(file, verbose=verbose)
    md <- .parquet_file_info(path(seed), seed@backend)$metadata
    if (is.null(md$gsva_assay) || is.null(md$gsva_shell)) {
        dpath <- .display_path(file)
        cli_abort(c("x"=paste("The file {.file {dpath}} does not contain GSVA",
                              "output saved with 'saveParquetGSVA()'.")))
    }

    gsvacontainer <- .restore_gsva_shell(.decode_r_object(md$gsva_shell),
                                         .display_path(file))
    if (verbose) {
        cls <- class(gsvacontainer)[1]
        sze <- format(structure(nchar(md$gsva_shell, type="bytes"),
                                class="object_size"), units="auto")
        cli_alert_info(paste("Restored a {cls} object with {nrow(seed)} rows",
                             "and {ncol(seed)} columns from {sze} of metadata"))
    }
    X <- DelayedArray(seed)
    dimnames(X) <- dimnames(gsvacontainer)
    assay(gsvacontainer, md$gsva_assay) <- X

    assay <- .check_assay_ranks_rnorm(assayNames(gsvacontainer), assay)

    .se_to_gsva_output(gsvacontainer, assay)
}


## load GSVA output given as a path in the arguments 'rowNormExprData' of
## gsvaColRanks() or 'rankExprData' of gsvaColScores(), named in 'argname',
## which can be either a file with GSVA output saved with saveHDF5GSVA(), or
## a file, an 's3://' or 'gs://' URI, or an 'http://' or 'https://' URL with
## GSVA output saved with saveParquetGSVA(); 'assay' is the name of the
## GSVA assay to load
#' @importFrom cli cli_abort cli_alert_info
.load_gsva_path <- function(path, assay, argname, verbose) {
    if (length(path) != 1L || is.na(path))
        cli_abort(c("x"="'{argname}' must be a single character string."))

    if (.is_uri(path))
        parquet <- TRUE
    else if (file.exists(path) && !dir.exists(path)) {
        if (.is_gsva_h5_file(path))
            parquet <- FALSE
        else if (.is_parquet_file(path))
            parquet <- TRUE
        else
            cli_abort(c("x"=paste("{.file {path}} is neither a file with GSVA",
                                  "output saved with 'saveHDF5GSVA()', nor",
                                  "with 'saveParquetGSVA()'.")))
    } else
        cli_abort(c("x"="{path} cannot be found in the filesystem"))

    ## both loading functions give their own messages
    if (parquet)
        loadParquetGSVA(path, assay=assay, verbose=verbose)
    else
        loadHDF5GSVA(path, assay=assay, verbose=verbose)
}
