## ----- DelayedArray seed for GSVA output stored in Apache Parquet files -----
##
## GSVA intermediate output of row-normalized values from gsvaRowNorm() and
## ranks from gsvaColRanks() are consumed as blocks of columns by, respectively,
## gsvaColRanks() and gsvaColScores(). For this reason, it makes sense to use
## Apache Parquet as storage format when serializing these two outputs, since
## it is a format designed for fast and efficient storage and retrieval of
## column-oriented data. After some benchmarking, we found out that we should
## use a different layout for dense and sparse matrices, and that the best
## performance was provided by the following two layouts:
##
## * "dense": a table with one row per matrix entry and the columns 'col'
##   (int32) and 'value' (int32 or double), sorted by column and row, so that
##   the values of each matrix column are contiguous.
##
## * "sparse": a table with one row per matrix column and the columns 'col'
##   (int32), 'row' (list<int32>) and 'value' (list<int32> or list<double>),
##   which store the non-zero values of each matrix column and their row
##   indices. Within each column, the first row index is stored as such, and
##   the next ones as differences with respect to the previous row index,
##   which makes them compress much better.
##
## In both layouts, each row group stores whole matrix columns, and the
## file-level key-value metadata stores the following fields:
##
##   gsva_format          "GSVA-Parquet"
##   gsva_format_version  "1"
##   gsva_layout          "dense" or "sparse"
##   gsva_dim             "<nrow>,<ncol>"
##   gsva_type            "integer" or "double"
##   gsva_firsts          first column (1-based) of each row group, separated
##                        by commas
##
## The index in 'gsva_firsts' is necessary because arrow splits row groups
## larger than 2^20 rows, and its R interface does not give access to the
## metadata of each row group. Row and column names, and any other metadata,
## are stored separately by the functions that write these files.

.GSVA_PARQUET_FORMAT <- "GSVA-Parquet"
.GSVA_PARQUET_FORMAT_VERSIONS <- "1"

## maximum number of values (rows x columns) read with a single call to the
## arrow Parquet reader; reading many row groups at once needs memory for all
## their columns and, in the sparse layout, can overflow the 32-bit offsets
## of arrow list arrays beyond 2^31 values, while reading a single row group
## at a time adds an overhead that dominates when row groups are small
.GSVA_PARQUET_MAXBATCHVALS <- 2^27


#' `GsvaParquetSeed` class
#'
#' Internal class of `DelayedArray` seeds used by GSVA to read its output
#' stored in Apache Parquet format. It is not intended to be used directly.
#'
#' @return `dim()` returns the dimensions of the matrix; `type()` its type,
#' either `"integer"` or `"double"`; `path()` the path or URI to the Parquet
#' file; `is_sparse()` whether it uses the sparse layout; `chunkdim()` the
#' dimensions of the matrix chunks stored in each row group, or `NULL` when
#' they are irregular; `extract_array()` an ordinary matrix; and
#' `extract_sparse_array()` an `SVT_SparseMatrix` object.
#'
#' @aliases GsvaParquetSeed-class
#' @aliases dim,GsvaParquetSeed-method
#' @aliases type,GsvaParquetSeed-method
#' @aliases path,GsvaParquetSeed-method
#' @aliases is_sparse,GsvaParquetSeed-method
#' @aliases chunkdim,GsvaParquetSeed-method
#' @aliases extract_array,GsvaParquetSeed-method
#' @aliases extract_sparse_array,GsvaParquetSeed-method
#' @name GsvaParquetSeed-class
#' @rdname GsvaParquetSeed-class
#' @keywords internal
#'
#' @importClassesFrom S4Arrays Array
setClass("GsvaParquetSeed",
         contains="Array",
         slots=c(path="character",
                 dim="integer",
                 type="character",
                 layout="character",
                 firsts="integer",
                 backend="character"))

setValidity("GsvaParquetSeed", function(object) {
    msg <- NULL
    if (length(object@path) != 1L || is.na(object@path))
        msg <- c(msg, "'path' must be a single character string")
    if (length(object@dim) != 2L || anyNA(object@dim) || any(object@dim < 0L))
        msg <- c(msg, "'dim' must be two non-negative integer values")
    if (length(object@type) != 1L ||
        !object@type %in% c("integer", "double"))
        msg <- c(msg, "'type' must be either \"integer\" or \"double\"")
    if (length(object@layout) != 1L ||
        !object@layout %in% c("dense", "sparse"))
        msg <- c(msg, "'layout' must be either \"dense\" or \"sparse\"")
    if (length(object@backend) != 1L ||
        !object@backend %in% c("arrow", "duckdb"))
        msg <- c(msg, "'backend' must be either \"arrow\" or \"duckdb\"")
    ncol <- object@dim[2]
    firsts <- object@firsts
    if (length(ncol) == 1L && !is.na(ncol) && ncol > 0L &&
        (length(firsts) == 0L || anyNA(firsts) || firsts[1] != 1L ||
         any(diff(firsts) <= 0L) || firsts[length(firsts)] > ncol))
        msg <- c(msg, paste("'firsts' must be an increasing sequence of",
                            "column indices starting at 1"))
    if (is.null(msg)) TRUE else msg
})


## ----- arrow availability and Parquet readers -----

## because directly importing the arrow package into GSVA would add a heavy
## dependency for a functionality that is not always needed, we instead add
## it as a suggested package and import it only when the user calls the
## functions that read or write GSVA output in Parquet format.

.is_uri <- function(path) grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", path)

## credentials that a URI may carry, such as 'KEY:SECRET' in
## 's3://KEY:SECRET@bucket/path', or NULL if it carries none. Because AWS
## secret keys may contain '/', the user information before the first '@' is
## considered credentials when it contains ':', or no '/', except for the user
## name 'anonymous'
.uri_credentials <- function(path) {
    if (!.is_uri(path))
        return(NULL)
    rest <- sub("^[a-zA-Z][a-zA-Z0-9+.-]*://", "", path)
    at <- regexpr("@", rest, fixed=TRUE)
    if (at < 1L)
        return(NULL)
    userinfo <- substring(rest, 1L, at - 1L)
    if (userinfo == "anonymous" ||
        (grepl("/", userinfo, fixed=TRUE) && !grepl(":", userinfo, fixed=TRUE)))
        return(NULL)
    userinfo
}

## hide the credentials that a URI, or a message showing it, may carry
.hide_credentials <- function(text, path=text) {
    cred <- .uri_credentials(path)
    if (is.null(cred) || !nzchar(cred))
        return(text)
    gsub(cred, "<credentials>", text, fixed=TRUE)
}

## URI or path to show in messages
.display_path <- function(path) .hide_credentials(path)

#' @importFrom cli cli_abort
.require_arrow <- function(path=NULL) {
    if (!requireNamespace("arrow", quietly=TRUE)) {
        msg <- paste("The 'arrow' package is required to read and write GSVA",
                     "output in Apache Parquet format. Please install it",
                     "with 'install.packages(\"arrow\")'.")
        cli_abort(c("x"=msg))
    }

    if (!is.null(path) && .is_uri(path)) {
        dpath <- .display_path(path)
        scheme <- tolower(sub("://.*$", "", path))
        if (scheme == "s3" && !arrow::arrow_with_s3())
            cli_abort(c("x"=paste("The installed 'arrow' package was built",
                                  "without support for Amazon S3, which is",
                                  "necessary to access {.val {dpath}}.")))
        else if (scheme %in% c("gs", "gcs") && !arrow::arrow_with_gcs())
            cli_abort(c("x"=paste("The installed 'arrow' package was built",
                                  "without support for Google Cloud Storage,",
                                  "which is necessary to access",
                                  "{.val {dpath}}.")))
        else if (!scheme %in% c("s3", "gs", "gcs", "file"))
            cli_abort(c("x"=paste("URIs with scheme {.val {scheme}} are not",
                                  "supported; only 's3://', 'gs://', 'gcs://'",
                                  "'http://' and 'https://' URIs can be used",
                                  "to access remote GSVA output in Apache",
                                  "Parquet format. For local files use either a",
                                  "filesystem path or the 'file://' scheme.")))
    }

    invisible(TRUE)
}

## default limit, in seconds, of the retries of Google Cloud Storage requests
## that arrow uses for 'gs://' URIs without the 'retry_limit_seconds' query
## parameter; without it, when no Google Cloud credentials are available,
## e.g., outside Google Cloud, the Google Cloud C++ client may keep retrying
## to reach the Google Cloud metadata server for many minutes without any
## message, and with it, opening the file fails after about 20 seconds
.GSVA_GCS_RETRY_LIMIT_SECONDS <- 15

.is_gcs_uri <- function(path) grepl("^(gs|gcs)://", path, ignore.case=TRUE)

.is_s3_uri <- function(path) grepl("^s3://", path, ignore.case=TRUE)

## add the 'retry_limit_seconds' query parameter to a 'gs://' URI without it
.gcs_uri_with_retry_limit <- function(uri,
                                      seconds=.GSVA_GCS_RETRY_LIMIT_SECONDS) {
    if (!.is_gcs_uri(uri) || grepl("[?&]retry_limit_seconds=", uri))
        return(uri)
    paste0(uri, if (grepl("?", uri, fixed=TRUE)) "&" else "?",
           "retry_limit_seconds=", seconds)
}

#' @importFrom cli cli_abort
.open_parquet_reader <- function(path) {
    if (!.is_uri(path))
        return(arrow::ParquetFileReader$create(path))

    uri <- .gcs_uri_with_retry_limit(path)
    tryCatch({
        fsp <- arrow::FileSystem$from_uri(uri)
        f <- fsp$fs$OpenInputFile(fsp$path)
        arrow::ParquetFileReader$create(f)
    }, error=function(e) {
        ## arrow error messages may contain braces, which cli interpolates
        ## unless they are given through a variable
        errmsg <- sub("\\s*error_info=\\{.*$", "", conditionMessage(e))
        errmsg <- .hide_credentials(errmsg, path)
        dpath <- .display_path(path)
        msg <- c("x"="Cannot open {.val {dpath}}.",
                 "i"="{errmsg}")
        ## accessing a bucket without anonymous access requires credentials
        if (.is_gcs_uri(path) && !grepl("^[a-zA-Z]+://anonymous@", path))
            msg <- c(msg,
                     "i"=paste("If the bucket is public, use anonymous access",
                               "with a URI of the form",
                               "'gs://anonymous@<bucket>/<path>'. Otherwise,",
                               "set up Google Cloud credentials, e.g., with",
                               "'gcloud auth application-default login' or",
                               "the environment variable",
                               "'GOOGLE_APPLICATION_CREDENTIALS'."))
        else if (.is_s3_uri(path) &&
                 grepl(paste("ACCESS_DENIED|AccessDenied|InvalidAccessKeyId",
                             "SignatureDoesNotMatch", sep="|"), errmsg))
            msg <- c(msg,
                     "i"=paste("The bucket may be private, or your AWS",
                               "credentials may not give access to it. The",
                               "arrow package reads AWS credentials from the",
                               "environment variables 'AWS_ACCESS_KEY_ID' and",
                               "'AWS_SECRET_ACCESS_KEY' (and",
                               "'AWS_SESSION_TOKEN' for temporary",
                               "credentials), or from the file",
                               "'~/.aws/credentials', whose profile can be",
                               "selected with 'AWS_PROFILE'."))
        cli_abort(msg, call=NULL)
    })
}

## Parquet readers hold pointers to external (C++) objects that cannot be
## serialized to BiocParallel workers, and should not be shared with forked
## processes. For this reason, seeds only store the path to the file, and
## readers are opened on demand and cached by process. For local files, the
## cache key also includes the modification time and size of the file, so
## that a file overwritten in the same path is opened again.
.gsva_parquet_readers <- new.env(parent=emptyenv())

.parquet_reader <- function(path) {
    key <- paste(Sys.getpid(), path, sep="|")
    if (!.is_uri(path)) {
        finfo <- file.info(path, extra_cols=FALSE)
        if (is.na(finfo$size))
            cli_abort(c("x"="Cannot find the Parquet file {.file {path}}."))
        key <- paste(key, as.numeric(finfo$mtime), finfo$size, sep="|")
    }

    r <- .gsva_parquet_readers[[key]]
    if (is.null(r)) {
        ## drop readers of previous versions of the same file in this process
        stale <- startsWith(ls(.gsva_parquet_readers),
                            paste0(Sys.getpid(), "|", path, "|"))
        rm(list=ls(.gsva_parquet_readers)[stale], envir=.gsva_parquet_readers)
        r <- .open_parquet_reader(path)
        assign(key, r, envir=.gsva_parquet_readers)
    }

    r
}


## ----- constructor -----

#' @importFrom cli cli_abort
.parse_gsva_parquet_metadata <- function(md, path) {
    dpath <- .display_path(path)
    if (is.null(md$gsva_format) || md$gsva_format != .GSVA_PARQUET_FORMAT)
        cli_abort(c("x"=paste("The file {.file {dpath}} does not contain GSVA",
                              "output in Apache Parquet format.")))

    if (is.null(md$gsva_format_version) ||
        !md$gsva_format_version %in% .GSVA_PARQUET_FORMAT_VERSIONS)
        cli_abort(c("x"=paste("The file {.file {dpath}} stores GSVA output in",
                              "a version of the GSVA Parquet format",
                              "({md$gsva_format_version}) that this version",
                              "of GSVA cannot read.")))

    fields <- c("gsva_layout", "gsva_dim", "gsva_type", "gsva_firsts")
    missingfields <- fields[!fields %in% names(md)]
    if (length(missingfields) > 0)
        cli_abort(c("x"=paste("The file {.file {dpath}} lacks the metadata",
                              "field{?s} {.val {missingfields}}.")))

    dim <- suppressWarnings(as.integer(strsplit(md$gsva_dim, ",")[[1]]))
    firsts <- suppressWarnings(as.integer(strsplit(md$gsva_firsts, ",")[[1]]))

    list(layout=md$gsva_layout, dim=dim, type=md$gsva_type, firsts=firsts)
}

## TRUE when the metadata of the remote Parquet file in 'path' has already
## been downloaded in this process
.remote_parquet_info_cached <- function(path) {
    if (.is_http_url(path))
        return(.duckdb_file_info_cached(path))
    !is.null(.gsva_parquet_readers[[paste(Sys.getpid(), path, sep="|")]])
}

## key-value metadata and number of rows, (leaf) columns and row groups of the
## Parquet file in 'path', read with the given backend
.parquet_file_info <- function(path, backend, verbose=FALSE) {
    if (backend == "duckdb")
        return(.duckdb_file_info(path, verbose=verbose))

    r <- .parquet_reader(path)
    list(metadata=r$GetSchema()$metadata, num_rows=r$num_rows,
         num_columns=r$num_columns, num_row_groups=r$num_row_groups)
}

## 'path' is either the path to a local file, an 's3://' or 'gs://' URI, read
## with arrow, or an 'http://' or 'https://' URL, read with DuckDB; 'backend'
## allows one to choose the backend, internally used only for testing
#' @importFrom cli cli_abort
GsvaParquetSeed <- function(path, backend=NULL, verbose=FALSE) {
    if (!is.character(path) || length(path) != 1L || is.na(path))
        cli_abort(c("x"="'path' must be a single character string."))

    if (is.null(backend))
        backend <- if (.is_http_url(path)) "duckdb" else "arrow"
    if (backend == "duckdb")
        .require_duckdb()
    else
        .require_arrow(path)

    if (!.is_uri(path)) {
        if (!file.exists(path))
            cli_abort(c("x"="Cannot find the Parquet file {.file {path}}."))
        path <- normalizePath(path, mustWork=TRUE)
    }

    info <- .parquet_file_info(path, backend, verbose=verbose)
    dpath <- .display_path(path)
    md <- .parse_gsva_parquet_metadata(info$metadata, path)

    seed <- new("GsvaParquetSeed", path=path, dim=md$dim, type=md$type,
                layout=md$layout, firsts=md$firsts, backend=backend)

    ## check that the index of row groups corresponds to the file contents
    nrows <- if (seed@layout == "dense") prod(as.numeric(seed@dim))
             else seed@dim[2]
    ncols <- if (seed@layout == "dense") 2L else 3L
    if (info$num_rows != nrows || info$num_columns != ncols ||
        (seed@dim[2] > 0L && info$num_row_groups != length(seed@firsts)))
        cli_abort(c("x"=paste("The contents of the file {.file {dpath}} do",
                              "not correspond to its GSVA metadata.")))

    seed
}


## ----- reading columns -----

## columns stored in the row groups 'gs' (1-based), in increasing order
.parquet_seed_cols <- function(x, gs) {
    lasts <- c(x@firsts[-1] - 1L, x@dim[2])
    unlist(mapply(`:`, x@firsts[gs], lasts[gs], SIMPLIFY=FALSE),
           use.names=FALSE)
}

## read the matrix columns 'js', a strictly increasing vector of column
## indices, in batches of consecutive row groups with at most
## '.GSVA_PARQUET_MAXBATCHVALS' values (rows x columns), keeping only the
## requested columns of each batch before reading the next one. The function
## 'readbatch(gs, cols)' reads the row groups 'gs' (1-based) holding the
## columns 'cols', and returns them as an ordinary matrix
## or as an 'SVT_SparseMatrix' object.
#' @importFrom BiocGenerics cbind
.read_parquet_seed_cols <- function(x, js, readbatch) {
    g <- findInterval(js, x@firsts)
    ug <- unique(g)
    nvals <- as.numeric(x@dim[1]) * (c(x@firsts[-1] - 1L, x@dim[2])[ug] -
                                     x@firsts[ug] + 1)
    batch <- integer(length(ug))
    b <- 1L
    acc <- 0
    for (k in seq_along(ug)) {
        if (acc > 0 && acc + nvals[k] > .GSVA_PARQUET_MAXBATCHVALS) {
            b <- b + 1L
            acc <- 0
        }
        batch[k] <- b
        acc <- acc + nvals[k]
    }

    if (b == 1L) { ## single batch, avoid copying it into another object
        cols <- .parquet_seed_cols(x, ug)
        m <- readbatch(ug, cols)
        if (!identical(js, cols))
            m <- m[, match(js, cols), drop=FALSE]
        return(m)
    }

    sparse <- x@layout == "sparse"
    if (sparse)
        res <- vector("list", b)
    else {
        res <- matrix(vector(x@type, 1L), nrow=x@dim[1], ncol=length(js))
        nfilled <- 0L
    }
    for (bi in seq_len(b)) {
        gs <- ug[batch == bi]
        cols <- .parquet_seed_cols(x, gs)
        m <- readbatch(gs, cols)
        m <- m[, match(js[g %in% gs], cols), drop=FALSE]
        if (sparse)
            res[[bi]] <- m
        else {
            res[, nfilled + seq_len(ncol(m))] <- m
            nfilled <- nfilled + ncol(m)
        }
    }
    if (sparse)
        res <- do.call(cbind, res)

    res
}

## values of the row groups 'gs' (1-based), holding the matrix columns 'cols',
## of a seed with dense layout
.fetch_dense_values <- function(x, gs, cols) {
    if (x@backend == "duckdb")
        return(.duckdb_read_cols(x, gs, "value")$value)

    tb <- .parquet_reader(x@path)$ReadRowGroups(gs - 1L,
                                                column_indices=1L) ## no 'col'
    as.vector(tb$value)
}

## number of non-zero values of each column, and row indices (stored as
## differences) and non-zero values, of the row groups 'gs' (1-based),
## holding the matrix columns 'cols', of a seed with sparse layout
.fetch_sparse_values <- function(x, gs, cols) {
    if (x@backend == "duckdb") {
        res <- .duckdb_read_cols(x, gs, "col, row, value")
        if (is.unsorted(res$col))
            res <- res[order(res$col), ]
        return(list(lens=lengths(res$row),
                    rows=unlist(res$row, use.names=FALSE),
                    vals=unlist(res$value, use.names=FALSE)))
    }

    tb <- .parquet_reader(x@path)$ReadRowGroups(gs - 1L,
                                                column_indices=1:2) ## no 'col'
    list(lens=as.vector(arrow::call_function("list_value_length", tb$row)),
         rows=as.vector(arrow::call_function("list_flatten", tb$row)),
         vals=as.vector(arrow::call_function("list_flatten", tb$value)))
}

#' @importFrom cli cli_abort
.read_dense_batch <- function(x) {
    function(gs, cols) {
        v <- .fetch_dense_values(x, gs, cols)
        if (length(v) != as.numeric(x@dim[1]) * length(cols)) {
            dpath <- .display_path(x@path)
            cli_abort(c("x"=paste("Unexpected number of values read from",
                                  "the file {.file {dpath}}.")))
        }
        storage.mode(v) <- x@type
        matrix(v, nrow=x@dim[1])
    }
}

## decode row indices stored as differences within each column, where 'lens'
## gives the number of non-zero values in each column
.decode_row_deltas <- function(d, lens) {
    if (length(d) == 0L)
        return(integer(0))
    cs <- cumsum(if (sum(as.numeric(d)) < .Machine$integer.max) d
                 else as.numeric(d))
    prevend <- c(0L, cumsum(lens)[-length(lens)])
    offset <- numeric(length(lens))
    offset[prevend > 0L] <- cs[prevend[prevend > 0L]]
    as.integer(cs - rep(offset, lens))
}

#' @importFrom cli cli_abort
#' @importFrom BiocGenerics "type<-"
.read_sparse_batch <- function(x) {
    function(gs, cols) {
        sv <- .fetch_sparse_values(x, gs, cols)
        lens <- sv$lens
        rows <- sv$rows
        vals <- sv$vals
        if (length(lens) != length(cols) || length(rows) != sum(lens) ||
            length(vals) != length(rows)) {
            dpath <- .display_path(x@path)
            cli_abort(c("x"=paste("Unexpected number of values read from",
                                  "the file {.file {dpath}}.")))
        }
        m <- new("dgCMatrix", Dim=c(x@dim[1], length(cols)),
                 i=.decode_row_deltas(rows, lens) - 1L,
                 p=c(0L, cumsum(lens)), x=as.double(vals))
        m <- as(m, "SVT_SparseMatrix")
        if (x@type != "double")
            type(m) <- x@type
        m
    }
}

## read the columns in 'j', which can be NULL (all columns), unsorted and
## contain duplicates, and then the rows in 'i', which can be NULL (all rows)
.extract_parquet_seed <- function(x, i, j) {
    if (is.null(j))
        j <- seq_len(x@dim[2])
    nrows <- if (is.null(i)) x@dim[1] else length(i)
    sparse <- x@layout == "sparse"

    if (length(j) == 0L || x@dim[1] == 0L) {
        m <- matrix(vector(x@type, 0L), nrow=nrows, ncol=length(j))
        if (sparse)
            m <- as(m, "SVT_SparseMatrix")
        return(m)
    }

    uj <- sort(unique(as.integer(j)))
    readbatch <- if (sparse) .read_sparse_batch(x) else .read_dense_batch(x)
    m <- .read_parquet_seed_cols(x, uj, readbatch)
    if (!identical(as.integer(j), uj))
        m <- m[, match(j, uj), drop=FALSE]
    if (!is.null(i))
        m <- m[i, , drop=FALSE]

    m
}


## ----- methods -----

setMethod("dim", "GsvaParquetSeed", function(x) x@dim)

#' @importFrom BiocGenerics type
setMethod("type", "GsvaParquetSeed", function(x) x@type)

#' @importFrom BiocGenerics path
setMethod("path", "GsvaParquetSeed", function(object, ...) object@path)

#' @importFrom S4Arrays is_sparse
setMethod("is_sparse", "GsvaParquetSeed", function(x) x@layout == "sparse")

## the row groups define a regular grid of chunks when they all store the
## same number of columns, except for the last one, which may store fewer
#' @importFrom DelayedArray chunkdim
setMethod("chunkdim", "GsvaParquetSeed", function(x) {
    if (x@dim[2] == 0L)
        return(NULL)
    w <- diff(c(x@firsts, x@dim[2] + 1L))
    if (any(w[-length(w)] != w[1]) || w[length(w)] > w[1])
        return(NULL)
    c(x@dim[1], w[1])
})

#' @importFrom S4Arrays extract_array
setMethod("extract_array", "GsvaParquetSeed", function(x, index) {
    m <- .extract_parquet_seed(x, index[[1]], index[[2]])
    if (is(m, "SVT_SparseMatrix"))
        m <- as.matrix(m)
    m
})

#' @importFrom SparseArray extract_sparse_array
setMethod("extract_sparse_array", "GsvaParquetSeed", function(x, index) {
    m <- .extract_parquet_seed(x, index[[1]], index[[2]])
    if (!is(m, "SVT_SparseMatrix")) {
        m <- as(m, "SVT_SparseMatrix")
        if (x@type != "double")
            type(m) <- x@type
    }
    m
})

## wrap a GSVA Parquet file into a 'DelayedMatrix' object
#' @importFrom DelayedArray DelayedArray
.GsvaParquetMatrix <- function(path) DelayedArray(GsvaParquetSeed(path))

## TRUE when the local file in 'path' starts with the magic number of the
## Apache Parquet format
.is_parquet_file <- function(path) {
    con <- file(path, open="rb")
    on.exit(close(con))
    identical(readBin(con, what="raw", n=4L), charToRaw("PAR1"))
}


## ----- writing GSVA Parquet files -----

## maximum number of rows that arrow writes in a single row group
.GSVA_PARQUET_MAXRGROWS <- 2^20

## default number of columns per row group in the sparse layout, where each
## row stores one matrix column; with fewer columns, reading a block of
## columns needs more calls to the reader, and with more columns, reading a
## few columns reads many more of them than necessary
.GSVA_PARQUET_SPARSE_COLS_PER_RGROUP <- 500L

## maximum size of the encoded R objects stored in the key-value metadata of
## the file footer; arrow fails to read footers with metadata strings larger
## than about 64-95 MB
.GSVA_PARQUET_MAXMDSIZE <- 64e6

## base64 encoding of raw vectors, since the key-value metadata of Parquet
## files written with arrow from R can only hold character strings
.B64CHARS <- utf8ToInt(paste0(c(LETTERS, letters, 0:9, "+", "/"),
                              collapse=""))

.base64_encode <- function(r) {
    pad <- (3L - length(r) %% 3L) %% 3L
    x <- matrix(as.integer(c(r, as.raw(rep(0L, pad)))), nrow=3L)
    v <- x[1L, ] * 65536L + x[2L, ] * 256L + x[3L, ]
    codes <- .B64CHARS[rbind(v %/% 262144L, (v %/% 4096L) %% 64L,
                             (v %/% 64L) %% 64L, v %% 64L) + 1L]
    if (pad > 0L)
        codes[length(codes) - seq_len(pad) + 1L] <- utf8ToInt("=")
    intToUtf8(codes)
}

#' @importFrom cli cli_abort
.base64_decode <- function(s) {
    codes <- utf8ToInt(s)
    pad <- sum(codes[length(codes) - 0:1] == utf8ToInt("="))
    codes[codes == utf8ToInt("=")] <- utf8ToInt("A")
    dec <- rep(NA_integer_, 256L)
    dec[.B64CHARS + 1L] <- 0:63
    v <- dec[codes + 1L]
    if (anyNA(v) || length(v) %% 4L != 0L)
        cli_abort(c("x"="Invalid base64-encoded metadata."))
    v <- matrix(v, nrow=4L)
    w <- v[1L, ] * 262144L + v[2L, ] * 4096L + v[3L, ] * 64L + v[4L, ]
    r <- as.raw(rbind(w %/% 65536L, (w %/% 256L) %% 256L, w %% 256L))
    r[seq_len(length(r) - pad)]
}

## serialize and compress an R object into a character string, and back
## this approach ties the format used in Parquet files to R and therefore
## they can only be read back by R. in the future, we may use a different
## approach that enables interoperability with other programming languages
.encode_r_object <- function(x)
    .base64_encode(memCompress(serialize(x, connection=NULL), type="gzip"))

.decode_r_object <- function(s)
    unserialize(memDecompress(.base64_decode(s), type="gzip"))

## number of columns per row group
#' @importFrom cli cli_abort
.parquet_cols_per_rgroup <- function(colsPerRowGroup, nrow, ncol, sparse) {
    maxdense <- .GSVA_PARQUET_MAXRGROWS %/% nrow
    if (identical(colsPerRowGroup, "auto"))
        K <- if (sparse) .GSVA_PARQUET_SPARSE_COLS_PER_RGROUP else maxdense
    else {
        if (!.isNumLength1(colsPerRowGroup) || !is.finite(colsPerRowGroup) ||
            length(colsPerRowGroup) != 1L || colsPerRowGroup < 1 ||
            colsPerRowGroup != round(colsPerRowGroup))
            cli_abort(c("x"=paste("'colsPerRowGroup' must be either \"auto\"",
                                  "or a positive integer number.")))
        ## a row group cannot store more columns than the matrix has, and
        ## capping 'colsPerRowGroup' to them before converting it into an
        ## integer also keeps it within the range of R integers
        K <- as.integer(min(colsPerRowGroup, max(ncol, 1L)))
        if (!sparse && K > maxdense)
            cli_abort(c("x"=paste("With {nrow} rows, 'colsPerRowGroup' can be",
                                  "at most {maxdense} to store dense values,",
                                  "because row groups cannot have more than",
                                  "2^20 values.")))
    }
    if (!sparse && K < 1L)
        cli_abort(c("x"=paste("Dense matrices with more than 2^20 rows cannot",
                              "be stored in Apache Parquet format.")))
    as.integer(min(K, .Machine$integer.max, max(ncol, 1L)))
}

## table with the columns 'js' of the matrix block 'm', with column indices
## 'cols' in the whole matrix, in the dense or sparse layout
#' @importFrom cli cli_abort
.parquet_rgroup_table <- function(m, js, cols, type, sch, sparse) {
    vt <- if (type == "integer") arrow::int32() else arrow::float64()
    if (!sparse) {
        v <- as.vector(m[, js, drop=FALSE])
        if (type == "integer" && is.double(v) &&
            any(v != round(v), na.rm=TRUE))
            cli_abort(c("x"="Non-integer values cannot be stored as integers."))
        storage.mode(v) <- type
        return(arrow::Table$create(col=rep(cols, each=nrow(m)), value=v,
                                   schema=sch))
    }

    m <- m[, js, drop=FALSE]
    lens <- diff(m@p)
    ## row indices as differences within each column, except the first one
    i1 <- m@i + 1L
    d <- i1
    if (length(i1) > 1L)
        d[-1L] <- diff(i1)
    starts <- m@p[-length(m@p)][lens > 0L] + 1L
    d[starts] <- i1[starts]
    x <- m@x
    if (type == "integer") {
        if (any(x != round(x), na.rm=TRUE))
            cli_abort(c("x"="Non-integer values cannot be stored as integers."))
        x <- as.integer(x)
    }
    f <- factor(rep(seq_along(lens), lens), levels=seq_along(lens))
    arrow::Table$create(col=cols,
                        row=arrow::Array$create(unname(split(d, f)),
                                                type=arrow::list_of(arrow::int32())),
                        value=arrow::Array$create(unname(split(x, f)),
                                                  type=arrow::list_of(vt)),
                        schema=sch)
}

## write the matrix 'X' to the Parquet file 'file' as values of type 'type',
## in the sparse layout when 'X' is sparse, and in the dense one otherwise,
## storing the character strings in the named list 'md' as additional
## key-value metadata. The file is first written with a temporary name and
## then renamed, so that an existing file is only replaced when writing
## succeeds, and it can be read while it is being replaced.
#' @importFrom cli cli_abort cli_progress_bar cli_progress_update
#' @importFrom cli cli_progress_done
#' @importFrom S4Arrays is_sparse
#' @importFrom DelayedArray getAutoBlockSize
.write_gsva_parquet <- function(X, file, type, colsPerRowGroup="auto",
                                md=list(), verbose=FALSE) {
    nr <- nrow(X)
    nc <- ncol(X)
    if (nr == 0L || nc == 0L)
        cli_abort(c("x"=paste("Cannot store in Apache Parquet format a matrix",
                              "without rows or columns.")))
    sparse <- is_sparse(X)
    K <- .parquet_cols_per_rgroup(colsPerRowGroup, nr, nc, sparse)
    firsts <- seq(1L, nc, by=K)

    vt <- if (type == "integer") arrow::int32() else arrow::float64()
    if (sparse) {
        sch <- arrow::schema(col=arrow::int32(),
                             row=arrow::list_of(arrow::int32()),
                             value=arrow::list_of(vt))
        leaves <- c("col", "row.list.element", "value.list.element")
        dict <- c(TRUE, TRUE, type == "integer")
    } else {
        sch <- arrow::schema(col=arrow::int32(), value=vt)
        leaves <- c("col", "value")
        dict <- c(TRUE, type == "integer")
    }
    md <- c(list(gsva_format=.GSVA_PARQUET_FORMAT,
                 gsva_format_version=.GSVA_PARQUET_FORMAT_VERSIONS[1],
                 gsva_layout=if (sparse) "sparse" else "dense",
                 gsva_dim=paste(nr, nc, sep=","),
                 gsva_type=type,
                 gsva_firsts=paste(firsts, collapse=",")), md)
    if (any(nchar(unlist(md), type="bytes") > .GSVA_PARQUET_MAXMDSIZE))
        cli_abort(c("x"=paste("The metadata of the GSVA output is too large to",
                              "be stored in Apache Parquet format. Consider",
                              "removing large objects from it, such as images",
                              "in a 'SpatialExperiment' object, or using",
                              "'saveHDF5GSVA()' instead.")))
    sch <- sch$WithMetadata(md)
    ## per-column properties of list columns must refer to their leaf columns
    props <- arrow::ParquetWriterProperties$create(leaves, compression="zstd",
                                                   use_dictionary=dict)

    tmpfile <- tempfile(pattern=paste0(".", basename(file), "."),
                        tmpdir=dirname(file))
    sink <- arrow::FileOutputStream$create(tmpfile)
    done <- FALSE
    on.exit({
        if (!done) {
            try(sink$close(), silent=TRUE)
            unlink(tmpfile)
        }
    })
    writer <- arrow::ParquetFileWriter$create(sch, sink, properties=props)

    ## read the input by blocks of whole row groups of at most the automatic
    ## block size, or a single row group if it does not fit in that size
    rgbytes <- as.numeric(nr) * K * 8
    rgsperblock <- max(1L, floor(getAutoBlockSize() / rgbytes))
    blockfirsts <- firsts[seq(1L, length(firsts), by=rgsperblock)]
    if (verbose)
        idpb <- cli_progress_bar("Writing columns", total=nc)
    for (b in seq_along(blockfirsts)) {
        bcols <- blockfirsts[b]:(if (b < length(blockfirsts))
                                     blockfirsts[b + 1L] - 1L else nc)
        m <- X[, bcols, drop=FALSE]
        m <- if (sparse) as(m, "dgCMatrix") else as.matrix(m)
        for (f in firsts[firsts >= bcols[1] & firsts <= bcols[length(bcols)]]) {
            cols <- f:min(f + K - 1L, nc)
            tb <- .parquet_rgroup_table(m, cols - bcols[1] + 1L, cols, type,
                                        sch, sparse)
            writer$WriteTable(tb, chunk_size=nrow(tb)) ## one row group
        }
        if (verbose)
            cli_progress_update(id=idpb, inc=length(bcols))
    }
    writer$Close()
    sink$close()
    if (verbose)
        cli_progress_done(id=idpb)

    if (!.replace_file(tmpfile, file)) {
        if (file.exists(file))
            cli_abort(c("x"="Cannot replace the file {.file {file}}.",
                        "i"=paste("The file may be open in another program",
                                  "or R session.")))
        cli_abort(c("x"="Cannot write the file {.file {file}}."))
    }
    done <- TRUE

    invisible(file)
}

## close the Parquet readers of the local file in 'path' cached in this
## process, releasing the memory mapping of the file, which on Windows
## prevents replacing it
.close_parquet_readers <- function(path) {
    path <- normalizePath(path, mustWork=FALSE)
    keys <- ls(.gsva_parquet_readers)
    cached <- startsWith(keys, paste0(Sys.getpid(), "|", path, "|"))
    if (any(cached)) {
        rm(list=keys[cached], envir=.gsva_parquet_readers)
        invisible(gc(verbose=FALSE)) ## run the finalizers of the readers
    }

    invisible(any(cached))
}

## move the file 'tmpfile' to 'file', replacing it if it exists, and keeping
## it if the replacement fails; returns TRUE if 'file' has been replaced.
## file.rename() replaces existing files where permissions allow it, but on
## Windows it fails when the file is open, e.g., by a Parquet reader of this
## R session, which is closed first, or by another program. When renaming
## over the existing file fails, it is first moved aside and, if the new file
## cannot then take its place, moved back. 'rename' allows one to simulate
## failures, internally used only for testing
.replace_file <- function(tmpfile, file, rename=file.rename) {
    rename <- match.fun(rename)
    mv <- function(from, to) isTRUE(suppressWarnings(rename(from, to)))

    if (file.exists(file))
        .close_parquet_readers(file)
    if (mv(tmpfile, file))
        return(TRUE)
    if (!file.exists(file))
        return(FALSE)

    bakfile <- tempfile(pattern=paste0(".", basename(file), ".bak."),
                        tmpdir=dirname(file))
    if (!mv(file, bakfile))
        return(FALSE)
    if (mv(tmpfile, file)) {
        unlink(bakfile)
        return(TRUE)
    }
    mv(bakfile, file) ## restore the existing file

    FALSE
}


## TRUE when 'x' is a 'DelayedArray' object with data from a GSVA Parquet file,
## also after delayed operations such as subsetting
#' @importFrom DelayedArray seedApply
.is_parquet_backed <- function(x) {
    is(x, "DelayedArray") &&
        any(unlist(seedApply(x, is, "GsvaParquetSeed"), use.names=FALSE))
}
