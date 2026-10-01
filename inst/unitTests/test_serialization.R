test_ranksserialization <- function() {

    message("Running unit tests for ranks serialization")

    suppressPackageStartupMessages({
        library(Matrix)
        library(GSEABase)
        library(SummarizedExperiment)
    })

    p <- 10 ## number of genes
    n <- 30 ## number of samples

    ## consider three disjoint gene sets
    gsets <- list(gset1=paste0("g", 1:3),
                  gset2=paste0("g", 4:6),
                  gset3=paste0("g", 7:10))

    ## build a random sparse count matrix with 85% sparsity
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.15)) ## 85% sparsity
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- matrix(cnt, nrow=p, ncol=n,
                  dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    cnt <- Matrix(cnt, sparse=TRUE)

    se <- SummarizedExperiment(assays=list(counts=cnt))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(se, gsets, verbose=FALSE)

    ## calculate row-normalized expression values
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## calculate GSVA column ranks
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)

    ## calculate GSVA scores
    es <- gsvaColScores(gsvacolranks, verbose=FALSE)
 
    ## check out that 'saveHDF5GSVA()' throws an error when the input
    ## object is not one of the classes of the union of classes defined
    ## by 'GsvaExprData'
    rnormdir <- tempfile()
    checkException(saveHDF5GSVA(gsvapar, rnormdir))

    ## save the GSVA row-normalized expression values to disk
    savedrnormdir <- saveHDF5GSVA(gsvarownorm, rnormdir)
    checkTrue(file.exists(savedrnormdir) && identical(savedrnormdir, rnormdir))
                                               
    ## save the GSVA rank values to disk
    ranksdir <- tempfile()
    savedranksdir <- saveHDF5GSVA(gsvacolranks, ranksdir)
    checkTrue(file.exists(savedranksdir) && identical(savedranksdir, ranksdir))
                                               
    ## load the GSVA row-normalized values from disk               
    loaded_gsvarownorm <- loadHDF5GSVA(savedrnormdir)      
                                               
    ## load the GSVA ranks from disk               
    loaded_gsvacolranks <- loadHDF5GSVA(savedranksdir)      
                                               
    ## check that the loaded row-normalized values provide the
    ## same ranks as the ones calculated from the original values
    gsvacolranks_from_loaded_gsvarownorm <- gsvaColRanks(loaded_gsvarownorm, verbose=FALSE)
    checkEqualsNumeric(assay(gsvacolranks, "gsvaranks"),
		       assay(gsvacolranks_from_loaded_gsvarownorm, "gsvaranks"))

    ## check that the loaded ranks provide the
    ## same scores as the original ranks
    es_from_loaded_gsvaranks <- gsvaColScores(loaded_gsvacolranks, verbose=FALSE)
    checkTrue(identical(es, es_from_loaded_gsvaranks))

    ## check it again saving and loading the ranks stored in a
    ## non-'SummarizedExperiment' object
    gsvapar <- gsvaParam(cnt, gsets, verbose=FALSE)
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    savedrnormdir <- saveHDF5GSVA(gsvarownorm, rnormdir, replace=TRUE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    savedranksdir <- saveHDF5GSVA(gsvacolranks, ranksdir, replace=TRUE)

    loaded_gsvarownorm <- loadHDF5GSVA(savedrnormdir)
    gsvacolranks_from_loaded_gsvarownorm <- gsvaColRanks(loaded_gsvarownorm, verbose=FALSE)
    checkEqualsNumeric(gsvacolranks, gsvacolranks_from_loaded_gsvarownorm)

    loaded_gsvacolranks <- loadHDF5GSVA(savedranksdir)
    es_from_loaded_gsvaranks <- gsvaColScores(loaded_gsvacolranks, verbose=FALSE)
    checkEqualsNumeric(assay(es), es_from_loaded_gsvaranks)
}

test_parquetserialization <- function() {

    if (!requireNamespace("arrow", quietly=TRUE)) {
        message("Skipping unit tests for Parquet serialization (no 'arrow')")
        return(invisible(TRUE))
    }

    message("Running unit tests for Parquet serialization")

    suppressPackageStartupMessages({
        library(Matrix)
        library(GSEABase)
        library(SummarizedExperiment)
        library(S4Arrays)
    })

    p <- 50 ## number of genes
    n <- 120 ## number of samples

    gsets <- list(gset1=paste0("g", 1:10),
                  gset2=paste0("g", 11:25),
                  gset3=paste0("g", 26:50))

    ## build a random sparse count matrix with 85% sparsity
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.15))
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- matrix(cnt, nrow=p, ncol=n,
                  dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    cnt <- Matrix(cnt, sparse=TRUE)
    se <- SummarizedExperiment(assays=list(counts=cnt))

    ## dense expression values
    y <- matrix(rnorm(n*p), nrow=p, ncol=n, dimnames=dimnames(cnt))

    getvals <- function(x, a)
        as.matrix(if (is(x, "SummarizedExperiment")) assay(x, a) else x)

    for (input in list(se, cnt, y)) {
        gsvapar <- gsvaParam(input, gsets, verbose=FALSE)
        gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
        gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
        es <- gsvaColScores(gsvacolranks, verbose=FALSE)

        rnormfile <- tempfile(fileext=".parquet")
        ranksfile <- tempfile(fileext=".parquet")
        checkIdentical(rnormfile, saveParquetGSVA(gsvarownorm, rnormfile))
        ## small row groups, to read columns from several of them
        checkIdentical(ranksfile, saveParquetGSVA(gsvacolranks, ranksfile,
                                                  colsPerRowGroup=7))

        loaded_gsvarownorm <- loadParquetGSVA(rnormfile)
        loaded_gsvacolranks <- loadParquetGSVA(ranksfile)
        checkTrue(is(loaded_gsvacolranks, class(gsvacolranks)[1]) ||
                  is(loaded_gsvacolranks, "DelayedMatrix"))

        ## loaded values are read from the files, keeping their sparsity
        ## and storing ranks as integers
        lrnks <- if (is(loaded_gsvacolranks, "SummarizedExperiment"))
                     assay(loaded_gsvacolranks, "gsvaranks")
                 else loaded_gsvacolranks
        checkTrue(GSVA:::.is_parquet_backed(lrnks))
        insparse <- is_sparse(if (is(input, "SummarizedExperiment"))
                                  assay(input) else input)
        checkIdentical(insparse, is_sparse(lrnks))
        checkIdentical("integer", type(lrnks))
        checkIdentical(dimnames(getvals(gsvacolranks, "gsvaranks")),
                       dimnames(lrnks))

        checkEqualsNumeric(getvals(gsvarownorm, "gsvarnorm"),
                           getvals(loaded_gsvarownorm, "gsvarnorm"))

        ## the loaded row-normalized values provide the same ranks
        gsvacolranks2 <- gsvaColRanks(loaded_gsvarownorm, verbose=FALSE)
        checkEqualsNumeric(getvals(gsvacolranks, "gsvaranks"),
                           getvals(gsvacolranks2, "gsvaranks"))

        ## the loaded ranks provide the same scores, also for other gene sets
        es2 <- gsvaColScores(loaded_gsvacolranks, verbose=FALSE)
        checkEqualsNumeric(getvals(es, "es"), getvals(es2, "es"))
        gsets2 <- list(gset4=paste0("g", c(1:5, 40:50)),
                       gset5=paste0("g", 20:35))
        es3 <- gsvaColScores(gsvacolranks, geneSets=gsets2, verbose=FALSE)
        es4 <- gsvaColScores(loaded_gsvacolranks, geneSets=gsets2,
                             verbose=FALSE)
        checkEqualsNumeric(getvals(es3, "es"), getvals(es4, "es"))

        ## existing files are only replaced when 'replace=TRUE'
        checkException(saveParquetGSVA(gsvacolranks, ranksfile), silent=TRUE)
        saveParquetGSVA(gsvarownorm, ranksfile, replace=TRUE)
        checkEqualsNumeric(getvals(gsvarownorm, "gsvarnorm"),
                           getvals(loadParquetGSVA(ranksfile), "gsvarnorm"))

        unlink(c(rnormfile, ranksfile))
    }

    ## errors
    f <- tempfile(fileext=".parquet")
    checkException(saveParquetGSVA(gsvapar, f), silent=TRUE)
    checkException(saveParquetGSVA(gsvacolranks, f, colsPerRowGroup=0),
                   silent=TRUE)
    checkException(saveParquetGSVA(gsvacolranks, f,
                                   colsPerRowGroup=2^20), silent=TRUE)
    checkException(loadParquetGSVA(f), silent=TRUE)
    arrow::write_parquet(data.frame(a=1:3), f)
    checkException(loadParquetGSVA(f), silent=TRUE)
    unlink(f)
}

test_serializationpaths <- function() {

    message("Running unit tests for GSVA output given as paths")

    suppressPackageStartupMessages({
        library(Matrix)
        library(SummarizedExperiment)
    })

    p <- 40 ## number of genes
    n <- 60 ## number of samples
    gsets <- list(gset1=paste0("g", 1:10),
                  gset2=paste0("g", 11:25),
                  gset3=paste0("g", 26:40))
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))

    getvals <- function(x, a)
        as.matrix(if (is(x, "SummarizedExperiment")) assay(x, a) else x)

    formats <- "HDF5"
    if (requireNamespace("arrow", quietly=TRUE))
        formats <- c(formats, "Parquet")

    for (input in list(y, SummarizedExperiment(assays=list(exprs=y)))) {
        gsvarownorm <- gsvaRowNorm(gsvaParam(input, gsets, verbose=FALSE),
                                   verbose=FALSE)
        gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
        es <- gsvaColScores(gsvacolranks, verbose=FALSE)

        for (fmt in formats) {
            if (fmt == "HDF5") {
                rnormpath <- saveHDF5GSVA(gsvarownorm, tempfile())
                rankspath <- saveHDF5GSVA(gsvacolranks, tempfile())
            } else {
                rnormpath <- saveParquetGSVA(gsvarownorm,
                                             tempfile(fileext=".parquet"))
                rankspath <- saveParquetGSVA(gsvacolranks,
                                             tempfile(fileext=".parquet"))
            }

            gsvacolranks2 <- gsvaColRanks(rnormpath, verbose=FALSE)
            checkEqualsNumeric(getvals(gsvacolranks, "gsvaranks"),
                               getvals(gsvacolranks2, "gsvaranks"))
            es2 <- gsvaColScores(rankspath, verbose=FALSE)
            checkEqualsNumeric(getvals(es, "es"), getvals(es2, "es"))

            ## the saved data must be of the kind expected by each function
            checkException(gsvaColRanks(rankspath, verbose=FALSE),
                           silent=TRUE)
            checkException(gsvaColScores(rnormpath, verbose=FALSE),
                           silent=TRUE)

            unlink(c(rnormpath, rankspath), recursive=TRUE)
        }
    }

    ## paths that cannot be loaded
    checkException(gsvaColRanks(tempfile(), verbose=FALSE), silent=TRUE)
    checkException(gsvaColScores(c(tempfile(), tempfile()), verbose=FALSE),
                   silent=TRUE)
    f <- tempfile()
    writeLines("not GSVA output", f)
    checkException(gsvaColScores(f, verbose=FALSE), silent=TRUE)
    unlink(f)
}

test_serializationscoresandeset <- function() {

    message("Running unit tests for serialization of scores and ExpressionSet")

    suppressPackageStartupMessages({
        library(Biobase)
        library(SummarizedExperiment)
    })

    p <- 40 ## number of genes
    n <- 60 ## number of samples
    gsets <- list(gset1=paste0("g", 1:10),
                  gset2=paste0("g", 11:25),
                  gset3=paste0("g", 26:40))
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))

    getvals <- function(x, a) {
        if (is(x, "SummarizedExperiment"))
            x <- assay(x, a)
        else if (is(x, "ExpressionSet"))
            x <- exprs(x)
        as.matrix(x)
    }

    formats <- "HDF5"
    if (requireNamespace("arrow", quietly=TRUE))
        formats <- c(formats, "Parquet")

    inputs <- list(y, SummarizedExperiment(assays=list(exprs=y)),
                   ExpressionSet(y))
    for (input in inputs) {
        gsvarownorm <- gsvaRowNorm(gsvaParam(input, gsets, verbose=FALSE),
                                   verbose=FALSE)
        gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
        es <- gsvaColScores(gsvacolranks, verbose=FALSE)

        for (fmt in formats) {
            savefun <- function(x) {
                if (fmt == "HDF5")
                    saveHDF5GSVA(x, tempfile())
                else
                    saveParquetGSVA(x, tempfile(fileext=".parquet"))
            }
            loadfun <- if (fmt == "HDF5") loadHDF5GSVA else loadParquetGSVA

            ## GSVA scores, keeping their gene sets
            espath <- savefun(es)
            loaded_es <- loadfun(espath)
            checkEqualsNumeric(getvals(es, "es"), getvals(loaded_es, "es"))
            if (is(input, "SummarizedExperiment")) {
                checkTrue(is(loaded_es, "SummarizedExperiment"))
                checkIdentical(as.list(rowData(es)$gs),
                               as.list(rowData(loaded_es)$gs))
            } else {
                ## an 'ExpressionSet' is loaded as a matrix
                checkTrue(is(loaded_es, "DelayedMatrix"))
                checkIdentical("es", attr(loaded_es, "assay", exact=TRUE))
                checkIdentical(attr(es, "geneSets"),
                               attr(loaded_es, "geneSets"))
            }

            ## row-normalized values and ranks
            rnormpath <- savefun(gsvarownorm)
            rankspath <- savefun(gsvacolranks)
            gsvacolranks2 <- gsvaColRanks(loadfun(rnormpath), verbose=FALSE)
            checkEqualsNumeric(getvals(gsvacolranks, "gsvaranks"),
                               getvals(gsvacolranks2, "gsvaranks"))
            es2 <- gsvaColScores(loadfun(rankspath), verbose=FALSE)
            checkEqualsNumeric(getvals(es, "es"), getvals(es2, "es"))

            unlink(c(espath, rnormpath, rankspath), recursive=TRUE)
        }
    }
}

test_gcsuriretrylimit <- function() {

    message("Running unit tests for the retry limit of 'gs://' URIs")

    addlimit <- GSVA:::.gcs_uri_with_retry_limit

    ## 'gs://' URIs without a retry limit get the default one
    checkIdentical("gs://bucket/x.parquet?retry_limit_seconds=15",
                   addlimit("gs://bucket/x.parquet", 15))
    checkIdentical("gs://anonymous@bucket/x.parquet?retry_limit_seconds=15",
                   addlimit("gs://anonymous@bucket/x.parquet", 15))
    checkIdentical("gs://bucket/x.parquet?scheme=http&retry_limit_seconds=15",
                   addlimit("gs://bucket/x.parquet?scheme=http", 15))

    ## a retry limit given in the URI is kept, and other URIs are not changed
    checkIdentical("gs://bucket/x.parquet?retry_limit_seconds=60",
                   addlimit("gs://bucket/x.parquet?retry_limit_seconds=60", 15))
    checkIdentical("s3://bucket/x.parquet",
                   addlimit("s3://bucket/x.parquet", 15))
    checkIdentical("/tmp/x.parquet", addlimit("/tmp/x.parquet", 15))
}

test_uricredentials <- function() {

    message("Running unit tests for hiding credentials in URIs")

    cred <- GSVA:::.uri_credentials
    disp <- GSVA:::.display_path
    hide <- GSVA:::.hide_credentials

    ## credentials, also with secrets containing '/' or '+'
    checkIdentical("KEY:SECRET", cred("s3://KEY:SECRET@bucket/x.parquet"))
    checkIdentical("s3://<credentials>@bucket/x.parquet",
                   disp("s3://KEY:SECRET@bucket/x.parquet"))
    checkIdentical("s3://<credentials>@bucket/x.parquet",
                   disp("s3://KEY:SE/CR+ET@bucket/x.parquet"))
    checkIdentical("gs://<credentials>@bucket/x.parquet",
                   disp("gs://user@bucket/x.parquet"))

    ## no credentials
    checkTrue(is.null(cred("s3://bucket/x.parquet")))
    checkTrue(is.null(cred("gs://anonymous@bucket/x.parquet")))
    checkTrue(is.null(cred("s3://bucket/dir@x/y.parquet")))
    checkTrue(is.null(cred("/tmp/x@y.parquet")))
    checkIdentical("gs://anonymous@bucket/x.parquet",
                   disp("gs://anonymous@bucket/x.parquet"))
    checkIdentical("/tmp/x@y.parquet", disp("/tmp/x@y.parquet"))

    ## credentials within a message showing the URI
    checkIdentical("Cannot parse 's3://<credentials>@bucket/x.parquet'",
                   hide("Cannot parse 's3://KEY:SECRET@bucket/x.parquet'",
                        "s3://KEY:SECRET@bucket/x.parquet"))
}

## the DuckDB backend, used to read GSVA output in Parquet format through
## HTTP(S), reads the same values as the arrow backend, which is checked here
## with local files, also readable with DuckDB without its extension 'httpfs'
test_parquetduckdbbackend <- function() {

    if (!requireNamespace("arrow", quietly=TRUE) ||
        !requireNamespace("duckdb", quietly=TRUE) ||
        !requireNamespace("DBI", quietly=TRUE)) {
        message(paste("Skipping unit tests for the DuckDB backend",
                      "(no 'arrow' or 'duckdb')"))
        return(invisible(TRUE))
    }

    message("Running unit tests for the DuckDB backend")

    suppressPackageStartupMessages({
        library(Matrix)
        library(S4Arrays)
        library(SparseArray)
        library(DelayedArray)
    })

    p <- 40 ## number of genes
    n <- 90 ## number of samples
    gsets <- list(gset1=paste0("g", 1:10), gset2=paste0("g", 11:30))
    set.seed(123)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.15))
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- Matrix(matrix(cnt, nrow=p, ncol=n, dimnames=dimnames(y)),
                  sparse=TRUE)

    files <- character(0)
    for (input in list(y, cnt)) {
        gsvarownorm <- gsvaRowNorm(gsvaParam(input, gsets, verbose=FALSE),
                                   verbose=FALSE)
        gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
        ## small row groups, to read columns from several of them
        for (x in list(gsvarownorm, gsvacolranks))
            files <- c(files, saveParquetGSVA(x, tempfile(fileext=".parquet"),
                                              colsPerRowGroup=7))
    }

    colsets <- list(NULL, 10:30, c(80L, 3L, 3L, 45L, 90L, 1L))
    for (f in files) {
        sa <- GSVA:::GsvaParquetSeed(f)
        sd <- GSVA:::GsvaParquetSeed(f, backend="duckdb")
        checkIdentical("duckdb", sd@backend)
        checkIdentical(dim(sa), dim(sd))
        checkIdentical(type(sa), type(sd))
        checkIdentical(is_sparse(sa), is_sparse(sd))
        checkIdentical(chunkdim(sa), chunkdim(sd))

        ia <- GSVA:::.parquet_file_info(sa@path, "arrow")
        id <- GSVA:::.parquet_file_info(sd@path, "duckdb")
        gsvafields <- grep("^gsva_", names(ia$metadata), value=TRUE)
        checkIdentical(ia$metadata[gsvafields], id$metadata[gsvafields])
        checkEquals(ia[c("num_rows", "num_columns", "num_row_groups")],
                    id[c("num_rows", "num_columns", "num_row_groups")],
                    check.attributes=FALSE)

        for (j in colsets) {
            checkIdentical(extract_array(sa, list(NULL, j)),
                           extract_array(sd, list(NULL, j)))
            checkIdentical(extract_array(sa, list(c(5L, 1L, 40L), j)),
                           extract_array(sd, list(c(5L, 1L, 40L), j)))
            if (is_sparse(sa)) {
                spa <- extract_sparse_array(sa, list(NULL, j))
                spd <- extract_sparse_array(sd, list(NULL, j))
                checkIdentical(as.matrix(spa), as.matrix(spd))
            }
        }
    }

    unlink(files)
}
