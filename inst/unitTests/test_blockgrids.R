## the blocks of rows and columns in which GSVA processes data stored in
## chunks, such as in HDF5 or Parquet files, are aligned with those chunks,
## without growing beyond the blocks calculated without alignment, while the
## blocks of data not stored in chunks are the same as without alignment
test_blockgrids <- function() {

    message("Running unit tests for block grids aligned with chunks")

    suppressPackageStartupMessages({
        library(Matrix)
        library(DelayedArray)
        library(HDF5Array)
    })

    ## blocks calculated without alignment, as before aligning them: with a
    ## finite maximum memory, as many rows or columns as fit in it and,
    ## otherwise, the default block length, split among the workers
    unaligned <- function(X, nworkers, maxmem, whdim) {
        autogrid <- if (whdim == 1) rowAutoGrid else colAutoGrid
        if (is.finite(maxmem)) {
            n <- GSVA:::.units_per_block(X, whdim, nworkers, maxmem)
            return(if (whdim == 1) rowAutoGrid(X, nrow=n)
                   else colAutoGrid(X, ncol=n))
        }
        if (nworkers == 1 && !is(X, "DelayedMatrix"))
            return(DummyArrayGrid(dim(X)))
        mbl <- min(.Machine$integer.max, getAutoBlockLength(type(X)) / nworkers)
        ebl <- max(1, ceiling(dim(X)[whdim] / nworkers) *
                      as.numeric(dim(X)[-whdim]))
        autogrid(X, block.length=min(mbl, ebl))
    }

    set.seed(123)
    m <- matrix(rnorm(777*3000), nrow=777, ncol=3000)
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(3000*8*100) ## 100 rows of doubles
    on.exit(setAutoBlockSize(oldautoblocksize))

    blockdims <- function(g) vapply(seq_along(g), function(i)
                                    dim(g[[as.integer(i)]]), integer(2))

    for (X in list(m, rsparsematrix(777, 3000, 0.1), DelayedArray(m))) {
        for (nw in c(1, 3)) {
            for (mm in c(Inf, 5e6)) {
                checkIdentical(blockdims(unaligned(X, nw, mm, 1)),
                               blockdims(GSVA:::.rowgridsize(X, nw, mm)))
                checkIdentical(blockdims(unaligned(X, nw, mm, 2)),
                               blockdims(GSVA:::.colgridsize(X, nw, mm)))
            }
        }
    }

    ## chunks of 37 rows and 31 columns: blocks of 100 rows are aligned to 74
    ## rows, while blocks of 33 rows, with three workers, are kept because
    ## they are smaller than one chunk; blocks of 100000 / 777 = 128 columns
    ## are aligned to 124 columns and blocks of 42 columns, with three
    ## workers, to 31 columns, while blocks of 25 columns, with five workers,
    ## are kept
    H <- writeHDF5Array(m, chunkdim=c(37L, 31L))
    checkIdentical(74L, nrow(GSVA:::.rowgridsize(H)[[1L]]))
    checkIdentical(33L, nrow(GSVA:::.rowgridsize(H, 3)[[1L]]))
    setAutoBlockSize(100000*8)
    checkIdentical(124L, ncol(GSVA:::.colgridsize(H)[[1L]]))
    checkIdentical(31L, ncol(GSVA:::.colgridsize(H, 3)[[1L]]))
    checkIdentical(25L, ncol(GSVA:::.colgridsize(H, 5)[[1L]]))
}

## the number of rows or columns of a block is the one that fits in the memory
## left by the input in main memory and the output, among the workers, while
## blocks are not smaller than the automatic block size, nor more than one per
## worker
test_units_per_block <- function() {

    message("Running unit tests for the size of blocks within a memory budget")

    suppressPackageStartupMessages(library(DelayedArray))

    oldautoblocksize <- getAutoBlockSize()
    on.exit(setAutoBlockSize(oldautoblocksize))
    setAutoBlockSize(8 * 1000 * 10) ## 10 columns of 1000 doubles

    X <- matrix(0, nrow=1000, ncol=500)
    insize <- as.numeric(object.size(X))
    colbytes <- insize / 500
    ## the budgets below are for the allocations of R, a fraction of 'maxmem'
    frac <- GSVA:::.mem_fraction_R
    upb <- function(nworkers, maxmem, ...)
        GSVA:::.units_per_block(X, 2L, nworkers, maxmem / frac, ...)

    ## memory for 50 columns per block, with 'workfactor' times the size of a
    ## column of working memory and an output of the same size as the input
    maxmem <- insize + 2 * insize + 50 * (2 + 1) * colbytes
    checkIdentical(upb(1, maxmem), 50L)
    ## shared by two workers
    checkIdentical(upb(2, maxmem), 25L)
    ## more working memory per column
    checkIdentical(upb(1, maxmem, workfactor=4), 30L)
    ## an output of 8 bytes per gene set per column for 100 gene sets
    maxmem <- insize + 2 * 800 * 500 + 50 * (2 * colbytes + 800)
    checkIdentical(upb(1, maxmem, outfactor=0, outextra=800), 50L)
    ## not smaller than the automatic block size, of 10 columns
    checkIdentical(upb(1, insize), 10L)
    ## not more than one block per worker
    checkIdentical(upb(3, Inf), 167L)
    ## data on disk is not in main memory and has the size of its dense form
    D <- DelayedArray(X)
    maxmem <- 2 * insize + 50 * (2 + 1) * colbytes
    checkIdentical(GSVA:::.units_per_block(D, 2L, 1, maxmem / frac), 50L)
    ## blocks have at most .Machine$integer.max values, which large sparse
    ## data in main memory could otherwise exceed
    S <- Matrix::sparseMatrix(i=1, j=1, x=1, dims=c(60000, 50000))
    checkIdentical(GSVA:::.units_per_block(S, 1L, 1, 2^40),
                   as.integer(floor(.Machine$integer.max / 50000)))
}

test_blockprocessing_bpparam_unchanged <- function() {

    message("Running unit tests for block processing keeping BPPARAM unchanged")

    suppressPackageStartupMessages(library(BiocParallel))

    set.seed(123)
    X <- matrix(rnorm(200 * 150), nrow=200, ncol=150)
    newbpparam <- function() {
        if (.Platform$OS.type == "unix")
            MulticoreParam(workers=2, stop.on.error=TRUE, progressbar=FALSE)
        else
            SnowParam(workers=2, stop.on.error=TRUE, progressbar=FALSE)
    }
    identity_fun <- function(x, verbose) x
    failing_fun <- function(x, verbose) stop("simulated failure")

    ## the parallel execution in chunks sets stop.on.error=FALSE and, with
    ## verbose=TRUE, progressbar=TRUE, in BPPARAM, which is a reference class
    ## object, and these settings are restored when the calculations finish,
    ## also when they fail
    for (proc in list(GSVA:::.processMatrixRows, GSVA:::.processMatrixCols)) {
        bpparam <- newbpparam()
        res <- suppressMessages(proc(X, FUN=identity_fun, verbose=TRUE,
                                     BPPARAM=bpparam))
        checkEqualsNumeric(X, res)
        checkTrue(bpstopOnError(bpparam))
        checkTrue(!bpprogressbar(bpparam))

        bpparam <- newbpparam()
        err <- tryCatch(suppressMessages(proc(X, FUN=failing_fun, verbose=TRUE,
                                              BPPARAM=bpparam)),
                        error=conditionMessage)
        checkTrue(is.character(err) &&
                  grepl("Cancelling execution", err, fixed=TRUE))
        checkTrue(bpstopOnError(bpparam))
        checkTrue(!bpprogressbar(bpparam))
    }
}

test_blockprocessing_dead_worker <- function() {

    if (.Platform$OS.type != "unix") ## forked worker processes
        return(invisible(TRUE))

    message("Running unit tests for block processing with a dead worker")

    suppressPackageStartupMessages(library(BiocParallel))

    ## a forked worker process that ends without returning its result, as
    ## when it crashes, gives an error telling so
    X <- matrix(rnorm(200 * 150), nrow=200, ncol=150)
    killing_fun <- function(x, verbose) {
        tools::pskill(Sys.getpid(), tools::SIGKILL)
        x
    }
    err <- tryCatch(suppressWarnings(suppressMessages(
                        GSVA:::.processMatrixRows(X, FUN=killing_fun,
                                                  verbose=FALSE,
                                                  BPPARAM=MulticoreParam(2)))),
                    error=conditionMessage)
    checkTrue(is.character(err) &&
              grepl("parallel worker process ended", err, fixed=TRUE))
}
