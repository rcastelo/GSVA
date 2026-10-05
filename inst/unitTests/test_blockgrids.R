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

    ## blocks calculated without alignment, as before aligning them
    unaligned <- function(X, nworkers, maxmem, autogrid) {
        typesze <- c("integer"=4, "double"=8)
        if (is.infinite(maxmem) && nworkers == 1 && !is(X, "DelayedMatrix"))
            return(DummyArrayGrid(dim(X)))
        mbl <- getAutoBlockLength(type(X))
        if (!is.infinite(maxmem))
            mbl <- max(mbl, ceiling(maxmem / typesze[type(X)]))
        mbl <- min(.Machine$integer.max, mbl / nworkers)
        ebl <- max(1, ceiling(nrow(X) / nworkers) * as.numeric(ncol(X)))
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
                checkIdentical(blockdims(unaligned(X, nw, mm, rowAutoGrid)),
                               blockdims(GSVA:::.rowgridsize(X, nw, mm)))
                checkIdentical(blockdims(unaligned(X, nw, mm, colAutoGrid)),
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
