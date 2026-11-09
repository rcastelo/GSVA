## the blocks of rows and columns in which GSVA processes data stored in
## chunks, such as in HDF5 or Parquet files, are aligned with those chunks,
## without growing beyond the blocks calculated without alignment, while the
## blocks of data not stored in chunks are the same as without alignment,
## balanced among the workers, see test_balanced_blocks()
test_blockgrids <- function() {

    message("Running unit tests for block grids aligned with chunks")

    suppressPackageStartupMessages({
        library(Matrix)
        library(DelayedArray)
        library(HDF5Array)
    })

    ## blocks calculated without alignment, as before aligning them: with a
    ## finite maximum memory, as many rows or columns as fit in it and,
    ## otherwise, the default block length, split among the workers, and then
    ## balanced among them
    unaligned <- function(X, nworkers, maxmem, whdim) {
        if (is.finite(maxmem))
            n <- GSVA:::.units_per_block(X, whdim, nworkers, maxmem)
        else {
            if (nworkers == 1 && !is(X, "DelayedMatrix"))
                return(DummyArrayGrid(dim(X)))
            mbl <- min(.Machine$integer.max,
                       getAutoBlockLength(type(X)) / nworkers)
            ebl <- max(1, ceiling(dim(X)[whdim] / nworkers) *
                          as.numeric(dim(X)[-whdim]))
            n <- min(max(1, floor(min(mbl, ebl) / max(1, dim(X)[-whdim]))),
                     dim(X)[whdim])
        }
        n <- GSVA:::.balance_blocks(n, X, whdim, nworkers)
        if (whdim == 1) rowAutoGrid(X, nrow=n) else colAutoGrid(X, ncol=n)
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

    ## blocks narrower than half a chunk, within a maximum memory of 1e9 bytes
    ## that fits blocks of 120 columns for all the workers together, give a
    ## warning suggesting the smallest maximum memory, in whole GiB, whose
    ## blocks with the same workers take whole chunks: 40 columns for each of
    ## 8 workers fit in 8/3 * 1e9 bytes, 2.48 GiB
    H <- writeHDF5Array(m, chunkdim=c(37L, 40L))
    upb <- function(nw, mm) as.integer(floor(120 * mm / 1e9 / nw))
    checkEquals(GSVA:::.chunk_maxmem(40L, 8L, 1e9, upb), 3 * 1024^3)
    checkIdentical(GSVA:::.maxmem_text(3 * 1024^3), "3G")
    ## or MiB, below 1 GiB: 8/3 * 1e8 bytes, 254.3 MiB
    upbm <- function(nw, mm) as.integer(floor(120 * mm / 1e8 / nw))
    checkEquals(GSVA:::.chunk_maxmem(40L, 8L, 1e8, upbm), 255 * 1024^2)
    checkIdentical(GSVA:::.maxmem_text(255 * 1024^2), "255M")
    ## a text that 'maxmem' takes
    checkEquals(GSVA:::.memtext2bytes(GSVA:::.maxmem_text(3 * 1024^3)),
                3 * 1024^3)
    chunkwarning <- function(n, nworkers, upb, step=NULL)
        tryCatch(GSVA:::.check_chunk_blocks(H, 2L, n, nworkers, 1e9, upb,
                                            step),
                 warning=conditionMessage)
    w <- chunkwarning(15L, 8L, upb)
    checkTrue(is.character(w) && grepl("blocks of 15 columns", w) &&
              grepl("The calculations process", w) &&
              grepl('maxmem="3G"', w, fixed=TRUE) &&
              grepl("with\\s+the\\s+8\\s+parallel\\s+workers", w) &&
              grepl("Fewer parallel workers", w))
    ## naming the step of the calculations
    w <- chunkwarning(15L, 8L, upb, step="Ranking columns")
    checkTrue(is.character(w) && grepl("Ranking columns processes", w))
    ## or a larger maximum memory, when none gives whole chunks
    w <- chunkwarning(15L, 8L, function(nw, mm) 15L)
    checkTrue(is.character(w) && grepl("larger", w))
    checkTrue(is.na(GSVA:::.chunk_maxmem(40L, 8L, 1e9, function(nw, mm) 15L)))
    ## while blocks of at least half a chunk, or a single block, do not
    checkTrue(!chunkwarning(20L, 6L, upb))
    checkTrue(!chunkwarning(ncol(H), 1L, upb))
    checkTrue(!GSVA:::.check_chunk_blocks(m, 2L, 15L, 8L, 1e9, upb))
    ## the warnings of several steps collected, as gsva() does, give a single
    ## warning with the largest of their maximum memory, here for 12 workers
    ## with blocks of 10 columns, 4 GiB, or none when no step warns
    w <- tryCatch(GSVA:::.with_chunk_blocks({
                      GSVA:::.check_chunk_blocks(H, 2L, 15L, 8L, 1e9, upb,
                                                 "Normalizing rows")
                      GSVA:::.check_chunk_blocks(H, 2L, 10L, 12L, 1e9, upb,
                                                 "Ranking columns")
                  }), warning=conditionMessage)
    checkTrue(is.character(w) && grepl("Some steps", w) &&
              grepl("Normalizing rows: blocks of 15 columns", w) &&
              grepl("Ranking columns: blocks of 10 columns", w) &&
              grepl('maxmem="4G"', w, fixed=TRUE) &&
              grepl("blocks of all of them", w))
    checkIdentical(GSVA:::.with_chunk_blocks(1L), 1L)
    checkTrue(is.null(GSVA:::gsva_global$chunkblocks))
    ## and a single step collected gives its own warning
    w <- tryCatch(GSVA:::.with_chunk_blocks(
                      GSVA:::.check_chunk_blocks(H, 2L, 15L, 8L, 1e9, upb,
                                                 "Ranking columns")),
                  warning=conditionMessage)
    checkTrue(is.character(w) && grepl("Ranking columns processes", w))
    ## also for a subset of the on-disk data, with the chunks of its seed
    checkTrue(is.null(chunkdim(H[, 1:2000])))
    checkIdentical(as.integer(GSVA:::.ondisk_chunkdim(H[, 1:2000])),
                   c(37L, 40L))
    w <- tryCatch(GSVA:::.check_chunk_blocks(H[, 1:2000], 2L, 15L, 8L, 1e9,
                                             upb),
                  warning=conditionMessage)
    checkTrue(is.character(w) && grepl("chunks of 40 columns", w))
    ## whose blocks are aligned to the chunks when their boundaries are the
    ## ones of the subset: consecutive columns from the start of a chunk,
    ## also after operations on each value, or transposing it, while rows
    ## that are not consecutive leave the alignment of the columns
    checkIdentical(GSVA:::.align_to_chunks(100L, H[, 1:2000], 2L), 80L)
    checkIdentical(GSVA:::.align_to_chunks(100L, H[, 41:2000], 2L), 80L)
    checkIdentical(GSVA:::.align_to_chunks(100L, H[, 2:2000], 2L), 100L)
    checkIdentical(GSVA:::.align_to_chunks(100L, H[, c(1:40, 81:200)], 2L),
                   100L)
    checkIdentical(GSVA:::.align_to_chunks(100L, log1p(H)[, 1:2000], 2L), 80L)
    checkIdentical(GSVA:::.align_to_chunks(100L, t(H)[1:2000, ], 1L), 80L)
    checkIdentical(GSVA:::.align_to_chunks(100L, H[c(3, 1, 2), 1:2000], 2L),
                   80L)
    checkIdentical(GSVA:::.align_to_chunks(50L, H[c(3, 1, 2), ], 1L), 50L)
    checkTrue(is.na(GSVA:::.ondisk_chunkdim(H[, 2:2000], aligned=TRUE)[2]))
    ## the number of blocks of a subset is balanced in whole chunks, see
    ## .balance_blocks(): the 1960 columns, 49 chunks, in blocks of 2 chunks
    ## take three rounds of 10 workers, while blocks of 1 chunk take five
    ## rounds, five sixths of their time, and blocks of 20 chunks take two
    ## rounds of 2 workers, which blocks of 13 chunks also take, while blocks
    ## of 5 chunks would take five rounds to save only 1/26 of their time
    checkIdentical(GSVA:::.balance_blocks(80L, H[, 41:2000], 2L, 10L), 40L)
    checkIdentical(GSVA:::.balance_blocks(800L, H[, 41:2000], 2L, 2L), 520L)
}

## the number of rows or columns of a block is the one that fits in the memory
## left by the input in main memory and the output, among the workers, while
## blocks are not smaller than the automatic block size, nor more than one per
## worker
test_units_per_block <- function() {

    message("Running unit tests for the size of blocks within a memory budget")

    suppressPackageStartupMessages({
        library(DelayedArray)
        library(BiocParallel)
    })

    oldautoblocksize <- getAutoBlockSize()
    on.exit(setAutoBlockSize(oldautoblocksize))
    setAutoBlockSize(8 * 1000 * 10) ## 10 columns of 1000 doubles

    X <- matrix(0, nrow=1000, ncol=500)
    insize <- as.numeric(object.size(X))
    colbytes <- insize / 500
    ## the budgets below are for the allocations of R, a fraction of 'maxmem'
    frac <- GSVA:::.mem_fraction_R
    upb <- function(nworkers, maxmem, ..., workermem=0)
        GSVA:::.units_per_block(X, 2L, nworkers, maxmem / frac, ...,
                                workermem=workermem)

    ## memory for 50 columns per block, with 'workfactor' times the size of a
    ## column of working memory and an output of the same size as the input;
    ## budgets have half a column more, to be robust to rounding errors
    maxmem <- insize + 2 * insize + 50.5 * (2 + 1) * colbytes
    checkIdentical(upb(1, maxmem), 50L)
    ## shared by two workers
    checkIdentical(upb(2, maxmem), 25L)
    ## which take, with the main R process, the memory of their R processes
    wm <- 1e6
    checkIdentical(GSVA:::.units_per_block(X, 2L, 2, maxmem / frac + 3 * wm,
                                           workermem=wm), 25L)
    checkIdentical(GSVA:::.units_per_block(X, 2L, 1, maxmem / frac + wm,
                                           workermem=wm), 50L)
    ## which, by default, is the memory of this process besides the data of
    ## the calculations, see .worker_mem(), and forked workers take a copy of
    ## one third of sparse data, unless they are stored in a dgCMatrix object
    oldopt <- options(GSVA.workermem=NULL)
    on.exit(options(oldopt), add=TRUE)
    extra <- GSVA:::.proc_extra_mem()
    checkEqualsNumeric(GSVA:::.worker_mem(Inf), extra)
    checkTrue(GSVA:::.worker_mem() > GSVA:::.worker_mem(insize))
    checkEqualsNumeric(GSVA:::.worker_mem(session=FALSE),
                       0.6 * 1024^3 + extra)
    checkTrue(GSVA:::.worker_mem(full=TRUE) <= GSVA:::.worker_mem())
    ## measured with a full collection of garbage only for large data
    checkTrue(!GSVA:::.large_data(dim(X)))
    checkTrue(GSVA:::.large_data(c(1e4, 1e6), 4))
    checkEqualsNumeric(GSVA:::.fork_copied_mem(3e6, sparse=TRUE), 1e6)
    checkEqualsNumeric(GSVA:::.fork_copied_mem(3e6, sparse=TRUE, dgc=TRUE), 0)
    checkEqualsNumeric(GSVA:::.fork_copied_mem(3e6, sparse=FALSE), 0)
    ## the data shared with forked workers are the matrices of the assays of
    ## the input data, not the rest of its data, such as the names of cells
    se <- SummarizedExperiment::SummarizedExperiment(
              assays=list(counts=X, logcounts=DelayedArray(X)),
              colData=S4Vectors::DataFrame(barcode=sprintf("cell%04d", 1:500)))
    checkEqualsNumeric(GSVA:::.container_shared_mem(se, "counts"), insize)
    checkTrue(GSVA:::.container_held_mem(se, "counts") > insize)
    ## up to the small changes of the memory of this process between calls
    wm <- GSVA:::.worker_mem(insize)
    upbdef <- GSVA:::.units_per_block(X, 2L, 2, maxmem / frac + 3 * wm)
    checkTrue(upbdef >= 24L && upbdef <= 26L)
    options(oldopt)
    ## unless they leave no memory for the blocks, which are then as large as
    ## without them
    checkIdentical(GSVA:::.units_per_block(X, 2L, 2, maxmem / frac,
                                           workermem=1e12), 25L)
    ## blocks sent to the workers from this process, which takes a copy of
    ## the block of each worker, take one more column per column of a block
    maxmem <- insize + 2 * insize + 2 * 25.5 * (2 + 1 + 1) * colbytes
    checkIdentical(upb(2, maxmem, sendblocks=TRUE), 25L)
    checkIdentical(upb(2, maxmem), 34L)
    ## except data on disk, which the workers take from it themselves
    checkIdentical(GSVA:::.units_per_block(DelayedArray(X), 2L, 2,
                                           2 * 25.5 * (2 + 1) * colbytes / frac,
                                           workermem=0, sendblocks=TRUE), 25L)
    ## and, with them, the memory estimated for the step
    mf <- list(workfactor=2, outfactor=1, outextra=0)
    datamem <- function(sendblocks)
        GSVA:::.step_data_mem(dim(X), 2L, insize, 8, FALSE, TRUE, mf, 0,
                              sendblocks=sendblocks)
    checkEqualsNumeric(datamem(TRUE) - datamem(FALSE), 10 * colbytes)
    checkTrue(GSVA:::.sends_blocks(X, SnowParam(2), 2L))
    checkTrue(!GSVA:::.sends_blocks(X, MulticoreParam(2), 2L))
    checkTrue(!GSVA:::.sends_blocks(DelayedArray(X), SnowParam(2), 2L))
    checkTrue(!GSVA:::.sends_blocks(X, SnowParam(2), 0L))
    ## more working memory per column
    maxmem <- insize + 2 * insize + 50.5 * (2 + 1) * colbytes
    checkIdentical(upb(1, maxmem, workfactor=4), 30L)
    ## an output of 8 bytes per gene set per column for 100 gene sets
    maxmem <- insize + 2 * 800 * 500 + 50.5 * (2 * colbytes + 800)
    checkIdentical(upb(1, maxmem, outfactor=0, outextra=800), 50L)
    ## which is not assembled in memory when it is written to disk by blocks
    maxmem <- insize + 50.5 * (2 * colbytes + 800)
    checkIdentical(upb(1, maxmem, outfactor=0, outextra=800, sinkout=TRUE),
                   50L)
    ## not smaller than the automatic block size, of 10 columns, for all the
    ## workers together
    checkIdentical(upb(1, insize), 10L)
    checkIdentical(upb(2, insize), 5L)
    ## not more than one block per worker
    checkIdentical(upb(3, Inf), 167L)
    ## data on disk is not in main memory, has the size of its dense form,
    ## and its output, written to disk by blocks, is not assembled in memory
    D <- DelayedArray(X)
    maxmem <- 50.5 * (2 + 1) * 1000 * 8
    checkIdentical(GSVA:::.units_per_block(D, 2L, 1, maxmem / frac), 50L)
    ## binding blocks of columns of SVT_SparseMatrix objects reuses their
    ## memory, so that the output is assembled once, and not twice
    S <- as(Matrix::rsparsematrix(1000, 500, density=0.1), "SVT_SparseMatrix")
    ssize <- as.numeric(object.size(S))
    scolbytes <- ssize / 500
    maxmem <- ssize + ssize + 50.5 * (2 + 1) * scolbytes
    checkIdentical(GSVA:::.units_per_block(S, 2L, 1, maxmem / frac), 50L)
    ## while binding dgCMatrix objects takes twice the size of the output
    S <- Matrix::rsparsematrix(1000, 500, density=0.1)
    ssize <- as.numeric(object.size(S))
    scolbytes <- ssize / 500
    maxmem <- ssize + 2 * ssize + 50.5 * (2 + 1) * scolbytes
    checkIdentical(GSVA:::.units_per_block(S, 2L, 1, maxmem / frac), 50L)
    ## blocks of sparse data converted into dense matrices take the working
    ## memory and output of their dense form, while the input keeps its size
    dcolbytes <- 1000 * 8
    maxmem <- ssize + 2 * dcolbytes * 500 + 50.5 * (2 + 1) * dcolbytes
    checkIdentical(GSVA:::.units_per_block(S, 2L, 1, maxmem / frac,
                                           dense=TRUE), 50L)
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

## with workers not forked from this process, such as socket workers, blocks
## of data in main memory are sent to the workers one at a time, with the same
## results, retries of failed blocks and errors as with other workers
test_blockprocessing_socket_workers <- function() {

    message("Running unit tests for block processing with socket workers")

    suppressPackageStartupMessages(library(BiocParallel))

    set.seed(123)
    X <- matrix(rnorm(200 * 150), nrow=200, ncol=150)
    flag <- tempfile()
    on.exit(unlink(flag), add=TRUE)
    ## fails the first time it is called, in any worker
    flaky_fun <- function(x, verbose, flag) {
        if (!file.exists(flag)) {
            file.create(flag)
            stop("simulated transient failure")
        }
        x * 2
    }
    failing_fun <- function(x, verbose) stop("simulated failure")
    ## started once, so that workers load GSVA only once
    bpparam <- bpstart(SnowParam(workers=2, progressbar=FALSE))
    on.exit(bpstop(bpparam), add=TRUE)
    for (proc in list(GSVA:::.processMatrixRows, GSVA:::.processMatrixCols)) {
        res <- suppressMessages(proc(X, FUN=function(x, verbose) x * 2,
                                     verbose=FALSE, BPPARAM=bpparam))
        checkEqualsNumeric(X * 2, res)
        unlink(flag)
        res <- suppressMessages(proc(X, FUN=flaky_fun, flag=flag,
                                     verbose=FALSE, BPPARAM=bpparam))
        checkEqualsNumeric(X * 2, res)
        err <- tryCatch(suppressMessages(proc(X, FUN=failing_fun,
                                              verbose=FALSE,
                                              BPPARAM=bpparam)),
                        error=conditionMessage)
        checkTrue(is.character(err) &&
                  grepl("Cancelling execution", err, fixed=TRUE))
    }
}

## with a DoparParam back-end, which does not support bpiterate(), blocks of
## data in main memory are sent to the workers in rounds of as many blocks as
## workers, with the same results, retries of failed blocks and errors as with
## other workers. the foreach back-end of doParallel forks this process, which
## is only possible in unix, so that the workers do not need to load GSVA
test_blockprocessing_dopar_workers <- function() {

    if (.Platform$OS.type != "unix" ||
        !requireNamespace("doParallel", quietly=TRUE))
        return(invisible(TRUE))

    message("Running unit tests for block processing with a DoparParam back-end")

    suppressPackageStartupMessages(library(BiocParallel))

    doParallel::registerDoParallel(cores=2)
    on.exit(foreach::registerDoSEQ(), add=TRUE)
    bpparam <- DoparParam()

    set.seed(123)
    X <- matrix(rnorm(200 * 150), nrow=200, ncol=150)
    flag <- tempfile()
    on.exit(unlink(flag), add=TRUE)
    ## fails the first time it is called on a block whose first value is
    ## 'x11', in any worker
    flaky_fun <- function(x, verbose, flag, x11) {
        if (x[1, 1] == x11 && !file.exists(flag)) {
            file.create(flag)
            stop("simulated transient failure")
        }
        x * 2
    }
    failing_fun <- function(x, verbose) stop("simulated failure")
    for (proc in list(GSVA:::.processMatrixRows, GSVA:::.processMatrixCols)) {
        res <- suppressMessages(proc(X, FUN=function(x, verbose) x * 2,
                                     verbose=FALSE, BPPARAM=bpparam))
        checkEqualsNumeric(X * 2, res)
        err <- tryCatch(suppressMessages(proc(X, FUN=failing_fun,
                                              verbose=FALSE,
                                              BPPARAM=bpparam)),
                        error=conditionMessage)
        checkTrue(is.character(err) &&
                  grepl("Cancelling execution", err, fixed=TRUE))
    }
    checkTrue(bpstopOnError(bpparam))

    ## seven blocks of rows run in four rounds by two workers, where a failed
    ## block of the second round is executed again, recomputing only that
    ## round
    rngs <- lapply(seq(1, 200, by=30), function(s)
                   IRanges::IRanges(s, min(s + 29, nrow(X))))
    bpstopOnError(bpparam) <- FALSE
    res <- bptry(GSVA:::.bp_blocks(X, 1L, rngs, flaky_fun, NULL, flag=flag,
                                   x11=X[91, 1], BPPARAM=bpparam))
    checkIdentical(length(attr(res, "rounds")), 4L)
    checkIdentical(bpok(res), seq_along(rngs) != 4L)
    res <- GSVA:::.bp_blocks(X, 1L, rngs, flaky_fun, NULL, flag=flag,
                             x11=X[91, 1], BPREDO=res, BPPARAM=bpparam)
    checkTrue(all(bpok(res)))
    checkEqualsNumeric(X * 2, do.call("rbind", res))
}

## with more than one worker, the number of blocks is a multiple of the number
## of workers, so that no worker is idle in the last round of blocks, by
## narrowing the blocks while keeping their number of rounds, and blocks of
## data stored in chunks keep taking whole chunks
test_balanced_blocks <- function() {

    message("Running unit tests for blocks balanced among the workers")

    suppressPackageStartupMessages({
        library(DelayedArray)
        library(HDF5Array)
    })

    X <- matrix(0, nrow=2, ncol=1000000)
    bb <- function(n, nworkers) GSVA:::.balance_blocks(n, X, 2L, nworkers)
    ## blocks of the scores of 1,000,000 cells with 10 workers, which were 11,
    ## 23 and 55, in 2, 3 and 6 rounds, become 20, 30 and 60
    for (x in list(c(93982, 50000, 20), c(43546, 33334, 30),
                   c(18328, 16667, 60))) {
        w <- bb(as.integer(x[1]), 10L)
        checkIdentical(w, as.integer(x[2]))
        checkEqualsNumeric(ceiling(ncol(X) / w), x[3])
    }
    ## blocks already balanced, a single worker or a single block, unchanged
    checkIdentical(bb(100000L, 10L), 100000L)
    checkIdentical(bb(93982L, 1L), 93982L)
    checkIdentical(bb(1000000L, 10L), 1000000L)
    ## blocks that do not split evenly: 3000 columns in 120 blocks of 25 with
    ## 7 workers, 18 rounds, become 125 blocks of 24, still 18 rounds
    X <- matrix(0, nrow=2, ncol=3000)
    checkIdentical(bb(25L, 7L), 24L)

    ## through the grids of columns and rows: with an automatic block size of
    ## 100 columns of doubles, 3 workers take blocks of 33 columns or rows,
    ## 34 blocks, which become 36 blocks of 31
    oldautoblocksize <- getAutoBlockSize()
    on.exit(setAutoBlockSize(oldautoblocksize))
    setAutoBlockSize(8 * 1100 * 100)
    X <- matrix(0, nrow=1100, ncol=1100)
    checkIdentical(length(GSVA:::.colgridsize(X, 3L)), 36L)
    checkIdentical(length(GSVA:::.rowgridsize(X, 3L)), 36L)
    checkIdentical(ncol(GSVA:::.colgridsize(X, 3L)[[1L]]), 31L)
    setAutoBlockSize(oldautoblocksize)

    ## data in chunks of 31 columns, 3000 columns in 97 chunks: blocks of 6
    ## chunks with 10 workers, 17 blocks in 2 rounds, become 20 blocks of 5
    ## chunks, while blocks of 4 chunks with 3 workers, 25 blocks in 9 rounds,
    ## are kept, because narrower blocks of whole chunks would save less than
    ## a tenth of their time, e.g., 33 blocks of 3 chunks in 11 rounds
    H <- writeHDF5Array(matrix(0, nrow=2, ncol=3000), chunkdim=c(2L, 31L))
    checkIdentical(GSVA:::.balance_blocks(186L, H, 2L, 10L), 155L)
    checkIdentical(length(colAutoGrid(H, ncol=155L)), 20L)
    checkIdentical(GSVA:::.balance_blocks(124L, H, 2L, 3L), 124L)
    ## 3000 columns in 30 chunks of 100 columns: blocks of 2 chunks, 15 blocks
    ## in 2 rounds of 10 workers, become 30 blocks of 1 chunk in 3 rounds,
    ## which take three quarters of their time
    H <- writeHDF5Array(matrix(0, nrow=2, ncol=3000), chunkdim=c(2L, 100L))
    checkIdentical(GSVA:::.balance_blocks(200L, H, 2L, 10L), 100L)
}

test_ondisk_mem_factors <- function() {
    ## a sparse matrix with a quarter of nonzero values, in every chunk and in
    ## every column, sampled through its chunks or, without them, through
    ## blocks of columns
    x <- matrix(0, nrow=40, ncol=400)
    x[seq(1L, length(x), by=4L)] <- 1
    H <- writeHDF5Array(x, chunkdim=c(10L, 50L), as.sparse=TRUE)
    checkEqualsNumeric(GSVA:::.sample_nonzero_fraction(H), 1/4)
    checkEqualsNumeric(GSVA:::.sample_nonzero_fraction(DelayedArray(as(x, "dgCMatrix")),
                                                       maxvalues=400), 1/4)

    ## the factors of ranking columns relative to the sparse size of blocks
    ## read from disk are scaled by the ratio of that size to their dense
    ## size, 3/8 with doubles, value and row index, and 1/2 with integers,
    ## unlike the ones of dense or in-memory data, of steps without such
    ## factors or without a fraction of nonzero values
    mf <- GSVA:::.ondisk_mem_factors(GSVA:::.step_mem_factors("colranks", H),
                                     H, 1/4)
    checkEqualsNumeric(mf$workfactor, 14.3 * 3/8)
    checkEqualsNumeric(mf$outfactor, 0.7 * 3/8)
    xi <- x
    storage.mode(xi) <- "integer"
    Hi <- writeHDF5Array(xi, chunkdim=c(10L, 50L), as.sparse=TRUE)
    mf <- GSVA:::.ondisk_mem_factors(GSVA:::.step_mem_factors("scores", Hi,
                                                              ngs=10), Hi, 1/4)
    checkEqualsNumeric(mf$workfactor, 2.5 * 1/2)
    checkEqualsNumeric(mf$outextra, 80)
    mf <- GSVA:::.step_mem_factors("colranks", H)
    checkIdentical(GSVA:::.ondisk_mem_factors(mf, H)$workfactor, 7)
    D <- writeHDF5Array(x, chunkdim=c(10L, 50L))
    mf <- GSVA:::.step_mem_factors("colranks", D)
    checkIdentical(GSVA:::.ondisk_mem_factors(mf, D, 1/4)$workfactor, 2)
    S <- as(x, "dgCMatrix")
    mf <- GSVA:::.step_mem_factors("colranks", S)
    checkIdentical(GSVA:::.ondisk_mem_factors(mf, S, 1/4)$workfactor, 7)
    mf <- GSVA:::.step_mem_factors("rowstats", H)
    checkIdentical(GSVA:::.ondisk_mem_factors(mf, H, 1/4)$workfactor, 2)
}
