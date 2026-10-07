## This file contains the internal functions that are used to process the rows
## or the columns of a matrix in blocks, stored in main memory or in an on-disk
## data structure such as an HDF5Matrix object, opening parallelism through a
## BiocParallelParam object when different from NULL, and reporting progress
## using the 'cli' package when possible. The block size is automatically
## determined based on the number of workers and the maximum available memory,
## if specified. The functions in this file are not exported and are not
## intended to be used directly by the users of the package. They are only used
## internally by the functions that need to process the rows or the columns of
## a matrix in blocks.

## fraction of the maximum main memory available to the allocations of R, the
## rest being taken by the memory that R releases but that the process keeps,
## and by allocations outside R, e.g., by the HDF5 library; in the steps of
## GSVA on single-cell data, the peak of memory allocated by R was 70% to 80%
## of the peak of memory resident in the process
.mem_fraction_R <- 0.7

## memory taken by the steps of GSVA to process each row or column of the
## input data 'X', relative to its size in 'X' ('workfactor'), and by the
## output of each row or column, relative to that size ('outfactor') plus a
## number of bytes ('outextra'), with 'ngs' gene sets in the steps giving
## enrichment scores, which are dense. the
## factors were measured on single-cell data and include a margin; the output
## of normalizing rows is double, and of ranking columns integer. instead of
## 'X', whether its values are sparse or integer can be given directly. the
## optional factor 'assembly' gives the memory taken by the output in main
## memory, relative to its size, while it is produced, see .fixed_mem(). the
## row normalization with CLR ('clr=TRUE') reads its input once through blocks
## of columns, see .rowStats(), and gives its output in one piece, which takes
## twice its size when the input is a dgCMatrix object ('dgc=TRUE') converted
## into an SVT_SparseMatrix object, see .rownorm_clr(); when the input is
## stored on disk, its output is a delayed operation that takes no memory
#' @importFrom BiocGenerics type
#' @importFrom S4Arrays is_sparse
.step_mem_factors <- function(step=c("rownorm", "rowstats", "colranks",
                                     "scores", "average", "plage", "zscore",
                                     "ssgsea"), X=NULL, ngs=0,
                              sparse=is_sparse(X),
                              int=(type(X) == "integer"), clr=FALSE,
                              dgc=is(X, "dgCMatrix")) {
    step <- match.arg(step)
    if (step == "rownorm" && clr)
        return(list(workfactor=if (sparse) 2 else 3,
                    outfactor=if (!int) 1 else if (sparse) 1.5 else 2,
                    outextra=0, assembly=if (dgc) 2 else 1))
    switch(step,
           rownorm=list(workfactor=if (sparse) 6 else 4,
                        outfactor=if (!int) 1 else if (sparse) 1.5 else 2,
                        outextra=0),
           ## the row statistics of each block of columns, four or five
           ## doubles per row, are small and reduced at the end, see
           ## .rowStats(); dense blocks take the mask of missing values and
           ## the logarithm of their values
           rowstats=list(workfactor=if (sparse) 2 else 3, outfactor=0,
                         outextra=0),
           colranks=list(workfactor=if (sparse) 7 else 2,
                         outfactor=if (int) 1 else if (sparse) 0.7 else 0.5,
                         outextra=0),
           scores=list(workfactor=2, outfactor=0, outextra=8 * ngs),
           ## the average method gives its scores directly, while PLAGE,
           ## z-score and ssGSEA give them from an intermediate matrix as large
           ## as their input, of scaled rows or column ranks
           average=list(workfactor=2, outfactor=0, outextra=8 * ngs),
           plage=, zscore=,
           ssgsea=list(workfactor=2, outfactor=1, outextra=8 * ngs))
}

## memory that remains allocated while a matrix with 'nunits' rows (whdim=1)
## or columns (whdim=2) of 'unitbytes' bytes each, sparse or not, is processed
## in blocks: the matrix itself, when it is in main memory, and the output,
## see .step_mem_factors(). when the matrix is stored on disk, the output
## proportional to it is written to disk block by block. otherwise, the output
## is assembled from the blocks at the end, which takes twice its size, except
## when binding blocks of columns of sparse data, which reuses the memory of
## the blocks, or as given by 'assembly'. the output given by 'outextra' is
## always assembled in memory
.fixed_mem <- function(nunits, unitbytes, whdim, sparse, inmemory, outfactor,
                       outextra, assembly=NULL) {
    fixed <- 2 * outextra * nunits
    if (inmemory) {
        if (is.null(assembly))
            assembly <- if (whdim == 2L && sparse) 1 else 2
        fixed <- fixed + nunits * unitbytes +
                 assembly * outfactor * unitbytes * nunits
    }

    fixed
}

## number of rows (whdim=1) or columns (whdim=2) of the blocks in which the
## matrix 'X' can be processed by 'nworkers' workers within the maximum main
## memory 'maxmem', where 'workfactor', 'outfactor' and 'outextra' are given
## by .step_mem_factors(). the memory available to the blocks of the workers
## is what is left by the memory that remains allocated, see .fixed_mem(),
## from the fraction .mem_fraction_R, available to the allocations of R, of
## the memory that the R processes of more than one worker, of 'workermem'
## bytes each, see .worker_mem(), leave from 'maxmem', when they leave enough
## for the memory that remains allocated; 'workermem=0' when the workers do
## not share 'maxmem', such as the jobs of gsvaMap(). the size of a row or column of 'X' stored on disk is the one of its
## dense form, which overestimates the size of sparse data. the blocks of all
## workers together are not smaller than the automatic block size of the
## DelayedArray package, see getAutoBlockSize(), because the overhead of
## processing many smaller blocks makes calculations too slow, so that a
## smaller 'maxmem' is not honored, and
## they have at most as many rows or columns as needed to give one block to
## each worker, and as the DelayedArray package supports
#' @importFrom BiocGenerics type
#' @importFrom utils object.size
#' @importFrom DelayedArray getAutoBlockLength
#' @importFrom S4Arrays is_sparse
.units_per_block <- function(X, whdim, nworkers, maxmem, workfactor=2,
                             outfactor=1, outextra=0,
                             workermem=.worker_mem()) {
    nunits <- dim(X)[whdim]
    inmemory <- !is(X, "DelayedArray")
    if (inmemory)
        unitbytes <- as.numeric(object.size(X)) / max(1, nunits)
    else
        unitbytes <- as.numeric(dim(X)[-whdim]) *
                     if (type(X) == "integer") 4 else 8
    fixed <- .fixed_mem(nunits, unitbytes, whdim, is_sparse(X), inmemory,
                        outfactor, outextra)
    ## with more than one worker, their R processes take memory from 'maxmem',
    ## unless they leave no memory for the blocks, when 'maxmem' cannot be
    ## honored anyway, and smaller blocks would only make calculations slower
    if (nworkers > 1 &&
        .mem_fraction_R * (maxmem - nworkers * workermem) > fixed)
        maxmem <- maxmem - nworkers * workermem
    avail <- (.mem_fraction_R * maxmem - fixed) / nworkers
    nperblock <- floor(avail / ((workfactor + outfactor) * unitbytes +
                                outextra))
    otherdim <- max(1, as.numeric(dim(X)[-whdim]))
    nperblock <- max(nperblock, floor(getAutoBlockLength(type(X)) / nworkers /
                                      otherdim))
    ## the DelayedArray package does not support blocks with more than
    ## .Machine$integer.max values
    nperblock <- min(ceiling(nunits / nworkers),
                     floor(.Machine$integer.max / otherdim), nperblock)

    as.integer(max(1, nperblock))
}

## adapted from .define_multiworker_grid() in beachmat/R/colBlockApply.R; with a
## finite maximum main memory 'maxmem', blocks have as many rows as fit in it,
## see .units_per_block(), and otherwise the default block length of the
## DelayedArray package, when there are several workers or 'X' is stored on
## disk, or form a single block
#' @importFrom BiocGenerics type
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength
.rowgridsize <- function(X, nworkers=1, maxmem=Inf, workfactor=2, outfactor=1,
                         outextra=0, workermem=.worker_mem()) {
  grid <- DummyArrayGrid(dim(X))
  if (is.finite(maxmem)) {
      nrowblock <- .units_per_block(X, 1L, nworkers, maxmem, workfactor,
                                    outfactor, outextra, workermem)
      grid <- rowAutoGrid(X, nrow=.align_to_chunks(nrowblock, X, 1L))
  } else if (nworkers > 1 || is(X, "DelayedMatrix")) {
      ## assuming all workers share memory, the maximum block length has to reduce
      ## by the number of workers and, in any case, it cannot exceed
      ## .Machine$integer.max
      max.block.length <- min(.Machine$integer.max,
                              getAutoBlockLength(type(X)) / nworkers)
      expected.block.length <- max(1, ceiling(nrow(X) / nworkers) * as.numeric(ncol(X)))
      block.length <- min(max.block.length, expected.block.length)
      ## number of rows per block, as calculated by rowAutoGrid()
      nrowblock <- min(max(1, floor(block.length / max(1, ncol(X)))), nrow(X))
      grid <- rowAutoGrid(X, nrow=.align_to_chunks(nrowblock, X, 1L))
  }
  grid
}

## adapted from .define_multiworker_grid() in beachmat/R/colBlockApply.R; with a
## finite maximum main memory 'maxmem', blocks have as many columns as fit in
## it, see .units_per_block(), and otherwise the default block length of the
## DelayedArray package, when there are several workers or 'X' is stored on
## disk, or form a single block
#' @importFrom BiocGenerics type
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength
.colgridsize <- function(X, nworkers=1, maxmem=Inf, workfactor=2, outfactor=1,
                         outextra=0, workermem=.worker_mem()) {
  grid <- DummyArrayGrid(dim(X))
  if (is.finite(maxmem)) {
      ncolblock <- .units_per_block(X, 2L, nworkers, maxmem, workfactor,
                                    outfactor, outextra, workermem)
      grid <- colAutoGrid(X, ncol=.align_to_chunks(ncolblock, X, 2L))
  } else if (nworkers > 1 || is(X, "DelayedMatrix")) {
      ## assuming all workers share memory, the maximum block length has to reduce
      ## by the number of workers and, in any case, it cannot exceed
      ## .Machine$integer.max
      max.block.length <- min(.Machine$integer.max,
                              getAutoBlockLength(type(X)) / nworkers)
      expected.block.length <- max(1, ceiling(ncol(X) / nworkers) * as.numeric(nrow(X)))
      block.length <- min(max.block.length, expected.block.length)
      ## number of columns per block, as calculated by colAutoGrid()
      ncolblock <- min(max(1, floor(block.length / max(1, nrow(X)))), ncol(X))
      grid <- colAutoGrid(X, ncol=.align_to_chunks(ncolblock, X, 2L))
  }
  grid
}

## when the data in 'X' is stored in chunks, such as in HDF5 or Parquet files,
## round down the block width 'n' along the dimension 'whdim' to a multiple of
## the chunk width, so that blocks do not split chunks, which would then be
## read more than once; when a block is narrower than a chunk, its width is
## kept, to avoid exceeding the memory bound behind the block width
#' @importFrom DelayedArray chunkdim
.align_to_chunks <- function(n, X, whdim) {
    cd <- if (is(X, "DelayedArray")) chunkdim(X) else NULL
    if (!is.null(cd) && cd[whdim] > 0 && n >= cd[whdim])
        n <- (n %/% cd[whdim]) * cd[whdim]
    as.integer(n)
}

#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges ranges
.splitRowsInRanges <- function(grid) {
    rir <- lapply(grid, function(r) ranges(r)[1])
    rir
}

## which elements of the result of bplapply() or bpiterate() did not fail;
## bpiterate() leaves the elements that failed as NULL, and stores their
## errors in the attribute 'errors', named by the position of the element
#' @importFrom BiocParallel bpok
.bp_ok <- function(res) {
    ok <- bpok(res)
    errs <- attr(res, "errors")
    if (!is.null(errs))
        ok[as.integer(names(errs))] <- FALSE

    ok
}

## error of the first element of the result of bplapply() or bpiterate() that
## failed
.bp_first_error <- function(res) {
    errs <- attr(res, "errors")
    if (!is.null(errs))
        return(errs[[1L]])

    res[[which(!bpok(res))[1L]]]
}

#' @importFrom cli cli_alert_warning
.report_parallel_errors <- function(res) {
    bpokmask <- .bp_ok(res)
    msg <- paste("{sum(!bpokmask)} execution thread(s) gave an error,",
                 "reporting the first one.")
    cli_alert_warning(msg)
    msg <- attr(.bp_first_error(res), "traceback")
    msg <- gsub("\\}", "]", gsub("\\{", "[", msg))
    out <- lapply(msg, cli_alert_warning)
}

## evaluate 'expr', a call to bplapply() or bpiterate(), replacing the
## uninformative error that BiocParallel gives when a forked worker process
## ends without returning its result, such as "wrong args for environment
## subassignment", by an error that tells what happened. that error of
## BiocParallel follows the warning of parallel::mccollect() that a parallel
## job did not deliver its result, whose message is matched in English only, so
## that in other languages the error of BiocParallel is given as it is. errors
## of the workers, which BiocParallel gives as 'bperror' objects, are not
## replaced
#' @importFrom cli cli_abort
.with_dead_worker_check <- function(expr) {
    nodelivery <- FALSE
    withCallingHandlers(expr,
        warning=function(w) {
            if (grepl("parallel jobs? did not deliver", conditionMessage(w)))
                nodelivery <<- TRUE
        },
        error=function(e) {
            if (nodelivery && !inherits(e, "bperror"))
                cli_abort(c("x"=paste("A parallel worker process ended",
                                      "without returning its result."),
                            "i"=paste("This happens, for instance, when the",
                                      "process crashes, e.g., with a",
                                      "segmentation fault, or when the",
                                      "operating system or the workload",
                                      "manager kills it, e.g., for exceeding",
                                      "a memory limit."),
                            "i"=paste("The error output of the process,",
                                      "such as the log of a batch job, may",
                                      "show the cause. Running the same",
                                      "calculations with",
                                      "{.code BPPARAM=SerialParam()} may also",
                                      "show it, but a crash would then end",
                                      "the R session.")),
                          parent=e)
        })
}

## bplapply() and bpiterate() with the arguments in '...', see
## .with_dead_worker_check()
#' @importFrom BiocParallel bplapply
.gsva_bplapply <- function(...) .with_dead_worker_check(bplapply(...))

#' @importFrom BiocParallel bpiterate
.gsva_bpiterate <- function(...) .with_dead_worker_check(bpiterate(...))

## apply 'WRAPPED_FUN' to a block of rows or columns of a matrix sent to a
## worker; defined outside the functions processing matrices by blocks, so
## that sending it to a worker does not send their whole matrix
BLOCK_FUN_WRAPPER <- function(block, WRAPPED_FUN, ...) {
    WRAPPED_FUN(block, ..., verbose=FALSE)
}

## blocks of rows (whdim=1) or columns (whdim=2) of the matrix 'X' in main
## memory, with the ranges 'rngs', processed by the workers of 'BPPARAM'. workers
## not forked from this process, such as socket workers, receive the blocks
## one at a time, as they become free, while forked workers, and workers
## processing a matrix stored on disk, take their blocks from 'X' themselves
## with 'FUN_WRAPPER', which would send the whole 'X' to workers not forked
## from this process. 'BPREDO' gives the result of a previous call to
## recompute only its failed elements
#' @importFrom BiocParallel MulticoreParam
#' @importFrom IRanges start end
.bp_blocks <- function(X, whdim, rngs, FUN, FUN_WRAPPER, ..., BPREDO=list(),
                       BPPARAM) {
    if (is(X, "DelayedArray") || is(BPPARAM, "MulticoreParam"))
        return(.gsva_bplapply(rngs, FUN=FUN_WRAPPER, verbose=FALSE,
                              idpbe=NULL, WRAPPED_FUN=FUN, ..., BPREDO=BPREDO,
                              BPPARAM=BPPARAM))

    i <- 0L
    ITER <- function() {
        i <<- i + 1L
        if (i > length(rngs))
            return(NULL)
        rng <- start(rngs[[i]]):end(rngs[[i]])
        if (whdim == 1L) X[rng, , drop=FALSE] else X[, rng, drop=FALSE]
    }
    .gsva_bpiterate(ITER, FUN=BLOCK_FUN_WRAPPER, WRAPPED_FUN=FUN, ...,
                    BPREDO=BPREDO, BPPARAM=BPPARAM)
}

## process the rows of a matrix with a given function FUN, opening parallelism
## through a BiocParallelParam object BPPARAM, when different from NULL, and
## reporting progress using the 'cli' package when possible. the arguments
## 'workfactor', 'outfactor' and 'outextra' give the memory that FUN takes to
## process the rows of 'X' within the maximum main memory 'maxmem', see
## .units_per_block()

#' @importFrom BiocGenerics type
#' @importFrom cli cli_abort cli_progress_bar cli_alert_warning
#' @importFrom BiocParallel bplapply bpnworkers "bpprogressbar<-" bptry
#' @importFrom BiocParallel bpok bpstopOnError "bpstopOnError<-"
#' @importFrom memuse howbig
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixRows <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL,
                               maxmem=Inf, workfactor=2, outfactor=1,
                               outextra=0) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        stopifnot(is(BPPARAM, "BiocParallelParam"))
        nworkers <- bpnworkers(BPPARAM)
    }

    grid <- .rowgridsize(X, nworkers, maxmem, workfactor, outfactor, outextra)
    rir <- .splitRowsInRanges(grid)
    if (length(rir) > 1 && verbose) {
        sze <- howbig(as.numeric(width(rir[[1]])), as.numeric(ncol(X)),
                      representation="dense", type=type(X))
        msg <- sprintf("Splitting calculations in %d chunks of [%d, %d] and %s",
                       length(rir), width(rir[[1]]), ncol(X), as.character(sze))
        cli_alert_info(msg)
    } else if (length(rir) == 1)          ## serial execution in one single call
        return(FUN(X, ..., verbose=verbose))

    FUN_WRAPPER <- function(rowsrng, verbose, idpbe, WRAPPED_FUN, ...) {
        rng <- rowsrng
        if (!is(X, "DelayedMatrix"))
            rng <- start(rowsrng):end(rowsrng)
        res <- WRAPPED_FUN(X[rng, , drop=FALSE], ..., verbose=FALSE)
        if (verbose && is(idpbe, "environment"))
            cli_progress_update(id=get("idpb", envir=idpbe), width(rowsrng))
        return(res)
    }
        
    totalnrows <- nrow(X)
    res <- NULL
    if (is.null(BPPARAM) || nworkers <= 1L) { ## serial execution in chunks
        env <- NULL
        if (verbose) {
            env <- new.env(parent=globalenv())
            assign("idpb", cli_progress_bar(progressmsg, total=totalnrows),
                   envir=env)
        }
        res <- lapply(rir, FUN=FUN_WRAPPER, verbose=verbose,
                      idpbe=env, WRAPPED_FUN=FUN, ...)
        if (verbose)
            cli_progress_done(get("idpb", envir=env))
    } else {                                  ## parallel execution in chunks
        ## 'BPPARAM' is a reference class object, so the following changes
        ## reach the object of the caller, which is restored on exit
        oldprogressbar <- bpprogressbar(BPPARAM)
        oldstoponerror <- bpstopOnError(BPPARAM)
        on.exit({
            bpprogressbar(BPPARAM) <- oldprogressbar
            bpstopOnError(BPPARAM) <- oldstoponerror
        }, add=TRUE)
        if (verbose)
            bpprogressbar(BPPARAM) <- TRUE    ## reporting progress wo/ cli
        bpstopOnError(BPPARAM) <- FALSE
        res <- bptry(.bp_blocks(X, 1L, rir, FUN, FUN_WRAPPER, ...,
                                BPPARAM=BPPARAM))
        bpokmask <- .bp_ok(res)
        if (any(!bpokmask)) {
            .report_parallel_errors(res)
            cli_alert_warning("Trying to execute again the failing thread(s)")
            res <- bptry(.bp_blocks(X, 1L, rir, FUN, FUN_WRAPPER, ...,
                                    BPREDO=res, BPPARAM=BPPARAM))
            bpokmask <- .bp_ok(res)
            if (any(!bpokmask)) {
                .report_parallel_errors(res)
                cli_abort(c("x"="Cancelling execution"))
            }
        }
    }
    res <- do.call("rbind", res)

    return(res)
}

#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges ranges
.splitColsInRanges <- function(grid) {
    cir <- lapply(grid, function(r) ranges(r)[2])
    cir
}

## process the columns of a matrix with a given function FUN, opening parallelism
## through a BiocParallelParam object BPPARAM, when different from NULL, and
## reporting progress using the 'cli' package when possible. the arguments
## 'workfactor', 'outfactor' and 'outextra' give the memory that FUN takes to
## process the columns of 'X' within the maximum main memory 'maxmem', see
## .units_per_block()
## the results of FUN on each block of columns are bound by columns or,
## when 'combine' is a function of two results, reduced with it

#' @importFrom BiocGenerics type
#' @importFrom cli cli_abort
#' @importFrom BiocParallel bplapply bpnworkers "bpprogressbar<-" bptry
#' @importFrom BiocParallel bpok bpstopOnError "bpstopOnError<-"
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixCols <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL,
                               maxmem=Inf, workfactor=2, outfactor=1,
                               outextra=0, combine=NULL) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        stopifnot(is(BPPARAM, "BiocParallelParam"))
        nworkers <- bpnworkers(BPPARAM)
    }

    grid <- .colgridsize(X, nworkers, maxmem, workfactor, outfactor, outextra)
    cir <- .splitColsInRanges(grid)

    if (length(cir) > 1 && verbose) {
        sze <- howbig(as.numeric(nrow(X)), as.numeric(width(cir[[1]])),
                      representation="dense", type=type(X))
        cli_alert_info(sprintf("Splitting calculations in %d chunks of [%d, %d] and %s",
                               length(cir), nrow(X), width(cir[[1]]), as.character(sze)))
    } else if (length(cir) == 1)              ## serial execution in one single call
        return(FUN(X, ..., verbose=verbose))

    FUN_WRAPPER <- function(colsrng, verbose, idpbe, WRAPPED_FUN, ...) {
        rng <- colsrng
        if (!is(X, "DelayedMatrix"))
            rng <- start(colsrng):end(colsrng)
        res <- WRAPPED_FUN(X[, rng, drop=FALSE], ..., verbose=FALSE)
        if (verbose && is(idpbe, "environment"))
            cli_progress_update(id=get("idpb", envir=idpbe), width(colsrng))
        return(res)
    }
        
    totalncols <- ncol(X)
    res <- NULL
    if (is.null(BPPARAM) || nworkers <= 1L) { ## serial execution in chunks
        env <- new.env(parent=globalenv())
        assign("idpb", cli_progress_bar(progressmsg, total=totalncols), envir=env)
        res <- lapply(cir, FUN=FUN_WRAPPER, verbose=verbose,
                      idpbe=env, WRAPPED_FUN=FUN, ...)
        cli_progress_done(get("idpb", envir=env))
    } else {                                  ## parallel execution in chunks
        ## 'BPPARAM' is a reference class object, so the following changes
        ## reach the object of the caller, which is restored on exit
        oldprogressbar <- bpprogressbar(BPPARAM)
        oldstoponerror <- bpstopOnError(BPPARAM)
        on.exit({
            bpprogressbar(BPPARAM) <- oldprogressbar
            bpstopOnError(BPPARAM) <- oldstoponerror
        }, add=TRUE)
        if (verbose)
            bpprogressbar(BPPARAM) <- TRUE    ## reporting progress wo/ cli
        bpstopOnError(BPPARAM) <- FALSE
        res <- bptry(.bp_blocks(X, 2L, cir, FUN, FUN_WRAPPER, ...,
                                BPPARAM=BPPARAM))
        bpokmask <- .bp_ok(res)
        if (any(!bpokmask)) {
            .report_parallel_errors(res)
            cli_alert_warning("Trying to execute again the failing thread(s)")
            res <- bptry(.bp_blocks(X, 2L, cir, FUN, FUN_WRAPPER, ...,
                                    BPREDO=res, BPPARAM=BPPARAM))
            bpokmask <- .bp_ok(res)
            if (any(!bpokmask)) {
                .report_parallel_errors(res)
                cli_abort(c("x"="Cancelling execution"))
            }
        }
    }

    if (is.function(combine))
        return(Reduce(combine, res))

    mines <- maxes <- NULL
    if (!is.null(attr(res[[1]], "min"))) { ## min and max enrichment scores stored by ssGSEA
        mines <- min(vapply(X=res, FUN=function(x) attr(x, "min"), FUN.VALUE=numeric(1)))
        maxes <- max(vapply(X=res, FUN=function(x) attr(x, "max"), FUN.VALUE=numeric(1)))
    }

    res <- do.call("cbind", res)

    if (!is.null(mines)) { ## min and max enrichment scores stored by ssGSEA
        attr(res, "min") <- mines
        attr(res, "max") <- maxes
    }

    return(res)
}
