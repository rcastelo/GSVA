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

## adapted from .define_multiworker_grid() in beachmat/R/colBlockApply.R
#' @importFrom BiocGenerics type
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength
.rowgridsize <- function(X, nworkers=1, maxmem=Inf) {
  typesze <- c("integer"=4, "double"=8) ## 4 bytes for integers, 8 bytes for doubles
  grid <- DummyArrayGrid(dim(X))
  if (!is.infinite(maxmem) || nworkers > 1 || is(X, "DelayedMatrix")) {
      ## initially maximum block length is the maximum of the default auto block
      ## length and the maximum available memory divided by the size of stored number
      ## if no finite maximum available memory is specified, then it becomes the
      ## default auto block length
      max.block.length <- getAutoBlockLength(type(X))
      if (!is.infinite(maxmem))
          max.block.length <- max(max.block.length, ceiling(maxmem / typesze[type(X)]))
      ## assuming all workers share memory, the maximum block length has to reduce
      ## by the number of workers to avoid exceeding the maximum available memory
      ## and, in any case, it cannot exceed .Machine$integer.max
      max.block.length <- min(.Machine$integer.max, max.block.length / nworkers)
      expected.block.length <- max(1, ceiling(nrow(X) / nworkers) * as.numeric(ncol(X)))
      block.length <- min(max.block.length, expected.block.length)
      grid <- rowAutoGrid(X, block.length=block.length)
  }
  grid
}

## adapted from .define_multiworker_grid() in beachmat/R/colBlockApply.R
#' @importFrom BiocGenerics type
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength
.colgridsize <- function(X, nworkers=1, maxmem=Inf) {
  typesze <- c("integer"=4, "double"=8) ## 4 bytes for integers, 8 bytes for doubles
  grid <- DummyArrayGrid(dim(X))
  if (!is.infinite(maxmem) || nworkers > 1 || is(X, "DelayedMatrix")) {
      ## initially maximum block length is the maximum of the default auto block
      ## length and the maximum available memory divided by the size of stored number
      ## if no finite maximum available memory is specified, then it becomes the
      ## default auto block length
      max.block.length <- getAutoBlockLength(type(X))
      if (!is.infinite(maxmem))
          max.block.length <- max(max.block.length, ceiling(maxmem / typesze[type(X)]))
      ## assuming all workers share memory, the maximum block length has to reduce
      ## by the number of workers to avoid exceeding the maximum available memory
      ## and, in any case, it cannot exceed .Machine$integer.max
      max.block.length <- min(.Machine$integer.max, max.block.length / nworkers)
      expected.block.length <- max(1, ceiling(nrow(X) / nworkers) * as.numeric(ncol(X)))
      block.length <- min(max.block.length, expected.block.length)
      grid <- colAutoGrid(X, block.length=block.length)
  }
  grid
}

#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges ranges
.splitRowsInRanges <- function(grid) {
    rir <- lapply(grid, function(r) ranges(r)[1])
    rir
}

## process the rows of a matrix with a given function FUN, opening parallelism
## through a BiocParallelParam object BPPARAM, when different from NULL, and
## reporting progress using the 'cli' package when possible

#' @importFrom BiocGenerics type
#' @importFrom cli cli_abort cli_progress_bar cli_alert_warning
#' @importFrom BiocParallel bplapply bpnworkers bpprogressbar bptry bpok
#' @importFrom memuse howbig
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixRows <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL,
                               maxmem=Inf) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        stopifnot(is(BPPARAM, "BiocParallelParam"))
        nworkers <- bpnworkers(BPPARAM)
    }

    grid <- .rowgridsize(X, nworkers, maxmem)
    rir <- .splitRowsInRanges(grid)
    if (length(rir) > 1 && verbose) {
        sze <- howbig(as.numeric(width(rir[[1]])), as.numeric(ncol(X)),
                      representation="dense", type=type(X))
        cli_alert_info(sprintf("Splitting calculations in %d chunks of [%d, %d] and %s",
                               length(rir), width(rir[[1]]), ncol(X), as.character(sze)))
    } else if (length(rir) == 1)                     ## serial execution in one single call
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
        if (verbose)
            bpprogressbar(BPPARAM) <- TRUE    ## reporting progress wo/ cli
        bptry(res <- bplapply(rir, FUN=FUN_WRAPPER, verbose=FALSE,
                              idpbe=NULL, WRAPPED_FUN=FUN, ...,
                              BPPARAM=BPPARAM))
        bpokmask <- bpok(res)
        if (any(!bpokmask)) {
            msg <- paste("{sum(!bpokmask)} execution thread(s) give an error,",
                         "reporting the first one.")
            cli_alert_warning(msg)
            print(attr(res[[which(!bpokmask)]], "traceback"))
            cli_alert_warning("Trying to execute again the failing thread(s)")
            bptry(res <- bplapply(rir, FUN=FUN_WRAPPER, verbose=FALSE,
                                  idpbe=NULL, WRAPPED_FUN=FUN, ...,
                                  BPREDO=res, BPPARAM=BPPARAM))
            bpokmask <- bpok(res)
            if (any(!bpokmask)) {
                msg <- paste("{sum(!bpokmask)} execution thread(s) give an",
                             "error, reporting the first one.")
                cli_alert_warning(msg)
                print(attr(res[[which(!bpokmask)]], "traceback"))
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
## reporting progress using the 'cli' package when possible

#' @importFrom BiocGenerics type
#' @importFrom cli cli_abort
#' @importFrom BiocParallel bplapply bpnworkers bptry bpok
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixCols <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL,
                               maxmem=Inf) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        if (!is(BPPARAM, "BiocParallelParam")) {
            msg <- paste("'BPPARAM' must be a BiocParallelParam derivative.",
                         "Please consult the 'BiocParallel' package.")
            cli_abort(c("x"=msg))
        }
        nworkers <- bpnworkers(BPPARAM)
    }

    grid <- .colgridsize(X, nworkers, maxmem)
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
        if (verbose)
            bpprogressbar(BPPARAM) <- TRUE    ## reporting progress wo/ cli
        bptry(res <- bplapply(cir, FUN=FUN_WRAPPER, verbose=FALSE,
                              idpbe=NULL, WRAPPED_FUN=FUN, ...,
                              BPPARAM=BPPARAM))
        bpokmask <- bpok(res)
        if (any(!bpokmask)) {
            msg <- paste("{sum(!bpokmask)} execution thread(s) give an error,",
                         "reporting the first one.")
            cli_alert_warning(msg)
            print(attr(res[[which(!bpokmask)]], "traceback"))
            cli_alert_warning("Trying to execute again the failing thread(s)")
            bptry(res <- bplapply(cir, FUN=FUN_WRAPPER, verbose=FALSE,
                                  idpbe=NULL, WRAPPED_FUN=FUN, ...,
                                  BPREDO=res, BPPARAM=BPPARAM))
            bpokmask <- bpok(res)
            if (any(!bpokmask)) {
                msg <- paste("{sum(!bpokmask)} execution thread(s) give an",
                             "error, reporting the first one.")
                cli_alert_warning(msg)
                print(attr(res[[which(!bpokmask)]], "traceback"))
                cli_abort(c("x"="Cancelling execution"))
            }
        }
    }

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
