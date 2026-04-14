## ----- private methods for wrapping/unwrapping assay data in containers -----

## unwrapData: extract a data matrix from a container object
setMethod("unwrapData", signature("matrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("dgCMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("SVT_SparseMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("DelayedMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("ExpressionSet"),
          function(container, assay) {
              return(exprs(container))
          })

setMethod("unwrapData", signature("SummarizedExperiment"),
          function(container, assay) {
              if (length(assays(container)) == 0L)
                  stop("The input SummarizedExperiment object has no assay data.")

              if (missing(assay) || is.na(assay)) {
                  assay <- names(assays(container))[1]
              } else {
                  if (!is.character(assay))
                      stop("The 'assay' argument must contain a character string.")

                  assay <- assay[1]

                  if (!assay %in% names(assays(container)))
                      stop(sprintf("Assay %s not found in the input SummarizedExperiment object.", assay))
              }

              return(assays(container)[[assay]])
          })

setMethod("unwrapData", signature("SingleCellExperiment"),
          function(container, assay) {
              if (length(assays(container)) == 0L)
                  stop("The input SingleCellExperiment object has no assay data.")

              if (missing(assay) || is.na(assay)) {
                  assay <- names(assays(container))[1]
              } else {
                  if (!is.character(assay))
                      stop("The 'assay' argument must contain a character string.")

                  assay <- assay[1]

                  if (!assay %in% names(assays(container)))
                      stop(sprintf("Assay %s not found in the input SingleCellExperiment object.", assay))
              }

              return(assays(container)[[assay]])
          })

setMethod("unwrapData", signature("SpatialExperiment"),
          function(container, assay) {
            if (length(assays(container)) == 0L)
              stop("The input SpatialExperiment object has no assay data.")
            
            if (missing(assay) || is.na(assay)) {
              assay <- names(assays(container))[1]
            } else {
              if (!is.character(assay))
                stop("The 'assay' argument must contain a character string.")
              
              assay <- assay[1]
              
              if (!assay %in% names(assays(container)))
                stop(sprintf("Assay %s not found in the input SpatialExperiment object.", assay))
            }
            
            return(assays(container)[[assay]])
          })


## wrapData: put the resulting data and gene sets into the original data container type
setMethod("wrapData", signature(container="matrix"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="dgCMatrix"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="SVT_SparseMatrix"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="DelayedMatrix"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="ExpressionSet"),
          function(container, dataMatrix, geneSets) {
              rval <- new("ExpressionSet", exprs=dataMatrix,
                          phenoData=phenoData(container),
                          experimentData=experimentData(container),
                          annotation="")
              if (!missing(geneSets))
                  attr(rval, "geneSets") <- geneSets
              
              return(rval)
          })

setMethod("wrapData", signature(container="SummarizedExperiment"),
          function(container, dataMatrix, geneSets) {
              rdata <- adata <- NULL
              if (!missing(geneSets)) {
                  adata <- SimpleList(es=dataMatrix)
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## assume missing geneSets imples dataMatrix are ranks
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  adata <- c(assays(container[mask, ]),
                             SimpleList(gsvaranks=dataMatrix))
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SummarizedExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container))
              if (!missing(geneSets))
                  metadata(rval)$annotation <- NULL

              return(rval)
          })

setMethod("wrapData", signature(container="SingleCellExperiment"),
          function(container, dataMatrix, geneSets) {
              rdata <- adata <- NULL
              if (!missing(geneSets)) {
                  adata <- SimpleList(es=dataMatrix)
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## assume missing geneSets imples dataMatrix are ranks
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  adata <- c(assays(container[mask, ]),
                             SimpleList(gsvaranks=dataMatrix))
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SingleCellExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container))
              if (!missing(geneSets))
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })

setMethod("wrapData", signature(container="SpatialExperiment"),
          function(container, dataMatrix, geneSets) {
              rdata <- adata <- NULL
              if (!missing(geneSets)) {
                  adata <- SimpleList(es=dataMatrix)
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## assume missing geneSets imples dataMatrix are ranks
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  adata <- c(assays(container[mask, ]),
                             SimpleList(gsvaranks=dataMatrix))
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SpatialExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container),
                  imgData=imgData(container),
                  spatialCoords=spatialCoords(container))
              if (!missing(geneSets))
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })


## return direct subclasse of a union class, e.g., 'GsvaExprData'
## or 'GsvaGeneSets'

## obtain direct subclasses of a union class for error reporting purposes
#' @importFrom methods isClass getClass is
#' @importFrom cli cli_abort
.getDirectSubclasses <- function(unionClass) {
    if (!isClass(unionClass))
        cli_abort(c("x"=sprintf("'%s' is not a valid union class name",
                                unionClass)))

    subClasses <- getClass(unionClass)@subclasses
    dist <- vapply(subClasses, function(x) x@distance, numeric(1))
    directSubclasses <- names(dist)[dist == 1]
    directSubclasses
}

.check_input_expr_gene_sets <- function(exprData, geneSets) {
    if (!is(exprData, "GsvaExprData")) {
        subclassnames <- .getDirectSubclasses("GsvaExprData")
        msg <- sprintf(paste("argument 'exprData' must be an object of one of",
                            "the following classes: %s"),
                       paste(subclassnames, collapse=", "))
        cli_abort(c("x"=msg))
    }
    if (!is(geneSets, "GsvaGeneSets")) {
        subclassnames <- .getDirectSubclasses("GsvaGeneSets")
        msg <- sprintf(paste("argument 'geneSets' must be an object of one of",
                            "the following classes: %s"),
                       paste(subclassnames, collapse=", "))
        cli_abort(c("x"=msg))
    }
}

## generate dummy names, e.g. row/col names for object M that knows 'nrow()'
.dummyNames <- function(M, n=nrow(M), prefix="row") {
    fmt <- sprintf("%s%%0%dd", prefix, floor(log10(n)) + 1)
    sprintf(fmt, seq_len(n))
}

## check for presence of valid row/feature names
## and abort or generate dummy names
.check_rowNames <- function(expr, useDummyNames=TRUE, verbose) {
    ## CHECK: is this the right place to check this?
    ## 21/10/24: let's do it at parameter constructor
    if (is.null(rownames(expr))) {
        if (useDummyNames) {
            if (verbose) {
                cli_alert_info("Using dummy rownames for the input assay object.")
                rownames(expr) <- .dummyNames(expr)
            }
        } else {
            cli_abort(c("x"="The input assay object doesn't have rownames"))
        }
    } else if (anyDuplicated(rownames(expr))) {
        cli_abort(c("x"="The input assay object has duplicated rownames"))
    }

    return(expr)
}


## check if the expression data object at hand is a multi-assay container
## required by .check_assayNames() (see below) but may be re-used elsewhere
## *currently* only TRUE for SummarizedExperiment and its subclasses, else FALSE
.isMultiAssayContainer <- function(xd) {
    return(is(xd, "SummarizedExperiment"))
}


## check assay parameter and assay names:
##   abort if selected assay name not found in existing assay name list
##   alert if no assay name selected while there are assay names (and select
##     the first one)
##   alert if assay name selected while there are none
##
##   tricky: no assay names AND no assay name selected is the normal case for
##     non-container data objects, e.g. matrix,  or single-assay containers,
##     e.g. ExpressionSet, and should hence be silently accepted.  Multi-assay
##     containers, e.g. SummarizedExperiment, MAY contain more than one assays
##     BUT no assay names.  In this case we also abort because our parameter
##     objects can only store assay names and we prefer requiring assay names
##     for multi-assay containers (which is not unreasonable!) over making
##     things even more complicated for little practical gain.
##
## 2025-03-12  axel: as an afterthought, if we have a multi-assay container with
##   assay names AND no assay is selected AND one of the assay names happens to
##   be 'logcounts' --> use this one by default rather than the first in list.
#' @importFrom utils head
.check_assayNames <- function(a, xd, verbose) {
    an <- gsvaAssayNames(xd)

    if(length(a) != 1) {
        msg <- "argument 'assay' must be of length 1 (it is {length(a)})"
        cli_abort(msg)
    }
    
    if(.isCharNonEmpty(an)) {   # we have assay names
        an <- .omitEmptyChar(an)
        
        if(is.na(a)) {          # but none selected: by default, 
            ## select a common value by default if available -- see afterthought
            ## if unavailable, just select the first available assay name
            def <- grep("logcounts", an, fixed=TRUE, value=TRUE)
            assay <- if(length(def) > 0) head(def, 1) else head(an, 1)
            if (verbose) {
                msg <- "No assay name provided; using default assay '{assay}'"
                cli_alert_info(msg)
            }
        } else {                # check the provided assay name before using it
            if(a %in% an) {
                assay <- a      # found it: OK!
            } else {            # assay name provided but not found: ERROR
                msg <- paste("invalid argument assay='{a}': not part of",
                             "assay names in input argument 'exprData'.")
                cli_abort(msg)
            }
        }
    } else {                    # we don't have no assay names at all
        if(.isMultiAssayContainer(xd)) {  # these must have assay names: ERROR
            msg <- "exprData object of class '{class(xd)}' has no assay names."
            cli_abort(msg)
        } else {                       # i.e. there is exactly one unnamed assay
            if(verbose && !is.na(a)) { # and the provided name is useless but harmless
                msg <- paste("argument assay='{a}' ignored since input argument",
                             "'exprData' has not assay names.")
                cli_alert_info(msg)
            }

            assay <- NA_character_
        }
    }

    return(assay)
}



## converts a dgCMatrix into a list of its columns, based on
## https://rpubs.com/will_townes/sparse-apply
## it is only slightly more efficient than .sparseToList() below BUT simpler
## and does NOT offer converting to a list of rows which is far less efficient
## on a dgCMatrix object.  if you need lists of rows, simply transpose before
## calling this function, t() is reasonably fast as is calling vapply() on its
## result
#' @importFrom Matrix nnzero
.sparse2columnList <- function(m) {
    return(unname(split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))))
}

## actually, it's not just an apply() but also in-place modification
## ellipsis added for cases such as when FUN=rank where we may need
## to set the parameter 'ties.method' of the 'rank()' function
#' @importFrom BiocParallel SerialParam bplapply
.sparseColumnApplyAndReplace <- function(m, FUN, ...) {
    x <- m@x
    x <- lapply(.sparse2columnList(m), FUN=FUN, ...)
    m@x <- unlist(x, use.names=FALSE)
    if (is.integer(m@x)) ## rank(ties.method="first") returns integers
        mode(m@x) <- "numeric" ## dgCMatrix holds only doubles and logicals
    return(m)
}

#' @importFrom cli cli_abort cli_alert_warning
.check_for_na_values <- function(exprData, assay, checkNA, use) {
    autonaclasseswocheck <- c("matrix", "ExpressionSet",
                              "SummarizedExperiment",
                              "RangedSummarizedExperiment")
    mask <- class(exprData) %in% autonaclasseswocheck
    checkNAyesno <- switch(checkNA, yes="yes", no="no",
                           ifelse(any(mask), "yes", "no"))
    didCheckNA <- any_na <- FALSE
    if (checkNAyesno == "yes") {
        any_na <- anyNA(unwrapData(exprData, assay))
        didCheckNA <- TRUE
        if (any_na) {
            if (use == "all.obs")
                cli_abort(c("x"="Input expression data has NA values."))
            else if (use == "everything")
                cli_alert_warning(paste("Input expression data has NA values,",
                                        "which will be propagated through",
                                        "calculations"))
            else ## na.rm
                cli_alert_warning(paste("Input expression data has NA values,",
                                        "which will be discarded from",
                                        "calculations"))
        }
    }

    list(any_na=any_na, didCheckNA=didCheckNA)
}

.get_NAuse <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam"))
  return(object@use)
}

.get_checkNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam"))
  return(object@checkNA)
}

.get_didCheckNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam"))
  return(object@didCheckNA)
}

.get_ondisk <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam"))
  return(object@ondisk)
}

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

#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges ranges
.splitColsInRanges <- function(grid) {
    cir <- lapply(grid, function(r) ranges(r)[2])
    cir
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
        if (!is(BPPARAM, "BiocParallelParam")) {
            msg <- paste("'BPPARAM' must be a BiocParallelParam derivative.",
                         "Please consult the 'BiocParallel' package.")
            cli_abort(c("x"=msg))
        }
        nworkers <- bpnworkers(BPPARAM)
    }

    grid <- .rowgridsize(X, nworkers, maxmem)
    rir <- .splitRowsInRanges(grid)
    if (length(rir) > 1 && verbose) {
        typesze <- c("integer"=4, "double"=8) ## 4 bytes for integers, 8 bytes for doubles
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
        typesze <- c("integer"=4, "double"=8) ## 4 bytes for integers, 8 bytes for doubles
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

#' @importFrom IRanges IRanges
#' @importFrom S4Arrays is_sparse ArrayViewport
#' @importFrom DelayedArray chunkdim chunkGrid
#' @importFrom SparseArray nzcount
#' @importFrom cli cli_alert_info cli_abort
.estimate_nzcount <- function(exprData, assay, verbose) {
    X <- unwrapData(exprData, assay)
    ## coerce to double to ensure we can deal with numbers larger than 2^31
    nr <- as.numeric(nrow(X))
    nc <- as.numeric(ncol(X))
    nzc <- tot <- nr*nc
    if (is_sparse(X)) {
        estimated_flag <- FALSE
        if (is(X, "dgCMatrix") || is(X, "SVT_SparseMatrix"))
            nzc <- nzcount(X)
        else if (is(X, "DelayedMatrix")) {
            if (nc < 2000)
                nzc <- nzcount(as(X, "dgCMatrix"))
            else {
                block_dim <- chunkdim(X)
                if (is.null(block_dim)) {
                    grid <- defaultAutoGrid(X)
                    block_dim <- dim(grid[[1L]])
                }
                block_dim <- c(min(c(nr, block_dim[1])), min(c(nc, block_dim[2]))) ## just in case there's only one block
                vp <- ArrayViewport(dim(X), IRanges(c(1, 1), width=block_dim))     ## just use the first block
                block <- read_block(X, vp)
                nzc <- ceiling(tot * as.numeric(nzcount(block)) / prod(block_dim))
                estimated_flag <- TRUE
            }
        } else
            cli_abort(c("x"="{class(X)} sparse matrix class cannot be handled."))

        if (verbose) {
            estmsg <- ""
            if (estimated_flag)
                estmsg <- " (estimated)"
            cli_alert_info(sprintf("%.0f nonzeros (%s than 2^31) and %.2f%% sparsity%s",
                                   nzc, ifelse(nzc > .Machine$integer.max, "more", "less"),
                                   100 - (100 * nzc / tot), estmsg))
        }
    }

    return(nzc)
}

#' @importFrom cli cli_abort
.memtext2bytes <- function(x) {
  if (is.numeric(x))
      return(x)

  x <- gsub(",", ".", x)
  pat <- "(\\d*(.\\d+)*)(.*)"
  num  <- as.numeric(sub(pat, "\\1", x))
  unit <- sub(pat, "\\3", x)
  unit[unit==""] <- "1"

  fac <- c("1"=1, "K"=1024, "M"=1024^2, "G"=1024^3, "T"=1024^4)
  if (!toupper(unit) %in% names(fac)) {
      msg <- "Unknown memory unit '{unit}', please use either K, M, G or T."
      cli_abort(c("x"=msg))
  }

  num * unname(fac[toupper(unit)])
}

#' @importFrom cli cli_abort cli_alert_info
#' @importFrom memuse Sys.meminfo
.check_maxmem <- function(param, x, verbose) {
    if (length(x) > 1 || (!is.numeric(x) && !is.character(x))) {
        msg <- paste("'maxmem' should be a vector of length 1 of either a",
                     "number in bytes or a character string formed by a",
                     "number followed by the suffix K, M, G or T.")
        cli_abort(c("x"=msg))
    }

    if (is.character(x) && x == "auto") {
        totalram <- Sys.meminfo()$totalram
        maxmem <- as.numeric(totalram * 0.9) ## auto takes 90% of RAM
        X <- unwrapData(get_exprData(param), get_assay(param))
        if (verbose && is(X, "DelayedArray") &&
            gsva_global$show_start_and_end_messages)
            cli_alert_info(sprintf("Maximum available main memory (90%%): %s",
                                   as.character(totalram * 0.9)))
    } else if (is.numeric(x))
        maxmem <- x
    maxmem <- .memtext2bytes(maxmem)
    maxmem
}


## verifies that the 'ondisk' parameter is either 'auto', 'yes' or 'no' and, if
## 'auto', checks whether the input data fits in the maximum available main
## memory and sets 'ondisk' to 'yes' or 'no' accordingly. If the input data is
## a DelayedArray, also reports whether it fits in the maximum available main
## memory. If 'ondisk' is set to 'yes', then this function returns TRUE,
## otherwise it returns FALSE.

#' @importFrom cli cli_abort cli_alert_info
#' @importFrom BiocGenerics type
#' @importFrom S4Arrays is_sparse
#' @importFrom memuse howbig
.check_ondisk <- function(param, maxmem, verbose) {
    ondisk <- .get_ondisk(param)
    if (ondisk == "auto") {
        X <- unwrapData(get_exprData(param), get_assay(param))
        tot <- as.numeric(nrow(X)) * as.numeric(ncol(X))
        rep <- "dense"
        spa <- 1
        if (is_sparse(X)) {
            rep <- "sparse"
            spa <- nzcount(param) / tot
        }
        sze <- howbig(as.numeric(nrow(X)), as.numeric(ncol(X)),
                      representation=rep, sparsity=spa, type=type(X))
        ondisk <- "no"
        if (as.numeric(sze) > maxmem) {
            ondisk <- "yes"
            if (is(X, "DelayedArray") && verbose) {
                msg <- paste("On-disk input data does not fit in the maximum",
                             "available main memory")
                cli_alert_info(msg)
            }
        } else if (is(X, "DelayedArray") && verbose) {
            msg <- paste("On-disk input data fits in the maximum available",
                         "main memory")
            cli_alert_info(msg)
        }

    } else if (ondisk != "yes" && ondisk != "no")
        cli_abort(c("x"="'ondisk' should be either 'auto', 'yes' or 'no'"))

    ondisk == "yes"
}


## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.isPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
  (name %in% loadedNamespaces())
}

.objPkgClass <- function(obj) {
    oc <- class(obj)
    pkg <- attr(oc, "package", exact=TRUE)
    opc <- if(is.null(pkg)) {
               oc[1]
           } else {
               paste(pkg[1], oc[1], sep = "::")
           }
    return(opc)
}

#' @importFrom Biobase selectSome
.showSome <- function(x) {
    paste0(paste(selectSome(x, 4), collapse=", "),
           " (", length(x), " total)")
}

#' @importFrom utils capture.output
.catObj <- function(x, prefix = "  ") {
    if(is.null(x)) {
        cat(paste0(prefix, "none."))
    } else {
        cat(paste0(prefix, capture.output(gsvaShow(x))), sep="\n")
    }
}

.isCharNonEmpty <- function(x) {
    return((!is.null(x)) &&
           (length(x) > 0) &&
           (is.character(x)) &&
           (!all(is.na(x))) &&
           any(nzchar(x)))
}

.omitEmptyChar <- function(x) {
    if(.isCharNonEmpty(x)) {
        return(x[(nzchar(x)) & (!is.na(x))])
    } else {
        return(character(0))
    }
}

.isCharLength1 <- function(x) {
    return((.isCharNonEmpty(x)) && (length(x) == 1))
}

## annotation package checks
.isAnnoPkgValid <- function(ap) {
    return(.isCharLength1(ap))
}

#' @importFrom utils installed.packages
.isAnnoPkgInstalled <- function(ap) {
    ap <- c(ap, paste0(ap, ".db"))
    return(any(ap %in% rownames(installed.packages())))
}
