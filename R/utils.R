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

#' @importFrom cli cli_abort
setMethod("unwrapData", signature("SummarizedExperiment"),
          function(container, assay) {
              assay <- .check_unwrapping_assay(container,
                                               "SummarizedExperiment", assay)

              return(assays(container)[[assay]])
          })

#' @importFrom cli cli_abort
setMethod("unwrapData", signature("SingleCellExperiment"),
          function(container, assay) {
              assay <- .check_unwrapping_assay(container,
                                               "SingleCellExperiment", assay)

              return(assays(container)[[assay]])
          })

setMethod("unwrapData", signature("SpatialExperiment"),
          function(container, assay) {
              assay <- .check_unwrapping_assay(container,
                                               "SpatialExperiment", assay)

              return(assays(container)[[assay]])
          })


## wrapData: put the resulting data and gene sets into the original data container type
setMethod("wrapData", signature(container="matrix"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              attr(dataMatrix, "gsvaParam") <- .gsvaParam_as_list(param)
              attr(dataMatrix, "assay") <- assay
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="dgCMatrix"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              attr(dataMatrix, "gsvaParam") <- .gsvaParam_as_list(param)
              attr(dataMatrix, "assay") <- assay
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="SVT_SparseMatrix"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              attr(dataMatrix, "gsvaParam") <- .gsvaParam_as_list(param)
              attr(dataMatrix, "assay") <- assay
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="DelayedMatrix"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              attr(dataMatrix, "gsvaParam") <- .gsvaParam_as_list(param)
              attr(dataMatrix, "assay") <- assay
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="ExpressionSet"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              rval <- new("ExpressionSet", exprs=dataMatrix,
                          phenoData=phenoData(container),
                          experimentData=experimentData(container),
                          annotation="")
              attr(rval, "gsvaParam") <- .gsvaParam_as_list(param)
              attr(rval, "assay") <- assay
              if (!missing(geneSets))
                  attr(rval, "geneSets") <- geneSets
              
              return(rval)
          })

.check_existing_assay <- function(container, assay) {
    if (assay %in% assayNames(container))
        cli_abort(c("x"=paste("Assay {assay} already exists in the input",
                              "container object.")))
}

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
setMethod("wrapData", signature(container="SummarizedExperiment"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## missing geneSets implies adding an assay
                  .check_existing_assay(container, assay)
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays)
                      adata <- c(assays(container[mask, ]), adata)
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SummarizedExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              if (!missing(geneSets)) ## row data has been replaced
                  metadata(rval)$annotation <- NULL

              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("wrapData", signature(container="SingleCellExperiment"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## missing geneSets implies adding an assay
                  .check_existing_assay(container, assay)
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays)
                      adata <- c(assays(container[mask, ]), adata)
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SingleCellExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              if (!missing(geneSets)) ## row data has been replaced
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("wrapData", signature(container="SpatialExperiment"),
          function(container, dataMatrix, param, assay, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## missing geneSets implies adding an assay
                  .check_existing_assay(container, assay)
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays)
                      adata <- c(assays(container[mask, ]), adata)
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SpatialExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container),
                  imgData=imgData(container),
                  spatialCoords=spatialCoords(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              if (!missing(geneSets)) ## row data has been replaced
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

.check_bpparam <- function(BPPARAM) {
    if (!is(BPPARAM, "BiocParallelParam")) {
        msg <- paste("Argument 'BPPARAM' must be a",
                     "'BiocParallelParam' derivative. Please",
                     "consult the BiocParallel package.")
        cli_abort(c("x"=msg))
    }
}

.check_sparse_load_input_expr <- function(expr, method, ondisk, verbose) {
    if (method != "GSVA" && is_sparse(expr)) { 
        msg <- paste("Input expression data is sparse, but the {method}",
                     "algorithm does not deal with sparsity",
                     "in any specific way, and data will be",
                     "converted into a dense matrix format")
        cli_alert_warning(msg)
    }

    if (is(expr, "DelayedMatrix") && !ondisk) {
        if (verbose)
            cli_alert_info("Loading input expression data into main memory")

        if (method == "GSVA" && is_sparse(expr))
            expr <- as(expr, "SVT_SparseMatrix")
        else
            expr <- as.matrix(expr)

        if (ncol(expr) > 10000) ## free up ASAP memory we need not anymore
          out <- gc()           ## and was allocated when reading a big object
    } 
 
    expr
}

#' @importFrom BiocParallel bpnworkers
#' @importFrom cli cli_alert_info
.check_open_parallelism <- function(expr, BPPARAM, minparrows, minparcols,
                                    verbose) {
    if (bpnworkers(BPPARAM) > 1 && nrow(expr) > minparrows &&
        ncol(expr) > minparcols) {
        if (verbose) {
            msg <- sprintf("Using a %s parallel back-end with %d workers",
                           class(BPPARAM), bpnworkers(BPPARAM))
                      cli_alert_info(msg)
        }
    } else
        BPPARAM <- NULL

    BPPARAM
}

#' @importFrom memuse howbig
#' @importFrom cli cli_alert_warning
.check_es_memory_requirements <- function(expr, gsets, ondisk, maxmem) {
    esreqmem <- howbig(as.numeric(length(gsets)), as.numeric(ncol(expr)),
                       representation="dense", sparsity=1, type="double")

    if (esreqmem > maxmem) {
        msg <- paste("The resulting (dense) matrix of enrichment scores will",
                     "not fit in the given maximum main memory size, and it",
                     "will be returned using an on-disk data structure")
        cli_alert_warning(msg)
        ondisk <- TRUE
    }

    ondisk
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
        cli_abort(c("x"=msg))
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
                cli_abort(c("x"=msg))
            }
        }
    } else {                    # we don't have no assay names at all
        if(.isMultiAssayContainer(xd)) {  # these must have assay names: ERROR
            msg <- "exprData object of class '{class(xd)}' has no assay names."
            cli_abort(c("x"=msg))
        } else {                       # i.e. there is exactly one unnamed assay
            if(verbose && !is.na(a)) { # and the provided name is useless but harmless
                msg <- paste("argument assay='{a}' ignored since input argument",
                             "'exprData' has no assay names.")
                cli_alert_info(msg)
            }

            assay <- NA_character_
        }
    }

    return(assay)
}

#' @importFrom cli cli_abort
.check_unwrapping_assay <- function(container, container_class, assay) {

    if (length(assays(container)) == 0L) {
        msg <- "The input {container_class} object has no assay data."
        cli_abort(c("x"=msg))
    }

    if (missing(assay) || is.na(assay)) {
        assay <- names(assays(container))[1]
    } else {
        if (!is.character(assay)) {
            msg <- "The 'assay' argument must contain a character string."
            cli_abort(c("x"=msg))
        }
        assay <- assay[1]
        if (!assay %in% names(assays(container))) {
            msg <- paste("Assay {assay} not found in the input",
                         "{container_class} object.")
            cli_abort(c("x"=msg))
        }
    }

    assay
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
            inherits(object, "gsvaParam") ||
            inherits(object, "zscoreParam") ||
            inherits(object, "plageParam"))
  return(object@ondisk)
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
            nzc <- as.numeric(nzcount(X))
        else if (is(X, "DelayedMatrix")) {
            if (nc < 2000)
                nzc <- as.numeric(nzcount(as(X, "dgCMatrix")))
            else {
                block_dim <- chunkdim(X)
                if (is.null(block_dim)) {
                    grid <- defaultAutoGrid(X)
                    block_dim <- dim(grid[[1L]])
                }
                block_dim <- c(min(c(nr, block_dim[1])), ## just in case there
                               min(c(nc, block_dim[2]))) ## is only one block
                ## just use the first block
                vp <- ArrayViewport(dim(X), IRanges(c(1, 1), width=block_dim))
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
.check_maxmem <- function(param, assay=get_assay(param), maxmem, verbose) {
    if (length(maxmem) > 1 || (!is.numeric(maxmem) && !is.character(maxmem))) {
        msg <- paste("'maxmem' should be a vector of length 1 of either a",
                     "number in bytes or a character string formed by a",
                     "number followed by the suffix K, M, G or T.")
        cli_abort(c("x"=msg))
    }

    if (is.character(maxmem) && maxmem == "auto") {
        totalram <- Sys.meminfo()$totalram
        maxmem <- as.numeric(totalram * 0.9) ## auto takes 90% of RAM
        X <- unwrapData(get_exprData(param), assay)
        if (verbose && is(X, "DelayedArray") &&
            gsva_global$show_start_and_end_messages)
            cli_alert_info(sprintf("Maximum available main memory (90%%): %s",
                                   as.character(totalram * 0.9)))
    }

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
.check_ondisk <- function(param, assay=get_assay(param), maxmem, verbose) {
    ondisk <- .get_ondisk(param)
    if (ondisk == "auto") {
        X <- unwrapData(get_exprData(param), assay)
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
