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


## metadata recorded in every GSVA output object: the version of GSVA that
## produced it, 'gsvaVersion', to handle in the future objects produced by
## different versions of GSVA, e.g., when they are stored on disk, and for
## column ranks, the number of rows on which they were calculated,
## 'ranksNrow', since the ranks are only valid for those rows; for other
## outputs 'ranksNrow' is NULL, which removes it when inherited from the input
#' @importFrom utils packageDescription
.gsva_output_metadata <- function(dataMatrix, assay) {
    list(gsvaVersion=packageDescription("GSVA")[["Version"]],
         ranksNrow=if (assay == "gsvaranks") nrow(dataMatrix) else NULL)
}

.wrapdata_nonSE <- function(dataMatrix, param, assay, first, last, rem, whdim,
                            dropAssays, geneSets) {
    stopifnot(!missing(param))
    stopifnot(!missing(assay))
    stopifnot(!missing(first))
    stopifnot(!missing(last))
    stopifnot(!missing(rem))
    stopifnot(!missing(whdim))
    stopifnot(!missing(dropAssays))
    omd <- .gsva_output_metadata(dataMatrix, assay)
    attr(dataMatrix, "gsvaParam") <- .gsvaParam_as_list(param)
    attr(dataMatrix, "assay") <- assay
    attr(dataMatrix, "gsvaVersion") <- omd$gsvaVersion
    attr(dataMatrix, "ranksNrow") <- omd$ranksNrow
    if (!is.na(first) || !is.na(last))
        attr(dataMatrix, "restrict") <- list(first=first, last=last,
                                             rem=rem, whdim=whdim)
    if (!missing(geneSets))
        attr(dataMatrix, "geneSets") <- geneSets
    return(dataMatrix)
}

## wrapData: put the resulting data and gene sets into the original data
##           container type
## parameters:
##  - container is the original data container object
##  - dataMatrix is a matrix of values to be wrapped into the container
##  - param is the parameter object used for the calculation that produced
##    dataMatrix
##  - assay is the name of the assay to be used in the container object to
##    store dataMatrix
##  - first if not NA, it is the index of the first row/column in the original
##    data container that was used to produce dataMatrix
##  - last if not NA, it is the index of the last row/column in the original
##    data container that was used to produce dataMatrix
##  - rem if not NA, it is the number of rows/columns that were discarded from
##    the original data container, e.g., when filtering out constant rows
##  - whdim if not NA, it is the dimension along which first and last are
##    defined, i.e., 1 for rows and 2 for columns
##  - dropAssays if TRUE and container uses assays, the original assays in the
##    container are dropped and only the new assay is kept. If FALSE, the new
##    assay is added to the existing assays
##  - geneSets if not missing, it is a list of gene sets that were used to
##    produce dataMatrix.
setMethod("wrapData", signature(container="matrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="dgCMatrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="SVT_SparseMatrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="DelayedMatrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="ExpressionSet"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              pdata <- phenoData(container)
              if (!is.na(first) && !is.na(last) && whdim == 2)
                  pdata <- pdata[first:last, , drop=FALSE]
              rval <- new("ExpressionSet", exprs=dataMatrix,
                          phenoData=pdata,
                          experimentData=experimentData(container),
                          annotation="")
              omd <- .gsva_output_metadata(dataMatrix, assay)
              attr(rval, "gsvaParam") <- .gsvaParam_as_list(param)
              attr(rval, "assay") <- assay
              attr(rval, "gsvaVersion") <- omd$gsvaVersion
              attr(rval, "ranksNrow") <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  attr(rval, "restrict") <- list(first=first, last=last,
                                                 rem=rem, whdim=whdim)
              if (!missing(geneSets))
                  attr(rval, "geneSets") <- geneSets
              
              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment assayNames
setMethod("wrapData", signature(container="SummarizedExperiment"),
          function(container, dataMatrix, param, assay, first, last, rem,
                   whdim, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## missing geneSets implies adding an assay
                  if (assay %in% assayNames(container))
                      cli_abort(c("x"=paste("Assay {assay} already exists in",
                                            "the input container object.")))
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays) {
                      if (!is.na(first) && !is.na(last) && whdim == 2)
                          adata <- c(assays(container[mask, first:last]), adata)
                      else
                          adata <- c(assays(container[mask, ]), adata)
                  }
                  rdata <- rowData(container)[mask, ]
              }
              cdata <- colData(container)
              if (!is.na(first) && !is.na(last) && whdim == 2)
                  cdata <- cdata[first:last, , drop=FALSE]
              rval <- SummarizedExperiment(assays=adata, colData=cdata,
                                           rowData=rdata,
                                           metadata=metadata(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              omd <- .gsva_output_metadata(dataMatrix, assay)
              metadata(rval)$gsvaVersion <- omd$gsvaVersion
              metadata(rval)$ranksNrow <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  metadata(rval)$restrict <- list(first=first, last=last,
                                                  rem=rem, whdim=whdim)
              if (!missing(geneSets)) ## row data has been replaced
                  metadata(rval)$annotation <- NULL

              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims altExps
setMethod("wrapData", signature(container="SingleCellExperiment"),
          function(container, dataMatrix, param, assay, first, last, rem,
                   whdim, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              rdimdata <- reducedDims(container)
              aexpsdata <- altExps(container)
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
                  rdimdata <- aexpsdata <- SimpleList()
              } else { ## missing geneSets implies adding an assay
                  if (assay %in% assayNames(container))
                      cli_abort(c("x"=paste("Assay {assay} already exists in",
                                            "the input container object.")))
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays) {
                      if (!is.na(first) && !is.na(last) && whdim == 2) {
                          adata <- c(assays(container[mask, first:last]), adata)
                          rdimdata <- lapply(rdimdata,
                                             function(x)
                                                 x[first:last, , drop=FALSE])
                          aexpsdata <- lapply(aexpsdata,
                                              function(x)
                                                  x[, first:last, , drop=FALSE])
                      } else
                          adata <- c(assays(container[mask, ]), adata)
                  }
                  rdata <- rowData(container)[mask, ]
              }
              cdata <- colData(container)
              if (!is.na(first) && !is.na(last) && whdim == 2)
                  cdata <- cdata[first:last, , drop=FALSE]
              rval <- SingleCellExperiment(assays=adata, colData=cdata,
                                           rowData=rdata, reducedDims=rdimdata,
                                           altExps=aexpsdata,
                                           metadata=metadata(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              omd <- .gsva_output_metadata(dataMatrix, assay)
              metadata(rval)$gsvaVersion <- omd$gsvaVersion
              metadata(rval)$ranksNrow <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  metadata(rval)$restrict <- list(first=first, last=last,
                                                  rem=rem, whdim=whdim)
              if (!missing(geneSets)) ## row data has been replaced
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims altExps
#' @importFrom SpatialExperiment SpatialExperiment
setMethod("wrapData", signature(container="SpatialExperiment"),
          function(container, dataMatrix, param, assay, first, last, rem,
                   whdim, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## missing geneSets implies adding an assay
                  if (assay %in% assayNames(container))
                      cli_abort(c("x"=paste("Assay {assay} already exists in",
                                            "the input container object.")))
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays) {
                      if (!is.na(first) && !is.na(last) && whdim == 2)
                          adata <- c(assays(container[mask, first:last]), adata)
                      else
                          adata <- c(assays(container[mask, ]), adata)
                  }
                  rdata <- rowData(container)[mask, ]
              }
              cdata <- colData(container)
              if (!is.na(first) && !is.na(last) && whdim == 2)
                  cdata <- cdata[first:last, , drop=FALSE]
              rval <- SpatialExperiment(
                  assays=adata,
                  colData=cdata,
                  rowData=rdata,
                  reducedDims=reducedDims(container),
                  altExps=altExps(container),
                  metadata=metadata(container),
                  imgData=imgData(container),
                  spatialCoords=spatialCoords(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              omd <- .gsva_output_metadata(dataMatrix, assay)
              metadata(rval)$gsvaVersion <- omd$gsvaVersion
              metadata(rval)$ranksNrow <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  metadata(rval)$restrict <- list(first=first, last=last,
                                                  rem=rem, whdim=whdim)
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
        msg <- paste("Argument 'BPPARAM' must be a 'BiocParallelParam'",
                     "derivative. Please consult the BiocParallel package.")
        cli_abort(c("x"=msg))
    }
}

## load into main memory as an SVT_SparseMatrix object the sparse on-disk data
## in 'X', reading one block of columns at a time, because coercing 'X' as a
## whole with as() takes 5 to 8 times the memory of the result. reading a
## block as sparse with the HDF5Array package takes about 160 bytes per
## nonzero value while it lasts, as measured on single-cell data, and the
## blocks of columns already read form the result, which takes their memory
## when bound. blocks span the width of the chunks of 'X' on disk, when it has
## them, because blocks splitting chunks make their reading slower, while
## wider blocks do not make it faster: the smallest multiple of that width
## not narrower than the default block of the DelayedArray package, which
## they have otherwise. blocks are narrower when their memory does not
## fit in the fraction .mem_fraction_R of the maximum main memory 'maxmem'
## left by the result, whose estimated size is 'insize', when both are known
#' @importFrom DelayedArray colAutoGrid read_block chunkdim
#' @importFrom BiocGenerics type
.load_sparse_by_blocks <- function(X, maxmem=Inf, insize=NA_real_) {
    if (ncol(X) == 0L)
        return(as(X, "SVT_SparseMatrix"))

    ncolblock <- ncol(colAutoGrid(X)[[1L]]) ## default block size
    cd <- chunkdim(X)
    if (!is.null(cd) && cd[2] > 0) ## multiple of the width of the chunks
        ncolblock <- cd[2] * ceiling(ncolblock / cd[2])
    if (!is.null(maxmem) && !is.null(insize) && is.finite(maxmem) &&
        !is.na(insize)) {
        eltbytes <- if (type(X) == "integer") 4 else 8
        nzpercol <- insize / (eltbytes + 4) / ncol(X)
        avail <- .mem_fraction_R * maxmem - insize
        ncolblock <- min(ncolblock, floor(avail / (160 * max(1, nzpercol))))
    }
    ## the DelayedArray package does not support blocks with more than
    ## .Machine$integer.max values
    ncolblock <- min(ncolblock, floor(.Machine$integer.max / max(1, nrow(X))))
    ncolblock <- as.integer(max(1, min(ncol(X), ncolblock)))
    grid <- colAutoGrid(X, ncol=ncolblock)
    blocks <- lapply(seq_along(grid), function(i)
                         as(read_block(X, grid[[i]], as.sparse=TRUE),
                            "SVT_SparseMatrix"))
    res <- do.call(cbind, blocks)
    dimnames(res) <- dimnames(X)

    res
}

.check_sparse_load_input_expr <- function(expr, method, first, last, whdim,
                                          ondisk, verbose) {
    ## methods that deal with sparse data, which is loaded in main memory as
    ## sparse, while the other methods load it as a dense matrix
    sparsemethods <- c("GSVA", "average")
    if (!method %in% sparsemethods && is_sparse(expr)) { 
        msg <- paste("Input expression data is sparse, but the {method}",
                     "algorithm does not deal with sparsity in a specific way,",
                     "and data will be converted into a dense matrix format")
        cli_alert_warning(msg)
    }

    if (!is.na(first) || !is.na(last)) {
        if (is.na(first))
            first <- 1
        if (is.na(last))
            last <- if (whdim == 1) nrow(expr) else ncol(expr)
        if (whdim == 1) {
            if (verbose) {
                msg <- "Restricting to rows {first}:{last}"
                cli_alert_info(msg)
            }
            expr <- expr[first:last, , drop=FALSE]
        } else if (whdim == 2) {
            if (verbose) {
                msg <- "Restricting to columns {first}:{last}"
                cli_alert_info(msg)
            }
            expr <- expr[, first:last, drop=FALSE]
        } else
            cli_abort(c("x"="Invalid internal value for 'whdim' argument."))
    }

    if (is(expr, "DelayedMatrix") && !ondisk) {
        if (verbose)
            cli_alert_info("Loading input expression data into main memory")

        if (method %in% sparsemethods && is_sparse(expr))
            expr <- .load_sparse_by_blocks(expr, maxmem=attr(ondisk, "maxmem"),
                                           insize=attr(ondisk, "insize"))
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
        if (!is.character(assay))
            cli_abort(c("x"="The 'assay' argument must be a character string."))
        assay <- assay[1]
        if (!assay %in% names(assays(container)))
            cli_abort(c("x"=paste("assay {assay} not found in the input",
                                  "{container_class} object.")))
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
            inherits(object, "gsvaParam") ||
            inherits(object, "avgParam"))
  return(object@use)
}

.get_checkNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam") ||
            inherits(object, "avgParam"))
  return(object@checkNA)
}

.get_didCheckNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam") ||
            inherits(object, "avgParam"))
  return(object@didCheckNA)
}

.get_ondisk <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam") ||
            inherits(object, "zscoreParam") ||
            inherits(object, "plageParam") ||
            inherits(object, "avgParam"))
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

## memory limit in bytes of the job or container where this R process runs,
## or Inf when there is none. on Linux, workload managers such as SLURM, and
## container engines, enforce this limit through the control group (cgroup)
## of the process, which can be set in any of the cgroups along its path, e.g.,
## at the level of a SLURM job and not of its steps, so that the smallest
## limit along that path is taken, using either cgroup version 1 or 2. the
## limit given by SLURM in its environment variables is also taken into
## account. the paths of the files and the environment variables are
## arguments to allow testing this function
.job_memory_limit <- function(procfile="/proc/self/cgroup",
                              cgroupdir="/sys/fs/cgroup",
                              env=Sys.getenv(c("SLURM_MEM_PER_NODE",
                                               "SLURM_MEM_PER_CPU",
                                               "SLURM_CPUS_ON_NODE"),
                                             unset=NA)) {
    readnum <- function(fname) {
        if (!file.exists(fname))
            return(NA_real_)
        suppressWarnings(tryCatch(as.numeric(readLines(fname, n=1L,
                                                       warn=FALSE)),
                                  error=function(e) NA_real_))
    }
    limit <- Inf
    cglines <- character(0)
    if (file.exists(procfile)) ## only on Linux
        cglines <- tryCatch(readLines(procfile, warn=FALSE),
                            error=function(e) character(0))
    for (cgline in cglines) {
        ## each line has the format 'hierarchy-ID:controllers:path', where
        ## cgroup version 2 has hierarchy-ID 0 and no controllers
        fields <- strsplit(cgline, ":", fixed=TRUE)[[1]]
        if (length(fields) < 3L)
            next
        if (fields[1] == "0" && fields[2] == "") {
            dir <- cgroupdir
            limitfile <- "memory.max" ## 'max' when there is no limit
        } else if ("memory" %in% strsplit(fields[2], ",", fixed=TRUE)[[1]]) {
            dir <- file.path(cgroupdir, "memory")
            limitfile <- "memory.limit_in_bytes"
        } else
            next
        path <- strsplit(paste(fields[-(1:2)], collapse=":"), "/",
                         fixed=TRUE)[[1]]
        path <- path[path != ""]
        for (k in rev(seq_len(length(path) + 1L) - 1L)) {
            v <- readnum(do.call(file.path,
                                 as.list(c(dir, path[seq_len(k)],
                                           limitfile))))
            if (!is.na(v) && v > 0)
                limit <- min(limit, v)
        }
    }

    ## SLURM gives memory limits in megabytes, where 0 means no limit
    env <- suppressWarnings(as.numeric(env))
    names(env) <- c("SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU",
                    "SLURM_CPUS_ON_NODE")
    if (!is.na(env["SLURM_MEM_PER_NODE"]) && env["SLURM_MEM_PER_NODE"] > 0)
        limit <- min(limit, env["SLURM_MEM_PER_NODE"] * 1024^2)
    else if (!is.na(env["SLURM_MEM_PER_CPU"]) && env["SLURM_MEM_PER_CPU"] > 0 &&
             !is.na(env["SLURM_CPUS_ON_NODE"]))
        limit <- min(limit, env["SLURM_MEM_PER_CPU"] *
                            env["SLURM_CPUS_ON_NODE"] * 1024^2)

    unname(limit)
}

#' @importFrom cli cli_abort cli_alert_info
#' @importFrom memuse Sys.meminfo mu
.check_maxmem <- function(param, assay=get_assay(param), maxmem, verbose) {
    if (length(maxmem) > 1 || (!is.numeric(maxmem) && !is.character(maxmem))) {
        msg <- paste("'maxmem' should be a vector of length 1 of either a",
                     "number in bytes or a character string formed by a",
                     "number followed by the suffix K, M, G or T.")
        cli_abort(c("x"=msg))
    }

    if (is.character(maxmem) && maxmem == "auto") {
        ## auto takes 90% of RAM or, when smaller, of the memory limit of the
        ## job or container where this R process runs
        totalram <- Sys.meminfo()$totalram
        joblimit <- .job_memory_limit()
        avail <- totalram
        what <- "main memory"
        if (joblimit < as.numeric(totalram)) {
            avail <- mu(joblimit)
            what <- "memory of the job"
        }
        maxmem <- as.numeric(avail * 0.9)
        X <- unwrapData(get_exprData(param), assay)
        if (verbose && is(X, "DelayedArray") &&
            gsva_global$show_start_and_end_messages)
            cli_alert_info(sprintf("Maximum available %s (90%%): %s", what,
                                   as.character(avail * 0.9)))
    }

    maxmem <- .memtext2bytes(maxmem)
    maxmem
}


## estimated size in bytes of the input data of a step in main memory, i.e.,
## of the assay 'assay' of the expression data in 'param', restricted to the
## rows (whdim=1) or columns (whdim=2) from 'first' to 'last', when given,
## where sparse data, in an SVT_SparseMatrix object, store each nonzero value
## and its row index, unless 'dense=TRUE', when they are loaded in main memory
## as a dense matrix, as the methods other than GSVA and average do
#' @importFrom cli cli_abort
#' @importFrom BiocGenerics type
#' @importFrom S4Arrays is_sparse
.input_mem_size <- function(param, assay, first, last, whdim,
                            recompute_nzcount=FALSE, dense=FALSE) {
    X <- unwrapData(get_exprData(param), assay)
    nr <- as.numeric(nrow(X))
    nc <- as.numeric(ncol(X))
    spa <- 1
    if (is_sparse(X) && !dense) {
        tot <- nr * nc
        if (recompute_nzcount)
            spa <- .estimate_nzcount(get_exprData(param), assay, FALSE) / tot
        else {
            spa <- nzcount(param) / tot
            ## the number of nonzero values in the parameter object
            ## may correspond to a larger data set, e.g., when columns
            ## have been removed, and then it has to be recomputed
            if (spa > 1)
                spa <- .estimate_nzcount(get_exprData(param), assay,
                                         FALSE) / tot
        }
        spa <- min(spa, 1)
    }
    if (!is.na(first) || !is.na(last)) {
        if (whdim == 1)
            nr <- as.numeric(last - first + 1)
        else if (whdim == 2)
            nc <- as.numeric(last - first + 1)
        else
            cli_abort(c("x"="Invalid internal value for 'whdim' argument."))
    }
    eltbytes <- if (type(X) == "integer") 4 else 8
    if (is_sparse(X) && !dense)
        return(spa * nr * nc * (eltbytes + 4))

    nr * nc * eltbytes
}

## memory in bytes taken by the data of a step processing in blocks a matrix
## with dimensions 'dims', along its rows (whdim=1) or columns (whdim=2), of
## size 'insize' in main memory and 'eltbytes' bytes per value, sparse or not,
## and stored in main memory or not ('inmemory'): the memory that remains
## allocated, see .fixed_mem(), and the working memory of the rows or columns
## of a block of the default size of the DelayedArray package, shared by all
## the workers, see .units_per_block(), which take their size in main memory
## or, when read from disk, the one of their dense form. 'mf' gives the memory
## factors of the step, see .step_mem_factors()
#' @importFrom DelayedArray getAutoBlockSize
.step_data_mem <- function(dims, whdim, insize, eltbytes, sparse, inmemory,
                           mf) {
    nunits <- max(1, dims[whdim])
    denseunitbytes <- as.numeric(dims[-whdim]) * eltbytes
    unitbytes <- if (inmemory) insize / nunits else denseunitbytes
    minunits <- min(nunits, max(1, floor(getAutoBlockSize() / denseunitbytes)))
    .fixed_mem(nunits, insize / nunits, whdim, sparse, inmemory,
               mf$outfactor, mf$outextra, mf$assembly) +
        mf$workfactor * minunits * unitbytes
}

## verifies that the 'ondisk' parameter is either 'auto', 'yes' or 'no' and, if
## 'auto', checks whether the calculations fit in the maximum available main
## memory with the input data in it, and sets 'ondisk' to 'yes' or 'no'
## accordingly. the calculations fit when the memory taken by their data, see
## .step_data_mem(), fits in the fraction .mem_fraction_R of 'maxmem'
## available to the allocations of R, where 'mf' gives the memory factors of
## the step, see .step_mem_factors(), and 'dense=TRUE' indicates that the input
## data is loaded in main memory as a dense matrix, even when it is sparse. If
## the input data is a DelayedArray, also reports whether the calculations fit
## in the maximum available main memory.
## If 'ondisk' is set to 'yes', then this function returns TRUE, otherwise it
## returns FALSE, with the estimated size of the input data in main memory and
## 'maxmem' as attributes.

#' @importFrom cli cli_abort cli_alert_info cli_alert_warning
#' @importFrom S4Arrays is_sparse
#' @importFrom BiocGenerics type
.check_ondisk <- function(param, assay=get_assay(param), first, last, whdim,
                          recompute_nzcount=FALSE, maxmem, verbose,
                          mf=list(workfactor=2, outfactor=1, outextra=0),
                          dense=FALSE) {
    ondisk <- .get_ondisk(param)
    if (ondisk != "auto" && ondisk != "yes" && ondisk != "no")
        cli_abort(c("x"="'ondisk' should be either 'auto', 'yes' or 'no'"))

    insize <- .input_mem_size(param, assay, first, last, whdim,
                              recompute_nzcount, dense)
    if (ondisk == "auto") {
        X <- unwrapData(get_exprData(param), assay)
        dims <- dim(X)
        if (!is.na(first) || !is.na(last))
            dims[whdim] <- last - first + 1
        need <- .step_data_mem(dims, whdim, insize,
                               if (type(X) == "integer") 4 else 8,
                               is_sparse(X) && !dense, TRUE, mf)
        ondisk <- "no"
        if (need > .mem_fraction_R * maxmem) {
            ondisk <- "yes"
            if (is(X, "DelayedArray") && verbose)
                cli_alert_info(paste("Calculations with on-disk input data",
                                     "loaded in main memory do not fit in",
                                     "the maximum available main memory"))
            ## steps giving enrichment scores return them on disk, which
            ## changes the class of their output
            if (mf$outextra > 0)
                cli_alert_warning(paste("The calculations do not fit in the",
                                        "maximum available main memory, and",
                                        "the resulting (dense) matrix of",
                                        "enrichment scores will be returned",
                                        "using an on-disk data structure"))
        } else if (is(X, "DelayedArray") && verbose)
            cli_alert_info(paste("Calculations with on-disk input data",
                                 "loaded in main memory fit in the maximum",
                                 "available main memory"))
    }

    ## the estimated size of the input data and the maximum main memory allow
    ## .check_sparse_load_input_expr() to load the input data into main
    ## memory by blocks that fit in the memory left by it
    structure(ondisk == "yes", insize=insize, maxmem=maxmem)
}

## memory in bytes that each R process running GSVA takes by itself, the main
## one and each of its parallel workers: loading GSVA and the packages it
## depends on takes about 0.8 GB, which forked workers end up taking too,
## because the garbage collector of R writes to the memory pages that they
## share with the main process. it can be set with the option 'GSVA.workermem'
.worker_mem <- function() {
    as.numeric(getOption("GSVA.workermem", 0.8 * 1024^3))
}

## memory in bytes that each R process takes outside R while it reads blocks
## of input data from disk, such as the buffers of the HDF5 library to
## decompress the chunks read: about 0.2 GB, measured on single-cell data
## stored in chunks of about 1000 x 1000 values
.disk_read_mem <- function() {
    0.2 * 1024^3
}

## number of parallel workers of 'BPPARAM' that process a matrix with
## dimensions 'dims', which is 0 when its calculations are not parallelized,
## see .check_open_parallelism()
#' @importFrom BiocParallel bpnworkers
.n_par_workers <- function(BPPARAM, dims, minparrows=100, minparcols=100) {
    if (is.null(BPPARAM) || bpnworkers(BPPARAM) <= 1 ||
        dims[1] <= minparrows || dims[2] <= minparcols)
        return(0L)

    as.integer(bpnworkers(BPPARAM))
}

## estimated memory in bytes required by a step processing in blocks a matrix
## with dimensions 'dims', along its rows (whdim=1) or columns (whdim=2), of
## size 'insize' in main memory and 'eltbytes' bytes per value, sparse or not,
## and stored in main memory or not ('inmemory'), with 'nworkers' parallel
## workers, where 'mf' gives the memory factors of the step, see
## .step_mem_factors(). the main R process and each worker take .worker_mem(),
## and the data of the step take .step_data_mem() from the fraction
## .mem_fraction_R of the memory available to the allocations of R. each R
## process also takes .disk_read_mem() when the input data is read from disk,
## i.e., when it is not loaded in main memory ('inmemory=FALSE'). forked
## workers do not take a copy of the input data in main memory: the garbage
## collector of R writes only to the memory page with the header of each
## object, so that the data of large vectors remain shared, as measured on
## Linux, while the many small objects of the packages loaded are copied, which
## .worker_mem() includes
.step_mem_need <- function(dims, insize, whdim, eltbytes, sparse, inmemory,
                           mf, nworkers) {
    permem <- .worker_mem()
    if (!inmemory)
        permem <- permem + .disk_read_mem()

    (1 + nworkers) * permem +
        .step_data_mem(dims, whdim, insize, eltbytes, sparse, inmemory, mf) /
        .mem_fraction_R
}

## warn, before the calculations of a step start, when the estimated memory
## 'need' that they require, with 'nworkers' parallel workers, exceeds the
## maximum main memory 'maxmem', suggesting fewer workers or a larger 'maxmem',
## or what 'hint' says, when given
#' @importFrom cli cli_warn
#' @importFrom memuse mu
.check_mem_need <- function(need, maxmem, nworkers, hint=NULL) {
    if (need <= maxmem)
        return(invisible(FALSE))

    needtxt <- as.character(mu(need))
    maxmemtxt <- as.character(mu(maxmem))
    wmtxt <- as.character(mu(.worker_mem()))
    msg <- c("!"=paste("The memory estimated for these calculations,",
                       "{needtxt}, exceeds the maximum available main memory",
                       "of {maxmemtxt}."))
    if (nworkers > 0)
        msg <- c(msg,
                 "i"=paste("The main R process and each of its {nworkers}",
                           "parallel workers take about {wmtxt} by themselves,",
                           "for loading GSVA and the packages it depends on."))
    if (!is.null(hint))
        msg <- c(msg, "i"=hint)
    else if (nworkers > 0)
        msg <- c(msg, "i"=paste("Consider using fewer parallel workers in",
                                "{.arg BPPARAM}, or a larger {.arg maxmem}."))
    else
        msg <- c(msg, "i"="Consider using a larger {.arg maxmem}.")
    cli_warn(msg)

    invisible(TRUE)
}

## check the memory required by a step processing in blocks the input data 'X',
## along its rows (whdim=1) or columns (whdim=2), restricted to those from
## 'first' to 'last', when given, with the on-disk decision 'ondisk' given by
## .check_ondisk(), the memory factors 'mf' of the step, the parallel back-end
## 'BPPARAM' and the maximum main memory 'maxmem', see .check_mem_need(). It is
## skipped while gsva() runs the steps, because gsva() checks the memory
## required by all of them before starting, and with the option
## 'GSVA.check_memory=FALSE'
#' @importFrom S4Arrays is_sparse
#' @importFrom BiocGenerics type
.check_step_mem <- function(X, whdim, first, last, ondisk, mf, BPPARAM,
                            maxmem) {
    if (!gsva_global$check_memory || !getOption("GSVA.check_memory", TRUE))
        return(invisible(FALSE))

    dims <- dim(X)
    if (!is.na(first) || !is.na(last))
        dims[whdim] <- last - first + 1
    eltbytes <- if (type(X) == "integer") 4 else 8
    insize <- attr(ondisk, "insize")
    if (is.null(insize)) ## e.g., data in Parquet files, processed from disk
        insize <- prod(as.numeric(dims)) * eltbytes
    nworkers <- .n_par_workers(BPPARAM, dims)
    need <- .step_mem_need(dims, insize, whdim, eltbytes, is_sparse(X),
                           !ondisk, mf, nworkers)

    .check_mem_need(need, maxmem, nworkers)
}

## check the memory required by each of the three steps that gsva() runs on
## the parameter object 'param', see .check_mem_need(), before running them,
## estimating the size of the input data of each step from the one of the
## previous step, and the on-disk decision of each step as .check_ondisk()
## does
#' @importFrom BiocGenerics type
#' @importFrom S4Arrays is_sparse
.check_gsva_mem <- function(param, BPPARAM, maxmem) {
    if (!getOption("GSVA.check_memory", TRUE))
        return(invisible(FALSE))

    X <- unwrapData(get_exprData(param), get_assay(param))
    dims <- dim(X)
    sparse <- is_sparse(X)
    nworkers <- .n_par_workers(BPPARAM, dims)
    ondisk <- .get_ondisk(param)
    ngs <- length(get_geneSets(param))
    insize <- .input_mem_size(param, get_assay(param), NA, NA, 1L)
    clr <- .get_rowNorm(param) == "clr"
    dgc <- is(X, "dgCMatrix")
    steps <- list(list(step="rownorm", whdim=1L, int=(type(X) == "integer")),
                  list(step="colranks", whdim=2L, int=FALSE),
                  list(step="scores", whdim=2L, int=TRUE))
    need <- 0
    for (st in steps) {
        mf <- .step_mem_factors(st$step, ngs=ngs, sparse=sparse, int=st$int,
                                clr=clr, dgc=dgc)
        eltbytes <- if (st$int) 4 else 8
        inmemory <- ondisk == "no"
        if (ondisk == "auto")
            inmemory <- .step_data_mem(dims, st$whdim, insize, eltbytes,
                                       sparse, TRUE, mf) <=
                        .mem_fraction_R * maxmem
        need <- max(need, .step_mem_need(dims, insize, st$whdim, eltbytes,
                                         sparse, inmemory, mf, nworkers))
        insize <- insize * mf$outfactor ## size of the input of the next step
    }

    .check_mem_need(need, maxmem, nworkers)
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

.isNumNonEmpty <- function(x) {
    return((!is.null(x)) &&
           (length(x) > 0) &&
           (is.numeric(x)) &&
           (!all(is.na(x))) &&
           (!all(is.nan(x))))
}

.isNumLength1 <- function(x) {
    return((.isNumNonEmpty(x)) && (length(x) == 1))
}

#' @importFrom cli cli_abort
.check_first_last_values <- function(X, dimfun, dimname, first, last, rmdt) {
    dimfun <- match.fun(dimfun)

    if (all(is.na(first)) && all(is.na(last)))
        return(list(first=NA_real_, last=NA_real_))

    if (!is.null(rmdt))
        cli_abort(c("x"=paste("Input expression already has chunk boundary",
                              "metadata and 'first' and 'last' cannot be set.")))

    if (all(is.na(first)))
        first <- 1
    if (all(is.na(last)))
        last <- dimfun(X)

    if (.isNumLength1(first) && .isNumLength1(last)) {
        if (first < 1 || last < 1 || first != as.integer(first) ||
            last != as.integer(last))
            cli_abort(c("x"="'first' and 'last' must be positive integers."))
        if (first > last)
            cli_abort(c("x"="'first' must be smaller or equal than 'last'."))
    } else
        cli_abort(c("x"="'first' and 'last' must be a single number."))

    if (first > dimfun(X) || last > dimfun(X))
        cli_abort(c("x"=paste("'first' and 'last' must be smaller or equal than",
                              "the number of {dimname} in the input data")))

    return(list(first=first, last=last))
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
