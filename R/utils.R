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
        msg <- sprintf("argument 'assay' must be of length 1 (is %d)", length(a))
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
                msg <- sprintf("No assay name provided; using default assay '%s'",
                               assay)
                cli_alert_info(msg)
            }
        } else {                # check the provided assay name before using it
            if(a %in% an) {
                assay <- a      # found it: OK!
            } else {            # assay name provided but not found: ERROR
                msg <- sprintf(paste0("invalid argument assay='%s': not part of ",
                                      "exprData's assay name list."), a)
                cli_abort(msg)
            }
        }
    } else {                    # we don't have no assay names at all
        if(.isMultiAssayContainer(xd)) {  # these must have assay names: ERROR
            msg <- sprintf("exprData object of class '%s' has no assay names.",
                           class(xd))
            cli_abort(msg)
        } else {                       # i.e. there is exactly one unnamed assay
            if(verbose && !is.na(a)) { # and the provided name is useless but harmless
                msg <- sprintf(paste0("argument assay='%s' ignored since exprData ",
                                      "has no assay names."), a)
                cli_alert_info(msg)
            }

            assay <- NA_character_
        }
    }

    return(assay)
}



## 2024-02-06  axel: function .filterGenes() is intended to detect genes (rows)
##  with constant expression (and, hence, no information), warn about them and
##  optionally remove them (in particular, ssGSEA's choice is to keep them).
##  the original approach tried to identify genes with a standard deviation of
##  exactly 0 but failed in certain cases of all identical values due to the use
##  of floating point arithmetic, see issues:
## https://github.com/rcastelo/GSVA/issues/54
## https://github.com/HenrikBengtsson/matrixStats/issues/204
##  an improvement in matrixStats::rowSds() fixes the original issue but cannot
##  guarantee that there won't be other problematic cases.
##
## We propose to detect cases of constant gene expression by comparing genewise
##  min and max values rather than computing the SD, which *should* avoid using
##  floating point arithmetic in favour of comparisons and scale linearly with
##  the number of samples (columns) -- we'll of course have to check that. ;-)
##
## A related but different issue has recently surfaced when methods PLAGE and
##  z-scores are applied to sparse matrices and attempt to scale the non-zero
##  values of genes: genes that are constant in their non-zero values will have
##  an SD of 0 and therefore scaling them will result in division by 0.

.rowNzRanges_dgCMatrix <- function(X, verbose=FALSE) {
    res <- .Call("row_rngs_nzrngs_RsparseMatrix_R", as(X, "RsparseMatrix"),
                 verbose=verbose)
    res
}

.rowNzRanges_SVT_SparseArray_byrow <- function(X, verbose=FALSE) {
    res <- .Call("row_rngs_nzrngs_SVT_SparseMatrix_R", X, verbose=verbose)
    res
}

.rowNzRanges_SVT_SparseArray_transpose_C <- function(X, verbose=FALSE) {
    res <- .Call("col_rngs_nzrngs_SVT_SparseMatrix_R", t(X), verbose=verbose)
    res
}

## after discussions at from https://github.com/Bioconductor/SparseArray/issues/22

.rowNzRanges_SVT_SparseArray <- function(X, verbose=FALSE) {
    res <- .Call("rowbycols_rngs_nzrngs_SVT_SparseMatrix_R", X, verbose=verbose)
    res
}

#' @importFrom SparseArray NaArray
.fast_replace_zeros_with_NAs <- function(x) {
    stopifnot(is(x, "SparseArray"))
    naa <- NaArray(dim=dim(x), type=type(x), dimnames=dimnames(x))
    naa@NaSVT <- x@SVT ## ASSUMING x@SVT HAS NO NA VALUES!!
    naa
}

#' @importFrom SparseArray NaArray nzwhich
.safe_replace_zeros_with_NAs <- function(x) {
    naa <- NaArray(dim=dim(x), type=type(x), dimnames=dimnames(x))
    nzidx <- nzwhich(x)
    naa[nzidx] <- x[nzidx]
    naa
}

#' @importFrom MatrixGenerics rowMins rowMaxs
.rowNzRanges_SVT_SparseArray_rowbycols_R <- function(X, anyna=FALSE, verbose=FALSE) {
    naa <- NULL
    if (anyna)
        naa <- .safe_replace_zeros_with_NAs(X)
    else
        naa <- .fast_replace_zeros_with_NAs(X)  # only if 'X' is guaranteed to be NA-free!

    ranges1 <- cbind(rowMins(X, na.rm=TRUE), rowMaxs(X, na.rm=TRUE))
    ranges2 <- suppressWarnings(cbind(rowMins(naa, na.rm=TRUE), rowMaxs(naa, na.rm=TRUE)))
    allzeros <- ranges1[ , 1L] == 0L & ranges1[ , 2L] == 0L
    ranges2[allzeros] <- NA_integer_
    cbind(ranges1, ranges2)
}

#' @importFrom S4Arrays DummyArrayGrid read_block
.rowNzRanges <- function(X, anyna=FALSE, verbose=FALSE) {
    res <- NULL
    if (is.matrix(X))
        res <- rowRanges(X, na.rm=TRUE)
    else if (is(X, "dgCMatrix"))
        res <- .rowNzRanges_dgCMatrix(X, verbose=verbose)
    else if (is(X, "SVT_SparseArray"))
        res <- .rowNzRanges_SVT_SparseArray(X, verbose=verbose)
    else if (is(X, "DelayedArray")) {
        grid <- DummyArrayGrid(dim(X))
        block <- read_block(X, grid[[1L]])
        if (is_sparse(block)) ## input HDF5 may be sparse or not
            res <- .rowNzRanges_SVT_SparseArray(block, verbose=verbose)
        else
            res <- rowRanges(X)
    } else
        cli_abort(c("x"=sprintf(".rowNzRanges: input object class %s not handled yet",
                                class(X))))
    res
}



#' @importFrom S4Arrays is_sparse
#' @importFrom sparseMatrixStats rowRanges
#' @importFrom DelayedArray blockApply setAutoBPPARAM rowRanges
#' @importFrom cli cli_alert_warning cli_abort cli_alert_info
#' @importFrom cli cli_progress_bar cli_progress_done
#' @importFrom BiocParallel SerialParam bpnworkers bpiterate bpprogressbar
.filterGenes <- function(expr, anyna=FALSE, removeConstant=TRUE,
                         removeNzConstant=TRUE, verbose=TRUE, BPPARAM=NULL,
                         maxmem=Inf) {
    rowrngs <- NULL

    if (verbose) {
        if (!is_sparse(expr))
            cli_alert_info("Searching for rows with constant values")
        else
            cli_alert_info("Searching for rows with constant (nonzero) values")
    }

    ## returns a matrix with as many rows as 'expr' and 2 columns if 'expr'
    ## is dense, and 4 columns if it is sparse, where the first two columns
    ## correspond to the minimum and maximum values of each row, while the
    ## third and fourth columns, if they exist, they correspond to the
    ## minimum and maximum nonzero values of each row, which will be NAs if
    ## there are no nonzero values.
    rowrngs <- .processMatrixRows(expr, .rowNzRanges, anyna=anyna,
                                  verbose=verbose, BPPARAM=BPPARAM, maxmem=maxmem)

    constantRows <- (rowrngs[, 1] == rowrngs[, 2])
    mask <- is.na(constantRows)
    if (any(mask))
        constantRows[mask] <- TRUE

    constantNzRows <- invalidRows <- invalidNzRows <- rep(FALSE, nrow(expr))
    if (ncol(rowrngs) > 2) { ## sparse input
        constantNzRows <- (rowrngs[, 3] == rowrngs[, 4])
        mask <- is.na(constantNzRows)
        if (any(mask)) ## no nonzero values imply constant nonzero values
            constantNzRows[mask] <- TRUE
    }

    if (verbose && any(constantRows)) {
        msg <- sprintf("%d rows with constant values throughout the columns",
                       sum(constantRows))
        cli_alert_warning(msg)
        if (removeConstant)
           cli_alert_warning("Rows with constant values are discarded")
    }

    nzmask <- constantNzRows & !constantRows
    if (verbose && any(nzmask)) {
        msg <- sprintf("%d rows with constant nonzero values throughout the samples",
                       sum(nzmask))
        cli_alert_warning(msg)
        if (removeNzConstant)
           cli_alert_warning("Rows with constant nonzero values are discarded")
    }

    removemask <- rep(FALSE, nrow(expr))
    if (removeConstant)
         removemask <- constantRows

    if (removeNzConstant && any(nzmask))
         removemask <- removemask | nzmask

    if (any(removemask)) {
        if (nrow(expr) - sum(removemask) < 2)
            cli_abort(c("x"="Less than two rows left in the input assay object"))

        expr <- expr[!removemask, ]
    }

    return(expr)
}


## maps gene sets content in 'gsets' to 'features', where 'gsets'
## is a 'list' object with character string vectors as elements,
## and 'features' is a character string vector object. it assumes
## features in both input objects follow the same nomenclature,

#' @importFrom cli cli_abort
#' @importFrom IRanges CharacterList match
.mapGeneSetsToFeatures <- function(gsets, features) {

    ## Aaron Lun's suggestion at
    ## https://github.com/rcastelo/GSVA/issues/39#issuecomment-765549620
    gsets2 <- CharacterList(gsets)
    mt <- match(gsets2, features)
    mapdgenesets <- as.list(mt[!is.na(mt)])

    if (length(unlist(mapdgenesets, use.names=FALSE)) == 0) {
      msg <- paste("No identifiers in the gene sets could be matched to the",
                   "identifiers in the expression data.")
      cli_abort(c("x"=msg))
    }

    mapdgenesets
}

## it assumes that all arguments have been already checked for correctness
#' @importFrom cli cli_abort cli_alert_warning
.filterAndMapGeneSets <- function(param, wgset=NA, filteredDataMatrix, verbose) {

    minSize <- get_minSize(param)
    maxSize <- get_maxSize(param)

    geneSets <- get_geneSets(param)
    if (!is.na(wgset))
        geneSets <- geneSets[wgset]

    ## we'll try to handle index lists of numeric/integer vectors as gene sets
    if (is(geneSets, "list") && all(vapply(geneSets, is.numeric, logical(1)))) {
        mappedGeneSets <- lapply(geneSets, function(idx) {
            as.integer(idx[idx > 0 & idx <= nrow(filteredDataMatrix)])
        })

        ## check and alert if we had to drop out-of-range indices
        diffGs <- names(geneSets)[lengths(geneSets) != lengths(mappedGeneSets)]
        if(length(diffGs) > 0) {
            singular <- length(diffGs) == 1
            msg <- sprintf(
                paste0("Out-of-range indices from %d index gene %s (%s) ",
                       "have been dropped."),
                length(diffGs),
                if(singular) "set" else "sets",
                paste0(sQuote(diffGs, q=FALSE), collapse = ", "))
            cli_alert_warning(msg)
        }
    } else { # not a list of index vectors, i.e., as before
        ## note that the method for 'GeneSetCollection' calls geneIds(), i.e., 
        ## whatever the input, from here on we have a list of character vectors
        anno <- get_annotation(param)
        if (identical(anno, NullIdentifier()))
            anno <- NULL
        geneSets <- mapGeneSetsToAnno(geneSets=geneSets,
                                      anno=anno,
                                      verbose=verbose)
        
        ## map to the actual features for which expression data is available
        ## note that the result is a list of integer vectors (indices to
        ## rownames) and not a list of character vector any longer
        mappedGeneSets <- .mapGeneSetsToFeatures(geneSets,
                                                 rownames(filteredDataMatrix))
    }
    
    ## remove gene sets from the analysis for which no features are available
    ## and meet the minimum and maximum gene-set size specified by the user
    filteredMappedGeneSets <- filterGeneSets(mappedGeneSets,
                                             minSize=minSize,
                                             maxSize=maxSize)
    
    if (length(filteredMappedGeneSets) == 0) {
        msg <- "No gene set left after mapping and filtering."
        cli_abort(c("x"=msg))
    }

    ## this should NEVER happen -- just to make sure it doesn't...
    if (anyDuplicated(names(filteredMappedGeneSets)) > 0) {
        msg <- "The gene set list contains duplicated gene set names."
        cli_abort(c("x"=msg))
    }

    if (any(lengths(filteredMappedGeneSets) == 1)) {
        msg <- "Some gene sets have size one. Consider setting minSize > 1"
        cli_alert_warning(msg)
    }

    return(filteredMappedGeneSets)
}

#' @importFrom cli cli_alert_warning
#' @importFrom BiocParallel SerialParam
.filterAndMapGenesAndGeneSets <- function(param,
                                          removeConstant=TRUE,
                                          removeNzConstant=TRUE,
                                          verbose=FALSE,
                                          BPPARAM=SerialParam()) {
    exprData <- get_exprData(param)
    dataMatrix <- unwrapData(exprData, get_assay(param))
    
    ## filter genes according to various criteria,
    ## e.g., constant expression
    filtDataMatrix <- .filterGenes(dataMatrix, anyna=anyNA(param),
                                   removeConstant=removeConstant,
                                   removeNzConstant=removeNzConstant,
                                   verbose, BPPARAM=BPPARAM)

    filtMappedGeneSets <- .filterAndMapGeneSets(param=param,
                                                filteredDataMatrix=filtDataMatrix,
                                                verbose=verbose)

    return(list(filteredDataMatrix=filtDataMatrix,
                filteredMappedGeneSets=filtMappedGeneSets))
}


## (re-)extract a list of gene names from a list of indices
## (indices resulting from the matching above)
.geneSetsIndices2Names <- function(indices, names) {
    return(lapply(indices, function(i, n) n[i], n=names))
}


## access to gene set attribute without explicit use of attributes
.geneSets <- function(obj) {
    gs <- attr(obj, "geneSets", exact=TRUE)

    if (is.null(gs))
        stop("The object does not contain information about gene sets.")

    return(gs)
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

## adapted from .define_multiworker_grid() in beachmat/R/colBlockApply.R
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength type
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
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength type
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

#' @importFrom cli cli_abort cli_progress_bar cli_alert_warning
#' @importFrom BiocParallel bplapply bpnworkers bpprogressbar bptry bpok
#' @importFrom memuse howbig
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixRows <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL, maxmem=Inf) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        if (!is(BPPARAM, "BiocParallelParam"))
            cli_abort(c("x"="'BPPARAM' must be a BiocParallelParam derivative"))
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
        res <- WRAPPED_FUN(X[rng, ], ..., verbose=FALSE)
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
            cli_alert_warning(sprintf("%d execution thread(s) give an error, reporting the first one"))
            print(attr(res[[which(!bpokmask)]], "traceback"))
            cli_alert_warning("Trying to execute again the failing thread(s)")
            bptry(res <- bplapply(rir, FUN=FUN_WRAPPER, verbose=FALSE,
                                  idpbe=NULL, WRAPPED_FUN=FUN, ...,
                                  BPREDO=res, BPPARAM=BPPARAM))
            if (any(!bpok(res))) {
                cli_alert_warning(sprintf("%d execution thread(s) give an error, reporting the first one"))
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

#' @importFrom cli cli_abort
#' @importFrom BiocParallel bplapply bpnworkers bptry bpok
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixCols <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL, maxmem=Inf) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        if (!is(BPPARAM, "BiocParallelParam"))
            cli_abort(c("x"="'BPPARAM' must be a BiocParallelParam derivative"))
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
        res <- WRAPPED_FUN(X[, rng], ..., verbose=FALSE)
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
            cli_alert_warning(sprintf("%d execution thread(s) give an error, reporting the first one"))
            print(attr(res[[which(!bpokmask)]], "traceback"))
            cli_alert_warning("Trying to execute again the failing thread(s)")
            bptry(res <- bplapply(cir, FUN=FUN_WRAPPER, verbose=FALSE,
                                  idpbe=NULL, WRAPPED_FUN=FUN, ...,
                                  BPREDO=res, BPPARAM=BPPARAM))
            if (any(!bpok(res))) {
                cli_alert_warning(sprintf("%d execution thread(s) give an error, reporting the first one"))
                print(attr(res[[which(!bpokmask)]], "traceback"))
                cli_abort(c("x"="Cancelling execution"))
            }
        }
    }
    res <- do.call("cbind", res)

    return(res)
}

## calculate number of nonzero values in an on-disk DelayedArray
.nzcountDA <- function(X) {
    ## coerce to double to ensure we can deal with numbers larger than 2^31
    nr <- as.numeric(nrow(X))
    nc <- as.numeric(ncol(X))
    block_dim <- chunkdim(X)
    grid_dim <- dim(chunkGrid(X))
    nzc <- 0
    for (i in seq_len(grid_dim[1]))
        for (j in seq_len(grid_dim[2])) {
            icoord <- (i-1)*block_dim[1]+1
            jcoord <- (j-1)*block_dim[2]+1
            bdim <- c(min(c(nr, i*block_dim[1]))-icoord+1, min(c(nc, j*block_dim[2]))-jcoord+1)
            vp <- ArrayViewport(dim(X), IRanges(c(icoord, jcoord), width=bdim))
            block <- read_block(X, vp)
            nzc <- nzc + as.numeric(nzcount(block))
        }
    nzc
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
        if (is(X, "dgCMatrix") || is(X, "SVT_SparseArray"))
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
            cli_abort(c("x"=sprintf("%s sparse matrix class cannot be handled", class(X))))

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
  if (!toupper(unit) %in% names(fac))
      cli_abort(c("x"=sprintf("Unknown memory unit '%s', please use either K, M, G or T", unit)))

  num * unname(fac[toupper(unit)])
}

#' @importFrom cli cli_abort cli_alert_info
#' @importFrom memuse Sys.meminfo
.check_maxmem <- function(param, x, verbose) {
    if (length(x) > 1 || (!is.numeric(x) && !is.character(x)))
        cli_abort(c("x"="'maxmem' should be a vector of length 1 of either a number in bytes or a character string"))

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

#' @importFrom cli cli_abort cli_alert_info
#' @importFrom S4Arrays is_sparse
#' @importFrom memuse howbig
.check_ondisk <- function(param, maxmem, verbose) {
    ondisk <- get_ondisk(param)
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
            if (is(X, "DelayedArray") && verbose)
                cli_alert_info("On-disk input data does not fit in the maximum available main memory")
        } else if (is(X, "DelayedArray") && verbose)
            cli_alert_info("On-disk input data fits in the maximum available main memory")

    } else if (ondisk != "yes" && ondisk != "no")
        cli_abort(c("x"="'ondisk' should be either 'auto', 'yes' or 'no'"))

    ondisk
}


## transforms a dgCMatrix into a list of its
## non-zero values by MARGIN (1 for row, 2 for column)
##
## currently unused because replaced by .sparse2columnList() and kept for the
## time being, just in case...
## .sparseToList <-function(dgCMat, MARGIN){
##   MARGIN <- as.integer(MARGIN)
##   J <- rep(1:ncol(dgCMat), diff(dgCMat@p))
##   I <- dgCMat@i + 1
##   x <- dgCMat@x
##   if (MARGIN == 1L) {
##     result <- split(x, I)
##     names(result) <- rownames(dgCMat)[as.numeric(names(result))]
##   } else if (MARGIN == 2L) {
##     result <- split(x, J)
##     names(result) <- colnames(dgCMat)[as.numeric(names(result))]
##   }
##   else {
##     warning("invalid MARGIN; return NULL")
##     result <- NULL
##   }
##   result
## }

## .dgCapply<-function(m, MARGIN, FUN){
##   x <- lapply(.sparseToList(m, MARGIN), FUN)
##   m@x <- unlist(x, use.names=FALSE)
##   m
## }


.guessIfCountData <- function(x, tolerance = sqrt(.Machine$double.eps)) {
    return(typeof(x) == "integer" ||
           (all(x >= 0) && all(x - round(x) < tolerance)))
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

## utility function to make sure matrix to dgCMatrix coercion is uniform
## (since direct coercion to dgCMatrix is deprecated (!) by Matrix pkg)
.matrix2dgCMatrix <- function(m) {
    return(as(as(as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix"))
}


### 2024-08-02  axel: the following three functions have been copied from
### GSEABase/R/utilities.R (v. 1.66.0) as our implementation of readGMT()
### is mostly based on (a copy of) GSEABase::getGmt() which is making use
### of these utility functions.  Since we decided that GSVA::readGMT() may
### return a list of gene sets as well as a GeneSetCollection, it should work
### if a user doesn't have GSEABase installed at all.

## Placeholder 'till something appropriate decided
.uniqueIdentifier <- local({
    node <- NULL
    pid <- NULL
    uid <- 0L
    function() {
        if (is.null(node)) {
            node <<- Sys.info()['nodename']
            pid <<- Sys.getpid()
        }
        uid <<- uid + 1L
        base::paste(node, pid, date(), uid, sep=":")
    }
})

.stopf <- function(...) {
    call <- match.call(call=sys.call(sys.parent(1)))
    msg <- paste(sprintf(...), collapse="\n")
    stop(simpleError(msg, call=call))
}

.warningf <- function(...) {
    call <- match.call(call=sys.call(sys.parent(1)))
    msg <- paste(sprintf(...), collapse="\n")
    warning(simpleWarning(msg, call=call))
}
### end of copy from GSEABase/R/utilities.R
