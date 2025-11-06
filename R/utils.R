## generate dummy names, e.g. row/col names for object M that knows 'nrow()'
.dummyNames <- function(M, n=nrow(M), prefix="row") {
    fmt <- sprintf("%s%%0%dd", prefix, floor(log10(n)) + 1)
    sprintf(fmt, seq_len(n))
}

## check for presence of valid row/feature names
##   and abort or generate dummy names
## #' @importFrom Biobase featureNames
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
    } else if (any(duplicated(rownames(expr)))) {
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

#' @importFrom S4Arrays is_sparse
#' @importFrom sparseMatrixStats rowRanges
#' @importFrom DelayedArray blockApply setAutoBPPARAM rowRanges
#' @importFrom cli cli_alert_warning cli_abort cli_alert_info
#' @importFrom cli cli_progress_bar cli_progress_done
#' @importFrom BiocParallel SerialParam bpnworkers bpiterate bpprogressbar
.filterGenes <- function(expr, removeConstant=TRUE, removeNzConstant=TRUE,
                         verbose=TRUE, BPPARAM=NULL) {
    geneRanges <- NULL
    ## open parallelism only if filtering have to be done for
    ## more than 100 genes on more than 100 samples
    if (!is.null(BPPARAM) && bpnworkers(BPPARAM) > 1 &&
        nrow(expr) > 100 && ncol(expr) > 100) {

        if (is(expr, "DelayedArray")) {     ## DelayedArray::rowRanges() takes care of
                                            ## grid layout with respect to parallelization
            if (verbose) {
                msg <- sprintf("Searching for genes/features with constant values using %d cores",
                               as.integer(bpnworkers(BPPARAM)))
                cli_alert_info(msg)
                bpprogressbar(BPPARAM) <- TRUE
            }
            setAutoBPPARAM(BPPARAM=BPPARAM)
            geneRanges <- rowRanges(expr)
            setAutoBPPARAM(NULL)

        } else {

            n_chunks <- bpnworkers(BPPARAM)
            if (verbose) {
                msg <- sprintf("Searching for genes/features with constant values using %d cores",
                               as.integer(bpnworkers(BPPARAM)))
                cli_alert_info(msg)
                idpb <- cli_progress_bar("Computing row ranges", total=n_chunks)
            }
            bpprogressbar(BPPARAM) <- FALSE
            geneRanges <- bpiterate(.row_iter(expr, idpb, n_chunks),
                                    rowRanges, na.rm=TRUE, useNames=FALSE,
                                    REDUCE=rbind, reduce.in.order=TRUE, BPPARAM=BPPARAM)

            if (verbose)
                cli_progress_done(idpb)
        }
    } else {
        BPPARAM <- NULL ## disable any other BPPARAM configuration
        if (verbose)
            cli_alert_info("Searching for genes/features with constant values")
        geneRanges <- rowRanges(expr, na.rm=TRUE, useNames=FALSE)
    }

    constantGenes <- (geneRanges[, 1] == geneRanges[, 2])

    if (any(constantGenes) || anyNA(constantGenes)) {
        invalidGenes <- (constantGenes | is.na(constantGenes))
        msg <- sprintf("%d genes/features with constant values throughout the samples",
                       sum(invalidGenes))
        cli_alert_warning(msg)
        if (removeConstant) {
            cli_alert_warning("Genes/features with constant values are discarded")
            expr <- expr[!invalidGenes, ]
        }
    }

    if (is_sparse(expr)) {
        if (!is.null(BPPARAM) && bpnworkers(BPPARAM) > 1 &&
            nrow(expr) > 100 && ncol(expr) > 100) {

            if (verbose) {
                msg <- sprintf("Searching for genes/features with constant nonzero values using %d cores",
                               as.integer(bpnworkers(BPPARAM)))
                cli_alert_info(msg)
            }
        } else if (verbose) {
            BPPARAM <- NULL ## disable any other BPPARAM configuration
            cli_alert_info("Searching for genes/features with constant nonzero values")
        }

        nzGeneRanges <- NULL
        if (is(expr, "DelayedMatrix")) {
            rowNonzeroRanges_byBlock <- function(block) {
                t(sapply(t(block)@SVT, function(x) range(x[[1]])))
            }
            nworkers <- 1L
            if (verbose && !is.null(BPPARAM) && bpnworkers(BPPARAM) > 1 &&
                nrow(expr) > 100 && ncol(expr) > 100) {
                bpprogressbar(BPPARAM) <- TRUE
                nworkers <- bpnworkers(BPPARAM)
            }
            nzGeneRanges <- blockApply(expr, rowNonzeroRanges_byBlock,
                                       grid=rowgridsize(expr, nworkers=nworkers),
                                       as.sparse=is_sparse(expr), BPPARAM=BPPARAM)
            nzGeneRanges <- do.call("rbind", nzGeneRanges)
        } else if (is(expr, "dgCMatrix")) {
            if (verbose && !is.null(BPPARAM) && bpnworkers(BPPARAM) > 1 &&
                nrow(expr) > 100 && ncol(expr) > 100) {
                bpprogressbar(BPPARAM) <- TRUE
                nzGeneRanges <- do.call("rbind", bplapply(.sparse2columnList(t(expr)),
                                                          FUN=range, BPPARAM=BPPARAM))
            } else
                nzGeneRanges <- do.call("rbind", lapply(.sparse2columnList(t(expr)),
                                                        FUN=range))
        } else if (is(expr, "SVT_SparseArray")) {
            if (verbose && !is.null(BPPARAM) && bpnworkers(BPPARAM) > 1 &&
                nrow(expr) > 100 && ncol(expr) > 100) {
                bpprogressbar(BPPARAM) <- TRUE
                nzGeneRanges <- do.call("rbind", bplapply(t(expr)@SVT,
                                                          FUN=function(x) range[[1]],
                                                          BPPARAM=BPPARAM))
            } else
                nzGeneRanges <- do.call("rbind", lapply(t(expr)@SVT,
                                                        FUN=function(x) range(x[[1]])))
        } else
            cli_abort("x"="Uknown sparse matrix class")

        stopifnot(is.matrix(nzGeneRanges)) ## QC
        stopifnot(nrow(nzGeneRanges) == nrow(expr)) ## QC

        ## nzGeneRanges <- NULL
        ## if (bpnworkers(BPPARAM) > 1 && nrow(expr) > 100 && ncol(expr) > 100) {
        ##   nzGeneRanges <- do.call("cbind", bplapply(nzGeneList, FUN=range, BPPARAM=BPPARAM))
        ## } else
        ##   nzGeneRanges <- vapply(nzGeneList, FUN=range, FUN.VALUE=double(2))
        constantNzGenes <- (nzGeneRanges[, 1] == nzGeneRanges[, 2])

        if (any(constantNzGenes) || anyNA(constantNzGenes)) {
            invalidNzGenes <- (constantNzGenes | is.na(constantNzGenes))
            msg <- sprintf("%d genes/features with constant nonzero values throughout the samples",
                           sum(invalidNzGenes))
            cli_alert_warning(msg)
            if (removeNzConstant) {
                cli_alert_warning("Genes/features with constant nonzero values are discarded")
                expr <- expr[!invalidNzGenes, ]
            }
        }
    }

    if (nrow(expr) < 2)
        cli_abort(c("x"="Less than two genes in the input assay object"))
    
    return(expr)
}


## maps gene sets content in 'gsets' to 'features', where 'gsets'
## is a 'list' object with character string vectors as elements,
## and 'features' is a character string vector object. it assumes
## features in both input objects follow the same nomenclature,

#' @importFrom cli cli_abort
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
    if(is(geneSets, "list") && all(sapply(geneSets, is.numeric))) {
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
    filteredDataMatrix <- .filterGenes(dataMatrix,
                                       removeConstant=removeConstant,
                                       removeNzConstant=removeNzConstant,
                                       verbose,
                                       BPPARAM=BPPARAM)

    filteredMappedGeneSets <- .filterAndMapGeneSets(param=param,
                                                    filteredDataMatrix=filteredDataMatrix,
                                                    verbose=verbose)

    return(list(filteredDataMatrix=filteredDataMatrix,
                filteredMappedGeneSets=filteredMappedGeneSets))
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
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixRows <- function(X, FUN, ..., progressmsg="Progress",
                               verbose=TRUE, minparrows=100, minparcols=100,
                               BPPARAM=NULL) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        if (!is(BPPARAM, "BiocParallelParam"))
            cli_abort("x"="'BPPARAM' must be a BiocParallelParam derivative")
        nworkers <- bpnworkers(BPPARAM)
    }

    grid <- DummyArrayGrid(dim(X))
    if (nworkers > 1 || is(X, "DelayedMatrix"))
        grid <- rowgridsize(X, nworkers)
    rir <- .splitRowsInRanges(grid)
    if (length(rir) == 1)                     ## serial execution in one single call
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
                cli_abort("x"="Cancelling execution")
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
#' @importFrom BiocParallel bplapply bpnworkers
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixCols <- function(X, FUN, ..., progressmsg="Progress",
                               verbose=TRUE, minparrows=100, minparcols=100,
                               BPPARAM=NULL) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        if (!is(BPPARAM, "BiocParallelParam"))
            cli_abort("x"="'BPPARAM' must be a BiocParallelParam derivative")
        nworkers <- bpnworkers(BPPARAM)
    }
    grid <- DummyArrayGrid(dim(X))
    if (nworkers > 1 || is(X, "DelayedMatrix"))
        grid <- colgridsize(X, nworkers)
    cir <- .splitColsInRanges(grid)
    if (length(cir) == 1)                     ## serial execution in one single call
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
            bptry(res <- bplapply(rir, FUN=FUN_WRAPPER, verbose=FALSE,
                                  idpbe=NULL, WRAPPED_FUN=FUN, ...,
                                  BPREDO=res, BPPARAM=BPPARAM))
            if (any(!bpok(res))) {
                cli_alert_warning(sprintf("%d execution thread(s) give an error, reporting the first one"))
                print(attr(res[[which(!bpokmask)]], "traceback"))
                cli_abort("x"="Cancelling execution")
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
    for (i in 1:grid_dim[1])
        for (j in 1:grid_dim[2]) {
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
            cli_abort("x"=sprintf("%s sparse matrix class cannot be handled", class(X)))

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

.showSome <- function(x) {
    paste0(paste(Biobase::selectSome(x, 4), collapse=", "),
           " (", length(x), " total)")
}

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
