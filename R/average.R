##
## methods for the average method
##

#' @importFrom S4Arrays is_sparse
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @aliases gsva,avgParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="avgParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {

              if (verbose) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
              }

              .check_bpparam(BPPARAM)

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       errorOnTooFewRows=TRUE,
                                                       verbose=verbose,
                                                       BPPARAM=BPPARAM)
              filtDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filtMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              maxmem <- .check_maxmem(param, maxmem=maxmem, verbose=verbose)
              ondisk <- .check_ondisk(param, first=NA, last=NA, whdim=2,
                                      recompute_nzcount=FALSE, maxmem=maxmem,
                                      verbose=verbose)

              filtDataMatrix <- .check_sparse_load_input_expr(filtDataMatrix,
                                                              "Z-score", first=NA,
                                                              last=NA, whdim=2,
                                                              ondisk, verbose)

              BPPARAM <- .check_open_parallelism(filtDataMatrix, BPPARAM,
                                                 minparrows=100, minparcols=100,
                                                 verbose)

              ondisk <- .check_es_memory_requirements(filtDataMatrix,
                                                      filtMappedGeneSets,
                                                      ondisk, maxmem)

              if (verbose) {
                  n <- length(filtMappedGeneSets)
                  cli_alert_info("Calculating average scores for {n} gene sets")
              }

              avg_es <- average(X=filtDataMatrix,
                                geneSets=filtMappedGeneSets,
                                method=.get_avgmethod(param),
                                any_na=anyNA(param),
                                na_use=.get_NAuse(param),
                                minSize=get_minSize(param),
                                ondisk=ondisk, verbose=verbose,
                                BPPARAM=BPPARAM, maxmem=maxmem)

              gs <- .geneSetsIndices2Names(
                  indices=filtMappedGeneSets,
                  names=rownames(filtDataMatrix))
              ## dropAssays=TRUE for consistency but doesn't apply here
              rval <- wrapData(get_exprData(param), avg_es, param, "es",
                               first=NA, last=NA, rem=NA, whdim=2,
                               dropAssays=TRUE, gs)
              
              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @title The `avgParam` class
#'
#' @description Objects of class `avgParam` contain the parameters for
#' running the average method.
#'
#' @details The average method takes a number of parameters shared
#' with all methods implemented by package GSVA but does not take any
#' method-specific parameters.
#'
#' @param exprData The expression data set. Must be one of the classes
#' supported by [`GsvaExprData-class`]. For a list of these classes, see its
#' help page using `help(GsvaExprData)`.
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].  For a list of these classes, see its help page using
#' `help(GsvaGeneSets)`.
#' 
#' @param assay Character vector of length 1. The name of the assay to use in
#' case `exprData` is a multi-assay container, otherwise ignored. By default,
#' an assay called 'logcounts' will be used if present, otherwise the first
#' assay is used.
#' 
#' @param annotation An object of class `GeneIdentifierType` from
#' package `GSEABase` describing the gene identifiers used as the row names of
#' the expression data set. See `GeneIdentifierType` for help on available
#' gene identifier types and how to construct them. This
#' information can be used to map gene identifiers occurring in the gene sets.
#' 
#' If the default value `NULL` is provided, an attempt will be made to extract
#' the gene identifier type from the expression data set provided as `exprData`
#' (by calling [`gsvaAnnotation`] on it). If still not successful, the
#' `NullIdentifier()` will be used as the gene identifier type, gene identifier
#' mapping will be disabled and gene identifiers used in expression data set and
#' gene sets can only be matched directly.
#' 
#' @param minSize Numeric vector of length 1. Minimum size of the resulting gene
#' sets after gene identifier mapping. By default, the minimum size is 1.
#' 
#' @param maxSize Numeric vector of length 1. Maximum size of the resulting gene
#' sets after gene identifier mapping. By default, the maximum size is `Inf`.
#' 
#' @param method character vector of length 1. Specific average method employed
#' in the calculation of the average scores. By now, the default and only
#' available value is `method="mean"`, which calculates the average scores as
#' the arithmetic mean of the expression values of the genes in each gene set
#' at each column sample or cell.
#'
#' @param checkNA Character vector of length 1 specifying whether the input
#' expression data should be checked for the presence of missing values (`NA`
#' or `NaN`). This must be one of the strings `"auto"` (default), `"yes"`, or
#' `"no"`. The default value `"auto"` means that the software will perform that
#' check only when the input expression data is provided as a base `matrix`, an
#' `ExpressionSet` or a `SummarizedExperiment` object, while every other type
#' of input expression data container (e.g., `SingleCellExperiment`, etc.) will
#' not be checked. If `checkNA="yes"`, then the input expression data will be
#' checked for missing values irrespective of the object class of the data
#' container, and if `checkNA="no"`, then that check will not be performed.
#'
#' @param use Character vector of length 1 specifying a policy for dealing with
#' missing values (`NA` or `NaN`) in the input expression data argument
#' `exprData`. It only applies when either `checkNA="yes"`, or `checkNA="auto"`
#' (see the `checkNA` parameter. The argument value must be one of the strings
#' `"everything"` (default), `"all.obs"`, or `"na.rm"`. The policy of the
#' default value `"everything"` consists of propagating missing values so that
#' the resulting enrichment score will be `NA`, whenever one or more of its
#' contributing values is missing, giving a warning when that happens. When
#' `use="all.obs"`, the presence of `NA`s in the input expression data will
#' produce an error. Finally, when `use="na.rm"`, missing values in the input
#' expression data will be removed from calculations, giving a warning when that
#' happens, and giving an error if no values are left after removing the missing
#' values.
#'
#' @param ondisk Character vector of length 1 denoting whether an on-disk backend
#' should be used to reduce the memory footprint. The default value
#' `ondisk="auto"` will attempt to load all the data in main memory when the
#' input nonzero values fit in main memory, otherwise it will attempt working
#' with an on-disk data structure that reduces de memory footprint. When
#' `ondisk="yes"` it will attempt to work with an on-disk data structure, while
#' when `ondisk="no"` it will attempt to load all the data in main memory.
#'
#' @param verbose Logical vector of length 1. It gives information about some
#' decisions made by the software during parameter object construction when
#' `verbose=TRUE` (default) and remains silent otherwise.
#'
#' @return A new [`avgParam-class`] object.
#'
#' @seealso [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class]
#'
#' @examples
#' suppressPackageStartupMessages({
#' library(GSEABase)
#' library(GSVA)
#' library(GSVAdata)
#' })
#'
#' data(geneprotExpCostaEtAl2021)
#' data(c2BroadSets)
#' 
#' ## for simplicity, use only a subset of the sample data
#' se <- geneExpCostaEtAl2021[1:1000, ]
#' gsc <- c2BroadSets[1:100]
#' avgp1 <- avgParam(se, gsc)
#' avgp1
#'
#' @importFrom methods new
#' @importFrom utils capture.output
#' @rdname avgParam-class
#' 
#' @export
avgParam <- function(exprData, geneSets,
                     assay=NA_character_, annotation=NULL,
                     minSize=1, maxSize=Inf, method="mean",
                     checkNA=c("auto", "yes", "no"),
                     use=c("everything", "all.obs", "na.rm"),
                     ondisk=c("auto", "yes", "no"),
                     verbose=TRUE) {

    .check_input_expr_gene_sets(exprData, geneSets)

    method <- match.arg(method)
    checkNA <- match.arg(checkNA)
    use <- match.arg(use)
    ondisk <- match.arg(ondisk)

    ## check assay parameter and assay names
    assay <- .check_assayNames(assay, exprData, verbose)

    ## check for presence of valid row/feature names
    exprData <- .check_rowNames(expr=exprData, useDummyNames=TRUE,
                                verbose=verbose)

    xa <- gsvaAnnotation(exprData)
    if (is.null(xa)) {
        if (is.null(annotation)) {
            annotation <- NullIdentifier()
        }
    } else {
        if (is.null(annotation)) {
            annotation <- xa
        } else if (verbose) {
            msg <- sprintf(paste("using argument annotation='%s' and",
                                 "ignoring exprData annotation ('%s')"),
                           capture.output(annotation), capture.output(xa))
            cli_alert_info(msg)
        }
    }

    naparam <- .check_for_na_values(exprData=exprData, assay=assay,
                                    checkNA=checkNA, use=use)

    nzc <- .estimate_nzcount(exprData, assay, verbose)

    new("avgParam", exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize, method=method,
        checkNA=checkNA, didCheckNA=naparam$didCheckNA,
        anyNA=naparam$any_na, use=use, nzcount=nzc, ondisk=ondisk)
}


## ----- validator -----

setValidity("avgParam", function(object) {
    inv <- NULL
    xd <- object@exprData
    dd <- dim(xd)
    an <- gsvaAssayNames(xd)
    oa <- object@assay
    
    if(dd[1] == 0) {
        inv <- c(inv, "@exprData has 0 rows")
    }
    if(dd[2] == 0) {
        inv <- c(inv, "@exprData has 0 columns")
    }
    if(length(object@geneSets) == 0) {
        inv <- c(inv, "@geneSets has length 0")
    }
    if(length(oa) != 1) {
        inv <- c(inv, "@assay must be of length 1")
    }
    if(.isCharLength1(oa) && .isCharNonEmpty(an) && (!(oa %in% an))) {
        inv <- c(inv, "@assay must be one of assayNames(@exprData)")
    }
    if(length(object@annotation) != 1) {
        inv <- c(inv, "@annotation must be of length 1")
    }
    if(!inherits(object@annotation, "GeneIdentifierType")) {
        inv <- c(inv, "@annotation must be a subclass of 'GeneIdentifierType'")
    }
    if(length(object@minSize) != 1) {
        inv <- c(inv, "@minSize must be of length 1")
    }
    if(object@minSize < 1) {
        inv <- c(inv, "@minSize must be at least 1 or greater")
    }
    if(length(object@maxSize) != 1) {
        inv <- c(inv, "@maxSize must be of length 1")
    }
    if(object@maxSize < object@minSize) {
        inv <- c(inv, "@maxSize must be at least @minSize or greater")
    }
    if(!.isCharLength1(object@method)) {
        inv <- c(inv, "@method must be a single character string")
    }
    if(!.isCharLength1(object@checkNA)) {
        inv <- c(inv, "@checkNA must be a single character string")
    }
    if(length(object@didCheckNA) != 1) {
        inv <- c(inv, "@didCheckNA must be of length 1")
    }
    if(is.na(object@didCheckNA)) {
        inv <- c(inv, "@didCheckNA must not be NA")
    }
    if(length(object@anyNA) != 1) {
        inv <- c(inv, "@anyNA must be of length 1")
    }
    if(is.na(object@anyNA)) {
        inv <- c(inv, "@anyNA must not be NA")
    }
    if(!.isCharLength1(object@use)) {
        inv <- c(inv, "@use must be a single character string")
    }
    if(length(object@nzcount) != 1) {
        inv <- c(inv, "@nzcount must be of length 1")
    }   
    if(is.na(object@nzcount)) {
        inv <- c(inv, "@nzcount must not be NA")
    }       
    if(!.isCharLength1(object@ondisk)) {
        inv <- c(inv, "@ondisk must be a single character string")
    }
    return(if(length(inv) == 0) TRUE else inv)
})


#' @param x An object of class [`avgParam-class`].
#'
#' @param recursive Not used with `x` being an object of
#' class [`avgParam-class`].
#'
#' @aliases anyNA,avgParam-method
#' @rdname avgParam-class
setMethod("anyNA", signature=c("avgParam"),
          function(x, recursive=FALSE)
            return(x@anyNA))


## ----- details method -----

#' @importFrom GSEABase details
#' @aliases details,avgParam-method
#' @rdname GsvaMethodParam-class
#' @exportMethod details
setMethod("details",
          signature=signature(object="avgParam"),
          function(object) {
              callNextMethod(object)
              cat("method: ", .get_avgmethod(object), "\n",
                  "checkNA: ", .get_checkNA(object), "\n", sep="")
              if (.get_didCheckNA(object)) {
                  if (anyNA(object)) {
                      cat("missing data: yes\n",
                          "na_use: ", .get_NAuse(object), "\n", sep="")
                  } else
                      cat("missing data: no\n")
              } else
                  cat("missing data: didn't check\n")
          })


## ------ internal functions ------

.get_avgmethod <- function(param) {
    return(param@method)
}

## calculate enrichment scores as average scores for all given genes sets
## through the columns of the input matrix Z
#' @importFrom MatrixGenerics colSums
.compute_average_scores_block <- function(Z, geneSetsIdx, method, any_na,
                                          na_use, minSize, wna_env, verbose) {
    idpb <- NULL
    if (verbose)
        idpb <- cli_progress_bar("Calculating average scores",
                                 total=2*length(geneSetsIdx))

    es <- t(vapply(lapply(geneSetsIdx,
                          function(i) {
                              if (verbose)
                                  cli_progress_update(id=idpb)
                              Z[i, , drop=FALSE]
                          }),
                   function(z) {
                       if (verbose)
                           cli_progress_update(id=idpb)
                       avg <- colMeans(z, na.rm=(any_na && na_use=="na.rm"))
                       if (any_na && na_use=="na.rm") {
                           nnas <- colSums(!is.na(z))
                           avg[nnas < minSize] <- NA
                           if (any(is.na(avg)))
                               assign("w", TRUE, envir=wna_env)
                       }
                       avg
                   }, numeric(ncol(Z))))

    if (verbose)
        cli_progress_done(idpb)

    return(es)
}

## this function computes enrichment scores as average scores for all gene
## sets in geneSetsIdx for a given rank matrix R, taking care that if
## 'ondisk=TRUE' because, e.g., the resulting matrix of average scores does not
## fit in main memory, the scores are written into an on-disk data structure
## (HDF5) instead of being returned in main memory.
#' @importFrom S4Arrays DummyArrayGrid
.compute_average_scores <- function(Z, geneSetsIdx, method, any_na, na_use,
                                    minSize, wna_env, ondisk, verbose) {
    p <- nrow(Z)
    n <- ncol(Z)
    es <- NULL

    if (is(Z, "DelayedMatrix") || ondisk) {
        sink <- HDF5RealizationSink(c(length(geneSetsIdx), ncol(Z)),
                                    as.sparse=FALSE) ## enrichment scores are dense
        grid <- DummyArrayGrid(dim(Z))
        grid_es <- DummyArrayGrid(dim(sink))

        if (length(grid) != length(grid_es) ||
            refdim(grid)[2] != refdim(grid_es)[2] ||
            dim(grid)[2] != dim(grid_es)[2]) {
            msg <- paste("Grid column blocks for ranks should match grid column",
                         "blocks for enrichment scores")
            cli_abort(c("x"=msg))
        }

        ## avp - ArrayViewport for reaching the (possibly sparse) expr. matrix
        ## avp_es - ArrayViewport for writing the enrichment dense scores matrix
        colScores_byBlock <- function(avp, avp_es, sink) {
            block <- read_block(Z, avp)
            block <- .compute_average_scores_block(block, geneSetsIdx, method,
                                                   any_na, na_use, minSize,
                                                   wna_env, verbose=verbose)
            write_block(sink, avp_es, block)
        }

        nblock <- length(grid)
        for (bid in seq_len(nblock))
            sink <- colScores_byBlock(grid[[bid]], grid_es[[bid]], sink)
        close(sink)
        es <- as(sink, "DelayedArray")
    } else
        es <- .compute_average_scores_block(Z, geneSetsIdx, method, any_na,
                                            na_use, minSize, wna_env,
                                            verbose=verbose)

    if (any_na && na_use =="na.rm")
        if (get("w", envir=wna_env)) {
            msg <- sprintf(paste("NA enrichment scores in gene sets with less than",
                                 "%d genes after removing missing values"), minSize)
            cli_alert_warning(msg)
        }

    return(es)
}


#' @importFrom cli cli_alert_info
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom BiocParallel bpnworkers bplapply bpprogressbar
#' @importFrom MatrixGenerics colSums
average <- function(X, geneSets, method="mean",
                    any_na=FALSE, na_use=c("everything", "all.obs", "na.rm"),
                    minSize=1, ondisk=FALSE, verbose=TRUE,
                    BPPARAM=NULL, maxmem=Inf) {
    method <- match.arg(method)
    na_use <- match.arg(na_use)

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    es <- NULL
    if (ncol(X) >= length(geneSets) || is(X, "DelayedMatrix") || ondisk) {
        es <- .processMatrixCols(X, .compute_average_scores, geneSets,
                                 method=method, any_na=any_na,
                                 na_use=na_use, minSize=minSize,
                                 wna_env=wna_env, ondisk=ondisk,
                                 verbose=verbose,
                                 minparrows=100, minparcols=100,
                                 progressmsg="Calculating average scores per gene set",
                                 BPPARAM=BPPARAM, maxmem=maxmem)
    } else {
        if (is.null(BPPARAM) || bpnworkers(BPPARAM) == 1L) {
            env <- NULL
            if (verbose) {
                env <- new.env(parent=globalenv())
                msg <- "Calculating average scores per gene set"
                assign("idpb", cli_progress_bar(msg, total=length(geneSets)),
                       envir=env)
            }
            es <- lapply(geneSets, function(gSetIdx, verbose, idpbe) {
                             if (verbose)
                                 cli_progress_update(id=get("idpb", envir=idpbe))
                             Xgset <- X[gSetIdx, , drop=FALSE]
                             avg <- colMeans(Xgset,
                                             na.rm=(any_na && na_use=="na.rm"))
                             if (any_na && na_use=="na.rm") {
                               nnas <- colSums(!is.na(Xgset))
                               avg[nnas < minSize] <- NA
                               if (any(is.na(avg)))
                                   assign("w", TRUE, envir=wna_env)
                             }
                             avg
                         }, verbose=verbose, idpbe=env)
            if (verbose)
                cli_progress_done(get("idpb", envir=env))
        } else {
            if (verbose)
                bpprogressbar(BPPARAM) <- TRUE ## reporting progress wo/ cli

            es <- bplapply(geneSets, function(gSetIdx) {
                               Xgset <- X[gSetIdx, , drop=FALSE]
                               avg <- colMeans(Xgset,
                                               na.rm=(any_na && na_use=="na.rm"))
                               if (any_na && na_use=="na.rm") {
                                 nnas <- colSums(!is.na(Xgset))
                                 avg[nnas < minSize] <- NA
                                 if (any(is.na(avg)))
                                     assign("w", TRUE, envir=wna_env)
                               }
                               avg
                           }, BPPARAM=BPPARAM)
        }
        es <- do.call(rbind, es)
    }

    if (any_na && na_use =="na.rm")
        if (get("w", envir=wna_env)) {
            msg <- sprintf(paste("NA enrichment scores in gene sets with less than",
                                 "%d genes after removing missing values"), minSize)
            cli_alert_warning(msg)
        }

    if (length(geneSets) == 1)
        es <- matrix(es, nrow=1)

    rownames(es) <- names(geneSets)
    colnames(es) <- colnames(X)

    es
}
