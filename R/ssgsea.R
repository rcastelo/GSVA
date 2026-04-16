##
## methods for the ssGSEA method from Barbie et al. (2009)
##

#' @importFrom S4Arrays is_sparse
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @aliases gsva,ssgseaParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="ssgseaParam"),
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
                                                       removeConstant=FALSE,
                                                       removeNzConstant=FALSE,
                                                       verbose=verbose,
                                                       BPPARAM=BPPARAM)
              filtDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filtMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              maxmem <- .check_maxmem(param, maxmem, verbose)
              ondisk <- .check_ondisk(param, maxmem, verbose)

              filtDataMatrix <- .check_sparse_load_input_expr(filtDataMatrix,
                                                              "ssGSEA",
                                                              ondisk, verbose)

              BPPARAM <- .check_open_parallelism(filtDataMatrix, BPPARAM,
                                                 minparrows=100, minparcols=100,
                                                 verbose)

              ondisk <- .check_es_memory_requirements(filtDataMatrix,
                                                      filtMappedGeneSets,
                                                      ondisk, maxmem)

              if (verbose) {
                  n <- length(filtMappedGeneSets)
                  cli_alert_info("Calculating ssGSEA scores for {n} gene sets")
              }

              ssgsea_es <- ssgsea(X=filtDataMatrix,
                                  geneSetsIdx=filtMappedGeneSets,
                                  alpha=.get_alpha(param), 
                                  normalization=.get_normalize(param),
                                  any_na=anyNA(param),
                                  na_use=.get_NAuse(param),
                                  minSize=get_minSize(param),
                                  ondisk=ondisk, verbose=verbose,
                                  BPPARAM=BPPARAM, maxmem=maxmem)

              gs <- .geneSetsIndices2Names(
                  indices=filtMappedGeneSets,
                  names=rownames(filtDataMatrix))
              rval <- wrapData(get_exprData(param), ssgsea_es, gs)
              
              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @title The `ssgseaParam` class
#'
#' @description Objects of class `ssgseaParam` contain the parameters for
#' running the `ssGSEA` method.
#'
#' @details In addition to a number of parameters shared with all methods
#' implemented by package GSVA, `ssGSEA` takes two method-specific parameters as
#' well as two more parameters for implementing a missing value policy. All of
#' these parameters are described in detail below.
#'
#' @param exprData The expression data set.  Must be one of the classes
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
#' @param alpha Numeric vector of length 1. The exponent defining the
#' weight of the tail in the random walk performed by the `ssGSEA` (Barbie et
#' al., 2009) method.  The default value is 0.25 as described in the paper.
#' 
#' @param normalize Logical vector of length 1; if `TRUE` runs the `ssGSEA`
#' method from Barbie et al. (2009) normalizing the scores by the absolute
#' difference between the minimum and the maximum, as described in their paper.
#' Otherwise this final normalization step is skipped.
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
#' @return A new [`ssgseaParam-class`] object.
#'
#' @seealso [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class],
#' [`matrix`],
#' \code{\link[Biobase]{ExpressionSet}},
### we are using the plain Rd above because
###  #' [`ExpressionSet`][Biobase::ExpressionSet-class],
### results in the following R CMD check NOTE:
### Non-topic package-anchored link(s) in Rd file 'ssgseaParam-class.Rd':
###  ‘[Biobase:class.ExpressionSet]{ExpressionSet}’
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment-class],
#' [`SingleCellExperiment`][SingleCellExperiment::SingleCellExperiment-class]
#'
#' @references Barbie, D.A. et al. Systematic RNA interference reveals that
#' oncogenic KRAS-driven cancers require TBK1.
#' *Nature*, 462(5):108-112, 2009.
#' \doi{10.1038/nature08460}
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
#' sp1 <- ssgseaParam(se, gsc)
#' sp1
#'
#' @importFrom methods new
#' @importFrom utils capture.output
#' @rdname ssgseaParam-class
#' 
#' @export
ssgseaParam <- function(exprData, geneSets,
                        assay=NA_character_, annotation=NULL,
                        minSize=1,maxSize=Inf,
                        alpha=0.25, normalize=TRUE,
                        checkNA=c("auto", "yes", "no"),
                        use=c("everything", "all.obs", "na.rm"),
                        ondisk=c("auto", "yes", "no"),
                        verbose=TRUE) {

    .check_input_expr_gene_sets(exprData, geneSets)

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
    
    new("ssgseaParam",
        exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize,
        alpha=alpha, normalize=normalize,
        checkNA=checkNA, didCheckNA=naparam$didCheckNA,
        anyNA=naparam$any_na, use=use, nzcount=nzc, ondisk=ondisk)
}


## ----- validator -----

setValidity("ssgseaParam", function(object) {
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
    if(length(object@alpha) != 1) {
        inv <- c(inv, "@alpha must be of length 1")
    }
    if(is.na(object@alpha)) {
        inv <- c(inv, "@alpha must not be NA")
    }
    if(length(object@normalize) != 1) {
        inv <- c(inv, "@normalize must be of length 1")
    }
    if(is.na(object@normalize)) {
        inv <- c(inv, "@normalize must not be NA")
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


#' @param x An object of class [`ssgseaParam-class`].
#'
#' @param recursive Not used with `x` being an object of
#' class [`ssgseaParam-class`].
#'
#' @aliases anyNA,ssgseaParam-method
#' @rdname ssgseaParam-class
setMethod("anyNA", signature=c("ssgseaParam"),
          function(x, recursive=FALSE)
            return(x@anyNA))


## ----- show -----

setMethod("show",
          signature=signature(object="ssgseaParam"),
          function(object) {
              callNextMethod(object)
              cat("alpha: ", .get_alpha(object), "\n", 
                  "normalize: ", .get_normalize(object), "\n",
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

.get_alpha <- function(object) {
  stopifnot(inherits(object, "ssgseaParam"))
  return(object@alpha)
}

.get_normalize <- function(object) {
  stopifnot(inherits(object, "ssgseaParam"))
  return(object@normalize)
}

.rndWalk <- function(gSetIdx, geneRanking, j, R, alpha) {
  indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  stepCDFinGeneSet <- cumsum((abs(R[geneRanking, j])^alpha * 
                      indicatorFunInsideGeneSet)) /
                      sum((abs(R[geneRanking, j])^alpha *
                      indicatorFunInsideGeneSet))
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet) /
                       sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

  sum(walkStat) 
}

## optimized version of the function .rndWalk by Alexey Sergushichev
## https://github.com/rcastelo/GSVA/pull/15
## based on his paper DOI:10.1101/060012
## with further optimizations with Bob Policastro discussed in
## https://github.com/rcastelo/GSVA/issues/71
.fastRndWalk <- function(gSetIdx, geneRanking, j, Ra) {
    n <- length(geneRanking)
    k <- length(gSetIdx)
    
    stepCDFinGeneSet <- 
        sum(Ra[geneRanking[gSetIdx], j] * (n - gSetIdx + 1)) /
        sum((Ra[geneRanking[gSetIdx], j]))    
    
    stepCDFoutGeneSet <- (n * (n + 1) / 2 - sum(n - gSetIdx + 1)) / (n - k)
    
    walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

    walkStat
}

.fastRndWalkNArm <- function(gSetIdx, geneRanking, j, Ra,
                             anyna=FALSE, na_use="everything", minSize=1,
                             wna_env=new.env()) {
    n <- length(geneRanking)
    if (anyna && na_use == "na.rm")
        gSetIdx <- na.omit(gSetIdx)
    k <- length(gSetIdx)
    
    walkStat <- NA_real_
    if (k >= minSize) {
      stepCDFinGeneSet <- 
          sum(Ra[geneRanking[gSetIdx], j] * (n - gSetIdx + 1)) /
          sum((Ra[geneRanking[gSetIdx], j]))    
    
      stepCDFoutGeneSet <- (n * (n + 1) / 2 - sum(n - gSetIdx + 1)) / (n - k)
    
      walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
    } else if (!get("w", envir=wna_env)) ## warn only once. it can only happen
      assign("w", TRUE, envir=wna_env)   ## with anyna=TRUE and na_use="na.rm"

    walkStat
}

#' @importFrom IRanges IntegerList match
#' @importFrom BiocParallel bpnworkers
#' @importFrom BiocGenerics "type<-"
#' @importFrom MatrixGenerics colRanks
#' @importFrom cli cli_alert_info cli_alert_warning cli_abort
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
ssgsea <- function(X, geneSetsIdx, alpha=0.25,
                   normalization=TRUE,
                   any_na=FALSE,
                   na_use=c("everything", "all.obs", "na.rm"),
                   minSize=1, ondisk=FALSE, verbose=TRUE,
                   BPPARAM=NULL,
                   maxmem=Inf) {
    na_use <- match.arg(na_use)

    R <- .processMatrixCols(X, FUN=compute.col.ranks, ties.method="average",
                            drop.sparsity=TRUE, verbose=verbose, minparrows=100,
                            minparcols=100, progressmsg="Calculating ranks",
                            BPPARAM=BPPARAM, maxmem=Inf)
    if (!is(R, "dgCMatrix")) ## dgCMatrix cannot be coerced to integer
      type(R) <- "integer"

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    ssgsea_es <- .processMatrixCols(R, FUN=.compute_ssgsea_scores,
                                    geneSetsIdx=geneSetsIdx, alpha=alpha,
                                    normalization=normalization, any_na=any_na,
                                    na_use=na_use, minSize=minSize,
                                    wna_env=wna_env, ondisk=ondisk,
                                    verbose=verbose,
                                    minparrows=100, minparcols=100,
                                    progressmsg="Calculating scores",
                                    BPPARAM=BPPARAM, maxmem=maxmem)

    if (any_na && na_use =="na.rm")
        if (get("w", envir=wna_env)) {
            msg <- sprintf(paste("NA enrichment scores in gene sets with less than",
                                 "%d genes after removing missing values"), minSize)
        cli_alert_warning(msg)
        }
  
    if (normalization) {
        if (verbose)
            cli_alert_info("Normalizing ssGSEA scores")
        ## normalize enrichment scores by using the entire data set, as indicated
        ## by Barbie et al., 2009, online methods, pg. 2

        rng <- c(attr(ssgsea_es, "min"), attr(ssgsea_es, "max"))

        if (length(rng) == 0 || any(is.na(rng) | !is.finite(rng))) {
            msg <- paste("Cannot calculate normalizing factor for the enrichment",
                         "scores, most likely due to NA values in the input data.")
            cli_abort(c("x"=msg))
        }
        attr(ssgsea_es, "min") <- attr(ssgsea_es, "max") <- NULL
        n <- ncol(ssgsea_es)
        ssgsea_es <- ssgsea_es[, seq_len(n), drop=FALSE] / (rng[2] - rng[1])
    }
  
    if (length(geneSetsIdx) == 1)
        ssgsea_es <- matrix(ssgsea_es, nrow=1)
  
    rownames(ssgsea_es) <- names(geneSetsIdx)
    colnames(ssgsea_es) <- colnames(X)
  
    ssgsea_es
}

## this function computes the ssGSEA scores for all gene sets in geneSetsIdx for
## a given rank matrix R, taking care that if 'ondisk=TRUE' because, e.g., the
## resulting matrix of ssGSEA scores does not fit in main memory, the scores are
## written into an on-disk data structure (HDF5) instead of being returned in
## main memory.
#' @importFrom IRanges IntegerList
#' @importFrom S4Arrays DummyArrayGrid
.compute_ssgsea_scores <- function(R, geneSetsIdx, alpha, normalization,
                                   any_na, na_use, minSize, wna_env, ondisk,
                                   verbose) {
    p <- nrow(R)
    n <- ncol(R)
    es <- NULL

    geneSetsIdx <- IntegerList(geneSetsIdx)

    if (is(R, "DelayedMatrix") || ondisk) {
        sink <- HDF5RealizationSink(c(length(geneSetsIdx), ncol(R)),
                                    as.sparse=FALSE) ## enrichment scores are dense
        grid <- DummyArrayGrid(dim(R))
        grid_es <- DummyArrayGrid(dim(sink))

        if (length(grid) != length(grid_es) ||
            refdim(grid)[2] != refdim(grid_es)[2] ||
            dim(grid)[2] != dim(grid_es)[2]) {
            msg <- paste("Grid column blocks for ranks should match grid column",
                         "blocks for enrichment scores")
            cli_abort(c("x"=msg))
        }

        ## avp - ArrayViewport for reaching the (possibly sparse) rank matrix
        ## avp_es - ArrayViewport for writing the enrichment dense scores matrix
        colScores_byBlock <- function(avp, avp_es, sink) {
            block <- read_block(R, avp)
            block <- .compute_ssgsea_scores_block(block, geneSetsIdx, alpha,
                                                  normalization, any_na,
                                                  na_use, minSize, wna_env,
                                                  verbose=verbose)
            if (normalization) {
                if (any_na) {
                    assign("mines", min(c(mines, attr(block, "min")), na.rm=TRUE),
                           envir=parent.frame(1))
                    assign("maxes", max(c(maxes, attr(block, "max")), na.rm=TRUE),
                           envir=parent.frame(1))
                } else {
                    assign("mines", min(c(mines, attr(block, "min")), na.rm=FALSE),
                           envir=parent.frame(1))
                    assign("maxes", max(c(maxes, attr(block, "max")), na.rm=FALSE),
                           envir=parent.frame(1))
                }
                attr(block, "min") <- attr(block, "max") <- NULL
            }
            write_block(sink, avp_es, block)
        }

        mines <- Inf
        maxes <- -Inf
        nblock <- length(grid)
        for (bid in seq_len(nblock))
            sink <- colScores_byBlock(grid[[bid]], grid_es[[bid]], sink)
        close(sink)
        es <- as(sink, "DelayedArray")
        if (normalization) {
            attr(es, "min") <- mines
            attr(es, "max") <- maxes
        }

    } else {
        es <- .compute_ssgsea_scores_block(R, geneSetsIdx, alpha,
                                           normalization, any_na, na_use,
                                           minSize, wna_env, verbose=verbose)
    }

    if (any_na && na_use =="na.rm")
        if (get("w", envir=wna_env)) {
            msg <- sprintf(paste("NA enrichment scores in gene sets with less than",
                                 "%d genes after removing missing values"), minSize)
        cli_alert_warning(msg)
        }

    return(es)
}

## here geneSetsIdx should be an 'IntegerList' object
#' @importFrom IRanges match
.compute_ssgsea_scores_block <- function(R, geneSetsIdx, alpha, normalization,
                                         any_na, na_use, minSize, wna_env, verbose) {
    stopifnot(is(geneSetsIdx, "IntegerList")) ## QC
    n <- ncol(R)
    idpb <- NULL
    if (verbose)
      idpb <- cli_progress_bar("Calculating ssGSEA scores", total=n)

    Ra <- R
    if (alpha != 1)
        Ra <- R^alpha

    mines <- Inf
    maxes <- -Inf

    es <- lapply(as.list(seq_len(n)), function(j) {
        if (any_na && na_use == "na.rm") {
            geneRanking <- order(R[, j], decreasing=TRUE, na.last=NA)
            geneSetsRankIdx <- match(geneSetsIdx, geneRanking)
            es_sample <- vapply(X=geneSetsRankIdx, FUN=.fastRndWalkNArm,
                                FUN.VALUE=numeric(1),
                                geneRanking, j, Ra, any_na, na_use,
                                minSize, wna_env, USE.NAMES=FALSE)
        } else {
            geneRanking <- order(R[, j], decreasing=TRUE)
            geneSetsRankIdx <- match(geneSetsIdx, geneRanking)
            es_sample <- vapply(X=geneSetsRankIdx, FUN=.fastRndWalk,
                                FUN.VALUE=numeric(1),
                                geneRanking, j, Ra)
        }
        if (verbose)
            cli_progress_update(id=idpb)
        if (normalization) {
            if (any_na) {
                assign("mines", min(c(mines, es_sample), na.rm=TRUE),
                       envir=parent.frame(2))
                assign("maxes", max(c(maxes, es_sample), na.rm=TRUE),
                       envir=parent.frame(2))
            } else {
                assign("mines", min(c(mines, es_sample), na.rm=FALSE),
                       envir=parent.frame(2))
                assign("maxes", max(c(maxes, es_sample), na.rm=FALSE),
                       envir=parent.frame(2))
            }
        }
        es_sample
    })
    if (verbose)
        cli_progress_done(idpb)
    es <- do.call("cbind", es)
    if (normalization) {
        attr(es, "min") <- mines
        attr(es, "max") <- maxes
    }
  
    return(es)
}
