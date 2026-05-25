##
## methods for the z-score method from Lee et al. (2008)
##

#' @importFrom S4Arrays is_sparse
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @aliases gsva,zscoreParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="zscoreParam"),
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
                                                       verbose=verbose,
                                                       BPPARAM=BPPARAM)
              filtDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filtMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              maxmem <- .check_maxmem(param, maxmem=maxmem, verbose=verbose)
              ondisk <- .check_ondisk(param, maxmem=maxmem, first=NA, last=NA,
                                      whdim=2, verbose=verbose)

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
                  cli_alert_info("Calculating Z-scores for {n} gene sets")
              }

              zscore_es <- zscore(X=filtDataMatrix,
                                  geneSets=filtMappedGeneSets,
                                  ondisk=ondisk, verbose=verbose,
                                  BPPARAM=BPPARAM, maxmem=maxmem)

              gs <- .geneSetsIndices2Names(
                  indices=filtMappedGeneSets,
                  names=rownames(filtDataMatrix))
              ## dropAssays=TRUE for consistency but doesn't apply here
              rval <- wrapData(get_exprData(param), zscore_es, param, "es",
                               first=NA, last=NA, whdim=2, dropAssays=TRUE, gs)
              
              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @title The `zscoreParam` class
#'
#' @description Objects of class `zscoreParam` contain the parameters for
#' running the combined z-scores method.
#'
#' @details The combined z-scores method takes a number of parameters shared
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
#' @return A new [`zscoreParam-class`] object.
#'
#' @seealso [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class]
#'
#' @references Lee, E. et al. Inferring pathway activity toward precise
#' disease classification.
#' *PLoS Comp Biol*, 4(11):e1000217, 2008.
#' \doi{10.1371/journal.pcbi.1000217}
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
#' zp1 <- zscoreParam(se, gsc)
#' zp1
#'
#' @importFrom methods new
#' @importFrom utils capture.output
#' @rdname zscoreParam-class
#' 
#' @export
zscoreParam <- function(exprData, geneSets,
                        assay=NA_character_, annotation=NULL,
                        minSize=1, maxSize=Inf, ondisk=c("auto", "yes", "no"),
                        verbose=TRUE) {

    .check_input_expr_gene_sets(exprData, geneSets)

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

    nzc <- .estimate_nzcount(exprData, assay, verbose)

    new("zscoreParam", exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize, nzcount=nzc, ondisk=ondisk)
}


## ----- validator -----

setValidity("zscoreParam", function(object) {
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


## ------ internal functions ------

#' @importFrom MatrixGenerics rowMeans rowSds
#' @importFrom SparseArray rowMeans rowSds
#' @importFrom DelayedMatrixStats rowSds
.scale_rows <- function(X, verbose) {
    ## scaled <- t(scale(t(X)))
    rmns <- rowMeans(X)
    rsds <- rowSds(X) ## produces tiny differences 10^-16 wrt scale(), but it
                      ## is more performant
    scaled <- (X - rmns) / rsds

    return(scaled)
}

## calculate enrichment scores as combined z-scores for all given genes sets
## through the columns of the input matrix Z
## of the input matrix
#' @importFrom MatrixGenerics colSums
.compute_z_scores_block <- function(Z, geneSetsIdx, verbose) {
    idpb <- NULL
    if (verbose)
        idpb <- cli_progress_bar("Calculating Z-scores",
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
                       colSums(z) / sqrt(nrow(z))
                   }, numeric(ncol(Z))))

    if (verbose)
        cli_progress_done(idpb)

    return(es)
}

## this function computes enrichment scores as combined z-scores for all gene
## sets in geneSetsIdx for a given rank matrix R, taking care that if
## 'ondisk=TRUE' because, e.g., the resulting matrix of ssGSEA scores does not
## fit in main memory, the scores are written into an on-disk data structure
## (HDF5) instead of being returned in main memory.
#' @importFrom S4Arrays DummyArrayGrid
.compute_z_scores <- function(Z, geneSetsIdx, ondisk, verbose) {
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

        ## avp - ArrayViewport for reaching the (possibly sparse) rank matrix
        ## avp_es - ArrayViewport for writing the enrichment dense scores matrix
        colScores_byBlock <- function(avp, avp_es, sink) {
            block <- read_block(Z, avp)
            block <- .compute_z_scores_block(block, geneSetsIdx, verbose)
            write_block(sink, avp_es, block)
        }

        nblock <- length(grid)
        for (bid in seq_len(nblock))
            sink <- colScores_byBlock(grid[[bid]], grid_es[[bid]], sink)
        close(sink)
        es <- as(sink, "DelayedArray")
    } else
        es <- .compute_z_scores_block(Z, geneSetsIdx, verbose)

    return(es)
}


#' @importFrom cli cli_alert_info
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom BiocParallel bpnworkers bplapply bpprogressbar
#' @importFrom MatrixGenerics colSums
zscore <- function(X, geneSets, ondisk=FALSE, verbose=TRUE,
                   BPPARAM=NULL, maxmem=Inf) {

    Z <- .processMatrixRows(X, .scale_rows, verbose=verbose,
                            minparrows=100, minparcols=100,
                            progressmsg="Centering and scaling rows",
                            BPPARAM=BPPARAM, maxmem=maxmem)

    es <- NULL
    if (ncol(Z) >= length(geneSets) || is(Z, "DelayedMatrix") || ondisk) {
        es <- .processMatrixCols(Z, .compute_z_scores, geneSets,
                                 ondisk=ondisk, verbose=verbose,
                                 minparrows=100, minparcols=100,
                                 progressmsg="Calculating Z-scores per gene set",
                                 BPPARAM=BPPARAM, maxmem=maxmem)
    } else {
        if (is.null(BPPARAM) || bpnworkers(BPPARAM) == 1L) {
            env <- NULL
            if (verbose) {
                env <- new.env(parent=globalenv())
                msg <- "Calculating Z-scores per gene set"
                assign("idpb", cli_progress_bar(msg, total=length(geneSets)),
                       envir=env)
            }
            es <- lapply(geneSets, function(gSetIdx, verbose, idpbe) {
                             if (verbose)
                                 cli_progress_update(id=get("idpb", envir=idpbe))
                             colSums(Z[gSetIdx, , drop=FALSE]) / sqrt(length(gSetIdx))
                         }, verbose=verbose, idpbe=env)
            if (verbose)
                cli_progress_done(get("idpb", envir=env))
        } else {
            if (verbose)
                bpprogressbar(BPPARAM) <- TRUE ## reporting progress wo/ cli

            es <- bplapply(geneSets, function(gSetIdx) {
                               colSums(Z[gSetIdx, , drop=FALSE]) / sqrt(length(gSetIdx))
                           }, BPPARAM=BPPARAM)
        }
        es <- do.call(rbind, es)
    }

    es
}
