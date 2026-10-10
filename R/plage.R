##
## methods for the PLAGE method from Tomfohr et al. (2005)
##

#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @aliases gsva,plageParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="plageParam"),
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
              ## the input data is loaded in main memory as a dense matrix,
              ## and the memory required includes the enrichment scores
              mf <- .step_mem_factors("plage", sparse=FALSE, int=FALSE,
                                      ngs=length(filtMappedGeneSets))
              ondisk <- .check_ondisk(param, first=NA, last=NA, whdim=2,
                                      recompute_nzcount=FALSE, maxmem=maxmem,
                                      verbose=verbose, mf=mf, dense=TRUE)
              .check_step_mem(filtDataMatrix, 2L, NA, NA, ondisk, mf, BPPARAM,
                              maxmem, dense=TRUE)

              filtDataMatrix <- .check_sparse_load_input_expr(filtDataMatrix,
                                                              "PLAGE", first=NA,
                                                              last=NA, whdim=2,
                                                              ondisk, verbose)

              BPPARAM <- .check_open_parallelism(filtDataMatrix, BPPARAM,
                                                 minparrows=100, minparcols=100,
                                                 verbose)


              if (verbose) {
                  n <- length(filtMappedGeneSets)
                  cli_alert_info("Calculating PLAGE scores for {n} gene sets")
              }

              plage_es <- plage(X=filtDataMatrix,
                                geneSets=filtMappedGeneSets,
                                ondisk=ondisk, verbose=verbose,
                                BPPARAM=BPPARAM, maxmem=maxmem)

              gs <- .geneSetsIndices2Names(
                  indices=filtMappedGeneSets,
                  names=rownames(filtDataMatrix))
              ## dropAssays=TRUE for consistency but doesn't apply here
              rval <- wrapData(get_exprData(param), plage_es, param, "es",
                               first=NA, last=NA, rem=NA, whdim=2,
                               dropAssays=TRUE, gs)

              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @title The `plageParam` class
#'
#' @description Objects of class `plageParam` contain the parameters for running
#' the `PLAGE` method.
#'
#' @details `PLAGE` takes a number of parameters shared with all methods
#' implemented by package GSVA but does not take any method-specific parameters.
#' These parameters are described in detail below.
#'
#' The `PLAGE` score of a gene set is the first right singular vector of the
#' standardized expression values of its genes, whose sign is arbitrary
#' (Bro et al., 2008). GSVA orients it so that the scores correlate positively
#' with the mean standardized expression of the genes in the gene set, in the
#' same way as module eigengenes are aligned with the average expression of
#' their module in the WGCNA package (Langfelder and Horvath, 2008). As a
#' result, higher `PLAGE` scores tend to correspond to higher expression of
#' the gene set, and the sign of the scores is the same across runs.
#'
#' @param exprData The expression data set.  Must be one of the classes
#' supported by [`GsvaExprData-class`].  For a list of these classes, see its
#' help page using `help(GsvaExprData)`.
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].  For a list of these classes, see its help page using
#' `help(GsvaGeneSets)`.
#' 
#' @param assay Character vector of length 1.  The name of the assay to use in
#' case `exprData` is a multi-assay container, otherwise ignored.  By default,
#' an assay called 'logcounts' will be used if present, otherwise the first
#' assay is used.
#' 
#' @param annotation An object of class `GeneIdentifierType` from
#' package `GSEABase` describing the gene identifiers used as the row names of
#' the expression data set.  See `GeneIdentifierType` for help on available
#' gene identifier types and how to construct them.  This
#' information can be used to map gene identifiers occurring in the gene sets.
#' 
#' If the default value `NULL` is provided, an attempt will be made to extract
#' the gene identifier type from the expression data set provided as `exprData`
#' (by calling [`gsvaAnnotation`] on it).  If still not successful, the
#' `NullIdentifier()` will be used as the gene identifier type, gene identifier
#' mapping will be disabled and gene identifiers used in expression data set and
#' gene sets can only be matched directly.
#' 
#' @param minSize Numeric vector of length 1.  Minimum size of the resulting gene
#' sets after gene identifier mapping. By default, the minimum size is 1.
#' 
#' @param maxSize Numeric vector of length 1.  Maximum size of the resulting gene
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
#' @return A new [`plageParam-class`] object.
#'
#' @seealso [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class]
#'
#' @references Tomfohr, J. et al. Pathway level analysis of gene expression
#' using singular value decomposition.
#' *BMC Bioinformatics*, 6:225, 2005.
#' \doi{10.1186/1471-2105-6-225}
#'
#' @references Bro, R., Acar, E. and Kolda, T.G. Resolving the sign ambiguity
#' in the singular value decomposition.
#' *Journal of Chemometrics*, 22:135-140, 2008.
#' \doi{10.1002/cem.1122}
#'
#' @references Langfelder, P. and Horvath, S. WGCNA: an R package for weighted
#' correlation network analysis.
#' *BMC Bioinformatics*, 9:559, 2008.
#' \doi{10.1186/1471-2105-9-559}
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
#' pp1 <- plageParam(se, gsc)
#' pp1
#'
#' @importFrom methods new
#' @importFrom utils capture.output
#' @rdname plageParam-class
#' 
#' @export
plageParam <- function(exprData, geneSets,
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
    if(is.null(xa)) {
        if(is.null(annotation)) {
            annotation <- NullIdentifier()
        }
    } else {
        if(is.null(annotation)) {
            annotation <- xa
        } else if (verbose) {
            msg <- sprintf(paste0("using argument annotation='%s' and ",
                                  "ignoring exprData annotation ('%s')"),
                           capture.output(annotation), capture.output(xa))
            cli_alert_info(msg)
        }
    }

    nzc <- .estimate_nzcount(exprData, assay, verbose)

    new("plageParam", exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize, nzcount=nzc, ondisk=ondisk)
}


## ----- validator -----

setValidity("plageParam", function(object) {
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

## ----- anyNA method -----

#' @param x An object of class [`plageParam-class`]. Since the PLAGE method does
#' not handle missing (`NA`) values, `anyNA()` always returns `FALSE` for
#' this class, without checking the input expression data.
#'
#' @param recursive Not used with `x` being an object of
#' class [`plageParam-class`].
#'
#' @aliases anyNA,plageParam-method
#' @rdname plageParam-class
setMethod("anyNA", signature=c("plageParam"),
          function(x, recursive=FALSE)
            return(FALSE))

## ------ internal functions ------

## the PLAGE score of a gene set in each column is the first right singular
## vector v of the scaled rows Z_gs of the gene set, Z_gs = U D V^T, which is
## calculated by blocks of columns, without taking the whole rows of the gene
## set: Z_gs Z_gs^T = U D^2 U^T, a small matrix of as many rows and columns as
## genes in the gene set, is the sum of the products of the blocks of columns
## of Z_gs, and v = Z_gs^T u / d, for the leading eigenvector u of that matrix
## and the square root d of its leading eigenvalue, is given by the blocks of
## columns of Z_gs. the sign of a singular vector is arbitrary, and v is
## oriented to correlate positively with the mean of the scaled rows of the
## gene set, (1 / k) 1^T Z_gs v = (d / k) sum(u), i.e., to have sum(u) >= 0,
## as the alignment of module eigengenes with their average expression in the
## WGCNA package

## products Z_gs Z_gs^T of the rows of each gene set, with the indices given by
## 'geneSetsIdx', in a block of columns 'Z' of the scaled rows
#' @importFrom S4Arrays DummyArrayGrid read_block
.plage_gram_block <- function(Z, geneSetsIdx, verbose=FALSE) {
    if (is(Z, "DelayedArray"))
        Z <- read_block(Z, DummyArrayGrid(dim(Z))[[1L]])
    Z <- as.matrix(Z)

    lapply(geneSetsIdx, function(i) tcrossprod(Z[i, , drop=FALSE]))
}

## PLAGE scores of each gene set, with the indices given by 'geneSetsIdx', in
## a block of columns 'Z' of the scaled rows, from the weights 'rowweights' of
## the rows of each gene set, u / d; the name of this argument, passed through
## '...' to .ondisk_blocks(), must not match partially its argument 'whdim'
#' @importFrom S4Arrays DummyArrayGrid read_block
.plage_scores_block <- function(Z, geneSetsIdx, rowweights, verbose=FALSE) {
    if (is(Z, "DelayedArray"))
        Z <- read_block(Z, DummyArrayGrid(dim(Z))[[1L]])
    Z <- as.matrix(Z)

    es <- vapply(seq_along(geneSetsIdx), function(g)
                     drop(crossprod(Z[geneSetsIdx[[g]], , drop=FALSE],
                                    rowweights[[g]])),
                 numeric(ncol(Z)))
    es <- t(matrix(es, ncol=length(geneSetsIdx)))
    dimnames(es) <- list(names(geneSetsIdx), colnames(Z))

    es
}

## groups of gene sets whose products Z_gs Z_gs^T, of 'bytes' bytes each, take
## together at most 'maxbytes' bytes, keeping each gene set in a group
.plage_groups <- function(bytes, maxbytes) {
    grp <- integer(length(bytes))
    g <- 1L
    acc <- 0
    for (i in seq_along(bytes)) {
        if (acc > 0 && acc + bytes[i] > maxbytes) {
            g <- g + 1L
            acc <- 0
        }
        grp[i] <- g
        acc <- acc + bytes[i]
    }

    grp
}

#' @importFrom cli cli_alert_info
#' @importFrom BiocParallel bpnworkers
plage <- function(X, geneSets, ondisk=FALSE, verbose=TRUE,
                  BPPARAM=NULL, maxmem=Inf) {

    ## the scaled rows, as large as the input, are delayed operations on it
    ## when it is stored on disk, and are stored on disk when it is in main
    ## memory but they do not fit in it; each block of sparse input gives
    ## dense scaled rows ('dense=TRUE')
    Z <- .processMatrixRows(X, .scale_rows, verbose=verbose,
                            minparrows=100, minparcols=100,
                            progressmsg="Centering and scaling rows",
                            BPPARAM=BPPARAM, maxmem=maxmem, dense=TRUE,
                            sinkout=(ondisk && !is(X, "DelayedMatrix")))

    ## the input in main memory remains allocated while the scaled rows are
    ## processed
    heldmem <- .held_mem() + .inmem_size(X)
    ## when the input is stored on disk, the scaled rows are delayed operations
    ## on it, and realizing each of their blocks of columns takes about three
    ## times the size of its dense form, as measured on single-cell data
    delayedextra <- if (is(X, "DelayedMatrix")) 2 else 0
    nworkers <- if (is.null(BPPARAM)) 1L else bpnworkers(BPPARAM)

    ## the products Z_gs Z_gs^T of the gene sets are accumulated over the
    ## blocks of columns of the scaled rows, by each worker and in the main
    ## process; with many large gene sets, these products would not fit in
    ## main memory, and the gene sets are processed in groups whose products
    ## take at most a quarter of the memory available, each group in a pass
    ## through the scaled rows
    grambytes <- 8 * as.numeric(lengths(geneSets))^2
    maxgram <- .mem_fraction_R * maxmem / 4 / (nworkers + 1)
    grp <- .plage_groups(grambytes, maxgram)
    mf <- .step_mem_factors("plagegram")
    w <- vector("list", length(geneSets))
    for (g in unique(grp)) {
        idx <- which(grp == g)
        grams <- .processMatrixCols(Z, .plage_gram_block,
                                    geneSetsIdx=geneSets[idx],
                                    verbose=verbose, minparrows=100,
                                    minparcols=100,
                                    progressmsg="Calculating PLAGE singular vectors",
                                    BPPARAM=BPPARAM, maxmem=maxmem,
                                    workfactor=mf$workfactor + delayedextra,
                                    outfactor=mf$outfactor,
                                    outextra=mf$outextra,
                                    heldmem=heldmem +
                                            (nworkers + 1) * sum(grambytes[idx]),
                                    combine=function(a, b) Map(`+`, a, b))
        w[idx] <- lapply(grams, function(gram) {
                             e <- eigen(gram, symmetric=TRUE)
                             u <- e$vectors[, 1L]
                             if (sum(u) < 0) ## see above for the sign of v
                                 u <- -u
                             u / sqrt(e$values[1L])
                         })
    }

    mf <- .step_mem_factors("plagescores", ngs=length(geneSets))
    es <- .processMatrixCols(Z, .plage_scores_block, geneSetsIdx=geneSets,
                             rowweights=w, verbose=verbose, minparrows=100,
                             minparcols=100,
                             progressmsg="Calculating PLAGE scores",
                             BPPARAM=BPPARAM, maxmem=maxmem,
                             workfactor=mf$workfactor + delayedextra,
                             outfactor=mf$outfactor, outextra=mf$outextra,
                             heldmem=heldmem,
                             sinkout=(ondisk || is(Z, "DelayedMatrix")))

    rownames(es) <- names(geneSets)
    colnames(es) <- colnames(X)

    es
}
