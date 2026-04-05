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
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       verbose=verbose,
                                                       BPPARAM=BPPARAM)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if(verbose)
                  cli_alert_info(sprintf("Calculating PLAGE scores for %d gene sets",
                                         length(filteredMappedGeneSets)))

              plageScores <- plage(X=filteredDataMatrix,
                                   geneSets=filteredMappedGeneSets,
                                   verbose=verbose,
                                   BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), plageScores, gs)

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
                       minSize=1, maxSize=Inf, verbose=TRUE) {
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

    new("plageParam", exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize)
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
    return(if(length(inv) == 0) TRUE else inv)
})


## ------ internal functions ------

#' @importFrom BiocSingular runExactSVD
#' @importFrom cli cli_progress_update
rightsingularsvdvectorgset <- function(gSetIdx, Z, verbose, idpb) {
  if(is(Z, "dgCMatrix")){
    s <- runExactSVD(Z[gSetIdx, ])
  } else {
    s <- svd(Z[gSetIdx, ])
  }
  if (verbose)
      cli_progress_update(id=idpb)
  s$v[, 1]
}

#' @importFrom cli cli_alert_info
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom BiocParallel bpnworkers SerialParam bplapply
plage <- function(X, geneSets, verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {
    Z <- NULL
    if (is(X, "dgCMatrix")){
        if (verbose)
            cli_alert_info("Centering and scaling non-zero values")

        Z <- t(.sparseColumnApplyAndReplace(t(X), FUN=scale))
    } else if (is.matrix(X)) {
        if (verbose)
            cli_alert_info("Centering and scaling values")

        Z <- t(scale(t(X)))
    } else
        stop(sprintf("Matrix class %s cannot be handled yet.", class(X)))

    if (bpnworkers(BPPARAM) > 1) {
        if (verbose)
            cli_progress_bar("Calculating PLAGE scores")
        es <- bplapply(geneSets, rightsingularsvdvectorgset, Z,
                       verbose=FALSE, idpb=character(0), BPPARAM=BPPARAM)
    } else {
        idpb <- NULL
        if (verbose)
            idpb <- cli_progress_bar("Calculating PLAGE scores",
                                     total=length(geneSets))
        es <- lapply(geneSets, rightsingularsvdvectorgset, Z,
                     verbose, idpb)
        if (verbose)
            cli_progress_done(idpb)
    }
        
    es <- do.call(rbind, es)
        
    ## why these extra steps for dense matrices?
    ## svd() removes all dimnames, while BiocSingular::runExactSVD()
    ## retains them.  colnames must be set; however, I don't see a reason
    ## to set rownames nor recreating the matrix if only one gene set.
    ## hmmm...
    if (length(geneSets) == 1)
        es <- matrix(es, nrow=1)
        
    rownames(es) <- names(geneSets)
    colnames(es) <- colnames(X)
    
    es
}
