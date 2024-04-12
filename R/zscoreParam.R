
#' @title The `zscoreParam` class
#'
#' @description Objects of class `zscoreParam` contain the parameters for running
#' the combined z-scores method.
#'
#' @details The combined z-scores method does not take any method-specific
#' parameters in addition to an expression data set and a collection of gene
#' sets.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`].  Type `help(GsvaExprData)` to consult
#' the available classes.
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].
#' 
#' @param assay The name of the assay to use in case `exprData` is a multi-assay
#' container, otherwise ignored.  By default, the first assay is used.
#' 
#' @param annotation The name of a Bioconductor annotation package for the gene
#' identifiers occurring in the row names of the expression data matrix.  This
#' can be used to map gene identifiers occurring in the gene sets if those are
#' provided in a [`GeneSetCollection`].  By default gene identifiers used in
#' expression data matrix and gene sets are matched directly.
#' 
#' @param minSize Minimum size of the resulting gene sets after gene identifier
#' mapping. By default, the minimum size is 1.
#' 
#' @param maxSize Maximum size of the resulting gene sets after gene identifier
#' mapping. By default, the maximum size is `Inf`.
#' 
#' @return A new [`zscoreParam-class`] object.
#'
#' @references Lee, E. et al. Inferring pathway activity toward precise
#' disease classification.
#' *PLoS Comp Biol*, 4(11):e1000217, 2008.
#' [DOI](https://doi.org/10.1371/journal.pcbi.1000217)
#'
#' @examples
#' library(GSVA)
#' library(GSVAdata)
#'
#' data(leukemia)
#' data(c2BroadSets)
#'
#' ## for simplicity, use only a subset of the sample data
#' ses <- leukemia_eset[1:1000, ]
#' gsc <- c2BroadSets[1:100]
#' zp1 <- zscoreParam(ses, gsc)
#' zp1
#'
#'
#' @importFrom methods new
#' @rdname zscoreParam-class
#' 
#' @export
zscoreParam <- function(exprData, geneSets,
                        assay=NA_character_, annotation=NA_character_,
                        minSize=1,maxSize=Inf) {
    an <- gsvaAssayNames(exprData)
    if((!is.na(assay)) && (!.isCharNonEmpty(an)))
        warning("argument assay='", assay,
                "' ignored since exprData has no assayNames()")
    if(is.na(assay) && .isCharNonEmpty(an))
        assay <- na.omit(an)[1]
    
    new("zscoreParam", exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize)
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
        inv <- c(inv, "@assay should be of length 1")
    }
    if(.isCharLength1(oa) && .isCharNonEmpty(an) && (!(oa %in% an))) {
        inv <- c(inv, "@assay should be one of assayNames(@exprData)")
    }
    if(length(object@annotation) != 1) {
        inv <- c(inv, "@annotation should be of length 1")
    }
    if(length(object@minSize) != 1) {
        inv <- c(inv, "@minSize should be of length 1")
    }
    if(object@minSize < 1) {
        inv <- c(inv, "@minSize should be at least 1 or greater")
    }
    if(length(object@maxSize) != 1) {
        inv <- c(inv, "@maxSize should be of length 1")
    }
    if(object@maxSize < object@minSize) {
        inv <- c(inv, "@maxSize should be at least @minSize or greater")
    }
    return(if(length(inv) == 0) TRUE else inv)
})

