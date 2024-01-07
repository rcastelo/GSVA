
#' @title The `ssgseaParam` class
#'
#' @description Objects of class `ssgseaParam` contain the parameters for running
#' the `ssGSEA` method.
#'
#' @details In addition to an expression data set and a collection of
#' gene sets, `ssGSEA` takes two method-specific parameters as described below.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`].
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
#' @param alpha Numeric vector of length 1.  The exponent defining the
#' weight of the tail in the random walk performed by the `ssGSEA` (Barbie et
#' al., 2009) method.  The default value is 0.25 as described in the paper.
#' 
#' @param normalize Logical vector of length 1; if `TRUE`  runs the `ssGSEA` method
#' from Barbie et al. (2009) normalizing the scores by the absolute difference
#' between the minimum and the maximum, as described in their paper. Otherwise
#' this last normalization step is skipped.
#' 
#' @return A new [`ssgseaParam-class`] object.
#'
#' @references Barbie, D.A. et al. Systematic RNA interference reveals that
#' oncogenic KRAS-driven cancers require TBK1.
#' *Nature*, 462(5):108-112, 2009.
#' [DOI](https://doi.org/10.1038/nature08460)
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
#' sp1 <- ssgseaParam(ses, gsc)
#' sp1
#'
#' @importFrom methods new
#' @rdname ssgseaParam-class
#' 
#' @export
ssgseaParam <- function(exprData, geneSets,
                        assay=NA_character_, annotation=NA_character_,
                        minSize=1,maxSize=Inf,
                        alpha=0.25, normalize=TRUE) {
    an <- gsvaAssayNames(exprData)
    if((!is.na(assay)) && (!.isCharNonEmpty(an)))
        warning("argument assay='", assay,
                "' ignored since exprData has no assayNames()")
    if(is.na(assay) && .isCharNonEmpty(an))
        assay <- na.omit(an)[1]
    
    new("ssgseaParam",
        exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize,
        alpha=alpha, normalize=normalize)
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
    if(length(object@alpha) != 1) {
        inv <- c(inv, "@alpha should be of length 1")
    }
    if(is.na(object@alpha)) {
        inv <- c(inv, "@alpha should not be NA")
    }
    if(length(object@normalize) != 1) {
        inv <- c(inv, "@normalize should be of length 1")
    }
    if(is.na(object@normalize)) {
        inv <- c(inv, "@normalize should not be NA")
    }
    return(if(length(inv) == 0) TRUE else inv)
})


## ----- getters -----

#' @noRd
get_alpha <- function(object) {
  stopifnot(inherits(object, "ssgseaParam"))
  return(object@alpha)
}

#' @noRd
do_normalize <- function(object) {
  stopifnot(inherits(object, "ssgseaParam"))
  return(object@normalize)
}


## ----- show -----

setMethod("show",
          signature=signature(object="ssgseaParam"),
          function(object) {
              callNextMethod(object)
              cat("alpha: ", get_alpha(object), "\n", 
                  "normalize: ", do_normalize(object), "\n", sep="")
          })
