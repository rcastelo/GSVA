
#' @title The `ssgseaParam` class
#'
#' @description Objects of class `ssgseaParam` contain the parameters for
#' running the `ssGSEA` method.
#'
#' @details In addition to a number of parameters shared with all methods
#' implemented by package GSVA, `ssGSEA` takes two method-specific parameters as
#' well as two more parameters for implementing a missing value policy.  All of
#' these parameters are described in detail below.
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
#' the first assay is used.
#' 
#' @param annotation An object of class [`GeneIdentifierType-class`] from
#' package `GSEABase` describing the gene identifiers used as the row names of
#' the expression data set.  See [`GeneIdentifierType`] for help on available
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
#' @param alpha Numeric vector of length 1.  The exponent defining the
#' weight of the tail in the random walk performed by the `ssGSEA` (Barbie et
#' al., 2009) method.  The default value is 0.25 as described in the paper.
#' 
#' @param normalize Logical vector of length 1; if `TRUE` runs the `ssGSEA`
#' method from Barbie et al. (2009) normalizing the scores by the absolute
#' difference between the minimum and the maximum, as described in their paper.
#' Otherwise this final normalization step is skipped.
#' 
#' @param checkNA Character vector of length 1 specifying whether the input
#' expression data should be checked for the presence of missing (`NA`) values.
#' This must be
#' one of the strings `"auto"` (default), `"yes"`, or `"no"`. The default value
#' `"auto"` means that the software will perform that check only when the input
#' expression data is provided as a base [`matrix`], an [`ExpressionSet`] or a
#' [`SummarizedExperiment`] object, while every other type of input expression
#' data container (e.g., [`SingleCellExperiment`], etc.) will not be checked.
#' If `checkNA="yes"`, then the input expression data will be checked for
#' missing values irrespective of the object class of the data container, and
#' if `checkNA="no"`, then that check will not be performed.
#'
#' @param use Character vector of length 1 specifying a policy for dealing with
#' missing values (`NA`s) in the input expression data argument `exprData`. It
#' only applies when either `checkNA="yes"`, or `checkNA="auto"` (see the
#' `checkNA` parameter. The argument value must be one of the strings
#' `"everything"` (default), `"all.obs"`, or `"na.rm"`. The policy of the
#' default value `"everything"` consists of propagating `NA`s so that the
#' resulting enrichment score will be `NA`, whenever one or more of its
#' contributing values is `NA`, giving a warning when that happens. When
#' `use="all.obs"`, the presence of `NA`s in the input expression data will
#' produce an error. Finally, when `use="na.rm"`, `NA` values in the input
#' expression data will be removed from calculations, giving a warning when that
#' happens, and giving an error if no values are left after removing the `NA`
#' values.
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
                        assay=NA_character_, annotation=NULL,
                        minSize=1,maxSize=Inf,
                        alpha=0.25, normalize=TRUE,
                        checkNA=c("auto", "yes", "no"),
                        use=c("everything", "all.obs", "na.rm")) {
    checkNA <- match.arg(checkNA)
    use <- match.arg(use)
    an <- gsvaAssayNames(exprData)
    if((!is.na(assay)) && (!.isCharNonEmpty(an))) {
        msg <- sprintf(paste0("argument assay='%s' ignored since exprData has ",
                              "no assayNames()"), assay)
        cli_alert_info(msg)
    }
    if(is.na(assay) && .isCharNonEmpty(an))
        assay <- na.omit(an)[1]

    ## check for presence of valid row/feature names
    .check_rownames(exprData)

    xa <- gsvaAnnotation(exprData)
    if(is.null(xa)) {
        if(is.null(annotation)) {
            annotation <- NullIdentifier()
        }
    } else {
        if(is.null(annotation)) {
            annotation <- xa
        } else {
            msg <- sprintf(paste0("using argument annotation='%s' and ",
                                  "ignoring exprData annotation ('%s')"),
                           capture.output(annotation), capture.output(xa))
            cli_alert_info(msg)
        }
    }

    naparam <- .check_for_na_values(exprData=exprData, assay=assay,
                                    checkNA=checkNA, use=use)
    
    new("ssgseaParam",
        exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize,
        alpha=alpha, normalize=normalize,
        checkNA=checkNA, didCheckNA=naparam$didCheckNA,
        anyNA=naparam$any_na, use=use)
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
        inv <- c(inv, "@use must be a single character string")
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
    return(if(length(inv) == 0) TRUE else inv)
})


## ----- getters -----

#' @noRd
get_alpha <- function(object) {
  stopifnot(inherits(object, "ssgseaParam"))
  return(object@alpha)
}

#' @noRd
get_normalize <- function(object) {
  stopifnot(inherits(object, "ssgseaParam"))
  return(object@normalize)
}

#' @noRd
get_checkNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam"))
  return(object@checkNA)
}

#' @noRd
get_didCheckNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam"))
  return(object@didCheckNA)
}

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

#' @noRd
get_NAuse <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam"))
  return(object@use)
}

## ----- show -----

setMethod("show",
          signature=signature(object="ssgseaParam"),
          function(object) {
              callNextMethod(object)
              cat("alpha: ", get_alpha(object), "\n", 
                  "normalize: ", get_normalize(object), "\n",
                  "checkNA: ", get_checkNA(object), "\n", sep="")
              if (get_didCheckNA(object)) {
                  if (anyNA(object)) {
                      cat("missing data: yes\n",
                          "na_use: ", get_NAuse(object), "\n", sep="")
                  } else
                      cat("missing data: no\n")
              } else
                  cat("missing data: didn't check\n")
          })
