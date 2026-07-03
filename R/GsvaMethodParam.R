
#' @title The `GsvaMethodParam` class
#'
#' @description A virtual superclass of the `GSVA` packages' method-specific
#' parameter classes. The method 'details()' provides a detailed summary of the
#' parameter values stored in this class and its subclasses.
#'
#' @param object An object of class `GsvaMethodParam` or one of its subclasses.
#'
#' @details The `GSVA` package implements four single-sample gene set analysis
#' methods (PLAGE, combined z-scores, ssGSEA, and GSVA) and a respective
#' method-specific parameter class that is used to invoke each of them with a
#' matching set of parameters.
#'
#' @seealso [`plageParam`], [`zscoreParam`], [`ssgseaParam`], [`gsvaParam`]
#' 
#' @name GsvaMethodParam-class
#' @rdname GsvaMethodParam-class
NULL


## ----- show -----

setMethod("show",
          signature=signature(object="GsvaMethodParam"),
          function(object) {
              cat("class: ", .objPkgClass(object), "\n", sep="")
              expr <- get_exprData(object)
              cat(sprintf("expression data dim: %d %d\n",
                          nrow(expr), ncol(expr)))
              cat(sprintf("number of gene sets: %d\n",
                          length(get_geneSets(object))))
              cat("details: use 'details(object)'\n")
          })

#' @aliases details,GsvaMethodParam-method
#' @rdname GsvaMethodParam-class
#' @exportMethod details
#' @importFrom GSEABase details
setMethod("details",
          signature=signature(object="GsvaMethodParam"),
          function(object) {
              cat("class: ", .objPkgClass(object), "\n",
                  "expression data:\n", sep="")
              .catObj(get_exprData(object))
              oa <- if(is.na(get_assay(object))) "none" else get_assay(object)
              cat("using assay: ", oa, "\n", sep="")
              cat("using annotation:\n")
              .catObj(get_annotation(object))
              cat("gene sets:\n")
              .catObj(get_geneSets(object))
              cat("gene set size: [", get_minSize(object), ", ",
                  get_maxSize(object),  "]\n", sep="")
              nzcmsg <- sprintf("nonzero values: %s than 2^31 (INT_MAX)\n",
                                ifelse(nzcount(object) > .Machine$integer.max,
                                       "more", "less"))
              cat(nzcmsg)
              cat("ondisk: ", .get_ondisk(object), "\n")
          })




## ----- getters -----

## somehow a method get_exprData() cannot dispatch on the members of the
## SummarizedExperiment class and the union GsvaExprData class hierarchy,
## and we would have to define one for each of them separately; see
## https://github.com/rcastelo/GSVA/issues/156
## https://stat.ethz.ch/pipermail/bioc-devel/2025-January/020788.html
## for this reason we define this getter using a function

## @importFrom cli cli_abort
get_exprData <- function(object) {
    if (is(object, "GsvaMethodParam"))
        return(object@exprData)
    else if (is(object, "GsvaExprData"))
        return(object)
    else
        cli_abort("object must be either 'GsvaMethodParam' or 'GsvaExprData'")
}

setMethod("get_geneSets", signature("GsvaMethodParam"),
          function(object) {
              return(object@geneSets)
          })

setMethod("get_assay", signature("GsvaMethodParam"),
          function(object) {
              return(object@assay)
          })

setMethod("get_annotation", signature("GsvaMethodParam"),
          function(object) {
              return(object@annotation)
          })

setMethod("get_minSize", signature("GsvaMethodParam"),
          function(object) {
              return(object@minSize)
          })

setMethod("get_maxSize", signature("GsvaMethodParam"),
          function(object) {
              return(object@maxSize)
          })

setMethod("nzcount", signature=c("GsvaMethodParam"),
          function(x) {
            return(x@nzcount)
          })



## ----- show component objects without overriding their show() method -----

setMethod("gsvaShow",
          signature=signature(object="GsvaExprData"),
          function(object) {
              show(object)
          })

setMethod("gsvaShow",
          signature=signature(object="matrix"),
          function(object) {
              cat("matrix [", nrow(object), ", ", ncol(object), "]\n",
                  "  rows: ", .showSome(rownames(object)), "\n",
                  "  cols: ", .showSome(colnames(object)), "\n", sep="")
          })

#' @importFrom Matrix nnzero
setMethod("gsvaShow",
          signature=signature(object="dgCMatrix"),
          function(object) {
              cat("dgCMatrix [", nrow(object), ", ", ncol(object), "]\n",
                  "  rows: ", .showSome(rownames(object)), "\n",
                  "  cols: ", .showSome(colnames(object)), "\n",
                  "  sparsity: ", 1 - nnzero(object) / length(object), "\n", sep="")
          })

setMethod("gsvaShow",
          signature=signature(object="GsvaGeneSets"),
          function(object) {
              show(object)
          })

setMethod("gsvaShow",
          signature=signature(object="list"),
          function(object) {
              cat("list\n",
                  "  names: ", .showSome(names(object)), "\n",
                  "  unique identifiers: ",
                  .showSome(unique(unname(unlist(object)))), "\n", sep="")
          })

## as it turns out, GeneSetCollection is() a list as well as a GsvaGeneSets
## and list is 'older' and hence wins when dispatching gsvaShow() :-|
setMethod("gsvaShow",
          signature=signature(object="GeneSetCollection"),
          function(object) {
              show(object)
          })

setMethod("gsvaShow",
          signature=signature(object="GeneIdentifierType"),
          function(object) {
              show(object)
          })


## ----- uniform access to assay names -----
## somehow a method gsvaAssayNames() cannot dispatch on the members of the
## SummarizedExperiment class and the union GsvaExprData class hierarchy,
## and we would have to define one for each of them separately; see
## https://github.com/rcastelo/GSVA/issues/156
## https://stat.ethz.ch/pipermail/bioc-devel/2025-January/020788.html
## for this reason we define this getter using a function

#' @importFrom cli cli_abort
#' @importFrom SummarizedExperiment assayNames
gsvaAssayNames <- function(object) {
    if (is(object, "SummarizedExperiment")) {
        a <- assayNames(object)
        return(if(.isCharNonEmpty(a)) a else NA_character_)
    } else if (is(object, "GsvaExprData")) {
        if (!is.null(attr(object, "assay")))
            return(attr(object, "assay"))
        return(NA_character_)
    } else
        cli_abort("object must be either 'SummarizedExperiment' or 'GsvaExprData'")
}
