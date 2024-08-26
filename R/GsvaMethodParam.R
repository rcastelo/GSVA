
#' @title The `GsvaMethodParam` class
#'
#' @description A virtual superclass of the `GSVA` packages' method-specific
#' parameter classes.
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
              cat("A ", .objPkgClass(object), " object\n",
                  "expression data:\n", sep="")
              .catObj(get_exprData(object))
              oa <- if(is.na(get_assay(object))) "none" else get_assay(object)
              cat("using assay: ", oa, "\n", sep="")
              oan <- if(is.na(get_annotation(object))) "none" else get_annotation(object)
              cat("using annotation: ", oan, "\n", sep="")
              cat("gene sets:\n")
              .catObj(get_geneSets(object))
              cat("gene set size: [", get_minSize(object), ", ",
                  get_maxSize(object),  "]\n", sep="")
          })


## ----- getters -----

setMethod("get_exprData", signature("GsvaMethodParam"),
          function(object) {
              return(object@exprData)
          })

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


## ----- uniform access to annotation -----

setMethod("gsvaAnnotation",
          signature=signature(object="GsvaExprData"),
          function(object) {
              ## in general
              return(NULL)
          })

setMethod("gsvaAnnotation",
          signature=signature(object="ExpressionSet"),
          function(object) {
              ## always a character giving the db pkg, potentially empty ("")
              return(annotation(object))
          })

setMethod("gsvaAnnotation", signature("SummarizedExperiment"),
          function(object) {
              ## NULL if unset; otherwise anything but we *expect* and handle
              ## a GSEABase::GeneIdentifierType with or without annotation(),
              ## i.e., db pkg, available.  Same for subclasses below.
              return(metadata(object)$annotation)
          })

setMethod("gsvaAnnotation", signature("SingleCellExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })

setMethod("gsvaAnnotation", signature("SpatialExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })



## ----- uniform access to assay names -----

setMethod("gsvaAssayNames",
          signature=signature(object="GsvaExprData"),
          function(object) {
              return(NA_character_)
          })

setMethod("gsvaAssayNames", signature("SummarizedExperiment"),
          function(object) {
              a <- assayNames(object)
              return(if(.isCharNonEmpty(a)) a else NA_character_)
          })

setMethod("gsvaAssayNames", signature("SingleCellExperiment"),
          function(object) {
              a <- assayNames(object)
              return(if(.isCharNonEmpty(a)) a else NA_character_)
          })
