
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
#' @seealso [`plageParam`][`zscoreParam`][`ssgseaParam`][`gsvaParam`]
#' 
#' @name GsvaMethodParam-class
#' @rdname GsvaMethodParam-class
NULL

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
