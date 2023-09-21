
#' @title The `GsvaMethodParam` class
#'
#' @description A virtual superclass of the `GSVA` packages' method-specific
#' parameter classes.
#'
#' @details The `GSVA` packages implements four single-sample gene set analysis
#' methods (`PLAGE`, combined z-scores, `ssGSEA`, and `GSVA`) and a respective
#' method-specific parameter class that is used to invoke each of them with a
#' matching set of parameters.
#'
#' @seealso [`plageParam`][`zscoreParam`][`ssgseaParam`][`gsvaParam`]
#' 
#' @aliases
#' get_exprData,GsvaMethodParam-method
#' get_geneSets,GsvaMethodParam-method
#' 
#' @rdname GsvaMethodParam-class


## ----- getters -----

#' @exportMethod get_exprData
setMethod("get_exprData", signature("GsvaMethodParam"),
          function(object) {
              return(object@exprData)
          })

#' @exportMethod get_geneSets
setMethod("get_geneSets", signature("GsvaMethodParam"),
          function(object) {
              return(object@geneSets)
          })
