
# ssGSEA Parameter Constructor --------------------------------------------

#' Construct a ssGSEA parameter object
#'
#' Construct and return a new object of class \linkS4class{ssgseaParam}.
#'
#' @param alpha Numeric vector of length 1; the exponent defining the
#'  weight of the tail in the random walk performed by the `ssGSEA` (Barbie et
#'  al., 2009) method.  The default value is 0.25 as described in the paper.
#' 
#' @param normalize Logical vector of length 1; if `TRUE`  runs the `ssGSEA` method
#'  from Barbie et al. (2009) normalizing the scores by the absolute difference
#'  between the minimum and the maximum, as described in their paper. Otherwise
#'  this last normalization step is skipped.
#' 
#' @return A new \linkS4class{ssgseaParam} object.
#'
#' @examples
#' sp <- ssgseaParam()
#'
#' @importFrom methods new
#' @rdname ssgseaParam
#' @export
ssgseaParam <- function(alpha = 0.25, normalize = TRUE) {
  new("ssgseaParam", alpha = alpha, normalize = normalize)
}


# Getters -----------------------------------------------------------------

#' Return the `alpha` parameter from a `ssgseaParam` object
#'
#' Returns the `alpha` parameter from a \linkS4class{ssgseaParam} object.
#'
#' @param obj A \linkS4class{ssgseaParam} object.
#' 
#' @return The requested `alpha` parameter from `obj` or an error if `obj` does
#' not inherit from \linkS4class{ssgseaParam}.
#'
#' @examples
#' sp <- ssgseaParam()
#' getAlpha(sp)
#' sp <- ssgseaParam(alpha = 0.42)
#' getAlpha(sp)
#'
#' @noRd
#' 
# #' @rdname ssgseaParam
# #' @export
getAlpha <- function(obj) {
  stopifnot(inherits(obj, "ssgseaParam"))
  return(obj@alpha)
}


#' Return the `normalize` flag from a `ssgseaParam` object
#'
#' Returns the `normalize` flag from a \linkS4class{ssgseaParam} object.
#'
#' @param obj A \linkS4class{ssgseaParam} object.
#' 
#' @return The requested `normalize` flag from `obj` or an error if `obj` does
#' not inherit from \linkS4class{ssgseaParam}.
#'
#' @examples
#' sp <- ssgseaParam()
#' doNormalize(sp)
#' sp <- ssgseaParam(doNormalize = FALSE)
#' doNormalize(sp)
#'
#' @noRd
#' 
# #' @rdname ssgseaParam
# #' @export
doNormalize <- function(obj) {
  stopifnot(inherits(obj, "ssgseaParam"))
  return(obj@normalize)
}
