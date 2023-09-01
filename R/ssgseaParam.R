
# ssGSEA Parameter Constructor --------------------------------------------

#' Construct a ssGSEA parameter object
#'
#' Construct and return a new object of class \linkS4class{SsGseaParam}.
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
#' @return A new \linkS4class{SsGseaParam} object.
#'
#' @examples
#' sp <- ssgseaParam()
#'
#' @importFrom methods new
#' @rdname ssgseaParam
#' @export
ssgseaParam <- function(alpha = 0.25, normalize = TRUE) {
  new("SsGseaParam", alpha = alpha, normalize = normalize)
}


# Getters -----------------------------------------------------------------

#' Return the `alpha` parameter from a `SsGseaParam` object
#'
#' Returns the `alpha` parameter from a \linkS4class{SsGseaParam} object.
#'
#' @param obj A \linkS4class{SsGseaParam} object.
#' 
#' @return The requested `alpha` parameter from `obj` or an error if `obj` does
#' not inherit from \linkS4class{SsGseaParam}.
#'
#' @examples
#' sp <- ssgseaParam()
#' get_alpha(sp)
#' sp <- ssgseaParam(alpha = 0.42)
#' get_alpha(sp)
#'
#' @noRd
get_alpha <- function(obj) {
  stopifnot(inherits(obj, "SsGseaParam"))
  return(obj@alpha)
}


#' Return the `normalize` flag from a `SsGseaParam` object
#'
#' Returns the `normalize` flag from a \linkS4class{SsGseaParam} object.
#'
#' @param obj A \linkS4class{SsGseaParam} object.
#' 
#' @return The requested `normalize` flag from `obj` or an error if `obj` does
#' not inherit from \linkS4class{SsGseaParam}.
#'
#' @examples
#' sp <- ssgseaParam()
#' do_normalize(sp)
#' sp <- ssgseaParam(doNormalize = FALSE)
#' do_normalize(sp)
#'
#' @noRd
do_normalize <- function(obj) {
  stopifnot(inherits(obj, "SsGseaParam"))
  return(obj@normalize)
}
