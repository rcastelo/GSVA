
#
# 2023-04-25  axel: methods for class ssgseaParam
#

# Constructor -----------------------------------------------------------

#' Build a \code{ssGSEA} parameter object
#'
#' Build an object of the class \code{ssgseaParam}
#'
#' @return A \linkS4class{ssgseaParam} object.
#'
#' @examples
#' pp <- ssgseaParam()
#'
#' @importFrom methods new
#' @rdname ssgseaParam
#' @export
ssgseaParam <- function(alpha = 0.25, normalize = TRUE) {
  new("ssgseaParam", alpha = alpha, normalize = normalize)
}



# Getters -----------------------------------------------------------------

setMethod("getAlpha",
          signature = c(obj = "ssgseaParam"),
          definition = function(obj) return(obj@alpha))

setMethod("doNormalize",
          signature = c(obj = "ssgseaParam"),
          definition = function(obj) return(obj@normalize))
