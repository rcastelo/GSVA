
#
# 2023-04-25  axel: methods for class plageParam
#

# Constructor -----------------------------------------------------------

#' Build a \code{plage} parameter object
#'
#' Build an object of the class \code{plageParam}
#'
#' @return A \linkS4class{plageParam} object.
#'
#' @examples
#' pp <- plageParam()
#'
#' @importFrom methods new
#' @rdname plageParam
#' @export
plageParam <- function() {
  new("plageParam")
}

