
# PLAGE Parameter Constructor ---------------------------------------------

#' Construct a PLAGE parameter object
#'
#' Construct and return a new object of class \linkS4class{plageParam}.
#'
#' @return A new \linkS4class{plageParam} object.
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

