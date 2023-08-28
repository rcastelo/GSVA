
# PLAGE Parameter Constructor ---------------------------------------------

#' Construct a PLAGE parameter object
#'
#' Construct and return a new object of class \linkS4class{PlageParam}.
#'
#' @return A new \linkS4class{PlageParam} object.
#'
#' @examples
#' pp <- plageParam()
#'
#' @importFrom methods new
#' @rdname plageParam
#' @export
plageParam <- function() {
  new("PlageParam")
}

