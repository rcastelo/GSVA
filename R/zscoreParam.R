
# z-Score Parameter Constructor -------------------------------------------

#' Construct a z-score parameter object
#'
#' Construct and return a new object of class \linkS4class{zscoreParam}.
#'
#' @return A new \linkS4class{zscoreParam} object.
#'
#' @examples
#' zp <- zscoreParam()
#'
#' @importFrom methods new
#' @rdname zscoreParam
#' @export
zscoreParam <- function() {
  new("zscoreParam")
}
