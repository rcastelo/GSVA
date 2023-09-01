
# z-Score Parameter Constructor -------------------------------------------

#' Construct a z-score parameter object
#'
#' Construct and return a new object of class \linkS4class{ZScoreParam}.
#'
#' @return A new \linkS4class{ZScoreParam} object.
#'
#' @examples
#' zp <- zscoreParam()
#'
#' @importFrom methods new
#' @rdname zscoreParam
#' @export
zscoreParam <- function() {
  new("ZScoreParam")
}
