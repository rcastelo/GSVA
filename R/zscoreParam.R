
#
# 2023-04-25  axel: methods for class zscoreParam
#

# Constructor -------------------------------------------------------------

#' Build a \code{zscore} parameter object
#'
#' Build an object of the class \code{zscoreParam}
#'
#' @return A \linkS4class{zscoreParam} object.
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
