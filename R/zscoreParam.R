
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
zscoreParam <- function(dataSet, geneSets) {
  new("zscoreParam", dataSet = dataSet, geneSets = geneSets)
}


### z-Score Parameter Validator ---------------------------------------------

setValidity("zscoreParam", function(object) {
    invalid <- character()
    dd <- dim(object@dataSet)
    if(dd[1] == 0) invalid <- c(invalid, "@dataSet has 0 rows")
    if(dd[2] == 0) invalid <- c(invalid, "@dataSet has 0 columns")
    if(length(object@geneSets) == 0) invalid <- c(invalid, "@geneSets has length 0")
    return(if(length(invalid) == 0) TRUE else invalid)
})
