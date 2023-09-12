
### PLAGE Parameter Constructor ---------------------------------------------

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
plageParam <- function(dataSet, geneSets) {
  new("plageParam", dataSet = dataSet, geneSets = geneSets)
}


### PLAGE Parameter Validator ---------------------------------------------

setValidity("plageParam", function(object) {
    invalid <- character()
    dd <- dim(object@dataSet)
    if(dd[1] == 0) invalid <- c(invalid, "@dataSet has 0 rows")
    if(dd[2] == 0) invalid <- c(invalid, "@dataSet has 0 columns")
    if(length(object@geneSets) == 0) invalid <- c(invalid, "@geneSets has length 0")
    return(if(length(invalid) == 0) TRUE else invalid)
})

