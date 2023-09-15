
## ----- constructor -----

#' Construct a combined z-scores parameter object
#'
#' Construct and return a new object of class \linkS4class{zscoreParam}.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`].
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].
#' 
#' @return A new \linkS4class{zscoreParam} object.
#'
## #' @examples
## #' zp <- zscoreParam()
## #'
## #' @importFrom methods new
#' @rdname zscoreParam
#' @export
zscoreParam <- function(exprData, geneSets) {
  new("zscoreParam", exprData=exprData, geneSets=geneSets)
}


## ----- validator -----

setValidity("zscoreParam", function(object) {
    invalid <- character()
    dd <- dim(object@exprData)
    if(dd[1] == 0) invalid <- c(invalid, "@exprData has 0 rows")
    if(dd[2] == 0) invalid <- c(invalid, "@exprData has 0 columns")
    if(length(object@geneSets) == 0) invalid <- c(invalid, "@geneSets has length 0")
    return(if(length(invalid) == 0) TRUE else invalid)
})


## ----- show -----

setMethod("show",
          signature=signature(object="zscoreParam"),
          function(object) {
              some <- function(x)
                  paste0(paste(Biobase::selectSome(x, 4), collapse=", "),
                        " (", length(x), " total)")
              ds <- get_exprData(object)
              dsDim <- sprintf(" [%s, %d]", nrow(ds), ncol(ds))
              gs <- get_geneSets(object)
              gsDim <- sprintf(" [%d, %d]", nrow(gs), ncol(gs))
              gs
              cat("Combined z-Scores Parameter object\n",
                  "  data set: ", class(ds)[1], dsDim, "\n",
                  "    rows: ", some(rownames(ds)), "\n",
                  "      (annotation: ", annotation(ds), ")", "\n",
                  "    columns: ", some(colnames(ds)), "\n",
                  "  gene sets: ", class(gs)[1], gsDim, "\n",
                  "    names: ", some(names(gs)), "\n",
                  sep="")
          })
