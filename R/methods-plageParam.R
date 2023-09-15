
## ----- constructor -----

#' Construct a `PLAGE` parameter object
#'
#' Construct and return a new object of class \linkS4class{plageParam}.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`].
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].
#' 
#' @return A new [`plageParam`] object.
#'
## #' @examples
## #' pp <- plageParam()
## #'
## #' @importFrom methods new
#' @rdname plageParam
#' @export
plageParam <- function(exprData, geneSets) {
  new("plageParam", exprData=exprData, geneSets=geneSets)
}


## ----- validator -----

setValidity("plageParam", function(object) {
    invalid <- character()
    dd <- dim(object@exprData)
    if(dd[1] == 0) invalid <- c(invalid, "@exprData has 0 rows")
    if(dd[2] == 0) invalid <- c(invalid, "@exprData has 0 columns")
    if(length(object@geneSets) == 0) invalid <- c(invalid, "@geneSets has length 0")
    return(if(length(invalid) == 0) TRUE else invalid)
})


## ----- show -----

setMethod("show",
          signature=signature(object="plageParam"),
          function(object) {
              some <- function(x)
                  paste0(paste(Biobase::selectSome(x, 4), collapse=", "),
                        " (", length(x), " total)")
              ds <- get_exprData(object)
              dsDim <- sprintf(" [%s, %d]", nrow(ds), ncol(ds))
              gs <- get_geneSets(object)
              gsDim <- sprintf(" [%d, %d]", nrow(gs), ncol(gs))
              gs
              cat("PLAGE Parameter object\n",
                  "  data set: ", class(ds)[1], dsDim, "\n",
                  "    rows: ", some(rownames(ds)), "\n",
                  "      (annotation: ", annotation(ds), ")", "\n",
                  "    columns: ", some(colnames(ds)), "\n",
                  "  gene sets: ", class(gs)[1], gsDim, "\n",
                  "    names: ", some(names(gs)), "\n",
                  sep="")
          })
