
#' @title The `plageParam` class
#'
#' @description Objects of class `plageParam` contain the parameters for running
#' the `PLAGE` method.
#'
#' @details `PLAGE` does not take any method-specific parameters in addition to
#' an expression data set and a collection of gene sets.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`].
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].
#' 
#' @param assay The name of the assay to use in case `exprData` is a multi-assay
#' container, otherwise ignored.  By default, the first assay is used.
#' 
#' @return A new [`plageParam-class`] object.
#'
#' @references Tomfohr, J. et al. Pathway level analysis of gene expression
#' using singular value decomposition.
#' *BMC Bioinformatics*, 6:225, 2005.
#' [DOI](https://doi.org/10.1186/1471-2105-6-225)
#'
#' @examples
#' library(GSVA)
#' library(GSVAdata)
#'
#' data(leukemia)
#' data(c2BroadSets)
#'
#' ## for simplicity, use only a subset of the sample data
#' ses <- leukemia_eset[1:1000, ]
#' gsc <- c2BroadSets[1:100]
#' pp1 <- plageParam(ses, gsc)
#' pp1
#'
#' @importFrom methods new
#' @rdname plageParam-class
#' 
#' @export
plageParam <- function(exprData, geneSets, assay=NA_character_) {
  new("plageParam", exprData=exprData, geneSets=geneSets, assay=assay)
}


## ----- validator -----

setValidity("plageParam", function(object) {
    inv <- NULL
    dd <- dim(object@exprData)
    if(dd[1] == 0) {
        inv <- c(inv, "@exprData has 0 rows")
    }
    if(dd[2] == 0) {
        inv <- c(inv, "@exprData has 0 columns")
    }
    if(length(object@geneSets) == 0) {
        inv <- c(inv, "@geneSets has length 0")
    }
    if(length(object@assay) != 1) {
        inv <- c(inv, "@assay should be of length 1")
    }
    return(if(length(inv) == 0) TRUE else inv)
})


## ----- show -----

#' @exportMethod show
setMethod("show",
          signature=signature(object="plageParam"),
          function(object) {
              some <- function(x)
                  paste0(paste(Biobase::selectSome(x, 4), collapse=", "),
                        " (", length(x), " total)")
              ds <- get_exprData(object)
              dsDim <- sprintf(" [%d, %d]", nrow(ds), ncol(ds))
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
