
#' @title The `zscoreParam` class
#'
#' @description Objects of class `zscoreParam` contain the parameters for running
#' the combined z-scores method.
#'
#' @details The combined z-scores method does not take any method-specific
#' parameters in addition to an expression data set and a collection of gene
#' sets.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`].
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].
#' 
#' @return A new [`zscoreParam-class`] object.
#'
#' @references Lee, E. et al. Inferring pathway activity toward precise
#' disease classification.
#' *PLoS Comp Biol*, 4(11):e1000217, 2008.
#' [DOI](https://doi.org/10.1371/journal.pcbi.1000217)
#'
#' @examples
#' library("GSEABase")
#' data(sample.ExpressionSet, package="Biobase")
#' 
#' ses <- sample.ExpressionSet[201:242,]
#' gsc <- GeneSetCollection(ses, setType=GOCollection())
#' zp1 <- zscoreParam(ses, gsc)
#' zp1
#' 
#' xes <- exprs(ses)
#' lgs <- geneIds(gsc)
#' zp2 <- zscoreParam(exprData=xes, geneSets=lgs)
#' zp2
#'
#'
#' @importFrom methods new
#' @rdname zscoreParam-class
#' 
#' @export
zscoreParam <- function(exprData, geneSets) {
  new("zscoreParam", exprData=exprData, geneSets=geneSets)
}


## ----- validator -----

setValidity("zscoreParam", function(object) {
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
    return(if(length(inv) == 0) TRUE else inv)
})


## ----- show -----

#' @exportMethod show
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
