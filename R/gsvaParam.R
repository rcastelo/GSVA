
#' @title The `gsvaParam` class
#'
#' @description Objects of class `gsvaParam` contain the parameters for running
#' the `GSVA` method.
#'
#' @details In addition to an expression data set and a collection of
#' gene sets, `GSVA` takes four method-specific parameters as described below.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`].
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].
#' 
#' @param kcdf Character vector of length 1 denoting the kernel to use during
#' the non-parametric estimation of the cumulative distribution function of
#' expression levels across samples.  By default, `kcdf="Gaussian"` which is
#' suitable when input expression values are continuous, such as microarray
#' fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or
#' log-TPMs.  When input expression values are integer counts, such as those
#' derived from RNA-seq experiments, then this argument should be set to
#' `kcdf="Poisson"`.
#'
#' @param tau Numeric vector of length 1.  The exponent defining the weight of
#' the tail in the random walk performed by the `GSVA` (Hänzelmann et al.,
#' 2013) method.  The default value is 1 as described in the paper.
#'
#' @param maxDiff Logical vector of length 1 which offers two approaches to
#' calculate the enrichment statistic (ES) from the KS random walk statistic.
#' * `FALSE`: ES is calculated as the maximum distance of the random walk
#' from 0.
#' * `TRUE` (the default): ES is calculated as the magnitude difference between
#' the largest positive and negative random walk deviations.
#'
#' @param absRanking Logical vector of length 1 used only when `maxDiff=TRUE`.
#' When `absRanking=FALSE` (default) a modified Kuiper statistic is used to
#' calculate enrichment scores, taking the magnitude difference between the
#' largest positive and negative random walk deviations. When
#' `absRanking=TRUE` the original Kuiper statistic that sums the largest
#' positive and negative random walk deviations, is used. In this latter case,
#' gene sets with genes enriched on either extreme (high or low) will be
#' regarded as ’highly’ activated.
#' 
#' @return A new [`gsvaParam-class`] object.
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' [DOI](https://doi.org/10.1186/1471-2105-14-7)
#'
#' @examples
#' library("GSEABase")
#' data(sample.ExpressionSet, package="Biobase")
#' 
#' ses <- sample.ExpressionSet[201:242,]
#' gsc <- GeneSetCollection(ses, setType=GOCollection())
#' gp1 <- gsvaParam(ses, gsc)
#' gp1
#' 
#' xes <- exprs(ses)
#' lgs <- geneIds(gsc)
#' gp2 <- gsvaParam(exprData=xes, geneSets=lgs, kcdf="Poisson", tau=0.42, maxDiff=FALSE)
#' gp2
#'
#' @importFrom methods new
#' @rdname gsvaParam-class
#' 
#' @export
gsvaParam <- function(exprData, geneSets, kcdf=c("Gaussian", "Poisson", "none"),
                      tau=1, maxDiff=TRUE, absRanking=FALSE) {
  kcdf <- match.arg(kcdf)
  new("gsvaParam",
      exprData=exprData, geneSets=geneSets,
      kcdf=kcdf, tau=tau, maxDiff=maxDiff, absRanking=absRanking)
}


## ----- validator -----

setValidity("gsvaParam", function(object) {
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
    if(length(object@tau) != 1) {
        inv <- c(inv, "@tau should be of length 1")
    }
    if(is.na(object@tau)) {
        inv <- c(inv, "@tau should not be NA")
    }
    if(length(object@maxDiff) != 1) {
        inv <- c(inv, "@maxDiff should be of length 1")
    }
    if(is.na(object@maxDiff)) {
        inv <- c(inv, "@maxDiff should not be NA")
    }
    if(length(object@absRanking) != 1) {
        inv <- c(inv, "@absRanking should be of length 1")
    }
    if(is.na(object@absRanking)) {
        inv <- c(inv, "@absRanking should not be NA")
    }
    return(if(length(inv) == 0) TRUE else inv)
})


## ----- getters -----

#' @noRd
get_kcdf <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@kcdf)
}

#' @noRd
get_tau <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@tau)
}

#' @noRd
get_maxDiff <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@maxDiff)
}

#' @noRd
get_absRanking <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@absRanking)
}


## ----- show -----

#' @exportMethod show
setMethod("show",
          signature=signature(object="gsvaParam"),
          function(object) {
              some <- function(x)
                  paste0(paste(Biobase::selectSome(x, 4), collapse=", "),
                        " (", length(x), " total)")
              ds <- get_exprData(object)
              dsDim <- sprintf(" [%s, %d]", nrow(ds), ncol(ds))
              gs <- get_geneSets(object)
              gsDim <- sprintf(" [%d, %d]", nrow(gs), ncol(gs))
              gs
              cat("GSVA Parameter object\n",
                  "  data set: ", class(ds)[1], dsDim, "\n",
                  "    rows: ", some(rownames(ds)), "\n",
                  "      (annotation: ", annotation(ds), ")", "\n",
                  "    columns: ", some(colnames(ds)), "\n",
                  "  gene sets: ", class(gs)[1], gsDim, "\n",
                  "    names: ", some(names(gs)), "\n",
                  sep="")
          })
