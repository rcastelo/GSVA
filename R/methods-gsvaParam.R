
## ----- constructor -----

#' Construct a `GSVA` parameter object
#'
#' Construct and return a new object of class \linkS4class{gsvaParam}.
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
#' @return A new \linkS4class{gsvaParam} object.
#'
## #' @examples
## #' gp <- gsvaParam()
## #'
## #' @importFrom methods new
#' @rdname gsvaParam
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
    invalid <- character()
    dd <- dim(object@exprData)
    if(dd[1] == 0) invalid <- c(invalid, "@exprData has 0 rows")
    if(dd[2] == 0) invalid <- c(invalid, "@exprData has 0 columns")
    if(length(object@geneSets) == 0)
        invalid <- c(invalid, "@geneSets has length 0")
    if(length(object@tau) != 1)
        invalid <- c(invalid, "@tau should be of length 1")
    if(is.na(object@tau)) invalid <- c(invalid, "@tau should not be NA")
    if(length(object@maxDiff) != 1)
        invalid <- c(invalid, "@maxDiff should be of length 1")
    if(is.na(object@maxDiff)) invalid <- c(invalid, "@maxDiff should not be NA")
    if(length(object@absRanking) != 1)
        invalid <- c(invalid, "@absRanking should be of length 1")
    if(is.na(object@absRanking))
        invalid <- c(invalid, "@absRanking should not be NA")
    return(if(length(invalid) == 0) TRUE else invalid)
})


## ----- getters -----

#' Return the `kcdf` parameter from a `gsvaParam` object
#'
#' Returns the `kcdf` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param object A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `kcdf` parameter from `object` or an error if `object`
#' does not inherit from \linkS4class{gsvaParam}.
#'
## #' @examples
## #' gp <- gsvaParam()
## #' get_kcdf(gp)
## #' gp <- gsvaParam(kcdf="Poisson")
## #' get_kcdf(gp)
## #'
#' @noRd
get_kcdf <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@kcdf)
}


#' Return the `tau` parameter from a `gsvaParam` object
#'
#' Returns the `tau` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param object A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `tau` parameter from `object` or an error if `object`
#' does not inherit from \linkS4class{gsvaParam}.
#'
## #' @examples
## #' gp <- gsvaParam()
## #' get_tau(gp)
## #' gp <- gsvaParam(tau=0.42)
## #' get_tau(gp)
## #'
#' @noRd
get_tau <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@tau)
}


#' Return the `maxDiff` parameter from a `gsvaParam` object
#'
#' Returns the `maxDiff` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param object A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `maxDiff` parameter from `object` or an error if
#' `object` does not inherit from \linkS4class{gsvaParam}.
#'
## #' @examples
## #' gp <- gsvaParam()
## #' get_maxDiff(gp)
## #' gp <- gsvaParam(maxDiff=FALSE)
## #' get_maxDiff(gp)
## #'
#' @noRd
get_maxDiff <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@maxDiff)
}


#' Return the `absRanking` parameter from a `gsvaParam` object
#'
#' Returns the `absRanking` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param object A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `absRanking` parameter from `object` or an error if
#' `object` does not inherit from \linkS4class{gsvaParam}.
#'
## #' @examples
## #' gp <- gsvaParam()
## #' get_absRanking(gp)
## #' gp <- gsvaParam(absRanking=TRUE)
## #' get_absRanking(gp)
## #'
#' @noRd
get_absRanking <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@absRanking)
}


## ----- show -----

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
