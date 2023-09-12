
# GSVA Parameter Constructor --------------------------------------------

#' Construct a GSVA parameter object
#'
#' Construct and return a new object of class \linkS4class{gsvaParam}.
#'
#' @param kcdf Character vector of length 1 denoting the kernel to use during
#'   the non-parametric estimation of the cumulative distribution function of
#'   expression levels across samples. By default, `kcdf="Gaussian"` which is
#'   suitable when input expression values are continuous, such as microarray
#'   fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or
#'   log-TPMs. When input expression values are integer counts, such as those
#'   derived from RNA-seq experiments, then this argument should be set to
#'   `kcdf="Poisson"`.
#'
#' @param tau Numeric vector of length 1; the exponent defining the weight of
#'   the tail in the random walk performed by the `GSVA` (Hänzelmann et al.,
#'   2013) method.  The default value is 1 as described in the paper.
#'
#' @param mx.diff Logical vector of length 1 which offers two approaches to
#'   calculate the enrichment statistic (ES) from the KS random walk statistic.
#'  * `FALSE`: ES is calculated as the maximum distance of the random walk
#'   from 0.
#'  * `TRUE` (the default): ES is calculated as the magnitude difference between
#'   the largest positive and negative random walk deviations.
#'
#' @param abs.ranking Logical vector of length 1 used only when `mx.diff=TRUE`.
#'   When `abs.ranking=FALSE` (default) a modified Kuiper statistic is used to
#'   calculate enrichment scores, taking the magnitude difference between the
#'   largest positive and negative random walk deviations. When
#'   `abs.ranking=TRUE` the original Kuiper statistic that sums the largest
#'   positive and negative random walk deviations, is used. In this latter case,
#'   gene sets with genes enriched on either extreme (high or low) will be
#'   regarded as ’highly’ activated.
#' 
#' @return A new \linkS4class{gsvaParam} object.
#'
#' @examples
#' gp <- gsvaParam()
#'
#' @importFrom methods new
#' @rdname gsvaParam
#' @export
gsvaParam <- function(dataSet, geneSets, kcdf = c("Gaussian", "Poisson", "none"), 
                      tau = 1, mx.diff = TRUE, abs.ranking = FALSE) {
  kcdf <- match.arg(kcdf)
  new("gsvaParam",
      dataSet = dataSet, geneSets = geneSets,
      kcdf = kcdf, tau = tau, mx.diff = mx.diff, abs.ranking = abs.ranking)
}


### GSVA Parameter Validator ---------------------------------------------

setValidity("gsvaParam", function(object) {
    invalid <- character()
    dd <- dim(object@dataSet)
    if(dd[1] == 0) invalid <- c(invalid, "@dataSet has 0 rows")
    if(dd[2] == 0) invalid <- c(invalid, "@dataSet has 0 columns")
    if(length(object@geneSets) == 0) invalid <- c(invalid, "@geneSets has length 0")
    if(length(object@tau) != 1) invalid <- c(invalid, "@tau should be of length 1")
    if(is.na(object@tau)) invalid <- c(invalid, "@tau should not be NA")
    if(length(object@mx.diff) != 1) invalid <- c(invalid, "@mx.diff should be of length 1")
    if(is.na(object@mx.diff)) invalid <- c(invalid, "@mx.diff should not be NA")
    if(length(object@abs.ranking) != 1) invalid <- c(invalid, "@abs.ranking should be of length 1")
    if(is.na(object@abs.ranking)) invalid <- c(invalid, "@abs.ranking should not be NA")
    return(if(length(invalid) == 0) TRUE else invalid)
})


# Getters -----------------------------------------------------------------

#' Return the `kcdf` parameter from a `gsvaParam` object
#'
#' Returns the `kcdf` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param obj A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `kcdf` parameter from `obj` or an error if `obj` does
#' not inherit from \linkS4class{gsvaParam}.
#'
#' @examples
#' gp <- gsvaParam()
#' get_kcdf(gp)
#' gp <- gsvaParam(kcdf = "Poisson")
#' get_kcdf(gp)
#'
#' @noRd
get_kcdf <- function(obj) {
  stopifnot(inherits(obj, "gsvaParam"))
  return(obj@kcdf)
}


#' Return the `tau` parameter from a `gsvaParam` object
#'
#' Returns the `tau` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param obj A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `tau` parameter from `obj` or an error if `obj` does
#' not inherit from \linkS4class{gsvaParam}.
#'
#' @examples
#' gp <- gsvaParam()
#' get_tau(gp)
#' gp <- gsvaParam(tau = 0.42)
#' get_tau(gp)
#'
#' @noRd
get_tau <- function(obj) {
  stopifnot(inherits(obj, "gsvaParam"))
  return(obj@tau)
}


#' Return the `mx.diff` parameter from a `gsvaParam` object
#'
#' Returns the `mx.diff` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param obj A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `mx.diff` parameter from `obj` or an error if `obj` does
#' not inherit from \linkS4class{gsvaParam}.
#'
#' @examples
#' gp <- gsvaParam()
#' get_mx.diff(gp)
#' gp <- gsvaParam(mx.diff = FALSE)
#' get_mx.diff(gp)
#'
#' @noRd
get_mx.diff <- function(obj) {
  stopifnot(inherits(obj, "gsvaParam"))
  return(obj@mx.diff)
}


#' Return the `abs.ranking` parameter from a `gsvaParam` object
#'
#' Returns the `abs.ranking` parameter from a \linkS4class{gsvaParam} object.
#'
#' @param obj A \linkS4class{gsvaParam} object.
#' 
#' @return The requested `abs.ranking` parameter from `obj` or an error if `obj`
#'   does not inherit from \linkS4class{gsvaParam}.
#'
#' @examples
#' gp <- gsvaParam()
#' get_abs.ranking(gp)
#' gp <- gsvaParam(abs.ranking = TRUE)
#' get_abs.ranking(gp)
#'
#' @noRd
get_abs.ranking <- function(obj) {
  stopifnot(inherits(obj, "gsvaParam"))
  return(obj@abs.ranking)
}

