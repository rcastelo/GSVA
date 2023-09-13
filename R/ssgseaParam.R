
### ssGSEA Parameter Constructor --------------------------------------------

#' Construct a ssGSEA parameter object
#'
#' Construct and return a new object of class \linkS4class{ssgseaParam}.
#'
#' @param alpha Numeric vector of length 1; the exponent defining the
#'  weight of the tail in the random walk performed by the `ssGSEA` (Barbie et
#'  al., 2009) method.  The default value is 0.25 as described in the paper.
#' 
#' @param normalize Logical vector of length 1; if `TRUE`  runs the `ssGSEA` method
#'  from Barbie et al. (2009) normalizing the scores by the absolute difference
#'  between the minimum and the maximum, as described in their paper. Otherwise
#'  this last normalization step is skipped.
#' 
#' @return A new \linkS4class{ssgseaParam} object.
#'
#' @examples
#' sp <- ssgseaParam()
#'
#' @importFrom methods new
#' @rdname ssgseaParam
#' @export
ssgseaParam <- function(dataSet, geneSets, alpha = 0.25, normalize = TRUE) {
    new("ssgseaParam",
        dataSet = dataSet, geneSets = geneSets, alpha = alpha, normalize = normalize)
}


### ssGSEA Parameter Validator ---------------------------------------------

setValidity("ssgseaParam", function(object) {
    invalid <- character()
    dd <- dim(object@dataSet)
    if(dd[1] == 0) invalid <- c(invalid, "@dataSet has 0 rows")
    if(dd[2] == 0) invalid <- c(invalid, "@dataSet has 0 columns")
    if(length(object@geneSets) == 0) invalid <- c(invalid, "@geneSets has length 0")
    if(length(object@alpha) != 1) invalid <- c(invalid, "@alpha should be of length 1")
    if(is.na(object@alpha)) invalid <- c(invalid, "@alpha should not be NA")
    if(length(object@normalize) != 1) invalid <- c(invalid, "@normalize should be of length 1")
    if(is.na(object@normalize)) invalid <- c(invalid, "@normalize should not be NA")
    return(if(length(invalid) == 0) TRUE else invalid)
})


### Getters -----------------------------------------------------------------

#' Return the `alpha` parameter from a `ssgseaParam` object
#'
#' Returns the `alpha` parameter from a \linkS4class{ssgseaParam} object.
#'
#' @param obj A \linkS4class{ssgseaParam} object.
#' 
#' @return The requested `alpha` parameter from `obj` or an error if `obj` does
#' not inherit from \linkS4class{ssgseaParam}.
#'
#' @examples
#' sp <- ssgseaParam()
#' get_alpha(sp)
#' sp <- ssgseaParam(alpha = 0.42)
#' get_alpha(sp)
#'
#' @noRd
get_alpha <- function(obj) {
  stopifnot(inherits(obj, "ssgseaParam"))
  return(obj@alpha)
}


#' Return the `normalize` flag from a `ssgseaParam` object
#'
#' Returns the `normalize` flag from a \linkS4class{ssgseaParam} object.
#'
#' @param obj A \linkS4class{ssgseaParam} object.
#' 
#' @return The requested `normalize` flag from `obj` or an error if `obj` does
#' not inherit from \linkS4class{ssgseaParam}.
#'
#' @examples
#' sp <- ssgseaParam()
#' do_normalize(sp)
#' sp <- ssgseaParam(doNormalize = FALSE)
#' do_normalize(sp)
#'
#' @noRd
do_normalize <- function(obj) {
  stopifnot(inherits(obj, "ssgseaParam"))
  return(obj@normalize)
}


### ssGSEA Parameter Show -------------------------------------------------

setMethod("show",
          signature=signature(
            object="ssgseaParam"),
          function(object) {
              some <- function(x)
                  paste0(paste(Biobase::selectSome(x, 4), collapse=", "),
                        " (", length(x), " total)")
              ds <- get_dataSet(object)
              dsDim <- sprintf(" [%s, %d]", nrow(ds), ncol(ds))
              gs <- get_geneSets(object)
              gsDim <- sprintf(" [%d, %d]", nrow(gs), ncol(gs))
              gs
              cat("ssGSEA Parameter object\n",
                  "  data set: ", class(ds)[1], dsDim, "\n",
                  "    rows: ", some(rownames(ds)), "\n",
                  "      (annotation: ", annotation(ds), ")", "\n",
                  "    columns: ", some(colnames(ds)), "\n",
                  "  gene sets: ", class(gs)[1], gsDim, "\n",
                  "    names: ", some(names(gs)), "\n",
                  sep="")
          })
