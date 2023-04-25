
#
# 2023-04-25  axel: S4 class definitions
#

# Virtual Superclass ------------------------------------------------------

#' empty parameter class
#'
#' Virtual empty superclass for deriving non-virtual empty method-specific
#' parameter classes for methods without method-specific parameters.  These
#' classes exist merely for method dispatch.
#'
## #' @seealso
## #' \code{\link{zscoreParam-class}}
## #' \code{\link{plageParam-class}}
## #' \code{\link{ssgseaParam-class}}
## #' \code{\link{gsvaParam-class}}
#'
#' @name emptyParam-class
#' @rdname emptyParam-class
#' @exportClass emptyParam
setClass("emptyParam",
         slots = character(),
         contains = "VIRTUAL")


# z-Score Parameter Class -------------------------------------------------

#' zscore parameter class
#'
#' Method-specific parameters for the z-score method.  Since this method does not
#' need any parameters, the class does not have any slots and exists merely for
#' method dispatch.
#'
## #' @seealso
## #' \code{\link{plageParam-class}}
## #' \code{\link{ssgseaParam-class}}
## #' \code{\link{gsvaParam-class}}
#'
#' @name zscoreParam-class
#' @rdname zscoreParam-class
#' @exportClass zscoreParam
setClass("zscoreParam",
         slots = character(),
         contains = "emptyParam")



# PLAGE Parameter Class -------------------------------------------------

#' PLAGE parameter class
#'
#' Method-specific parameters for the PLAGE method.  Since this method does not
#' need any parameters, the class does not have any slots and exists merely for
#' method dispatch.
#'
## #' @seealso
## #' \code{\link{zscoreParam-class}}
## #' \code{\link{ssgseaParam-class}}
## #' \code{\link{gsvaParam-class}}
#'
#' @name plageParam-class
#' @rdname plageParam-class
#' @exportClass plageParam
setClass("plageParam",
         slots = character(),
         contains = "emptyParam")


# ssGSEA Parameter Class -------------------------------------------------

#' ssGSEA parameter class
#'
#' Method-specific parameters for the ssGSEA method.  
#'
## #' @seealso
## #' \code{\link{zscoreParam-class}}
## #' \code{\link{plageParam-class}}
## #' \code{\link{gsvaParam-class}}
#'
#' @name ssgseaParam-class
#' @rdname ssgseaParam-class
#' @exportClass ssgseaParam
setClass("ssgseaParam",
         slots = c(alpha = "numeric", normalize = "logical"),
         contains = "emptyParam")

