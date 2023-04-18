
# Class Definitions -------------------------------------------------------


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
## #' \code{\link{ssGseaParam-class}}
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
## #' \code{\link{ssGseaParam-class}}
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
## #' \code{\link{ssGseaParam-class}}
## #' \code{\link{gsvaParam-class}}
#'
#' @name plageParam-class
#' @rdname plageParam-class
#' @exportClass plageParam
setClass("plageParam",
         slots = character(),
         contains = "emptyParam")


# Constructor Functions ---------------------------------------------------

# z-Score Parameter Constructor -------------------------------------------

#' Build a \code{zscore} parameter object
#'
#' Build an object of the class \code{zscoreParam}
#'
#' @return A \linkS4class{zscoreParam} object.
#'
#' @examples
#' zp <- zscoreParam()
#'
#' @importFrom methods new
#' @rdname zscoreParam
#' @export
zscoreParam <- function() {
  new("zscoreParam")
}


# PLAGE Parameter Constructor -------------------------------------------

#' Build a \code{plage} parameter object
#'
#' Build an object of the class \code{plageParam}
#'
#' @return A \linkS4class{plageParam} object.
#'
#' @examples
#' pp <- plageParam()
#'
#' @importFrom methods new
#' @rdname plageParam
#' @export
plageParam <- function() {
  new("plageParam")
}

