
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
## #' 
## #' @examples
## #' bamfiles <- list.files(system.file("extdata", package="atena"),
## #'                        pattern="*.bam", full.names=TRUE)
## #' TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
## #'                     package="atena"))
## #' ttpar <- TEtranscriptsParam(bamfiles, teFeatures=TE_annot, singleEnd=TRUE,
## #'                             ignoreStrand=TRUE, aggregateby = c("repName"))
## #' path(ttpar)
#'
#' @name zscoreParam-class
#' @rdname zscoreParam-class
#' @exportClass zscoreParam
setClass("zscoreParam",
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

