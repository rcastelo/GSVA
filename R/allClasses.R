
# Virtual Superclass ------------------------------------------------------

#' `emptyParam` class
#'
#' Virtual superclass of method-specific parameter classes.
#'
#' To ensure that the various methods implemented in `GSVA` are invoked with
#' appropriate parameters, S4 classes are introduced as containers for
#' method-specific parameters as well as for dispatching them to the
#' appropriate method -- which requires them to be instantiated.
#' 
#' Some simple methods, however, do not take any method-specific parameters at
#' all and hence their parameter classes have no slots and will be created as
#' virtual classes unless they have a superclass (see [`methods::setClass`]).
#' 
#' This is the purpose of this class.
#'
#' @seealso
#' [`zscoreParam-class`], 
#' [`plageParam-class`], 
#' [`ssgseaParam-class`], 
#' [`gsvaParam-class`]
#'
#' @name emptyParam-class
#' @rdname emptyParam-class
#' @exportClass emptyParam
setClass("emptyParam",
         slots = character(),
         contains = "VIRTUAL")


# z-Score Parameter Class -------------------------------------------------

#' `zscoreParam` class
#'
#' Method-specific parameters for the `zscore` method.
#' 
#' Since this method does not take any parameters, the parameter class does not
#' have any slots and exists merely for method dispatch.  It is derived from the
#' virtual superclass [`emptyParam-class`].
#'
#' @seealso
#' [`plageParam-class`],
#' [`ssgseaParam-class`],
#' [`gsvaParam-class`]
#'
#' @name zscoreParam-class
#' @rdname zscoreParam-class
#' @exportClass zscoreParam
setClass("zscoreParam",
         slots = character(),
         contains = "emptyParam")



# PLAGE Parameter Class -------------------------------------------------

#' `plageParam` class
#'
#' Method-specific parameters for the `plage` method.
#' 
#' Since this method does not take any parameters, the parameter class does not
#' have any slots and exists merely for method dispatch.  It is derived from the
#' virtual superclass [`emptyParam-class`].
#'
#' @seealso
#' [`zscoreParam-class`],
#' [`ssgseaParam-class`],
#' [`gsvaParam-class`]
#'
#' @name plageParam-class
#' @rdname plageParam-class
#' @exportClass plageParam
setClass("plageParam",
         slots = character(),
         contains = "emptyParam")


# ssGSEA Parameter Class -------------------------------------------------

#' `ssgseaParam` class
#'
#' Method-specific parameters for the `ssgsea` method.
#' 
#' This class has slots for storing the two parameters to the `ssgsea` method
#' described below.  It is derived from the virtual superclass [`emptyParam-class`].
#' 
#' @slot alpha Numeric vector of length 1; the exponent defining the
#'  weight of the tail in the random walk performed by the `ssGSEA` (Barbie et
#'  al., 2009) method.
#'
#' @slot normalize Logical vector of length 1; if `TRUE`  runs the `ssGSEA` method
#'  from Barbie et al. (2009) normalizing the scores by the absolute difference
#'  between the minimum and the maximum, as described in their paper. Otherwise
#'  this last normalization step is skipped.
#'
#' @seealso
#' [`zscoreParam-class`],
#' [`plageParam-class`],
#' [`gsvaParam-class`]
#'
#' @name ssgseaParam-class
#' @rdname ssgseaParam-class
#' @exportClass ssgseaParam
setClass("ssgseaParam",
         slots = c(alpha = "numeric", normalize = "logical"),
         contains = "emptyParam")


# GSVA Parameter Class ----------------------------------------------------

#' `gsvaParam` class
#'
#' Method-specific parameters for the `gsva` method.
#'
#' This class has slots for storing the two parameters to the `gsva` method
#' described below.  It is derived from the virtual superclass
#' [`emptyParam-class`].
#'
#' @slot kcdf Character vector of length 1 denoting the kernel to use during the
#'   non-parametric estimation of the cumulative distribution function of
#'   expression levels across samples. `kcdf="Gaussian"` is suitable when input
#'   expression values are continuous, such as microarray fluorescent units in
#'   logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input
#'   expression values are integer counts, such as those derived from RNA-seq
#'   experiments, then this argument should be set to `kcdf="Poisson"`.
#'
#' @slot tau Numeric vector of length 1; the exponent defining the weight of the
#'   tail in the random walk performed by the `GSVA` (Hänzelmann et al., 2013)
#'   method.
#'
#' @slot mx.diff Logical vector of length 1 which offers two approaches to
#'   calculate the enrichment statistic (ES) from the KS random walk statistic.
#'  * `FALSE`: ES is calculated as the maximum distance of the random walk
#'   from 0.
#'  * `TRUE`: ES is calculated as the magnitude difference between
#'   the largest positive and negative random walk deviations.
#'
#' @slot abs.ranking Logical vector of length 1 used only when `mx.diff=TRUE`.
#'   When `abs.ranking=FALSE` a modified Kuiper statistic is used to calculate
#'   enrichment scores, taking the magnitude difference between the largest
#'   positive and negative random walk deviations. When `abs.ranking=TRUE` the
#'   original Kuiper statistic that sums the largest positive and negative
#'   random walk deviations, is used. In this latter case, gene sets with genes
#'   enriched on either extreme (high or low) will be regarded as ’highly’
#'   activated.
#'
#' @seealso
#' [`zscoreParam-class`],
#' [`plageParam-class`],
#' [`ssgseaParam-class`]
#'
#' @name gsvaParam-class
#' @rdname gsvaParam-class
#' @exportClass gsvaParam
setClass("gsvaParam",
         slots = c(kcdf = "character", tau = "numeric", 
                   mx.diff = "logical", abs.ranking = "logical"),
         contains = "emptyParam")



