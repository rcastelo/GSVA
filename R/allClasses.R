
# Virtual Superclass ------------------------------------------------------

#' `EmptyParam` class
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
#' [`ZScoreParam-class`], 
#' [`PlageParam-class`], 
#' [`SsGseaParam-class`], 
#' [`GsvaParam-class`]
#'
#' @name EmptyParam-class
#' @rdname EmptyParam-class
#' @exportClass EmptyParam
setClass("EmptyParam",
         slots = character(),
         contains = "VIRTUAL")

setClassUnion("GsvaDataSet",
              c("matrix", "ExpressionSet", "SummarizedExperiment",
                "SingleCellExperiment", "dgCMatrix", "HDF5Array"))

setClassUnion("GsvaGeneSet",
              c("list", "GeneSetCollection"))


# z-Score Parameter Class -------------------------------------------------

#' `ZScoreParam` class
#'
#' Method-specific parameters for the `zscore` method.
#' 
#' Since this method does not take any parameters, the parameter class does not
#' have any slots and exists merely for method dispatch.  It is derived from the
#' virtual superclass [`EmptyParam-class`].
#'
#' @seealso
#' [`PlageParam-class`],
#' [`SsGseaParam-class`],
#' [`GsvaParam-class`]
#'
#' @name ZScoreParam-class
#' @rdname ZScoreParam-class
#' @exportClass ZScoreParam
setClass("ZScoreParam",
         slots = character(),
         contains = "EmptyParam")



# PLAGE Parameter Class -------------------------------------------------

#' `PlageParam` class
#'
#' Method-specific parameters for the `plage` method.
#' 
#' Since this method does not take any parameters, the parameter class does not
#' have any slots and exists merely for method dispatch.  It is derived from the
#' virtual superclass [`EmptyParam-class`].
#'
#' @seealso
#' [`ZScoreParam-class`],
#' [`SsGseaParam-class`],
#' [`GsvaParam-class`]
#'
#' @name PlageParam-class
#' @rdname PlageParam-class
#' @exportClass PlageParam
setClass("PlageParam",
         slots = character(),
         contains = "EmptyParam")


# ssGSEA Parameter Class -------------------------------------------------

#' `SsGseaParam` class
#'
#' Method-specific parameters for the `ssgsea` method.
#' 
#' This class has slots for storing the two parameters to the `ssgsea` method
#' described below.  It is derived from the virtual superclass [`EmptyParam-class`].
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
#' [`ZScoreParam-class`],
#' [`PlageParam-class`],
#' [`GsvaParam-class`]
#'
#' @name SsGseaParam-class
#' @rdname SsGseaParam-class
#' @exportClass SsGseaParam
setClass("SsGseaParam",
         slots = c(alpha = "numeric", normalize = "logical"),
         contains = "EmptyParam")


# GSVA Parameter Class ----------------------------------------------------

#' `GsvaParam` class
#'
#' Method-specific parameters for the `gsva` method.
#'
#' This class has slots for storing the two parameters to the `gsva` method
#' described below.  It is derived from the virtual superclass
#' [`EmptyParam-class`].
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
#' [`ZScoreParam-class`],
#' [`PlageParam-class`],
#' [`SsGseaParam-class`]
#'
#' @name GsvaParam-class
#' @rdname GsvaParam-class
#' @exportClass GsvaParam
setClass("GsvaParam",
         slots = c(kcdf = "character", tau = "numeric", 
                   mx.diff = "logical", abs.ranking = "logical"),
         contains = "EmptyParam")



