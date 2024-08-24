
#' @title The `gsvaParam` class
#'
#' @description Objects of class `gsvaParam` contain the parameters for running
#' the `GSVA` method.
#'
#' @details In addition to an expression data set and a collection of
#' gene sets, `GSVA` takes four method-specific parameters as described below.
#'
#' @param exprData The expression data.  Must be one of the classes
#' supported by [`GsvaExprData-class`]. Type `help(GsvaExprData)` to consult
#' the available classes.
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].
#'
#' @param assay The name of the assay to use in case `exprData` is a multi-assay
#' container, otherwise ignored.  By default, the first assay is used.
#' 
#' @param annotation The name of a Bioconductor annotation package for the gene
#' identifiers occurring in the row names of the expression data matrix.  This
#' can be used to map gene identifiers occurring in the gene sets if those are
#' provided in a [`GeneSetCollection`].  By default gene identifiers used in
#' expression data matrix and gene sets are matched directly.
#' 
#' @param minSize Minimum size of the resulting gene sets after gene identifier
#' mapping. By default, the minimum size is 1.
#' 
#' @param maxSize Maximum size of the resulting gene sets after gene identifier
#' mapping. By default, the maximum size is `Inf`.
#' 
#' @param kcdf Character vector of length 1 denoting the kernel to use during
#' the non-parametric estimation of the empirical cumulative distribution
#' function (ECDF) of expression levels across samples. By default,
#' `kcdf="auto"`, which will make GSVA to decide automatically any of the other
#' possible values. The value `kcdf="Gaussian"` is suitable when input
#' expression values are continuous, such as microarray fluorescent units in
#' logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input
#' expression values are integer counts, such as those derived from RNA-seq
#' experiments, then this argument should be set to `kcdf="Poisson"`. When we
#' do not want to use a kernel approach for the estimation of the ECDF, then
#' we should set `kcdf="none"`.
#'
#' @param kcdfNoneMinSampleSize Integer vector of length 1. When `kcdf="auto"`,
#' this parameter decides at what minimum sample size `kcdf="none"`, i.e., the
#' estimation of the empirical cumulative distribution function (ECDF) of
#' expression levels across samples is performed directly without using a
#' kernel; see the `kcdf` slot.
#'
#' @param tau Numeric vector of length 1.  The exponent defining the weight of
#' the tail in the random walk performed by the `GSVA` (Hänzelmann et al.,
#' 2013) method.  The default value is 1 as described in the paper.
#'
#' @param maxDiff Logical vector of length 1 which offers two approaches to
#' calculate the enrichment statistic (ES) from the KS random walk statistic.
#' * `FALSE`: ES is calculated as the maximum distance of the random walk
#' from 0. This approach produces a distribution of enrichment scores that is
#' bimodal, but it can give large enrichment scores to gene sets whose genes
#' are not concordantly activated in one direction only.
#' * `TRUE` (the default): ES is calculated as the magnitude difference between
#' the largest positive and negative random walk deviations. This default value
#' gives larger enrichment scores to gene sets whose genes are concordantly
#' activated in one direction only.
#'
#' @param absRanking Logical vector of length 1 used only when `maxDiff=TRUE`.
#' When `absRanking=FALSE` (default) a modified Kuiper statistic is used to
#' calculate enrichment scores, taking the magnitude difference between the
#' largest positive and negative random walk deviations. When
#' `absRanking=TRUE` the original Kuiper statistic that sums the largest
#' positive and negative random walk deviations is used.
#' 
#' @param sparse Logical vector of length 1 used only when the input expression
#' data in `exprData` is stored in a sparse matrix (e.g., a `dgCMatrix` or a
#' `singleCellExperiment` object storing the expression data in a `dgCMatrix`).
#' In such a case, when `sparse=TRUE`, a sparse version of the GSVA algorithm
#' will be applied. Otherwise, when `sparse=FALSE`, the classical version of
#' the GSVA algorithm will be used.
#'
#' @return A new [`gsvaParam-class`] object.
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' [DOI](https://doi.org/10.1186/1471-2105-14-7)
#'
#' @examples
#' library(GSVA)
#' library(GSVAdata)
#'
#' data(leukemia)
#' data(c2BroadSets)
#' 
#' ## for simplicity, use only a subset of the sample data
#' ses <- leukemia_eset[1:1000, ]
#' gsc <- c2BroadSets[1:100]
#' gp1 <- gsvaParam(ses, gsc)
#' gp1
#'
#'
#' @importFrom methods new
#' @rdname gsvaParam-class
#' 
#' @export
gsvaParam <- function(exprData, geneSets,
                      assay=NA_character_, annotation=NA_character_,
                      minSize=1,maxSize=Inf,
                      kcdf=c("auto", "Gaussian", "Poisson", "none"),
                      kcdfNoneMinSampleSize=50, tau=1, maxDiff=TRUE,
                      absRanking=FALSE, sparse=FALSE) {
    kcdf <- match.arg(kcdf)
    kcdfNoneMinSampleSize=as.integer(kcdfNoneMinSampleSize)

    an <- gsvaAssayNames(exprData)
    if((!is.na(assay)) && (!.isCharNonEmpty(an)))
        warning("argument assay='", assay,
                "' ignored since exprData has no assayNames()")
    if(is.na(assay) && .isCharNonEmpty(an))
        assay <- na.omit(an)[1]
    
    new("gsvaParam",
        exprData=exprData, geneSets=geneSets,
        assay=assay, annotation=annotation,
        minSize=minSize, maxSize=maxSize,
        kcdf=kcdf, kcdfNoneMinSampleSize=kcdfNoneMinSampleSize,
        tau=tau, maxDiff=maxDiff, absRanking=absRanking, sparse=sparse)
}


## ----- validator -----

setValidity("gsvaParam", function(object) {
    inv <- NULL
    xd <- object@exprData
    dd <- dim(xd)
    an <- gsvaAssayNames(xd)
    oa <- object@assay
    
    if(dd[1] == 0) {
        inv <- c(inv, "@exprData has 0 rows")
    }
    if(dd[2] == 0) {
        inv <- c(inv, "@exprData has 0 columns")
    }
    if(length(object@geneSets) == 0) {
        inv <- c(inv, "@geneSets has length 0")
    }
    if(length(oa) != 1) {
        inv <- c(inv, "@assay should be of length 1")
    }
    if(.isCharLength1(oa) && .isCharNonEmpty(an) && (!(oa %in% an))) {
        inv <- c(inv, "@assay should be one of assayNames(@exprData)")
    }
    if(length(object@annotation) != 1) {
        inv <- c(inv, "@annotation should be of length 1")
    }
    if(length(object@minSize) != 1) {
        inv <- c(inv, "@minSize should be of length 1")
    }
    if(object@minSize < 1) {
        inv <- c(inv, "@minSize should be at least 1 or greater")
    }
    if(length(object@maxSize) != 1) {
        inv <- c(inv, "@maxSize should be of length 1")
    }
    if(object@maxSize < object@minSize) {
        inv <- c(inv, "@maxSize should be at least @minSize or greater")
    }
    if(length(object@kcdfNoneMinSampleSize) != 1) {
        inv <- c(inv, "@kcdfNoneMinSampleSize should be of length 1")
    }
    if(object@kcdfNoneMinSampleSize < 0) {
        inv <- c(inv, "@kcdfNoneMinSampleSize should be a non-negative integer")
    }
    if(is.na(object@kcdfNoneMinSampleSize)) {
        inv <- c(inv, "@kcdfNoneMinSampleSize should not be NA")
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
    if(length(object@sparse) != 1) {
        inv <- c(inv, "@sparse should be of length 1")
    }
    if(is.na(object@sparse)) {
        inv <- c(inv, "@sparse should not be NA")
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
get_kcdfNoneMinSampleSize <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@kcdfNoneMinSampleSize)
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

#' @noRd
get_sparse <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@sparse)
}


## ----- show -----

setMethod("show",
          signature=signature(object="gsvaParam"),
          function(object) {
              callNextMethod(object)
              cat("kcdf: ", get_kcdf(object), "\n",
                  "kcdfNoneMinSampleSize: ", get_kcdfNoneMinSampleSize(object), "\n",
                  "tau: ", get_tau(object), "\n",
                  "maxDiff: ", get_maxDiff(object), "\n",
                  "absRanking: ", get_absRanking(object), "\n",
                  "sparse: ", get_sparse(object), "\n",
                  sep="")
          })
