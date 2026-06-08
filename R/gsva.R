#' @title Gene Set Variation Analysis
#' 
#' @description Estimates GSVA enrichment scores.
#' 
#' @param param A parameter object of one of the following classes:
#' * A [`gsvaParam`] object built using the constructor function
#' [`gsvaParam`].
#'   This object will trigger `gsva()` to use the GSVA algorithm by
#'   Hänzelmann et al. (2013).
#' * A [`plageParam`] object built using the constructor function
#' [`plageParam`].
#'   This object will trigger `gsva()` to use the PLAGE algorithm by
#'   Tomfohr et al. (2005).
#' * A [`zscoreParam`] object built using the constructor function
#' [`zscoreParam`].
#'   This object will trigger `gsva()` to use the combined z-score algorithm by
#'   Lee et al. (2008).
#' * A [`ssgseaParam`] object built using the constructor function
#' [`ssgseaParam`].
#'   This object will trigger `gsva()` to use the ssGSEA algorithm by
#'   Barbie et al. (2009).
#'
#' @param verbose Gives information about each calculation step. Default: `TRUE`.
#' 
#' @param BPPARAM An object of class `BiocParallelParam` specifying parameters
#'   related to the parallel execution of some of the tasks and calculations
#'   within this function.
#' 
#' @param maxmem A vector of length 1 either specifying a number in bytes, or
#' a character string with either the word `auto` (default), or a number
#' followed by a suffix indicating kilobytes (K), megabytes (M), gigabytes (G)
#' or terabytes (T), which GSVA will use to attempt bounding the maximum amount
#' of main memory used across all threads of execution to that given quantity.
#' By default `maxmem="auto"`, indicating that the maximum memory will be the
#' 90% of the total main memory, as calculated by
#' [`Sys.meminfo()`][memuse::Sys.meminfo]. To avoid setting any bound on the
#' maximum memory, use `maxmem=Inf`. Note that the amount of main memory used
#' in an R session or script may depend on other commands and packages used in
#' that same session or script.
#'
#' @return A gene-set by sample matrix of GSVA enrichment scores stored in a
#' container object of the same type as the input expression data container,
#' except for the fact that enrichment scores are always dense, irrespective of
#' whether the input is sparse, such as in single-cell data. If the input was a
#' base matrix, a [`dgCMatrix`][Matrix::dgCMatrix-class], a
#' [`SVT_SparseMatrix`][SparseArray::SVT_SparseMatrix-class], or a 
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix-class] object, then the output
#' will be either a base matrix object or a
#' [`DelayedMatrix`][DelayedArray::DelayedMatrix-class], with the gene sets
#' employed in the calculations stored in an attribute called `geneSets` of that
#' object. If the input was an `ExpressionSet` object, then the output will be
#' also an `ExpressionSet` object with the gene sets employed in the
#' calculations stored in an attribute called `geneSets`. If the input was an
#' object of either class
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment-class],
#' [`SingleCellExperiment`][SingleCellExperiment::SingleCellExperiment-class],
#' or [`SpatialExperiment`][SpatialExperiment::SpatialExperiment-class],
#' then the output will be of the same class, where enrichment scores will be
#' stored in an assay called `es` and the gene sets employed in the
#' calculations will be stored in the `rowData` slot of the object under the
#' column name `gs`.
#' 
#' @seealso [`gsvaParam`], [`plageParam`], [`zscoreParam`], [`ssgseaParam`],
#' [`BiocParallelParam`][BiocParallel::BiocParallelParam-class],
#' [`gsvaRowNorm`], [`gsvaColRanks`], [`gsvaColScores`]
#'
#' @aliases gsva
#' @name gsva
#' @rdname gsva
#' 
#' @references Barbie, D.A. et al. Systematic RNA interference reveals that
#' oncogenic KRAS-driven cancers require TBK1.
#' *Nature*, 462(5):108-112, 2009.
#' \doi{10.1038/nature08460}
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' \doi{10.1186/1471-2105-14-7}
#'
#' @references Lee, E. et al. Inferring pathway activity toward precise
#' disease classification.
#' *PLoS Comp Biol*, 4(11):e1000217, 2008.
#' \doi{10.1371/journal.pcbi.1000217}
#'
#' @references Tomfohr, J. et al. Pathway level analysis of gene expression
#' using singular value decomposition.
#' *BMC Bioinformatics*, 6:225, 2005.
#' \doi{10.1186/1471-2105-6-225}
#'
#' @examples
#' library(GSVA)
#' library(limma)
#' 
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#' nGrp1 <- 15 ## number of samples in group 1
#' nGrp2 <- n - nGrp1 ## number of samples in group 2
#' 
#' ## consider three disjoint gene sets
#' geneSets <- list(set1=paste("g", 1:3, sep=""),
#'                  set2=paste("g", 4:6, sep=""),
#'                  set3=paste("g", 7:10, sep=""))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
#' y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2
#' 
#' ## build design matrix
#' design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nGrp1), rep(1, nGrp2)))
#' 
#' ## fit linear model
#' fit <- lmFit(y, design)
#' 
#' ## estimate moderated t-statistics
#' fit <- eBayes(fit)
#' 
#' ## genes in set1 are differentially expressed
#' topTable(fit, coef="sampleGroup2vs1")
#' 
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#' 
#' ## estimate GSVA enrichment scores for the three sets
#' gsva_es <- gsva(gsvapar)
#' 
#' ## fit the same linear model now to the GSVA enrichment scores
#' fit <- lmFit(gsva_es, design)
#' 
#' ## estimate moderated t-statistics
#' fit <- eBayes(fit)
#' 
#' ## set1 is differentially expressed
#' topTable(fit, coef="sampleGroup2vs1")
NULL

#' @aliases gsva,gsvaParam-method
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {

              if (verbose) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
                  gsva_global$show_start_and_end_messages <- FALSE
              }

              .check_bpparam(BPPARAM)

              gsvarownr <- gsvaRowNorm(param=param, verbose=verbose,
                                       dropExistingAssays=TRUE,
                                       BPPARAM=BPPARAM, maxmem=maxmem)

              gsvaranks <- gsvaColRanks(rowNormExprData=gsvarownr,
                                        verbose=verbose,
                                        dropExistingAssays=TRUE,
                                        BPPARAM=BPPARAM,
                                        maxmem=maxmem)

              es <- gsvaColScores(rankExprData=gsvaranks, verbose=verbose,
                                  BPPARAM=BPPARAM, maxmem=maxmem)

              if (verbose) {
                  cli_alert_success("Calculations finished")
                  gsva_global$show_start_and_end_messages <- TRUE
              }
              
              return(es)
          })


#' @title The `gsvaParam` class
#'
#' @description Objects of class `gsvaParam` contain the parameters for running
#' the `GSVA` method.
#'
#' @details In addition to a number of parameters shared with all methods
#' implemented by package GSVA, `GSVA` takes six method-specific parameters.
#' All of these parameters are described in detail below.
#'
#' @param exprData The expression data set.  Must be one of the classes
#' supported by [`GsvaExprData-class`].  For a list of these classes, see its
#' help page using `help(GsvaExprData)`.
#'
#' @param geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].  For a list of these classes, see its help page using
#' `help(GsvaGeneSets)`.
#' 
#' @param assay Character vector of length 1.  The name of the assay to use in
#' case `exprData` is a multi-assay container, otherwise ignored.  By default,
#' an assay called 'logcounts' will be used if present, otherwise the first
#' assay is used.
#' 
#' @param annotation An object of class `GeneIdentifierType` from
#' package `GSEABase` describing the gene identifiers used as the row names of
#' the expression data set.  See `GeneIdentifierType` for help on available
#' gene identifier types and how to construct them.  This
#' information can be used to map gene identifiers occurring in the gene sets.
#' 
#' If the default value `NULL` is provided, an attempt will be made to extract
#' the gene identifier type from the expression data set provided as `exprData`
#' (by calling [`gsvaAnnotation`] on it).  If still not successful, the
#' `NullIdentifier()` will be used as the gene identifier type, gene identifier
#' mapping will be disabled and gene identifiers used in expression data set and
#' gene sets can only be matched directly.
#' 
#' @param minSize Numeric vector of length 1.  Minimum size of the resulting gene
#' sets after gene identifier mapping. By default, the minimum size is 1.
#' 
#' @param maxSize Numeric vector of length 1.  Maximum size of the resulting gene
#' sets after gene identifier mapping. By default, the maximum size is `Inf`.
#' 
#' @param kcdf Character vector of length 1 denoting the kernel to use during
#' the non-parametric estimation of the empirical cumulative distribution
#' function (ECDF) of expression levels across samples. The value `kcdf="auto"`
#' will allow GSVA to automatically choose one of the possible values. The
#' value `kcdf="Gaussian"` is suitable when input expression values are
#' continuous, such as microarray fluorescent units in logarithmic scale,
#' RNA-seq log-CPMs, log-RPKMs, or log-TPMs. When input expression values are
#' integer counts, such as those derived from RNA-seq experiments, then this
#' argument should be set to `kcdf="Poisson"`. When we do not want to use a
#' kernel approach for the estimation of the ECDF, then we should set
#' `kcdf="none"`.
#'
#' @param kcdfNoneMinSampleSize Integer vector of length 1. When `kcdf="auto"`,
#' this parameter decides at what minimum sample size `kcdf="none"`, i.e., the
#' estimation of the empirical cumulative distribution function (ECDF) of
#' expression levels across samples is performed directly without using a
#' kernel. By default, this value is set to 200; see the `kcdf` slot.
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
#' `SingleCellExperiment` object storing the expression data in a `dgCMatrix`).
#' In such a case, when `sparse=TRUE` (default), a sparse version of the GSVA
#' algorithm will be applied. Otherwise, when `sparse=FALSE`, the classical
#' version of the GSVA algorithm will be used.
#'
#' @param checkNA Character vector of length 1 specifying whether the input
#' expression data should be checked for the presence of missing values (`NA`
#' or `NaN`). This must be one of the strings `"auto"` (default), `"yes"`, or
#' `"no"`. The default value `"auto"` means that the software will perform that
#' check only when the input expression data is provided as a base `matrix`, an
#' `ExpressionSet` or a `SummarizedExperiment` object, while every other type
#' of input expression data container (e.g., `SingleCellExperiment`, etc.) will
#' not be checked. If `checkNA="yes"`, then the input expression data will be
#' checked for missing values irrespective of the object class of the data
#' container, and if `checkNA="no"`, then that check will not be performed.
#'
#' @param use Character vector of length 1 specifying a policy for dealing with
#' missing values (`NA` or `NaN`) in the input expression data argument
#' `exprData`. It only applies when either `checkNA="yes"`, or `checkNA="auto"`
#' (see the `checkNA` parameter. The argument value must be one of the strings
#' `"everything"` (default), `"all.obs"`, or `"na.rm"`. The policy of the
#' default value `"everything"` consists of propagating missing values so that
#' the resulting enrichment score will be `NA`, whenever one or more of its
#' contributing values is missing, giving a warning when that happens. When
#' `use="all.obs"`, the presence of `NA`s in the input expression data will
#' produce an error. Finally, when `use="na.rm"`, missing values in the input
#' expression data will be removed from calculations, giving a warning when that
#' happens, and giving an error if no values are left after removing the missing
#' values.
#'
#' @param filterRows Logical vector of length 1, indicating whether the rows in,
#' the input expression data, typically corresponding to transcripts, genes or
#' proteins, should be filtered for constant expression across columns,
#' typically corresponding to samples or cells, with respect to all available
#' (nonmissing) values and to the non-zero values. By default, this slot is set
#' to `TRUE` and the user may set it to `FALSE` when there is absolute certainty
#' that no such rows exist in the input expression data, since this may save
#' running time, especially with data sets with hundreds of thousands or
#' millions of columns.
#'
#' @param ondisk Character vector of length 1 denoting whether an on-disk backend
#' should be used to reduce the memory footprint. The default value
#' `ondisk="auto"` will attempt to load all the data in main memory when the
#' input nonzero values fit in main memory, otherwise it will attempt working
#' with an on-disk data structure that reduces de memory footprint. When
#' `ondisk="yes"` it will attempt to work with an on-disk data structure, while
#' when `ondisk="no"` it will attempt to load all the data in main memory.
#'
#' @param verbose Logical vector of length 1. It gives information about some
#' decisions made by the software during parameter object construction when
#' `verbose=TRUE` (default) and remains silent otherwise.
#'
#' @return A new [`gsvaParam-class`] object.
#'
#' @seealso [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class],
#' [`matrix`],
#' \code{\link[Biobase]{ExpressionSet}},
### we are using the plain Rd above because
###  #' [`ExpressionSet`][Biobase::ExpressionSet-class],
### results in the following R CMD check NOTE:
### Non-topic package-anchored link(s) in Rd file 'gsvaParam-class.Rd':
###  ‘[Biobase:class.ExpressionSet]{ExpressionSet}’
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment-class],
#' [`SingleCellExperiment`][SingleCellExperiment::SingleCellExperiment-class],
#' [`SpatialExperiment`][SpatialExperiment::SpatialExperiment-class]
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' \doi{10.1186/1471-2105-14-7}
#'
#' @examples
#' suppressPackageStartupMessages({
#' library(GSEABase)
#' library(GSVA)
#' library(GSVAdata)
#' })
#'
#' data(geneprotExpCostaEtAl2021)
#' data(c2BroadSets)
#' 
#' ## for simplicity, use only a subset of the sample data
#' se <- geneExpCostaEtAl2021[1:1000, ]
#' gsc <- c2BroadSets[1:100]
#' gp1 <- gsvaParam(se, gsc)
#' gp1
#'
#'
#' @importFrom methods new
#' @importFrom cli cli_abort cli_alert_warning
#' @importFrom utils capture.output
#' @rdname gsvaParam-class
#' 
#' @export
gsvaParam <- function(exprData, geneSets,
                      assay=NA_character_, annotation=NULL,
                      minSize=1, maxSize=Inf,
                      kcdf=c("auto", "Gaussian", "Poisson", "none"),
                      kcdfNoneMinSampleSize=200, tau=1, maxDiff=TRUE,
                      absRanking=FALSE, sparse=TRUE,
                      checkNA=c("auto", "yes", "no"),
                      use=c("everything", "all.obs", "na.rm"),
                      filterRows=TRUE,
                      ondisk=c("auto", "yes", "no"),
                      verbose=TRUE) {

    .check_input_expr_gene_sets(exprData, geneSets)

    kcdf <- match.arg(kcdf)
    kcdfNoneMinSampleSize <- as.integer(kcdfNoneMinSampleSize)
    checkNA <- match.arg(checkNA)
    use <- match.arg(use)
    ondisk <- match.arg(ondisk)

    ## check assay parameter and assay names
    assay <- .check_assayNames(assay, exprData, verbose)

    ## check for presence of valid row/feature names
    exprData <- .check_rowNames(expr=exprData, useDummyNames=TRUE,
                                verbose=verbose)

    xa <- gsvaAnnotation(exprData)
    if(is.null(xa)) {
        if(is.null(annotation)) {
            annotation <- NullIdentifier()
        }
    } else {
        if(is.null(annotation)) {
            annotation <- xa
        } else {
            msg <- sprintf(paste0("using argument annotation='%s' and ",
                                  "ignoring exprData annotation ('%s')"),
                           capture.output(annotation), capture.output(xa))
            cli_alert_info(msg)
        }
    }

    naparam <- .check_for_na_values(exprData=exprData, assay=assay,
                                    checkNA=checkNA, use=use)

    nzc <- .estimate_nzcount(exprData, assay, verbose)

    if (!is_sparse(unwrapData(exprData, assay))) ## use sparse regime only
        sparse <- FALSE                          ## when input is sparse

    if (!filterRows) {
        cli_alert_warning("filterRows=FALSE and rows with constant values will not be filtered out")
        cli_alert_warning("Use it only if you are sure that such rows are not present in the input data")
    }

    param <- new("gsvaParam",
                 exprData=exprData, geneSets=geneSets,
                 assay=assay, annotation=annotation,
                 minSize=minSize, maxSize=maxSize,
                 kcdf=kcdf, kcdfNoneMinSampleSize=kcdfNoneMinSampleSize,
                 tau=as.double(tau), maxDiff=maxDiff, absRanking=absRanking,
                 sparse=sparse, checkNA=checkNA, didCheckNA=naparam$didCheckNA,
                 anyNA=naparam$any_na, use=use, filterRows=filterRows,
                 nzcount=nzc, ondisk=ondisk)

    maxmem <- .check_maxmem(param, maxmem="auto", verbose=verbose)
    .check_ondisk(param, maxmem=maxmem, first=NA, last=NA, whdim=1,
                  verbose=verbose)

    return(param)
}

## ----- setters for gsvaParam -----

#' @param object For the replacement method, an object of class
#' [`gsvaParam-class`].
#'
#' @param value For the replacement method, an object of the classes supported
#' by [`GsvaGeneSets-class`].
#'
#' @aliases geneSets<-
#' @aliases geneSets<-,gsvaParam,GsvaGeneSets-method
#' @rdname gsvaParam-class
#' @exportMethod geneSets
setReplaceMethod("geneSets", signature=signature(object="gsvaParam",
                                                 value="GsvaGeneSets"),
                 function(object, value) {
                   object@geneSets <- value
                   object
                 })

## ----- validator -----

setValidity("gsvaParam", function(object) {
    inv <- NULL
    xd <- object@exprData
    dd <- dim(xd)
    ## an <- gsvaAssayNames(xd)
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
        inv <- c(inv, "@assay must be of length 1")
    }
    ## this is incompatible with using dropExistingAssays=TRUE
    ## if(.isCharLength1(oa) && .isCharNonEmpty(an) && (!(oa %in% an))) {
    ##     inv <- c(inv, "@assay must be one of assayNames(@exprData)")
    ## }
    if(length(object@annotation) != 1) {
        inv <- c(inv, "@annotation must be of length 1")
    }
    if(!inherits(object@annotation, "GeneIdentifierType")) {
        inv <- c(inv, "@annotation must be a subclass of 'GeneIdentifierType'")
    }
    if(length(object@minSize) != 1) {
        inv <- c(inv, "@minSize must be of length 1")
    }
    if(object@minSize < 1) {
        inv <- c(inv, "@minSize must be at least 1 or greater")
    }
    if(length(object@maxSize) != 1) {
        inv <- c(inv, "@maxSize must be of length 1")
    }
    if(object@maxSize < object@minSize) {
        inv <- c(inv, "@maxSize must be at least @minSize or greater")
    }
    if(length(object@kcdfNoneMinSampleSize) != 1) {
        inv <- c(inv, "@kcdfNoneMinSampleSize must be of length 1")
    }
    if(object@kcdfNoneMinSampleSize < 0) {
        inv <- c(inv, "@kcdfNoneMinSampleSize must be a non-negative integer")
    }
    if(is.na(object@kcdfNoneMinSampleSize)) {
        inv <- c(inv, "@kcdfNoneMinSampleSize must not be NA")
    }
    if(length(object@tau) != 1) {
        inv <- c(inv, "@tau must be of length 1")
    }
    if(is.na(object@tau)) {
        inv <- c(inv, "@tau must not be NA")
    }
    if(length(object@maxDiff) != 1) {
        inv <- c(inv, "@maxDiff must be of length 1")
    }
    if(is.na(object@maxDiff)) {
        inv <- c(inv, "@maxDiff must not be NA")
    }
    if(length(object@absRanking) != 1) {
        inv <- c(inv, "@absRanking must be of length 1")
    }
    if(is.na(object@absRanking)) {
        inv <- c(inv, "@absRanking must not be NA")
    }
    if(length(object@sparse) != 1) {
        inv <- c(inv, "@sparse must be of length 1")
    }
    if(is.na(object@sparse)) {
        inv <- c(inv, "@sparse must not be NA")
    }
    if(!.isCharLength1(object@checkNA)) {
        inv <- c(inv, "@use must be a single character string")
    }
    if(length(object@didCheckNA) != 1) {
        inv <- c(inv, "@didCheckNA must be of length 1")
    }
    if(is.na(object@didCheckNA)) {
        inv <- c(inv, "@didCheckNA must not be NA")
    }
    if(length(object@anyNA) != 1) {
        inv <- c(inv, "@anyNA must be of length 1")
    }
    if(is.na(object@anyNA)) {
        inv <- c(inv, "@anyNA must not be NA")
    }
    if(!.isCharLength1(object@use)) {
        inv <- c(inv, "@use must be a single character string")
    }
    if(length(object@filterRows) != 1) {
        inv <- c(inv, "@filterRows must be of length 1")
    }
    if(is.na(object@filterRows)) {
        inv <- c(inv, "@filterRows must not be NA")
    }
    if(length(object@nzcount) != 1) {
        inv <- c(inv, "@nzcount must be of length 1")
    }
    if(is.na(object@nzcount)) {
        inv <- c(inv, "@nzcount must not be NA")
    }
    if(!.isCharLength1(object@ondisk)) {
        inv <- c(inv, "@use must be a single character string")
    }
    return(if(length(inv) == 0) TRUE else inv)
})


## ----- getters -----

#' @noRd
.get_kcdf <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@kcdf)
}

#' @noRd
.get_kcdfNoneMinSampleSize <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@kcdfNoneMinSampleSize)
}

#' @noRd
.get_tau <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@tau)
}

#' @noRd
.get_maxDiff <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@maxDiff)
}

#' @noRd
.get_absRanking <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@absRanking)
}

#' @noRd
.get_sparse <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@sparse)
}

#' @noRd
.get_filterRows <- function(object) {
  stopifnot(inherits(object, "gsvaParam"))
  return(object@filterRows)
}

## getters for 'checkNA', 'didCheckNA', 'use' and 'ondisk' are
## in utils.R as they are shared with ssGSEA

#' @param x An object of class [`gsvaParam-class`].
#'
#' @param recursive Not used with `x` being an object of
#' class [`gsvaParam-class`].
#'
#' @aliases anyNA,gsvaParam-method
#' @rdname gsvaParam-class
setMethod("anyNA", signature=c("gsvaParam"),
          function(x, recursive=FALSE)
            return(x@anyNA))


## ----- details method -----

#' @importFrom methods callNextMethod
#' @importFrom GSEABase details
#' @aliases details,gsvaParam-method
#' @rdname GsvaMethodParam-class
#' @exportMethod details
setMethod("details",
          signature=signature(object="gsvaParam"),
          function(object) {
              callNextMethod(object)
              cat("kcdf: ", .get_kcdf(object), "\n",
                  "kcdfNoneMinSampleSize: ", .get_kcdfNoneMinSampleSize(object), "\n",
                  "tau: ", .get_tau(object), "\n",
                  "maxDiff: ", .get_maxDiff(object), "\n",
                  "absRanking: ", .get_absRanking(object), "\n",
                  sep="")
              cat("sparse: ", .get_sparse(object), "\n")
              cat("checkNA: ", .get_checkNA(object), "\n", sep="")
              if (.get_didCheckNA(object)) {
                  if (anyNA(object)) {
                      cat("missing data: yes\n",
                          "na_use: ", .get_NAuse(object), "\n", sep="")
                  } else
                      cat("missing data: no\n")
              } else
                  cat("missing data: didn't check\n")
              cat("filterRows: ", .get_filterRows(object), "\n")
          })

.gsvaParam_as_list <- function(x) {
    lst <- list(originalClassWasSE=is(get_exprData(x),
                                      "SummarizedExperiment"),
                geneSets=get_geneSets(x),
                assay=get_assay(x),
                annotation=get_annotation(x),
                minSize=get_minSize(x),
                maxSize=get_maxSize(x),
                nzcount=nzcount(x),
                ondisk=.get_ondisk(x))

    if (is(x, "ssgseaParam"))
        lst <- c(lst,
                 alpha=.get_alpha(x),
                 normalize=.get_normalize(x))

    if (is(x, "ssgseaParam") || is(x, "gsvaParam"))
        lst <- c(lst,
                 checkNA=.get_checkNA(x),
                 didCheckNA=.get_didCheckNA(x),
                 anyNA=anyNA(x),
                 use=.get_NAuse(x))

    if (is(x, "gsvaParam"))
        lst <- c(lst,
                 kcdf=.get_kcdf(x),
                 kcdfNoneMinSampleSize=.get_kcdfNoneMinSampleSize(x),
                 tau=.get_tau(x),
                 maxDiff=.get_maxDiff(x),
                 absRanking=.get_absRanking(x),
                 sparse=.get_sparse(x),
                 filterRows=.get_filterRows(x))
    return(lst)
}

## by now this is only called from gsva*() functions, i.e., no need
## to care about other methods
#' @importFrom S4Vectors metadata
.pull_param <- function(exprData) {

    p <- NULL
    if (is(exprData, "matrix") || is(exprData, "dgCMatrix") ||
        is(exprData, "SVT_SparseMatrix") || is(exprData, "DelayedMatrix") ||
        is(exprData, "HDF5Matrix") || is(exprData, "ExpressionSet")) {
        mask <- is.null(attr(exprData, "gsvaParam")) ||
                is.null(attr(exprData, "assay"))
        if (any(mask))
            cli_abort(c("x"="Missing metadata in the input expression data."))
        p <- attr(exprData, "gsvaParam")
        a <- attr(exprData, "assay")
        if (!a %in% c("gsvarownr", "gsvaranks"))
            cli_abort(c("x"="Wrong metadata in the input expression data."))
    } else { ## a SummarizedExperiment derivative
        if (is.null(metadata(exprData)$gsvaParam))
            cli_abort(c("x"="Missing metadata in the input expression data"))
        p <- metadata(exprData)$gsvaParam
        if (!any(assayNames(exprData) %in% c("gsvarownr", "gsvaranks"))) 
            cli_abort(c("x"="Wrong metadata in the input expression data."))
    }

    param <- new("gsvaParam",
                 exprData=exprData, geneSets=p$geneSets,
                 assay=p$assay, annotation=p$annotation,
                 minSize=p$minSize, maxSize=p$maxSize,
                 kcdf=p$kcdf, kcdfNoneMinSampleSize=p$kcdfNoneMinSampleSize,
                 tau=p$tau, maxDiff=p$maxDiff, absRanking=p$absRanking,
                 sparse=p$sparse, checkNA=p$checkNA, didCheckNA=p$didCheckNA,
                 anyNA=p$anyNA, use=p$use, filterRows=p$filterRows,
                 nzcount=p$nzcount, ondisk=p$ondisk)

    return(param)
}


#' @title GSVA ranks and scores
#'
#' @description Calculate GSVA scores in three steps: (1) normalize values of
#' expression by row; (2) calculate GSVA ranks by column from the previous
#' row-normalized values; and (3) calculate GSVA scores by column from the
#' previously calculated column ranks.
#'
#' @param param A [`gsvaParam-class`] object built using the constructor
#' function [`gsvaParam`].
#'
#' @param verbose Gives information about each calculation step. Default:
#' `TRUE`.
#'
#' @param dropExistingAssays Logical vector of length 1. It only applies when
#' the input expression data is stored using a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment]
#' derivative, which allows one to store more than one matrix of expression
#' values in different assay slots. By default `dropExistingAssays=FALSE` and
#' the new assay with the row-normalized expression values or the column ranks
#' will be stored as a new assay in the same input object. When
#' `dropExistingAssays=TRUE`, any existing assay will be dropped before adding
#' the new assay with the row-normalized expression values or the column ranks.
#'
#' @param first Numeric vector of length 1. First row, in the case of
#' `gsvaRowNorm()`, or first column, in the case of `gsvaColRanks()` and
#' `gsvaColScores()`, to which calculations should be restricted. By default,
#' `first=NA_real_`, which implies that calculations start at the first row or
#' column of the input expression data.
#'
#' @param last Numeric vector of length 1. Last row, in the case of
#' `gsvaRowNorm()`, or last column, in the case of `gsvaColRanks()` and
#' `gsvaColScores()`, to which calculations should be restricted. By default,
#' `last=NA_real_`, which implies that calculations end at the last row or
#' column of the input expression data.
#'
#' @param BPPARAM An object of class `BiocParallelParam` specifying parameters
#' related to the parallel execution of some of the tasks and calculations
#' within this function.
#'
#' @param maxmem A vector of length 1 either specifying a number in bytes, or
#' a character string with either the word `auto` (default), or a number
#' followed by a suffix indicating kilobytes (K), megabytes (M), gigabytes (G)
#' or terabytes (T), which GSVA will use to attempt bounding the maximum amount
#' of main memory used across all threads of execution to that given quantity.
#' By default `maxmem="auto"`, indicating that the maximum memory will be the
#' 90% of the total main memory, as calculated by
#' [`Sys.meminfo()`][memuse::Sys.meminfo]. To avoid setting any bound on the
#' maximum memory, use `maxmem=Inf`. Note that the amount of main memory used
#' in an R session or script may depend on other commands and packages used in
#' that same session or script.
#'
#' @seealso [`gsvaParam-class`], [`gsva`], [`gsvaEnrichment`],
#' [`BiocParallelParam`][BiocParallel::BiocParallelParam-class],
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' \doi{10.1186/1471-2105-14-7}
#'
#' @examples
#' library(GSVA)
#'
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#'
#' ## consider three disjoint gene sets
#' geneSets <- list(gset1=paste0("g", 1:3),
#'                  gset2=paste0("g", 4:6),
#'                  gset3=paste0("g", 7:10))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate row-normalized expression values
#' gsvarownormexpr <- gsvaRowNorm(gsvapar)
#'
#' ## calculate GSVA column ranks
#' gsvacolranks <- gsvaColRanks(gsvarownormexpr)
#'
#' ## calculate GSVA scores
#' gsva_es <- gsvaColScores(gsvacolranks)
#'
#' ## calculate now GSVA scores in a single step
#' gsva_es1 <- gsva(gsvapar)
#'
#' ## both approaches give the same result with the same input gene sets
#' all.equal(gsva_es1, gsva_es)
#'
#' ## however, results will be (obviously) different with different gene sets
#' geneSets2 <- list(gset1=paste0("g", 3:6),
#'                   gset2=paste0("g", c(1, 2, 7, 8)))
#'
#' ## note that there is no need to calculate the GSVA ranks again
#' ## geneSets(gsvarankspar) <- geneSets2
#' ## gsvaScores(gsvarankspar)
#'
#' @return In the case of 'gsvaRowNorm()', an object of the same class as the
#' input expresssion data given in the argument `exprData` of the `gsvaParam`
#' object, containing the row-normalized expression values. The resulting
#' object will have metadata with a copy of the input `gsvaParam` object,
#' except for the `exprData` slot, and in the case of being a derivative of a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment] object,
#' an additional assay called "gsvarownr" storing the row-normalized expression
#' values.
#'
#' @aliases gsvaRowNorm,gsvaParam-method
#' @name gsvaRowNorm
#' @rdname gsvaRanks
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @exportMethod gsvaRowNorm
setMethod("gsvaRowNorm", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   dropExistingAssays=FALSE,
                   first=NA_real_, last=NA_real_,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {

              if (verbose && gsva_global$show_start_and_end_messages) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
              }

              .check_bpparam(BPPARAM)

              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))

              checkedfl <- .check_first_last_values(dataMatrix, nrow, "rows",
                                                    first, last)
              first <- checkedfl$first
              last <- checkedfl$last

              maxmem <- .check_maxmem(param, maxmem=maxmem, verbose=verbose)
              ondisk <- .check_ondisk(param, first=first, last=last, whdim=1,
                                      maxmem=maxmem, verbose=verbose)

              dataMatrix <- .check_sparse_load_input_expr(dataMatrix, "GSVA",
                                                          first, last, whdim=1,
                                                          ondisk, verbose)

              filtDataMatrix <- dataMatrix
              BPPARAM <- .check_open_parallelism(filtDataMatrix, BPPARAM,
                                                 minparrows=100, minparcols=100,
                                                 verbose)

              if (.get_filterRows(param))
                  filtDataMatrix <- .filterGenes(dataMatrix, anyNA(param),
                                                 removeConstant=TRUE,
                                                 removeNzConstant=TRUE,
                                                 verbose, BPPARAM=BPPARAM,
                                                 maxmem=maxmem)
              else if (verbose) {
                  msg <- "Skipping filtering of constant rows (filterRows=FALSE)"
                  cli_alert_warning(msg)
              }
              
              if (verbose)
                  cli_alert_info(sprintf("Normalizing rows"))

              kcdfminssize <- .get_kcdfNoneMinSampleSize(param)
              gsvarownr <- .compute_row_norm(expr=filtDataMatrix,
                                             kcdf=.get_kcdf(param),
                                             kcdf.min.ssize=kcdfminssize,
                                             sparse=.get_sparse(param),
                                             any_na=anyNA(param),
                                             na_use=.get_NAuse(param),
                                             verbose=verbose,
                                             BPPARAM=BPPARAM,
                                             maxmem=maxmem)

              rownames(gsvarownr) <- rownames(filtDataMatrix)
              colnames(gsvarownr) <- colnames(filtDataMatrix)

              rval <- wrapData(get_exprData(param), gsvarownr, param,
                               "gsvarownr", first, last, whdim=1,
                               dropExistingAssays)

              if (verbose && gsva_global$show_start_and_end_messages)
                  cli_alert_success("Calculations finished")

              return(rval)
          })



#'
#' @param rowNormExprData A row-normalized expression data set obtained with
#' [`gsvaRowNorm`].  Must be one of the classes
#' supported by [`GsvaExprData-class`].  For a list of these classes, see its
#' help page using `help(GsvaExprData)`.
#'
#' @return In the case of 'gsvaColRanks()', an object of the same class as the
#' input expresssion data given in the argument `exprData` of the `gsvaParam`
#' object, containing the column rank values. The resulting object will have
#' metadata with a copy of the input `gsvaParam` object, except for the
#' `exprData` slot, and in the case of being a derivative of a
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment] object,
#' an additional assay called "gsvaranks" storing the column rank values.
#'
#' @aliases gsvaColRanks,GsvaExprData-method
#' @name gsvaColRanks
#' @rdname gsvaRanks
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @exportMethod gsvaColRanks
setMethod("gsvaColRanks", signature(rowNormExprData="GsvaExprData"),
          function(rowNormExprData,
                   verbose=TRUE,
                   dropExistingAssays=FALSE,
                   first=NA_real_, last=NA_real_,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {

              param <- .pull_param(rowNormExprData)

              if (verbose && gsva_global$show_start_and_end_messages) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
              }

              .check_bpparam(BPPARAM)

              dataMatrix <- unwrapData(rowNormExprData, "gsvarownr")

              checkedfl <- .check_first_last_values(dataMatrix, ncol, "columns",
                                                    first, last)
              first <- checkedfl$first
              last <- checkedfl$last

              maxmem <- .check_maxmem(param, assay="gsvarownr", maxmem=maxmem,
                                      verbose=verbose)
              ondisk <- .check_ondisk(param, assay="gsvarownr", first=first,
                                      last=last, whdim=2, maxmem=maxmem,
                                      verbose=verbose)

              dataMatrix <- .check_sparse_load_input_expr(dataMatrix, "GSVA",
                                                          first, last, whdim=2,
                                                          ondisk, verbose)
              sparse <- .get_sparse(param)
              if (sparse && !is_sparse(dataMatrix))
                  sparse <- FALSE

              gsvarnks <- .compute_gsva_ranks(Z=dataMatrix, sparse=sparse,
                                              verbose=verbose, BPPARAM=BPPARAM,
                                              maxmem=maxmem)

              rownames(gsvarnks) <- rownames(dataMatrix)
              colnames(gsvarnks) <- colnames(dataMatrix)

              rval <- wrapData(get_exprData(param), gsvarnks, param,
                               "gsvaranks", first, last, whdim=2,
                               dropExistingAssays)

              if (verbose && gsva_global$show_start_and_end_messages)
                  cli_alert_success("Calculations finished")

              return(rval)
          })


#' @param rankExprData A column-rank expression data set obtained with
#' [`gsvaColRanks`].  Must be one of the classes
#' supported by [`GsvaExprData-class`].  For a list of these classes, see its
#' help page using `help(GsvaExprData)`.
#'
#' @param geneSets An object of the classes supported by [`GsvaGeneSets-class`].
#' Currently, either a [`GeneSetCollection`][GSEABase::GeneSetCollection-class]
#' object or a `list` object.
#'
#' @return In the case of 'gsvaColScores()', an object of the same class as the
#' input expression data given in the argument `exprData` of the `gsvaParam`
#' object, containing the enrichment scores for the given gene sets. Note that
#' while it will have the same columns as the input expression data, the rows
#' will correspond to the gene sets for which the enrichment scores were
#' calculated.
#'
#' @aliases gsvaColScores,GsvaExprData-method
#' @name gsvaColScores
#' @rdname gsvaRanks
#'
#' @importFrom S4Arrays is_sparse
#' @importFrom cli cli_alert_info cli_alert_success
#' @exportMethod gsvaColScores
setMethod("gsvaColScores", signature(rankExprData="GsvaExprData"),
          function(rankExprData, geneSets, verbose=TRUE,
                   first=NA_real_, last=NA_real_,
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {

              param <- .pull_param(rankExprData)

              if (!missing(geneSets)) {
                  if (!is(geneSets, "GsvaGeneSets"))
                      cli_abort(c("x"=paste("'geneSets' must be a",
                                            "'GsvaGeneSets' object. See",
                                            "class ? GsvaGeneSets.")))
                  geneSets(param) <- geneSets
              }


              if (verbose && gsva_global$show_start_and_end_messages) {
                  pkgversion <- packageDescription("GSVA")[["Version"]]
                  cli_alert_info("GSVA version {pkgversion}")
              }

              .check_bpparam(BPPARAM)

              ## assuming rows in the rank data have been already filtered
              filtDataMatrix <- unwrapData(rankExprData, "gsvaranks")

              filtMappedGeneSets <- .filterAndMapGeneSets(param=param,
                                           filteredDataMatrix=filtDataMatrix,
                                           verbose=verbose)

              sparse <- .get_sparse(param)
              if (sparse && !is_sparse(filtDataMatrix))
                  sparse <- FALSE

              if (verbose) {
                if (sparse)
                    cli_alert_info("GSVA sparse algorithm")
                  else
                    cli_alert_info("GSVA dense (classical) algorithm")
              }

              checkedfl <- .check_first_last_values(filtDataMatrix, ncol,
                                                    "columns", first, last)
              first <- checkedfl$first
              last <- checkedfl$last

              maxmem <- .check_maxmem(param, assay="gsvaranks", maxmem=maxmem,
                                      verbose=verbose)
              ondisk <- .check_ondisk(param, assay="gsvaranks", first=first,
                                      last=last, whdim=2, maxmem=maxmem,
                                      verbose=verbose)


              filtDataMatrix <- .check_sparse_load_input_expr(filtDataMatrix,
                                                              "GSVA", first,
                                                              last, whdim=2,
                                                              ondisk, verbose)

              BPPARAM <- .check_open_parallelism(filtDataMatrix, BPPARAM,
                                                 minparrows=100, minparcols=100,
                                                 verbose)

              ondisk <- .check_es_memory_requirements(filtDataMatrix,
                                                      filtMappedGeneSets,
                                                      ondisk, maxmem)
              if (verbose) {
                  n <- length(filtMappedGeneSets)
                  cli_alert_info("Calculating GSVA scores for {n} gene sets")
              }

              gsva_es <- .processMatrixCols(filtDataMatrix,
                                            FUN=.compute_gsva_scores,
                                            geneSetsIdx=filtMappedGeneSets,
                                            tau=.get_tau(param),
                                            maxDiff=.get_maxDiff(param),
                                            absRanking=.get_absRanking(param),
                                            sparse=sparse, any_na=anyNA(param),
                                            na_use=.get_NAuse(param),
                                            minSize=get_minSize(param),
                                            ondisk=ondisk, verbose=verbose,
                                            minparrows=100, minparcols=100,
                                            BPPARAM=BPPARAM,
                                            maxmem=ceiling(maxmem/100)) ## use
                                            ## of memory increases here about
                                            ## 10-fold over block size memory

              rownames(gsva_es) <- names(filtMappedGeneSets)
              colnames(gsva_es) <- colnames(filtDataMatrix)

              gs <- .geneSetsIndices2Names(indices=filtMappedGeneSets,
                                           names=rownames(filtDataMatrix))

              ## dropAssays=TRUE for consistency but doesn't apply here
              rval <- wrapData(get_exprData(param), gsva_es, param, "es",
                               first, last, whdim=2, dropAssays=TRUE, gs)

              if (verbose && gsva_global$show_start_and_end_messages)
                  cli_alert_success("Calculations finished")

              return(rval)
          })

#' @title GSVA enrichment data and visualization
#'
#' @description Extract and plot enrichment data from GSVA scores.
#'
#' @param rankExprData A column-rank expression data set obtained with
#' [`gsvaColRanks`]. Must be one of the classes
#' supported by [`GsvaExprData-class`].  For a list of these classes, see its
#' help page using `help(GsvaExprData)`.
#'
#'
#' @param column The column for which we want to retrieve the enrichment data.
#' This parameter is only available in the `gsvaEnrichment()` method.
#'
#' @param geneSet Either a single positive integer number between 1 and the
#' number of available gene sets in parameter object stored in `rankExprData`,
#' or a single character string with the name of one of the gene sets available
#' in that object, or a vector of integers or character strings with the index
#' values or names of rows in `rankExprData` that should be considered as the
#' gene set for which the enrichment data should be retrieved.
#'
#' @param plot A character string indicating whether an enrichment plot should
#' be produced using either base R graphics (`plot="base"`) or the ggplot2 package
#' (`plot="ggplot"`), or not (`plot="no"`). In the latter case, the enrichment
#' data will be returned. By default `plot="auto"`, which implies that if this
#' method is called from an interactive session, a plot using base R graphics
#' will be produced and, otherwise, the enrichment data is returned.
#'
#' @param ... Further arguments passed to the `plot()` function when the
#' previous parameter `plot="base"`.
#'
#' @return When `plot="no"`, this method returns the enrichment data. When
#' `plot="ggplot"`, this method returns a `ggplot` object. When `plot="base"`
#' no value is returned.
#'
#' @seealso [`gsvaColRanks`], [`GsvaExprData-class`]
#'
#' @aliases gsvaEnrichment,GsvaExprData-method
#' @name gsvaEnrichment
#' @rdname gsvaEnrichment
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' \doi{10.1186/1471-2105-14-7}
#'
#' @examples
#' library(GSVA)
#'
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#' nGrp1 <- 15 ## number of samples in group 1
#' nGrp2 <- n - nGrp1 ## number of samples in group 2
#'
#' ## consider three disjoint gene sets
#' geneSets <- list(gset1=paste0("g", 1:3),
#'                  gset2=paste0("g", 4:6),
#'                  gset3=paste0("g", 7:10))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate GSVA ranks
#' gsvarownorm <- gsvaRowNorm(gsvapar)
#' gsvaranks <- gsvaColRanks(gsvarownorm)
#'
#' ## by default the enrichment data for the first column and the first
#' ## gene set in the input parameter object, are retrieved
#' gsvaEnrichment(gsvaranks)
#'
#' ## we can calculate the enrichment data for any of the gene sets given
#' ## in the input parameter object
#' gsvaEnrichment(gsvaranks, geneSet="gset2")
#'
#' ## we can calculate the enrichment data for a new gene set that did not
#' ## form part of the input parameter object
#' gsvaEnrichment(gsvaranks, geneSet=c("g1", "g4", "g7"))
#'
#' @importFrom cli cli_alert_info cli_abort cli_alert_danger
#' @importFrom utils installed.packages
#' @exportMethod gsvaEnrichment
setMethod("gsvaEnrichment", signature(rankExprData="GsvaExprData"),
          function(rankExprData, column=1, geneSet=1,
                   plot=c("auto", "base", "ggplot", "no"), ...) {

              if (length(column) != 1)
                  cli_abort(c("x"="'column' should be of length 1."))

              if (is.numeric(column) && !is.na(column)) {
                  if (column != as.integer(column) ||
                      column < 1 || column > ncol(rankExprData))
                      cli_abort(c("x"=paste("'column' should be a positive",
                                            "integer between 1 and the number",
                                            "of columns in the input data.")))
              } else
                  cli_abort(c("x"="'column' should be a positive integer."))
                    
              param <- .pull_param(rankExprData)

              plot <- match.arg(plot)

              geneSets <- get_geneSets(param)
              if (is.character(geneSet) && length(geneSet) == 1) {
                  if (!geneSet %in% names(geneSets)) {
                      msg <- paste("Gene set {geneSet} is missing from the input",
                                   "parameter object.")
                      cli_abort(c("x"=msg))
                  }
              } else if (is.numeric(geneSet) && length(geneSet) == 1 && !is.na(geneSet)) {
                  if (geneSet < 1 || geneSet > length(geneSets)) {
                       msg <- paste("When 'geneSet' is a single number, it",
                                    "should be a number between 1 and the",
                                    "number of gene sets ({length(geneSets)}).")
                       cli_abort(c("x"=msg))
                  }
              } else if (!is.character(geneSet) && !is.numeric(geneSet)) {
                  msg <- paste("input argument 'geneSet' should be either",
                               "numeric or character.")
                  cli_abort(c("x"=msg))
              }

              if (is.numeric(geneSet) && length(geneSet) > 1) {
                  if (any(geneSet != as.integer(geneSet)) || any(geneSet < 1) ||
                      any(geneSet > nrow(rankExprData))) {
                       msg <- paste("When 'geneSet' is a numeric vector,",
                                    "it should contain positive integers",
                                    "between 1 and the number of rows in the",
                                    "input data.")
                       cli_abort(c("x"=msg))
                  }
                  geneSets(param) <- list(geneSet)
                  geneSet <- 1
              }

              if (is.character(geneSet) && length(geneSet) > 1) {
                  if (!all(geneSet %in% rownames(rankExprData))) {
                      msg <- paste("When 'geneSet' is a character vector, all",
                                   "its values should be present in the row",
                                   "names of the input data.")
                      cli_abort(c("x"=msg))
                  }
                  geneSets(param) <- list(geneSet)
                  geneSet <- 1
              }

              tau <- .get_tau(param)
              maxDiff <- .get_maxDiff(param)
              absRanking <- .get_absRanking(param)
              sparse <- .get_sparse(param)
              any_na <- anyNA(param)
              na_use <- .get_NAuse(param)
              minsize <- get_minSize(param)

              exprData <- get_exprData(param)
              filtDataMatrix <- unwrapData(rankExprData, "gsvaranks")

              ## no need for verbosity when mapping a single gene set
              filtMappedGeneSets <- .filterAndMapGeneSets(param,
                                           wgset=geneSet, ## use that gene set
                                           filteredDataMatrix=filtDataMatrix,
                                           verbose=FALSE)

              geneSetIdx <- filtMappedGeneSets[[1]]
              edata <- .gsva_enrichment_data(R=filtDataMatrix,
                                             column=column,
                                             geneSetIdx=geneSetIdx,
                                             maxDiff=maxDiff,
                                             absRanking=absRanking,
                                             tau=tau,
                                             sparse=sparse,
                                             any_na=any_na,
                                             na_use=na_use,
                                             minSize=minsize)

              if (plot == "no" || (plot == "auto" && !interactive()))
                  return(edata)

              if (plot == "auto" || plot == "base")
                  .plot_enrichment_base(edata, ...) 
              else { ## plot == "ggplot"
                  instpkgs <- installed.packages(noCache=TRUE)[, "Package"]
                  if (!"ggplot2" %in% instpkgs)
                      cli_alert_danger("Please install the ggplot2 package")
                  else
                      .plot_enrichment_ggplot(edata)
              }
          })



#' @importFrom S4Arrays is_sparse
#' @importFrom DelayedArray seed
#' @importFrom cli cli_abort
compute.gene.cdf <- function(expr, Gaussk=TRUE, kernel=TRUE,
                             sparse=FALSE, any_na=FALSE,
                             na_use=c("everything", "all.obs", "na.rm"),
                             grid=NULL, verbose=TRUE, BPPARAM=NULL) {

    na_use <- match.arg(na_use)
    n.test.samples <- ncol(expr)
    n.genes <- nrow(expr)
    n.density.samples <- ncol(expr)

    if (any_na && na_use == "all.obs") {
        msg <- paste("missing values present in the input expression data and",
                     "'use=\"all.obs\".")
        cli_abort(c("x"=msg))
    }
    
    gene.cdf <- NA
    if (kernel) { ## kernel ECDF estimation
        if (is(expr, "dgCMatrix")) {
            if (sparse)
                gene.cdf <- .kcdfvals_sparse_to_sparse(expr, Gaussk, verbose)
            else
                gene.cdf <- .kcdfvals_sparse_to_dense(expr, Gaussk, verbose)
        } else if (is(expr, "SVT_SparseMatrix")) {
            if (sparse)
                gene.cdf <- .kcdfvals_svt_to_svt(expr, Gaussk, verbose)
            else
                gene.cdf <- .kcdfvals_svt_to_dense(expr, Gaussk, verbose)
        } else if (is(expr, "DelayedMatrix")) {
            if (sparse)
                gene.cdf <- .kcdfvals_sparseh5_to_sparseh5(expr, Gaussk=Gaussk,
                                                           grid=grid,
                                                           verbose)
            else {
                if (is_sparse(expr)) ## input HDF5 may be sparse or not
                    gene.cdf <- .kcdfvals_sparseh5_to_denseh5(expr,
                                                              Gaussk=Gaussk,
                                                              grid=grid,
                                                              verbose)
                else
                    gene.cdf <- .kcdfvals_denseh5_to_denseh5(expr,
                                                             Gaussk=Gaussk,
                                                             grid=grid,
                                                             verbose)
            }
        } else if (is.matrix(expr)) {
            A <- .Call("matrix_density_R",
                       as.double(t(expr)),
                       as.double(t(expr)),
                       n.density.samples,
                       n.test.samples,
                       n.genes,
                       as.integer(Gaussk),
                       any_na,
                       as.integer(factor(na_use,
                                         levels=c("everything", "all.obs",
                                                  "na.rm"))),
                       verbose)
            gene.cdf <- t(matrix(A, n.test.samples, n.genes))
        } else {
            msg <- "Matrix class {class(expr)} cannot be handled yet."
            cli_abort(c("x"=msg))
        }
    } else { ## direct ECDF estimation
        if (is(expr, "dgCMatrix")) {
            if (sparse)
                gene.cdf <- .ecdfvals_sparse_to_sparse(expr, verbose)
            else
                gene.cdf <- .ecdfvals_sparse_to_dense(expr, verbose)
        } else if (is(expr, "SVT_SparseMatrix")) {
            if (sparse)
                gene.cdf <- .ecdfvals_svt_to_svt(expr, verbose)
            else
                gene.cdf <- .ecdfvals_svt_to_dense(expr, verbose)
        } else if (is(expr, "DelayedMatrix")) {
            if (sparse)
                gene.cdf <- .ecdfvals_sparseh5_to_sparseh5(expr, grid=grid,
                                                           verbose=verbose)
            else {
                if (is_sparse(expr)) ## input HDF5 may be sparse or not
                    gene.cdf <- .ecdfvals_sparseh5_to_denseh5(expr, grid=grid,
                                                              verbose)
                else
                    gene.cdf <- .ecdfvals_denseh5_to_denseh5(expr, grid=grid,
                                                             verbose)
            }
        } else if (is.matrix(expr)) {
            if (any_na)
                gene.cdf <- .ecdfvals_dense_to_dense_nas(expr, verbose)
            else
                gene.cdf <- .ecdfvals_dense_to_dense(expr, verbose)
        } else {
            msg <- "Input container class {class(expr)} cannot be handled yet."
            cli_abort(c("x"=msg))
        }
    }

    if (ncol(expr) > 10000) ## free up ASAP memory we need not anymore and was
        out <- gc()         ## allocated during ECDF calculations on a big expr

    return(gene.cdf)	
}

#' @importFrom Matrix nnzero
.sufficient_ssize <- function(expr, kcdf.min.ssize) {
  ## in the sparse case stored in a 'dgCMatrix' or a 'SVT_SparseMatrix',
  ## by now, use the average nonzero values per row
  if (is_sparse(expr)) {
    nnz <- nnzero(expr)
    if (is.na(nnz)) {
        msg <- "The input sparse matrix of expression contains NA values."
        cli_abort(c("x"=msg))
    }
    return((nnz / nrow(expr)) >= kcdf.min.ssize)
  }

  ## in every other case, including the dense case, by now,
  ## just look at the number of columns
  return(ncol(expr) >= kcdf.min.ssize)
}

#' @importFrom S4Arrays is_sparse
.parse_kcdf_param <- function(expr, kcdf, kcdf.min.ssize, sparse, verbose) {
    kernel <- FALSE
    Gaussk <- TRUE  ## default (TRUE) is a Gaussian kernel, Poisson otherwise (FALSE)
    if (kcdf == "auto") {
        if (verbose)
            cli_alert_info("kcdf='auto' (default)")
        if (!.sufficient_ssize(expr, kcdf.min.ssize)) {
            kernel <- TRUE
            if (is(expr, "dgCMatrix")) { ## dgCMatrix does not store integers
                                         ## so we check them with x == floor(x)
                sam <- sample(expr@x, size=min(1000, length(expr@x)),
                              replace=FALSE)
                Gaussk <- any((sam < 0) | (sam != floor(sam)))
            } else if (is.integer(expr[1, 1]))
                Gaussk <- FALSE
        }
    } else {
        if (kcdf == "Gaussian") {
            kernel <- TRUE
            Gaussk <- TRUE
        } else if (kcdf == "Poisson") {
            kernel <- TRUE
            Gaussk <- FALSE
        } else
            kernel <- FALSE
    }

    if (verbose) {
        is_sparse_matrix <- is(expr, "dgCMatrix") ||
                            is(expr, "SVT_SparseMatrix") ||
                            (is(expr, "DelayedMatrix") && is_sparse(expr))
        if (is_sparse_matrix && sparse)
            cli_alert_info("GSVA sparse algorithm")
        else
            cli_alert_info("GSVA dense (classical) algorithm")
        if (kernel) {
            if (Gaussk)
                cli_alert_info("Row-wise ECDF estimation with Gaussian kernels")
            else
                cli_alert_info("Row-wise ECDF estimation with Poisson kernels")
        } else
            cli_alert_info("Direct row-wise ECDFs estimation")
    }

    list(kernel=kernel, Gaussk=Gaussk)
}



#' @importFrom cli cli_alert_info
.compute_row_norm <- function(expr, kcdf, kcdf.min.ssize,
                              sparse, any_na, na_use, verbose,
                              BPPARAM=NULL, maxmem=Inf) {

    kcdfparam <- .parse_kcdf_param(expr, kcdf, kcdf.min.ssize, sparse, verbose)
    kernel <- kcdfparam$kernel
    Gaussk <- kcdfparam$Gaussk

    if (verbose)
       cli_alert_info("Calculating row ECDFs")

    Z <- .processMatrixRows(expr, FUN=compute.gene.cdf, Gaussk=Gaussk,
                            kernel=kernel, sparse=sparse, any_na=any_na,
                            na_use=na_use, verbose=verbose, minparrows=100,
                            minparcols=100, BPPARAM=BPPARAM, maxmem=maxmem)

    return(Z)
}



## here 'ties.method="last"' allows one to obtain the result
## from 'order()' based on ranks
## pending how to propagate verbosity if necessary

#' @importFrom MatrixGenerics colRanks
#' @importFrom sparseMatrixStats colRanks
compute.col.ranks <- function(Z, ties.method="last", drop.sparsity=FALSE,
                              verbose=TRUE) {
    R <- NULL

    if (drop.sparsity && !is(Z, "DelayedMatrix"))
        Z <- as.matrix(Z)

    if (is(Z, "dgCMatrix")) { ## assumes expression values are positive
        R <- .sparseColumnApplyAndReplace(Z, rank, ties.method=ties.method)
    } else if (is(Z, "SVT_SparseMatrix")) {
        R <- .colRanks_SVT_SparseMatrix(Z, ties.method=ties.method)
    } else if (is(Z, "DelayedMatrix")) {
        R <- .colRanksHDF5(Z, ties.method=ties.method, drop.sparsity=drop.sparsity)
    } else {
        R <- colRanks(Z, ties.method=ties.method, preserveShape=TRUE)
    }

    if (ncol(Z) > 10000) ## free up ASAP memory we need not anymore and was
        out <- gc()      ## allocated during rank calculations on a big Z

    return(R)
}

#' @importFrom cli cli_alert_info
#' @importFrom cli cli_progress_done cli_abort
.compute_gsva_ranks <- function(Z, sparse, verbose, BPPARAM=NULL, maxmem=Inf) {
    if (verbose) {
        if (sparse)
            cli_alert_info("Calculating sparse column ranks")
        else
            cli_alert_info("Calculating column ranks")
    }
 
    ## here 'ties.method="last"' allows one to obtain the result
    ## from 'order()' based on ranks
    R <- .processMatrixCols(Z, FUN=compute.col.ranks, ties.method="last",
                            drop.sparsity=FALSE, verbose=verbose,
                            minparrows=100, minparcols=100,
                            BPPARAM=BPPARAM, maxmem=maxmem)

    return(R)
}

## here gSetIdx, decOrderStat and symRnkStat contain the positions with respect
## to the original order of genes in the data
.gsvaRndWalk <- function(gSetIdx, decOrderStat, symRnkStat, tau) {
    n <- length(decOrderStat)
    k <- length(gSetIdx)
    gSetRnk <- decOrderStat[gSetIdx]

    stepCDFinGeneSet <- integer(n)
    if (tau == 1)
      stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]
    else {
      stepCDFinGeneSet <- numeric(n)
      stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]^tau
    }

    stepCDFinGeneSet <- cumsum(stepCDFinGeneSet)
    stepCDFoutGeneSet <- rep(1L, n)
    stepCDFoutGeneSet[gSetRnk] <- 0L
    stepCDFoutGeneSet <- cumsum(stepCDFoutGeneSet)

    walkStat <- rep(NA_real_, n)
    if (stepCDFinGeneSet[n] > 0 && stepCDFoutGeneSet[n] > 0) {
        stepCDFinGeneSet <- stepCDFinGeneSet / stepCDFinGeneSet[n]
        stepCDFoutGeneSet <- stepCDFoutGeneSet / stepCDFoutGeneSet[n]

        walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
    }

    walkStat
}

.gsva_score_genesets_Rimp <- function(geneSetsIdx, decOrdStat, symRnkStat,
                                      maxDiff, absRanking, tau, any_na, na_use,
                                      minSize) {
   md <- lapply(geneSetsIdx, function(gSetIdx, decOrdStat, symRnkStat) {
             maxDev <- c(NA_real_, NA_real_)
             if (any_na) {
                 walkStat <- .gsvaRndWalk_nas(gSetIdx, decOrdStat, symRnkStat,
                                              tau, na_use, minSize)
                 if (any(!is.na(walkStat))) {
                     if (na_use == "na.rm")
                         maxDev <- c(max(c(0, max(walkStat, na.rm=TRUE))),
                                     min(c(0, min(walkStat, na.rm=TRUE))))
                     else
                         maxDev <- c(max(c(0, max(walkStat))),
                                     min(c(0, min(walkStat))))
                 }
             } else {
                 walkStat <- .gsvaRndWalk(gSetIdx, decOrdStat, symRnkStat, tau)
                 maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
             }
             maxDev
         }, decOrdStat, symRnkStat)
   md <- do.call("rbind", md)
   if (maxDiff && absRanking)
       md[, 2] <- -1 * md[, 2]
   sco <- rowSums(md)
   if (!maxDiff) {
       if (any_na) {
         mask <- is.na(sco)
         sco[!mask] <- md[cbind(seq_len(sum(!mask)), ifelse(sco[!mask] > 0, 1, 2))]
       } else
         sco <- md[cbind(seq_along(sco), ifelse(sco > 0, 1, 2))]
   }
   sco
}

## here gSetIdx, decOrderStat and symRnkStat contain the positions with respect
## to the original order of genes in the data
.gsvaRndWalk_nas <- function(gSetIdx, decOrderStat, symRnkStat, tau, na_use,
                             minSize=1L, wna_env) {
    n <- length(decOrderStat)
    gSetRnk <- decOrderStat[gSetIdx]

    if (anyNA(gSetRnk)) {
        if (na_use == "everything")
            return(rep(NA_real_, n))
        else if (na_use == "all.obs")
            cli_abort(c("x"="Input GSVA ranks have NA values."))
        else if (na_use == "na.rm") {
            gSetIdx <- gSetIdx[!is.na(gSetRnk)]
            gSetRnk <- gSetRnk[!is.na(gSetRnk)]
        }
    }
    k <- length(gSetIdx)

    walkStat <- rep(NA_real_, n)
    if (k >= minSize) {

        stepCDFinGeneSet <- integer(n)
        if (tau == 1)
          stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]
        else {
          stepCDFinGeneSet <- numeric(n)
          stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]^tau
        }

        stepCDFinGeneSet <- cumsum(stepCDFinGeneSet)
        stepCDFoutGeneSet <- rep(1L, n)
        stepCDFoutGeneSet[gSetRnk] <- 0L
        stepCDFoutGeneSet <- cumsum(stepCDFoutGeneSet)

        if (stepCDFinGeneSet[n] > 0 && stepCDFinGeneSet[n] > 0) {
            stepCDFinGeneSet <- stepCDFinGeneSet / stepCDFinGeneSet[n]
            stepCDFoutGeneSet <- stepCDFoutGeneSet / stepCDFoutGeneSet[n]

            walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
        }
    } else if (!get("w", envir=wna_env)) ## warn only once. it can only happen
        assign("w", TRUE, envir=wna_env) ## with na_use="na.rm"

    walkStat
}

## convert ranks into decreasing order statistics and symmetric rank statistics
.ranks2stats <- function(r, sparse) {
    mask <- r == 0
    p <- length(r)
    r_dense <- as.integer(r)          ## assume ranks are integer

    if (any(mask)) {                  ## sparse ranks into dense ranks
        nzs <- sum(mask)
        r_dense[!mask] <- r_dense[!mask] + nzs ## shift ranks of nonzero values
        r_dense[mask] <- seq.int(nzs)          ## zeros get increasing ranks
    }

    dos <- p - r_dense + 1L           ## dense ranks into decreasing order stats
    srs <- numeric(p)

    if (any(mask) && sparse) {
        r[!mask] <- r[!mask] + 1      ## shift ranks of nonzero values by one
        r[mask] <- 1                  ## all zeros get the same first rank
        srs <- abs(max(r)/2 - r)
    } else
        srs <- abs(p/2 - r_dense)

    list(dos=dos, srs=srs)
}

## convert ranks into decreasing order statistics and symmetric rank statistics
## skipping NA values
.ranks2stats_nas <- function(r, sparse) {
    na_mask <- is.na(r)

    if (all(na_mask))
        return(list(dos=rep(NA, length(r)), srs=rep(NA, length(r))))

    n_nas <- sum(na_mask)
    mode(n_nas) <- "integer"
    mask <- !na_mask & r == 0
    p <- length(r)
    r_dense <- as.integer(r)          ## assume ranks are integer

    if (any(mask)) {                  ## sparse ranks into dense ranks
        nzs <- sum(mask)
        mode(nzs) <- "integer"
        r_dense[!mask] <- r_dense[!mask] + nzs ## shift ranks of nonzero values
        r_dense[mask] <- seq.int(nzs)          ## zeros get increasing ranks
    }

    dos <- p - n_nas - r_dense + 1L   ## dense ranks into decreasing order stats
    srs <- numeric(p)

    if (any(mask) && sparse) {
        r[!mask] <- r[!mask] + 1L     ## shift ranks of nonzero values by one
        r[mask] <- 1L                 ## all zeros get the same first rank
        srs <- abs(max(r, na.rm=TRUE)/2 - r)
    } else
        srs <- abs((p - n_nas)/2 - r_dense)

    list(dos=dos, srs=srs)
}


## this function computes the GSVA scores for all gene sets in geneSetsIdx for
## a given rank matrix R, taking care that if 'ondisk=TRUE' because, e.g., the
## resulting matrix of GSVA scores does not fit in main memory, the scores are
## written into an on-disk data structure (HDF5) instead of being returned in
## main memory.
#' @importFrom cli cli_alert_info cli_alert_warning
#' @importFrom S4Arrays is_sparse refdim DummyArrayGrid read_block write_block
#' @importFrom DelayedArray close
.compute_gsva_scores <- function(R, geneSetsIdx, tau, maxDiff, absRanking,
                                 sparse, any_na, na_use, minSize, ondisk,
                                 verbose) {
    p <- nrow(R)
    n <- ncol(R)
    es <- NULL
    if (sparse && !is_sparse(R))
        sparse <- FALSE
    intrnks <- is.integer(R[1, 1])

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)
    es <- NULL
    if (is(R, "DelayedMatrix") || ondisk) {
        sink <- HDF5RealizationSink(c(length(geneSetsIdx), ncol(R)),
                                    as.sparse=FALSE) ## GSVA scores are dense
        grid <- DummyArrayGrid(dim(R))
        grid_es <- DummyArrayGrid(dim(sink))

        if (length(grid) != length(grid_es) ||
            refdim(grid)[2] != refdim(grid_es)[2] ||
            dim(grid)[2] != dim(grid_es)[2]) {
            msg <- paste("Grid column blocks for ranks should match grid",
                         "column blocks for enrichment scores")
            cli_abort(c("x"=msg))
        }

        ## avp - ArrayViewport for reaching the (possibly sparse) rank matrix
        ## avp_es - ArrayViewport for writing the enrichment dense scores matrix
        colScores_byBlock <- function(avp, avp_es, sink) {
            block <- read_block(R, avp)
            block <- .gsva_score_genesets(block, geneSetsIdx, intrnks, sparse,
                                          maxDiff, absRanking, tau, any_na,
                                          na_use, minSize, wna_env, verbose)
            write_block(sink, avp_es, block)
        }

        nblock <- length(grid)
        for (bid in seq_len(nblock))
            sink <- colScores_byBlock(grid[[bid]], grid_es[[bid]], sink)
        close(sink)
        es <- as(sink, "DelayedArray")

    } else {
        es <- .gsva_score_genesets(R, geneSetsIdx, intrnks, sparse, maxDiff,
                                   absRanking, tau, any_na, na_use, minSize,
                                   wna_env, verbose)
    }

    if (any_na && na_use == "na.rm")
        if (get("w", envir=wna_env)) {
            msg <- paste("NA enrichment scores in gene sets with less",
                         "than {minSize} genes after removing missing values")
            cli_alert_warning(msg)
        }

    if (ncol(R) > 10000) ## free up ASAP memory we need not anymore and was
        out <- gc()      ## allocated during score calculations on a big R

    return(es)
}

#' @importFrom S4Arrays is_sparse
.gsva_enrichment_data <- function(R, column, geneSetIdx, maxDiff,
                                  absRanking, tau, sparse, any_na,
                                  na_use, minSize) {
    n <- ncol(R)
    es <- NULL
    if (!is_sparse(R))
        sparse <- FALSE
    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    if (any_na) {
        rnkstats <- .ranks2stats_nas(R[, column], sparse)
        walkStat <- .gsvaRndWalk_nas(geneSetIdx, rnkstats$dos, rnkstats$srs,
                                     tau, na_use, minSize, wna_env=wna_env)
        maxDev <- whMaxDev <- c(NA, NA)
        if (any(!is.na(walkStat))) {
          if (na_use == "na.rm")
              maxDev <- c(max(c(0, max(walkStat, na.rm=TRUE))),
                          min(c(0, min(walkStat, na.rm=TRUE))))
          else
              maxDev <- c(max(c(0, max(walkStat))),
                          min(c(0, min(walkStat))))
        }
        if (length(which.max(walkStat)) > 0)
            whMaxDev[1] <- which.max(walkStat)
        if (length(which.min(walkStat)) > 0)
            whMaxDev[2] <- which.min(walkStat)
    } else {
        rnkstats <- .ranks2stats(R[, column], sparse)
        walkStat <- .gsvaRndWalk(geneSetIdx, rnkstats$dos, rnkstats$srs, tau)
        maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
        whMaxDev <- c(which.max(walkStat), which.min(walkStat))
        whMaxDev[maxDev == 0] <- NA
    }
    
    if (any_na && na_use == "na.rm")
        if (get("w", envir=wna_env)) {
            msg <- paste("Gene set has fewer than {minSize} genes after",
                         "removing missing values, no enrichment data",
                         "available")
            cli_alert_warning(msg)
            return(list())
        }

    if (maxDiff && absRanking)
        maxDev[2] <- -1 * maxDev[2]
    sco <- sum(maxDev)
    if (!maxDiff) {
        if (any_na) {
            if (!is.na(sco)) {
                sco <- maxDev[1]
                if (abs(maxDev[2]) > maxDev[1])
                    sco <- maxDev[2]
            }
        }
    }

    edat <- data.frame(rank=seq.int(nrow(R)),
                       stat=walkStat)
    rownames(edat)[na.omit(rnkstats$dos)] <- rownames(R)[!is.na(rnkstats$dos)]

    geneSetIdx <- geneSetIdx[!is.na(rnkstats$dos[geneSetIdx])]
    gsetrnk <- rnkstats$dos[geneSetIdx]
    lepos <- leneg <- NA
    if (!is.na(whMaxDev[1]))
        lepos <- geneSetIdx[gsetrnk <= whMaxDev[1]]
    if (!is.na(whMaxDev[2])) {
        if (!is.na(whMaxDev[1]) && whMaxDev[2] < whMaxDev[1]) {
            mask <- gsetrnk >= whMaxDev[2] & gsetrnk <= whMaxDev[1]
            lepos <- leneg <- geneSetIdx[mask]
        } else
            leneg <- geneSetIdx[gsetrnk >= whMaxDev[2]]
    }
    if (all(!is.na(lepos)))
        lepos <- rownames(R)[lepos]
    if (all(!is.na(leneg)))
        leneg <- rownames(R)[leneg]

    res <- list(stats=edat,
                gsetrnk=gsetrnk,
                maxPos=maxDev[1],
                whichMaxPos=whMaxDev[1],
                maxNeg=maxDev[2],
                whichMaxNeg=whMaxDev[2],
                leadingEdgePos=lepos,
                leadingEdgeNeg=leneg,
                score=sco,
                tau=tau,
                maxDiff=maxDiff,
                absRanking=absRanking,
                sparse=sparse)

    return(res)
}

#' @importFrom graphics plot abline grid lines segments
.plot_enrichment_base <- function(edata, ...) {
    ylim <- range(edata$stats$stat)
    hgsetticks <- (ylim[2] - ylim[1]) * 0.1
    plot(edata$stats, type="l", lwd=2, las=1, panel.first=grid(),
         xlab="Gene Ranking", ylab="Random Walk Statistic", col="green", ...)
    abline(h=0, lwd=2, lty=2, col="grey")
    lines(edata$stats, lwd=2, col="green")
    segments(edata$gsetrnk, -hgsetticks/2, edata$gsetrnk, hgsetticks/2, lwd=2)
    if (!is.na(edata$whichMaxPos) &&
        (edata$maxDiff || edata$maxPos >= abs(edata$maxNeg)))
        segments(edata$whichMaxPos, 0, edata$whichMaxPos, edata$maxPos,
                 lwd=2, lty=2, col="darkred")
    if (!is.na(edata$whichMaxNeg) &&
        (edata$maxDiff || edata$maxPos < abs(edata$maxNeg)))
        segments(edata$whichMaxNeg, 0, edata$whichMaxNeg, edata$maxNeg,
                 lwd=2, lty=2, col="darkred")
}

#' @importFrom cli cli_abort
#' @importFrom utils globalVariables
.plot_enrichment_ggplot <- function(edata, ...) {
    if (!.isPackageLoaded("ggplot2")) {
        loaded <- suppressPackageStartupMessages(requireNamespace("ggplot2"))
        if (!loaded)
            cli_abort(c("x"="ggplot2 could not be loaded"))
    }

    ylim <- range(edata$stats$stat)
    hgsetticks <- (ylim[2] - ylim[1]) * 0.1
    gsetticks <- data.frame(gsetrnk=edata$gsetrnk)
    ## from https://stackoverflow.com/a/39877048
    fintticks <- function(x) unique(floor(pretty(seq(min(x),
                                    (max(x) + 1) * 1.1))))
    .data <- ggplot2::.data
    ggplot2::ggplot(data=edata$stats) +
        ggplot2::scale_x_continuous(breaks=fintticks) +
        ggplot2::geom_line(ggplot2::aes(x=.data$rank, y=.data$stat),
                           color="green") +
        ggplot2::geom_segment(data=gsetticks,
                     mapping=ggplot2::aes(x=.data$gsetrnk, y=-hgsetticks/2,
                                 xend=.data$gsetrnk, yend=hgsetticks/2),
                     linewidth=1) +
        ggplot2::geom_hline(yintercept=0, colour="grey", linetype="dashed") +
        { if (!is.na(edata$whichMaxPos) &&
              (edata$maxDiff || edata$maxPos >= abs(edata$maxNeg)))
              ggplot2::geom_segment(data=data.frame(whichMaxPos=edata$whichMaxPos,
                                                    maxPos=edata$maxPos),
                           mapping=ggplot2::aes(x=.data$whichMaxPos, y=0,
                                       xend=.data$whichMaxPos,
                                       yend=.data$maxPos),
                           colour="darkred", linetype="dashed") } +
        { if (!is.na(edata$whichMaxPos) &&
              (edata$maxDiff || edata$maxPos < abs(edata$maxNeg)))
              ggplot2::geom_segment(data=data.frame(whichMaxNeg=edata$whichMaxNeg,
                                                    maxNeg=edata$maxNeg),
                           mapping=ggplot2::aes(x=.data$whichMaxNeg, y=0,
                                       xend=.data$whichMaxNeg,
                                       yend=.data$maxNeg),
                           colour="darkred", linetype="dashed") } +
        ggplot2::theme(panel.background=ggplot2::element_blank(),
              panel.grid.major=ggplot2::element_line(colour="grey",
                                                     linetype="dotted"),
              panel.grid.minor=ggplot2::element_line(colour=NA),
              axis.text=ggplot2::element_text(size=12),
              axis.title=ggplot2::element_text(size=14),
              panel.border=ggplot2::element_rect(colour="black", fill=NA)) +
        ggplot2::labs(x="Gene Ranking", y="Random Walk Statistic")
}

##
## functions interfacing C code
##

.fetch_row_nzvals <- function(X, i, whimin1=NULL) {
  stopifnot(is(X, "SVT_SparseMatrix")) ## QC
  stopifnot(is.numeric(i)) ## QC
  if (!is.null(whimin1)) {
      stopifnot(is.numeric(whimin1)) ## QC
      whimin1 <- as.integer(whimin1)
      stopifnot(length(whimin1) == ncol(X)) ## QC
  }
  .Call("fetch_row_nzvals_R", X, as.integer(i), whimin1)
}

.ecdfvals_svt_to_dense <- function(X, verbose) {
  stopifnot(is(X, "SVT_SparseMatrix")) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_svt_to_dense_R", X, verbose)
}

.ecdfvals_svt_to_svt <- function(X, verbose) {
  stopifnot(is(X, "SVT_SparseMatrix")) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_svt_to_svt_R", X, verbose)
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom DelayedArray seed gridReduce close
.ecdfvals_sparseh5_to_sparseh5 <- function(X, grid=NULL, verbose=FALSE) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=TRUE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowEcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .ecdfvals_svt_to_svt(block, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowEcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom DelayedArray seed gridReduce close
.ecdfvals_sparseh5_to_denseh5 <- function(X, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowEcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .ecdfvals_svt_to_dense(block, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowEcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce close
.ecdfvals_denseh5_to_denseh5 <- function(X, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowEcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .ecdfvals_dense_to_dense(block, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowEcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

.ecdfvals_sparse_to_sparse <- function(X, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_sparse_to_sparse_R", X, Xrsp, verbose)
}

.ecdfvals_sparse_to_dense <- function(X, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_sparse_to_dense_R", X, Xrsp, verbose)
}

.ecdfvals_dense_to_dense <- function(X, verbose) {
  stopifnot(is.matrix(X)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_dense_to_dense_R", X, verbose)
}

.ecdfvals_dense_to_dense_nas <- function(X, verbose) {
  stopifnot(is.matrix(X)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_dense_to_dense_nas_R", X, verbose)
}

.kcdfvals_svt_to_dense <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "SVT_SparseMatrix")) ## QC
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_svt_to_dense_R", X, Gaussk, verbose)
}

.kcdfvals_svt_to_svt <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "SVT_SparseMatrix")) ## QC
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_svt_to_svt_R", X, Gaussk, verbose)
}

.kcdfvals_sparse_to_sparse <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_sparse_to_sparse_R", X, Xrsp, Gaussk, verbose)
}

.kcdfvals_sparse_to_dense <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_sparse_to_dense_R", X, Xrsp, Gaussk, verbose)
}

#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce close
.kcdfvals_sparseh5_to_sparseh5 <- function(X, Gaussk, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=TRUE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowKcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .kcdfvals_svt_to_svt(block, Gaussk=Gaussk, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowKcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce close
.kcdfvals_sparseh5_to_denseh5 <- function(X, Gaussk, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowKcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .kcdfvals_svt_to_dense(block, Gaussk=Gaussk, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowKcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce close
.kcdfvals_denseh5_to_denseh5 <- function(X, Gaussk, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowKcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- t(matrix(.Call("matrix_density_R",
                            as.double(t(X)),
                            as.double(t(X)),
                            ncol(X),
                            ncol(X),
                            nrow(X),
                            as.integer(Gaussk),
                            FALSE, 1L,
                            verbose), ncol(X), nrow(X)))
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowKcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom cli cli_abort
.gsva_score_genesets <- function(R, geneSetsIdx, intrnks, sparse, maxDiff,
                                 absRanking, tau, any_na, na_use, minSize,
                                 wna_env, verbose) {
    minSize <- as.integer(minSize)
    stopifnot(is.list(geneSetsIdx)) ## QC
    stopifnot(length(geneSetsIdx) > 0) ## QC
    stopifnot(is.integer(geneSetsIdx[[1]])) ## QC
    stopifnot(is.logical(intrnks)) ## QC
    stopifnot(is.logical(sparse)) ## QC
    stopifnot(is.logical(maxDiff)) ## QC
    stopifnot(is.logical(absRanking)) ## QC
    stopifnot(is.numeric(tau)) ## QC but it still might be an integer!!
    stopifnot(is.logical(any_na)) ## QC
    stopifnot(is.character(na_use)) ## QC
    stopifnot(is.integer(minSize)) ## QC
    stopifnot(is.logical(verbose)) ## QC
    na_use <- as.integer(factor(na_use, levels=c("everything", "all.obs",
                                                 "na.rm")))
    sco <- .Call("gsva_score_genesets_R", R, geneSetsIdx, intrnks, sparse,
                 maxDiff, absRanking, as.double(tau), any_na, na_use, minSize,
                 verbose)

    if (any_na) {
      if (na_use == 2 && !is.null(attr(sco, "attrNAs")))
          cli_abort(c("x"="Input GSVA ranks have NA values."))

      if (na_use == 3 && !is.null(attr(sco, "attrNAs")))
          assign("w", TRUE, envir=wna_env)

      attr(sco, "attrNAs") <- NULL ## clean up the NA informing attribute
    }

    sco
}

## calculate ranks using on an SVT_SparseMatrix object
#' @importFrom BiocGenerics type
.colRanks_SVT_SparseMatrix <- function(X, ties.method="last") {
    R <- X
    whposlen <- which(lengths(X@SVT) > 0L)
    rnks <- lapply(lapply(X@SVT[whposlen], "[[", 1), rank,
                          ties.method=ties.method)
    R@SVT[whposlen] <- mapply(list, rnks, lapply(X@SVT[whposlen], "[[", 2),
                              SIMPLIFY=FALSE)
    if (ties.method == "last")
        R@type <- "integer" ## rank() w/ ties.method="last" returns integer

    R
}

## calculate ranks using an HDF5 backend

#' @importFrom BiocGenerics "type<-"
#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom MatrixGenerics colRanks
#' @importFrom BiocParallel SerialParam
#' @importFrom DelayedArray close
.colRanksHDF5 <- function(X, grid=NULL, ties.method="last",
                          drop.sparsity=FALSE) {
    stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

    sink <- HDF5RealizationSink(dim(X), H5type="H5T_STD_I32LE", ## integer ranks
                                as.sparse=is_sparse(X) && !drop.sparsity)
    if (is.null(grid))
        grid <- DummyArrayGrid(dim(X))

    colRanks_byBlock <- function(grid, sink) {
        block <- read_block(X, grid)
        if (is(block, "SVT_SparseMatrix") && drop.sparsity)
            block <- as.matrix(block)
        if (is(block, "SVT_SparseMatrix")) {
            block <- .colRanks_SVT_SparseMatrix(block, ties.method=ties.method)
        } else {
            block <- colRanks(block, ties.method=ties.method,
                              preserveShape=TRUE)
            if (ties.method == "last")
                type(block) <- "integer"
        }
        write_block(sink, grid, block)
    }

    sink <- gridReduce(colRanks_byBlock, grid, sink)
    close(sink)
    res <- as(sink, "DelayedArray")
    res
}
