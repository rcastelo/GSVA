
## ----- Virtual Superclasses -----

#' `GsvaExprData` class
#'
#' Virtual superclass of expression data classes supported by `GSVA`.
#'
#' `GSVA` supports expression data matrices in a growing number of containers
#' and representations.  This class union allows to store any of these in a slot
#' of another class as well as defining common methods for all of them.
#'
#' @seealso
#' [`matrix`],
#' [`dgCMatrix-class`],
#' [`ExpressionSet-class`],
#' [`SummarizedExperiment-class`],
#' [`SingleCellExperiment-class`],
#' [`SpatialExperiment-class`]
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom Biobase ExpressionSet
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importClassesFrom DelayedArray DelayedArray
#' @importClassesFrom HDF5Array HDF5Array
#'
#' @aliases GsvaExprData
#' @name GsvaExprData-class
#' @rdname GsvaExprData-class
#' @exportClass GsvaExprData
setClassUnion("GsvaExprData",
              c("matrix", "dgCMatrix", "ExpressionSet",
                "SummarizedExperiment", "SingleCellExperiment", "SpatialExperiment", "HDF5Array"))


#' `GsvaGeneSets` class
#'
#' Virtual superclass of gene set classes supported by `GSVA`.
#'
#' `GSVA` supports gene sets in either a list of character vectors or an object
#' of class `GSEABase::GeneSetCollection`.  This class union allows to store any
#' of these in a slot of another class as well as defining common methods for
#' them.
#'
#' @seealso
#' [`list`], 
#' [`GeneSetCollection`]
#'
#' @importClassesFrom GSEABase GeneSetCollection GeneSet GeneIdentifierType
#'
#' @name GsvaGeneSets-class
#' @rdname GsvaGeneSets-class
#' @exportClass GsvaGeneSets
setClassUnion("GsvaGeneSets",
              c("list", "GeneSetCollection"))


#' `GsvaMethodParam` class
#'
#' Virtual superclass of method parameter classes supported by `GSVA`.
#'
#' `GSVA` implements four single-sample gene set analysis methods: PLAGE,
#' combined z-scores, ssGSEA, and GSVA.  All of them take at least an expression
#' data matrix and one or more gene sets as input.  Further common parameters
#' include an assay name for use with multi-assay expression data containers,
#' the gene ID type used by the expression data set, and a minimum and maximum
#' size for gene sets to limit the range of gene set sizes used in an analysis.
#' This virtual class provides the necessary slots for this shared parameter set
#' and serves as the parent class for all `GSVA` method parameter classes.
#'
#' @slot exprData The expression data set.  Must be one of the classes
#' supported by [`GsvaExprData-class`].  For a list of these classes, see its
#' help page using `help(GsvaExprData)`.
#'
#' @slot geneSets The gene sets.  Must be one of the classes supported by
#' [`GsvaGeneSets-class`].  For a list of these classes, see its help page using
#' `help(GsvaGeneSets)`.
#' 
#' @slot assay Character vector of length 1.  The name of the assay to use in
#' case `exprData` is a multi-assay container, otherwise ignored.  By default,
#' the first assay is used.
#' 
#' @slot annotation An object of class [`GeneIdentifierType-class`] from package
#' `GSEABase` describing the gene identifiers used as the row names of the
#' expression data set.  See [`GeneIdentifierType`] for help on available gene
#' identifier types and how to construct them.  This
#' information can be used to map gene identifiers occurring in the gene sets.
#' By default, this slot has value `NullIdentifier` and gene identifiers used in
#' expression data set and gene sets are matched directly.
#' 
#' @slot minSize Numeric vector of length 1.  Minimum size of the resulting gene
#' sets after gene identifier mapping. By default, the minimum size is 1.
#' 
#' @slot maxSize Numeric vector of length 1.  Maximum size of the resulting gene
#' sets after gene identifier mapping. By default, the maximum size is `Inf`.
#' 

#' @seealso
#' [`GsvaExprData-class`],
#' [`GsvaGeneSets-class`],
#' [`zscoreParam-class`], 
#' [`plageParam-class`], 
#' [`ssgseaParam-class`], 
#' [`gsvaParam-class`],
#' [`GeneIdentifierType-class`]
#' [`GeneIdentifierType`]
#'
#' @name GsvaMethodParam-class
#' @rdname GsvaMethodParam-class
#' @exportClass GsvaMethodParam
setClass("GsvaMethodParam",
         slots=c(exprData="GsvaExprData",
                 geneSets="GsvaGeneSets",
                 assay="character",
                 annotation="GeneIdentifierType",
                 minSize="numeric",
                 maxSize="numeric"),
         contains="VIRTUAL")


## ----- PLAGE Parameter Class -----

#' `plageParam` class
#'
#' S4 class for PLAGE method parameter objects.
#' 
#' Since method PLAGE does not take any method-specific parameters, this
#' class does not add any slots to the common slots inherited from
#' [`GsvaMethodParam-class`].
#'
#' @seealso
#' [`GsvaExprData-class`],
#' [`GsvaGeneSets-class`],
#' [`GsvaMethodParam-class`],
#' [`zscoreParam-class`],
#' [`ssgseaParam-class`],
#' [`gsvaParam-class`]
#'
#' @name plageParam-class
#' @rdname plageParam-class
#' @exportClass plageParam
setClass("plageParam",
         contains="GsvaMethodParam",
         prototype=list(exprData=NULL,
                        geneSets=NULL,
                        assay=NA_character_,
                        annotation=NULL,
                        minSize=NA_integer_,
                        maxSize=NA_integer_))


## ----- Combined z-Scores Parameter Class -----

#' `zscoreParam` class
#'
#' S4 class for combined z-scores method parameter objects.
#' 
#' Since the combined z-scores method does not take any method-specific
#' parameters, this class does not add any slots to the common slots inherited
#' from#' [`GsvaMethodParam-class`].
#'
#' @seealso
#' [`GsvaExprData-class`],
#' [`GsvaGeneSets-class`],
#' [`GsvaMethodParam-class`],
#' [`plageParam-class`],
#' [`ssgseaParam-class`],
#' [`gsvaParam-class`]
#'
#' @name zscoreParam-class
#' @rdname zscoreParam-class
#' @exportClass zscoreParam
setClass("zscoreParam",
         contains="GsvaMethodParam",
         prototype=list(exprData=NULL,
                        geneSets=NULL,
                        assay=NA_character_,
                        annotation=NULL,
                        minSize=NA_integer_,
                        maxSize=NA_integer_))



## ----- ssGSEA Parameter Class -----

#' `ssgseaParam` class
#'
#' S4 class for ssGSEA method parameter objects.
#'
#' In addition to the common parameter slots inherited from
#' `[GsvaMethodParam]`, this class has slots for the two method-specific
#' parameters of the `ssGSEA` method described below as well as four more slots
#' for implementing a missing value policy.
#' 
#' @slot alpha Numeric vector of length 1.  The exponent defining the
#' weight of the tail in the random walk performed by the ssGSEA (Barbie et
#' al., 2009) method.
#'
#' @slot normalize Logical vector of length 1.  If `TRUE` runs the ssGSEA
#' method from Barbie et al. (2009) normalizing the scores by the absolute
#' difference between the minimum and the maximum, as described in their paper.
#' Otherwise this final normalization step is skipped.
#'
#' @slot checkNA Character vector of length 1. One of the strings `"auto"`
#' (default), `"yes"`, or `"no"`, which refer to whether the input expression
#' data should be checked for the presence of missing (`NA`) values.
#'
#' @slot didCheckNA Logical vector of length 1, indicating whether the input
#' expression data was checked for the presence of missing (`NA`) values.
#'
#' @slot anyNA Logical vector of length 1, indicating whether the input
#' expression data contains missing (`NA`) values.
#'
#' @slot use Character vector of length 1. One of the strings `"everything"`
#' (default), `"all.obs"`, or `"na.rm"`, which refer to three different policies
#' to apply in the presence of missing values in the input expression data; see
#' [`ssgseaParam`].
#'
#' @seealso
#' [`GsvaExprData-class`],
#' [`GsvaGeneSets-class`],
#' [`GsvaMethodParam-class`],
#' [`plageParam-class`],
#' [`zscoreParam-class`],
#' [`gsvaParam-class`]
#'
#' @name ssgseaParam-class
#' @rdname ssgseaParam-class
#' @exportClass ssgseaParam
setClass("ssgseaParam",
         slots=c(alpha="numeric",
                 normalize="logical",
                 checkNA="character",
                 didCheckNA="logical",
                 anyNA="logical",
                 use="character"),
         contains="GsvaMethodParam",
         prototype=list(exprData=NULL,
                        geneSets=NULL,
                        assay=NA_character_,
                        annotation=NULL,
                        minSize=NA_integer_,
                        maxSize=NA_integer_,
                        alpha=NA_real_,
                        normalize=NA,
                        checkNA=NA_character_,
                        didCheckNA=NA,
                        anyNA=NA,
                        use=NA_character_))


## ----- GSVA Parameter Class -----

#' `gsvaParam` class
#'
#' S4 class for GSVA method parameter objects.
#'
#' In addition to the common parameter slots inherited from `[GsvaMethodParam]`,
#' this class has slots for the six method-specific parameters of the GSVA
#' method described below.
#'
#' @slot kcdf Character vector of length 1 denoting the kernel to use during
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
#' @slot kcdfNoneMinSampleSize Integer vector of length 1. When `kcdf="auto"`,
#' this parameter decides at what minimum sample size `kcdf="none"`, i.e., the
#' estimation of the empirical cumulative distribution function (ECDF) of
#' expression levels across samples is performed directly without using a
#' kernel; see the `kcdf` slot.
#'
#' @slot tau Numeric vector of length 1.  The exponent defining the weight of
#' the tail in the random walk performed by the GSVA (Hänzelmann et al., 2013)
#' method.
#'
#' @slot maxDiff Logical vector of length 1 which offers two approaches to
#' calculate the enrichment statistic (ES) from the KS random walk statistic.
#' * `FALSE`: ES is calculated as the maximum distance of the random walk from 0.
#' * `TRUE`: ES is calculated as the magnitude difference between
#' the largest positive and negative random walk deviations.
#'
#' @slot absRanking Logical vector of length 1 used only when `maxDiff=TRUE`.
#' When `absRanking=FALSE` a modified Kuiper statistic is used to calculate
#' enrichment scores, taking the magnitude difference between the largest
#' positive and negative random walk deviations. When `absRanking=TRUE` the
#' original Kuiper statistic that sums the largest positive and negative
#' random walk deviations, is used. In this latter case, gene sets with genes
#' enriched on either extreme (high or low) will be regarded as ’highly’
#' activated.
#'
#' @slot sparse Logical vector of length 1 used only when the input expression
#' data in `exprData` is stored in a sparse matrix (e.g., a `dgCMatrix` or a
#' container object, such as a `SingleCellExperiment`, storing the expression
#' data in a `dgCMatrix`).
#' In such a case, when `sparse=TRUE`, a sparse version of the GSVA algorithm
#' will be applied. Otherwise, when `sparse=FALSE`, the classical version of
#' the GSVA algorithm will be used.
#'
#' @seealso
#' [`GsvaExprData-class`],
#' [`GsvaGeneSets-class`],
#' [`GsvaMethodParam-class`],
#' [`plageParam-class`],
#' [`zscoreParam-class`],
#' [`ssgseaParam-class`]
#'
#' @name gsvaParam-class
#' @rdname gsvaParam-class
#' @exportClass gsvaParam
setClass("gsvaParam",
         slots=c(kcdf="character",
                 kcdfNoneMinSampleSize="integer",
                 tau="numeric", 
                 maxDiff="logical",
                 absRanking="logical",
                 sparse="logical"),
         contains="GsvaMethodParam",
         prototype=list(exprData=NULL,
                        geneSets=NULL,
                        assay=NA_character_,
                        annotation=NULL,
                        minSize=NA_integer_,
                        maxSize=NA_integer_,
                        kcdf=NA_character_,
                        kcdfNoneMinSampleSize=NA_integer_,
                        tau=NA_real_,
                        maxDiff=NA,
                        absRanking=NA,
                        sparse=FALSE))
