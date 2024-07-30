
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
#' @importClassesFrom GSEABase GeneSetCollection
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
#' data matrix and one or many gene sets as input.  This virtual class provides
#' the necessary slots for this minimum parameter set and serves as all `GSVA`
#' method parameter classes,
#'
#' @seealso
#' [`GsvaExprData-class`],
#' [`GsvaGeneSets-class`],
#' [`zscoreParam-class`], 
#' [`plageParam-class`], 
#' [`ssgseaParam-class`], 
#' [`gsvaParam-class`]
#'
#' @name GsvaMethodParam-class
#' @rdname GsvaMethodParam-class
#' @exportClass GsvaMethodParam
setClass("GsvaMethodParam",
         slots=c(exprData="GsvaExprData",
                 geneSets="GsvaGeneSets",
                 assay="character",
                 annotation="character",
                 minSize="numeric",
                 maxSize="numeric"),
         contains="VIRTUAL")


## ----- PLAGE Parameter Class -----

#' `plageParam` class
#'
#' Method-specific parameters for the PLAGE method.
#' 
#' Since this method does not take any method-specific parameters, the parameter
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
                        annotation=NA_character_,
                        minSize=NA_integer_,
                        maxSize=NA_integer_))


## ----- Combined z-Scores Parameter Class -----

#' `zscoreParam` class
#'
#' Method-specific parameters for the combined z-scores method.
#' 
#' Since this method does not take any method-specific parameters, the parameter
#' class does not add any slots to the common slots inherited from
#' [`GsvaMethodParam-class`].
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
                        annotation=NA_character_,
                        minSize=NA_integer_,
                        maxSize=NA_integer_))



## ----- ssGSEA Parameter Class -----

#' `ssgseaParam` class
#'
#' Method-specific parameters for the ssGSEA method.
#'
#' In addition to the two common parameter slots inherited from
#' `[GsvaMethodParam]`, this class has slots for the two method-specific
#' parameters of the `ssGSEA` method described below.
#' 
#' @slot alpha Numeric vector of length 1.  The exponent defining the
#' weight of the tail in the random walk performed by the ssGSEA (Barbie et
#' al., 2009) method.
#'
#' @slot normalize Logical vector of length 1.  If `TRUE` runs the ssGSEA
#' method from Barbie et al. (2009) normalizing the scores by the absolute
#' difference between the minimum and the maximum, as described in their paper.
#' Otherwise this last normalization step is skipped.
#'
#' @slot checkNA Character string. One of the strings `"auto"` (default),
#' `"yes"`, or `"no"`, which refer to whether the input expression data should
#' be checked for the presence of missing (`NA`) values.
#'
#' @slot didCheckNA Logical vector of length 1, indicating whether the input
#' expression data was checked for the presence of missing (`NA`) values.
#'
#' @slot anyNA Logical vector of length 1, indicating whether the input
#' expression data contains missing (`NA`) values.
#'
#' @slot use Character string. One of the strings `"everything"` (default),
#' `"all.obs"`, or `"na.rm"`, which refer to three different policies to apply
#' in the presence of missing values in the input expression data; see
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
                        annotation=NA_character_,
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
#' Method-specific parameters for the GSVA method.
#'
#' In addition to the two common parameter slots inherited from
#' `[GsvaMethodParam]`, this class has slots for the two method-specific
#' parameters of the GSVA method described below.
#'
#' @slot kcdf Character vector of length 1 denoting the kernel to use during
#' the non-parametric estimation of the cumulative distribution function of
#' expression levels across samples. `kcdf="Gaussian"` is suitable when input
#' expression values are continuous, such as microarray fluorescent units in
#' logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input
#' expression values are integer counts, such as those derived from RNA-seq
#' experiments, then this argument should be set to `kcdf="Poisson"`.
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
#' @slot absRanking Logical vector of length 1 used only when `mx.diff=TRUE`.
#' When `abs.ranking=FALSE` a modified Kuiper statistic is used to calculate
#' enrichment scores, taking the magnitude difference between the largest
#' positive and negative random walk deviations. When `abs.ranking=TRUE` the
#' original Kuiper statistic that sums the largest positive and negative
#' random walk deviations, is used. In this latter case, gene sets with genes
#' enriched on either extreme (high or low) will be regarded as ’highly’
#' activated.
#'
#' @slot sparse Logical vector of length 1 used only when the input expression
#' data in `exprData` is stored in a sparse matrix (e.g., a `dgCMatrix` or a
#' `singleCellExperiment` object storing the expression data in a `dgCMatrix`).
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
                 tau="numeric", 
                 maxDiff="logical",
                 absRanking="logical",
                 sparse="logical"),
         contains="GsvaMethodParam",
         prototype=list(exprData=NULL,
                        geneSets=NULL,
                        assay=NA_character_,
                        annotation=NA_character_,
                        minSize=NA_integer_,
                        maxSize=NA_integer_,
                        kcdf=NA_character_,
                        tau=NA_real_,
                        maxDiff=NA,
                        absRanking=NA,
                        sparse=FALSE))
