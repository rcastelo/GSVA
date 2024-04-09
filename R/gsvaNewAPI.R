
#' @title Gene Set Variation Analysis
#' 
#' @description Estimates GSVA enrichment scores. The API of this function has
#' changed in the Bioconductor release 3.18 and this help page describes the
#' new API. The old API is defunct and will be removed in the next
#' Bioconductor release. If you are looking for the documentation of the old
#' API to the `gsva()` function, please consult [`GSVA-pkg-defunct`].
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
#' @param BPPARAM An object of class [`BiocParallelParam`] specifying parameters
#'   related to the parallel execution of some of the tasks and calculations
#'   within this function.
#' 
#' @return A gene-set by sample matrix (of `matrix` or [`dgCMatrix-class`] type, 
#'   depending on the input) of GSVA enrichment scores.
#' 
#' @seealso [`plageParam`], [`zscoreParam`], [`ssgseaParam`], [`gsvaParam`]
#'
#' @aliases gsva
#' @name gsva
#' @rdname gsva
#' 
#' @references Barbie, D.A. et al. Systematic RNA interference reveals that
#' oncogenic KRAS-driven cancers require TBK1.
#' *Nature*, 462(5):108-112, 2009.
#' [DOI](https://doi.org/10.1038/nature08460)
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' [DOI](https://doi.org/10.1186/1471-2105-14-7)
#'
#' @references Lee, E. et al. Inferring pathway activity toward precise
#' disease classification.
#' *PLoS Comp Biol*, 4(11):e1000217, 2008.
#' [DOI](https://doi.org/10.1371/journal.pcbi.1000217)
#'
#' @references Tomfohr, J. et al. Pathway level analysis of gene expression
#' using singular value decomposition.
#' *BMC Bioinformatics*, 6:225, 2005.
#' [DOI](https://doi.org/10.1186/1471-2105-6-225)
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
#' gsvapar <- gsvaParam(y, geneSets, maxDiff=TRUE)
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

#' @aliases gsva,plageParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="plageParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              famGaGS <- .filterAndMapGenesAndGeneSets(param)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose)
                  cat(sprintf("Setting parallel calculations through a %s back-end\n",
                              "with workers=%d and tasks=100.\n",
                              class(BPPARAM), bpnworkers(BPPARAM)))

              ## if(rnaseq)
              ##     stop("rnaseq=TRUE does not work with method='plage'.")

              if(verbose)
                  cat("Estimating PLAGE scores for",
                      length(filteredMappedGeneSets), "gene sets.\n")

              plageScores <- plage(X=filteredDataMatrix,
                                   geneSets=filteredMappedGeneSets,
                                   verbose=verbose,
                                   BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), plageScores, gs)

              return(rval)
          })


#' @aliases gsva,zscoreParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="zscoreParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              famGaGS <- .filterAndMapGenesAndGeneSets(param)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose)
                  cat(sprintf("Setting parallel calculations through a %s back-end\n",
                              "with workers=%d and tasks=100.\n",
                              class(BPPARAM), bpnworkers(BPPARAM)))

              ## if (rnaseq)
              ##     stop("rnaseq=TRUE does not work with method='zscore'.")

              if(verbose)
                  cat("Estimating combined z-scores for",
                      length(filteredMappedGeneSets), "gene sets.\n")

              zScores <- zscore(X=filteredDataMatrix,
                                geneSets=filteredMappedGeneSets,
                                verbose=verbose,
                                BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), zScores, gs)
              
              return(rval)
          })


#' @aliases gsva,ssgseaParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="ssgseaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=FALSE,
                                                       removeNzConstant=FALSE)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose)
                  cat(sprintf("Setting parallel calculations through a %s back-end\n",
                              "with workers=%d and tasks=100.\n",
                              class(BPPARAM), bpnworkers(BPPARAM)))

              if(verbose)
                  cat("Estimating ssGSEA scores for",
                      length(filteredMappedGeneSets), "gene sets.\n")

              ssgseaScores <- ssgsea(X=filteredDataMatrix,
                                     geneSets=filteredMappedGeneSets,
                                     alpha=get_alpha(param), 
                                     normalization=do_normalize(param),
                                     verbose=verbose,
                                     BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), ssgseaScores, gs)
              
              return(rval)
          })


#' @aliases gsva,gsvaParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              famGaGS <- .filterAndMapGenesAndGeneSets(param)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose)
                  cat(sprintf("Setting parallel calculations through a %s back-end\n",
                              "with workers=%d and tasks=100.\n",
                              class(BPPARAM), bpnworkers(BPPARAM)))

              if(verbose)
                  cat("Estimating GSVA scores for",
                      length(filteredMappedGeneSets),
                      "gene sets.\n")
              
              nSamples <- ncol(filteredDataMatrix)
              ## nGenes <- nrow(filteredDataMatrix)
              nGeneSets <- length(filteredMappedGeneSets)

              if (get_kcdf(param) == "Gaussian") {
                  rnaseq <- FALSE
                  kernel <- TRUE
              } else if (get_kcdf(param) == "Poisson") {
                  rnaseq <- TRUE
                  kernel <- TRUE
              } else
                  kernel <- FALSE
              
              gsvaScores <- compute.geneset.es(expr=filteredDataMatrix,
                                               gset.idx.list=filteredMappedGeneSets,
                                               sample.idxs=seq.int(nSamples),
                                               rnaseq=rnaseq,
                                               abs.ranking=get_absRanking(param),
                                               parallel.sz=if(inherits(BPPARAM, "SerialParam")) 1L else bpnworkers(BPPARAM),
                                               mx.diff=get_maxDiff(param),
                                               tau=get_tau(param),
                                               kernel=kernel,
                                               verbose=verbose,
                                               BPPARAM=BPPARAM)
              
              colnames(gsvaScores) <- colnames(filteredDataMatrix)
              rownames(gsvaScores) <- names(filteredMappedGeneSets)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), gsvaScores, gs)
              
              return(rval)
          })


### ----- methods for retrieving gene sets -----

#' @title Retrieve or Determine Gene Sets
#' 
#' @description Retrieves or determines the gene sets that have been used
#' or would be used in a `gsva()` gene set analysis.  These are not necessarily
#' the same as the input gene sets.  See Details.
#' 
#' @param obj An object of one of the following classes:
#' * An expression data object of one of the classes described in
#' [`GsvaExprData-class`] that is the return value of a call to `gsva()`.
#' * A parameter object of one of the classes described in
#' [`GsvaMethodParam-class`] that could be used in a call to `gsva()`.
#'
#' @return The `geneSets()` methods return a named list of character vectors
#' where each character vector contains the gene IDs of a gene set.
#' The `geneSetSizes()` methods return a named integer vector of gene set sizes.
#' 
#' @details The gene sets used in a `gsva()` gene set analysis, or just their
#' sizes, may be a valuable input to subsequent analyses.  However, they are not
#' necessarily the same as the original input gene sets, or their sizes: based
#' on user choices, the gene annotation used, or presence/absence of genes in
#' gene sets and expression data set, `gsva()` may have to modify them during
#' the preparation of an analysis run.
#' In order to make use of these gene sets or their sizes, you can either
#' * retrieve them from the object returned by `gsva()` by passing this object
#' to `geneSets()` or `geneSetSizes()`, or
#' * predict them by calling `geneSets()` or `geneSetSizes()` on the parameter
#' object that would also be passed to `gsva()`.  This is much slower and should
#' only be done if you do not intend to run an actual gene set analysis.
#'
#' `geneSetSizes()` is a convenience wrapper running `lengths()` on the list of
#' gene sets returned by `geneSets()`.
#'
#' @aliases geneSets geneSetSizes
#' @name geneSets
#' @rdname geneSets
#' 
NULL


#' @aliases geneSets,GsvaMethodParam-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("GsvaMethodParam"),
          function(obj) {
              famGaGS <- .filterAndMapGenesAndGeneSets(obj)

              return(.geneSetsIndices2Names(
                  indices=famGaGS[["filteredMappedGeneSets"]],
                  names=rownames(famGaGS[["filteredDataMatrix"]])
              ))
          })

#' @aliases geneSets,SummarizedExperiment-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("SummarizedExperiment"),
          function(obj) {
              return(as(rowData(obj)$gs, "list"))
          })

#' @aliases geneSets,SingleCellExperiment-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("SingleCellExperiment"),
          function(obj) {
              return(as(rowData(obj)$gs, "list"))
          })

#' @aliases geneSets,SpatialExperiment-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("SpatialExperiment"),
          function(obj) {
              return(as(rowData(obj)$gs, "list"))
          })

#' @aliases geneSets,GsvaExprData-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("GsvaExprData"),
          function(obj) {
              return(.geneSets(obj))
          })


#' @aliases geneSetSizes,GsvaMethodParam-method
#' @rdname geneSets
#' @exportMethod geneSetSizes
setMethod("geneSetSizes", signature("GsvaMethodParam"),
          function(obj) {
              return(lengths(geneSets(obj)))
          })

#' @aliases geneSetSizes,GsvaExprData-method
#' @rdname geneSets
#' @exportMethod geneSetSizes
setMethod("geneSetSizes", signature("GsvaExprData"),
          function(obj) {
              return(lengths(geneSets(obj)))
          })


### ----- methods for data pre-/post-processing -----

## unwrapData: extract a data matrix from a container object
setMethod("unwrapData", signature("matrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("dgCMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("ExpressionSet"),
          function(container, assay) {
              return(exprs(container))
          })

setMethod("unwrapData", signature("SummarizedExperiment"),
          function(container, assay) {
              if (length(assays(container)) == 0L)
                  stop("The input SummarizedExperiment object has no assay data.")

              if (missing(assay) || is.na(assay)) {
                  assay <- names(assays(container))[1]
              } else {
                  if (!is.character(assay))
                      stop("The 'assay' argument must contain a character string.")

                  assay <- assay[1]

                  if (!assay %in% names(assays(container)))
                      stop(sprintf("Assay %s not found in the input SummarizedExperiment object.", assay))
              }

              return(assays(container)[[assay]])
          })

setMethod("unwrapData", signature("SingleCellExperiment"),
          function(container, assay) {
              if (length(assays(container)) == 0L)
                  stop("The input SingleCellExperiment object has no assay data.")

              if (missing(assay) || is.na(assay)) {
                  assay <- names(assays(container))[1]
              } else {
                  if (!is.character(assay))
                      stop("The 'assay' argument must contain a character string.")

                  assay <- assay[1]

                  if (!assay %in% names(assays(container)))
                      stop(sprintf("Assay %s not found in the input SingleCellExperiment object.", assay))
              }

              return(assays(container)[[assay]])
          })

setMethod("unwrapData", signature("SpatialExperiment"),
          function(container, assay) {
            if (length(assays(container)) == 0L)
              stop("The input SpatialExperiment object has no assay data.")
            
            if (missing(assay) || is.na(assay)) {
              assay <- names(assays(container))[1]
            } else {
              if (!is.character(assay))
                stop("The 'assay' argument must contain a character string.")
              
              assay <- assay[1]
              
              if (!assay %in% names(assays(container)))
                stop(sprintf("Assay %s not found in the input SpatialExperiment object.", assay))
            }
            
            return(assays(container)[[assay]])
          })


## wrapData: put the resulting data and gene sets into the original data container type
setMethod("wrapData", signature(container="matrix"),
          function(container, dataMatrix, geneSets) {
              attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="dgCMatrix"),
          function(container, dataMatrix, geneSets) {
              attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="ExpressionSet"),
          function(container, dataMatrix, geneSets) {
              rval <- new("ExpressionSet", exprs=dataMatrix,
                          phenoData=phenoData(container),
                          experimentData=experimentData(container),
                          annotation="")
              attr(rval, "geneSets") <- geneSets
              
              return(rval)
          })

setMethod("wrapData", signature(container="SummarizedExperiment"),
          function(container, dataMatrix, geneSets) {
              rval <- SummarizedExperiment(
                  assays=SimpleList(es=dataMatrix),
                  colData=colData(container),
                  rowData=DataFrame(gs=CharacterList(geneSets)),
                  metadata=metadata(container))
              metadata(rval)$annotation <- NULL

              return(rval)
          })

setMethod("wrapData", signature(container="SingleCellExperiment"),
          function(container, dataMatrix, geneSets) {
              rval <- SingleCellExperiment(
                  assays=SimpleList(es=dataMatrix),
                  colData=colData(container),
                  rowData=DataFrame(gs=CharacterList(geneSets)),
                  metadata=metadata(container))
              metadata(rval)$annotation <- NULL
              
              return(rval)
          })

setMethod("wrapData", signature(container="SpatialExperiment"),
          function(container, dataMatrix, geneSets) {
              rval <- SpatialExperiment(
                  assays=SimpleList(es=dataMatrix),
                  colData=colData(container),
                  rowData=DataFrame(gs=CharacterList(geneSets)),
                  metadata=metadata(container),
                  imgData = imgData(container),
                  spatialCoords = spatialCoords(container))
              metadata(rval)$annotation <- NULL
              
              return(rval)
          })


## mapGeneSetsToAnno: translate feature IDs used in gene sets to specified
##                    annotation type (if any, and if possible)
setMethod("mapGeneSetsToAnno", signature("list"),
          function(geneSets, anno) {
              return(geneSets)
          })

setMethod("mapGeneSetsToAnno", signature("GeneSetCollection"),
          function(geneSets, anno) {
              if(.isAnnoPkgValid(anno)) {
                  if(!.isAnnoPkgInstalled(anno))
                      stop(sprintf("Please install the annotation package %s. If %s does not seem to exist as a package, please try to append the suffix .db to its name.", anno, anno))
                  ## TODO: provide a check for verbosity
                  cat("Mapping identifiers between gene sets and feature names\n")

                  ## map gene identifiers of the gene sets to the features in the chip
                  mappedGeneSets <- mapIdentifiers(geneSets,
                                                   AnnoOrEntrezIdentifier(anno))
                  rval <- geneIds(mappedGeneSets)

              } else {
                  ## TODO: provide a check for verbosity
                  cat("No annotation package name available in the input data object.",
                      "Attempting to directly match identifiers in data to gene sets.", sep="\n")

                  rval <- geneIds(geneSets)
              }

              return(rval)
          })

