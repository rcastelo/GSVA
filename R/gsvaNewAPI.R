
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
#' @return A gene-set by sample matrix of GSVA enrichment scores stored in a
#' container object of the same type as the input expression data container. If
#' the input was a base matrix or a [`dgCMatrix-class`] object, then the output will
#' be a base matrix object with the gene sets employed in the calculations
#' stored in an attribute called `geneSets`. If the input was an
#' [`ExpressionSet`] object, then the output will be also an [`ExpressionSet`]
#' object with the gene sets employed in the calculations stored in an
#' attributed called `geneSets`. If the input was an object of one of the
#' classes described in [`GsvaExprData`], such as a [`SingleCellExperiment`],
#' then the output will be of the same class, where enrichment scores will be
#' stored in an assay called `es` and the gene sets employed in the
#' calculations will be stored in the `rowData` slot of the object under the
#' column name `gs`.
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

#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @aliases gsva,plageParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="plageParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       verbose)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if(verbose)
                  cli_alert_info(sprintf("Calculating PLAGE scores for %d gene sets",
                                         length(filteredMappedGeneSets)))

              plageScores <- plage(X=filteredDataMatrix,
                                   geneSets=filteredMappedGeneSets,
                                   verbose=verbose,
                                   BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), plageScores, gs)

              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @aliases gsva,zscoreParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="zscoreParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       verbose)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              ## if (rnaseq)
              ##     stop("rnaseq=TRUE does not work with method='zscore'.")

              if(verbose)
                  cli_alert_info(sprintf("Calculating Z-scores for %d gene sets",
                                         length(filteredMappedGeneSets)))

              zScores <- zscore(X=filteredDataMatrix,
                                geneSets=filteredMappedGeneSets,
                                verbose=verbose,
                                BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), zScores, gs)
              
              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })

#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @aliases gsva,ssgseaParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="ssgseaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=FALSE,
                                                       removeNzConstant=FALSE,
                                                       verbose)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if(verbose)
                  cli_alert_info(sprintf("Calculating  ssGSEA scores for %d gene sets",
                                         length(filteredMappedGeneSets)))

              ssgseaScores <- ssgsea(X=filteredDataMatrix,
                                     geneSets=filteredMappedGeneSets,
                                     alpha=get_alpha(param), 
                                     normalization=do_normalize(param),
                                     any_na=anyNA(param),
                                     na_use=na_use(param),
                                     minSize=get_minSize(param),
                                     BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), ssgseaScores, gs)
              
              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @aliases gsva,gsvaParam-method
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       verbose)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if(verbose)
                  cli_alert_info(sprintf("Calculating GSVA scores for %d gene sets",
                                         length(filteredMappedGeneSets)))
              
              gsvaScores <- compute.geneset.es(expr=filteredDataMatrix,
                                               gset.idx.list=filteredMappedGeneSets,
                                               sample.idxs=seq.int(ncol(filteredDataMatrix)),
                                               kcdf=get_kcdf(param),
                                               kcdf.min.ssize=get_kcdfNoneMinSampleSize(param),
                                               abs.ranking=get_absRanking(param),
                                               parallel.sz=if(inherits(BPPARAM, "SerialParam")) 1L else bpnworkers(BPPARAM),
                                               mx.diff=get_maxDiff(param),
                                               tau=get_tau(param),
                                               sparse=get_sparse(param),
                                               verbose=verbose,
                                               BPPARAM=BPPARAM)
              
              colnames(gsvaScores) <- colnames(filteredDataMatrix)
              rownames(gsvaScores) <- names(filteredMappedGeneSets)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), gsvaScores, gs)

              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


### ----- methods for retrieving and setting annotation metadata -----

#' @title Store and Retrieve Annotation Metadata
#' 
#' @description Methods for storing and retrieving annotation metadata in
#' expression data objects that support it.  If gene sets and expression data
#' are using different but known gene identifier types and an appropriate
#' annotation database is available, gene set identifiers can be mapped to
#' expression data identifiers without manual user intervention, e.g. from
#' an MSigDb gene set using ENTREZ IDs or gene symbols to an expression data
#' set using ENSEMBL IDs.
#' 
#' @param object An expression data object of one of the classes described in
#' [`GsvaExprData-class`].  Simple `matrix` and `dgCMatrix` objects are not
#' capable of storing annotation metadata and will return `NULL`.
#'
#' @param value For the replacement methods, the annotation metadata to be
#' stored in the object.  For [`ExpressionSet-class`] objects, this must be a
#' character of length 1 specifying the name of the annotation database to be
#' used.  For [`SummarizedExperiment-class`] and its subclasses, this must be
#' a [`GeneIdentifierType`] created by one of the constructors from package
#' `GSEABase` where the `annotation` argument is typically the name of an
#' organism or annotation database, e.g. `org.Hs.eg.db`.  Simple `matrix` and
#' `dgCMatrix` objects are not capable of storing annotation metadata and the
#' attempt to do so will result in an error.
#'
#' @return For the retrieval methods, the annotation metadata stored in the
#' object of `NULL`.  For the replacement methods, the updated object.
#'
#' @aliases gsvaAnnotation gsvaAnnotation<-
#' @name gsvaAnnotation
#' @rdname gsvaAnnotation
#' 
NULL


#' @aliases gsvaAnnotation,GsvaExprData-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="GsvaExprData"),
          function(object) {
              return(attr(object, which="geneIdType", exact=TRUE))
          })

#' @aliases gsvaAnnotation<-,GsvaExprData,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                     object="GsvaExprData",
                     value="GeneIdentifierType"),
          function(object, value) {
              attr(object, which="geneIdType") <- value
              object
          })

#' @aliases gsvaAnnotation,ExpressionSet-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="ExpressionSet"),
          function(object) {
              ## always a character giving the db pkg, potentially empty ("")
              ## unfortunately, as it turns out, even character(0) sometimes
              ao <- annotation(object)
              if(.isCharLength1(ao)) {
                  return(AnnoOrEntrezIdentifier(ao))
              } else {
                  return(NULL)
              }
          })

#' @aliases gsvaAnnotation<-,ExpressionSet,character-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="ExpressionSet",
                   value="character"),
                 function(object, value) {
                     annotation(object) <- value
                     object
                 })

#' @aliases gsvaAnnotation<-,ExpressionSet,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="ExpressionSet",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     gsvaAnnotation(object) <- annotation(value)
                     object
                 })

#' @aliases gsvaAnnotation,SummarizedExperiment-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation", signature("SummarizedExperiment"),
          function(object) {
              ## NULL if unset; otherwise anything but we *expect* and handle
              ## a GSEABase::GeneIdentifierType with or without annotation(),
              ## i.e., db pkg, available.  Same for subclasses below.
              return(metadata(object)$annotation)
          })

#' @aliases gsvaAnnotation<-,SummarizedExperiment,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="SummarizedExperiment",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     metadata(object)$annotation <- value
                     object
                 })

#' @aliases gsvaAnnotation,SingleCellExperiment-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation", signature("SingleCellExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })

#' @aliases gsvaAnnotation<-,SingleCellExperiment,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="SingleCellExperiment",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     metadata(object)$annotation <- value
                     object
                 })

#' @aliases gsvaAnnotation,SpatialExperiment-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation", signature("SpatialExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })

#' @aliases gsvaAnnotation<-,SpatialExperiment,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="SpatialExperiment",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     metadata(object)$annotation <- value
                     object
                 })


#' @aliases gsvaAnnotation,list-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="list"),
          function(object) {
              return(attr(object, which="geneIdType", exact=TRUE))
          })

#' @aliases gsvaAnnotation<-,list,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                     object="list",
                     value="GeneIdentifierType"),
          function(object, value) {
              attr(object, which="geneIdType") <- value
              object
          })

#' @aliases gsvaAnnotation,GeneSetCollection-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="GeneSetCollection"),
          function(object) {
              lgit <- unique(lapply(object, geneIdType))
              return(if(length(lgit) == 1) lgit[[1]] else NULL)
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


### ----- helper functions for gene set I/O and preprocessing -----

#' @title Handling of Duplicated Gene Set Names
#' 
#' @description Offers a choice of ways for handling duplicated gene set names
#' that may not be suitable as input to other gene set analysis functions.
#' 
#' @param geneSets A named list of gene sets represented as character vectors
#' of gene IDs as e.g. returned by [`readGMT`].
#'
#' @param deduplUse A character vector of length 1 specifying one of several
#' methods to handle duplicated gene set names.
#' Duplicated gene set names are explicitly forbidden by the
#' [GMT file format specification](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
#' but can nevertheless be encountered in the wild.
#' The available choices are:
#' * `first` (the default): drops all gene sets whose names are [`duplicated`]
#' according to the base R function and retains only the first occurence of a
#' gene set name.
#' * `drop`:  removes *all* gene sets that have a duplicated name, including its
#' first occurrence.
#' * `union`: replaces gene sets with duplicated names by a single gene set
#' containing the union of all their gene IDs.
#' * `smallest`: drops gene sets with duplicated names and retains only the
#' smallest of them, i.e. the one with the fewest gene IDs.  If there are
#' several smallest gene sets, the first will be selected.
#' * `largest`: drops gene sets with duplicated names and retains only the
#' largest of them, i.e. the one with the most gene IDs.  If there are
#' several largest gene sets, the first will be selected.
#'
#' @return A named list of gene sets represented as character vectors of
#' gene IDs.
#' 
#' @aliases deduplicateGeneSets
#' @name deduplicateGeneSets
#' @rdname deduplicateGeneSets
#' @export
#' 
deduplicateGeneSets <- function(geneSets,
                                deduplUse = c("first", "drop", "union",
                                              "smallest", "largest")) {
    ddUse <- match.arg(deduplUse)
    isNameDuplicated <- duplicated(names(geneSets))
    duplicatedNames <- unique(names(geneSets[isNameDuplicated]))

    ## a nested list containing sublists of duplicated gene sets
    duplicatedGeneSets <- sapply(duplicatedNames,
                                 function(dn, gs) unname(gs[dn == names(gs)]),
                                 gs = geneSets, simplify=FALSE)

    ## transformation function operating on sublists of such nested lists,
    ## returning a single deduplicated gene set, i.e. character vector
    ddFunc <- switch(ddUse,
                     union=function(dgs) Reduce(union, dgs),
                     smallest=function(dgs) dgs[[which.min(lengths(dgs))]],
                     largest=function(dgs) dgs[[which.max(lengths(dgs))]])

    ## apply transformation function to deduplicate gene sets (if requested)
    if(!is.null(ddFunc))
        dedupl <- sapply(duplicatedGeneSets, FUN=ddFunc, simplify=FALSE)

    ## drop all duplicate gene sets (sufficient for default of "first")
    geneSets[isNameDuplicated] <- NULL

    ## remove or replace non-duplicated with deduplicated gene sets
    if(ddUse == "drop") {
        geneSets[duplicatedNames] <- NULL
    } else if(!is.null(ddFunc)) {
        geneSets[duplicatedNames] <- dedupl
    }

    return(geneSets)
}


deduplicateGmtLines <- function(geneSets,
                                deduplUse = c("first", "drop", "union",
                                              "smallest", "largest")) {
    ddUse <- match.arg(deduplUse)
    gsName <- sapply(geneSets, head, 1)
    isNameDuplicated <- which(duplicated(gsName))

    if(length(isNameDuplicated) > 0) {
        warning("GMT contains duplicated gene set names; deduplicated",
                " using method: ", ddUse)
        duplicatedNames <- unique(gsName[isNameDuplicated])
        lIdxDuplGS <- lapply(duplicatedNames,
                             function(DN, GSN) which(DN==GSN),
                             GSN = gsName)
        idxReplace <- sapply(lIdxDuplGS, head, 1)
        idxRemove <- unique(unlist(lapply(lIdxDuplGS, tail, -1)))

        ddFunc <- switch(ddUse,
                         union = function(idxGS, lGS) {
                             ld <- lGS[idxGS]
                             gsn <- lapply(ld, head, 1)[[1]]
                             gsd <- do.call("paste", c(lapply(ld, "[", 2), sep = " | "))
                             gsg <- Reduce(union, lapply(ld, tail, -2))
                             c(gsn, gsd, gsg)
                         },
                         smallest = function(idxGS, lGS) {
                             lGS[idxGS][[which.min(lengths(lGS[idxGS]))]]
                         },
                         largest = function(idxGS, lGS) {
                             lGS[idxGS][[which.max(lengths(lGS[idxGS]))]]
                         })

        if(!is.null(ddFunc)) {
            dedupl <- lapply(lIdxDuplGS, FUN = ddFunc, lGS = geneSets)
            geneSets[idxReplace] <- dedupl
        } else if(ddUse == "drop") {
            idxRemove <- union(idxRemove, idxReplace)
        }
        
        geneSets[idxRemove] <- NULL
    }

    return(geneSets)
}


#' @title Guess the gene identifier type from a list of character vectors
#' 
#' @description This function tries to derive the type of gene IDs used in a
#' named list of `character` vectors provided as input.
#' 
#' @param geneIdsList A named list of character vectors like the ones returned
#' by `geneIds()`.
#'
#' @return An object of a subclass of [`GeneIdentifierType`] derived from the
#' input.
#'
#' @details In order to make this function useful and keep it as simple as
#' possible, we limit ourselves to the most common types of gene identifiers:
#' "Gene IDs" consisting of digits only are considered ENTREZ IDs, anything
#' starting with 'ENS' an ENSEMBL identifier and anything else a HuGO gene
#' symbol.
#' 
#' @seealso [`GeneIdentifierType`]
#'
#' @aliases guessGeneIdType
#' @name guessGeneIdType
#' @rdname guessGeneIdType
#' @export
#' 
guessGeneIdType <- function(geneIdsList) {
    allIds <- unlist(geneIdsList)
    
    if(all(grepl("^[[:digit:]]+$", allIds))) {
        retVal <- EntrezIdentifier()
    } else if(all(grepl("^ENS", allIds))) {
        retVal <- ENSEMBLIdentifier()
    } else {
        retVal <- SymbolIdentifier()
    }

    return(retVal)
}


#' @title Construct a GeneSetCollection object from a list of character vectors
#' 
#' @description This function is essentially the reverse of
#' `GSEABase::geneIds()`, i.e., it takes as input a named list of `character`
#' vectors representing gene sets and returns the corresponding
#' GeneSetCollection object.
#' 
#' @param geneIdsList A named list of character vectors like the ones returned
#' by `geneIds()`.  Names must be unique; otherwise see `deduplicateGeneSets()`
#' for a number of strategies to resolve this issue.  
#'
#' @param geneIdType By default a character vector of length 1 with the special
#' value `"auto"` or an object of a subclass of [`GeneIdentifierType`].  If set
#' to `"auto"`, the function will try to derive the gene ID type from argument
#' `geneIdsList` using [`guessGeneIdType`].
#' The gene ID type of all `GeneSet` objects in the resulting
#' `GeneSetCollection` will be set to this value.
#' 
#' @param collectionType An object of class [`CollectionType`].  The collection
#' type of all `GeneSet` objects in the resulting `GeneSetCollection` will be
#' set to this value but can afterwards be modified for individual `GeneSet`s
#' if necessary.
#'
#' @return An object of class [`GeneSetCollection`] with all its [`GeneSet`]
#' objects using the gene ID and collection types specified by the corresponding
#' arguments.  Applying function `geneIds()` to this object should return a list
#' identical to the `geneIdsList` argument.
#' 
#' @seealso [`GeneSetCollection`], [`geneIds`], [`deduplicateGeneSets`],
#' [`guessGeneIdType`], [`GeneSet`]
#'
#' @aliases geneIdsToGeneSetCollection
#' @name geneIdsToGeneSetCollection
#' @rdname geneIdsToGeneSetCollection
#' @export
#' 
geneIdsToGeneSetCollection <- function(geneIdsList,
                                       geneIdType="auto",
                                       collectionType=NullCollection()) {
    if(geneIdType == "auto") {
        if(is.null(git <- gsvaAnnotation(geneIdsList))) {
            git <- guessGeneIdType(geneIdsList)
        }
    } else {
        git <- geneIdType
    }
    
    return(GeneSetCollection(mapply(function(gn, gs) {
        if(anyDuplicated(gs) > 0) {
            gs <- unique(gs)
            msg <- sprintf("Duplicated gene IDs removed from gene set %s", gn)
            cli_alert_warning(msg)
        }
        
        GeneSet(gs,
                geneIdType=git,
                collectionType=collectionType,
                setName=gn)
    }, gn=names(geneIdsList), gs=geneIdsList)))
}


#' @title Import Gene Sets from a GMT File
#' 
#' @description Imports a list of gene sets from a GMT (Gene Matrix Transposed)
#' format file, offering a choice of ways to handle duplicated gene set names.
#' 
#' @param con A connection object or character string containing e.g.
#' the file name or URL of a GTM file.  This is directly passed to [`readLines`]
#' and hence may contain anything that `readLines()` can handle.
#'
#' @param sep The character string separating members of each gene set in the
#' GMT file.
#'
#' @param geneIdType By default a character vector of length 1 with the special
#' value `"auto"` or an object of a subclass of [`GeneIdentifierType`].  If set
#' to `"auto"`, the function will try to derive the gene ID type from argument
#' `geneIdsList` using [`guessGeneIdType`].
#' Depending on the value of argument `valueType`, the gene ID type of the
#' resulting list or of all `GeneSet` objects in the resulting
#' `GeneSetCollection` will be set to this value.
#' 
#' @param collectionType Only used when `valueType == "GeneSetCollection"`. See
#' [`getGmt`] for more information.
#'
#' @param valueType A character vector of length 1 specifying the desired type
#' of return value.  It must be one of:
#' * `GeneSetCollection` (the default): a [`GeneSetCollection`] object as defined
#' and described by package [`GSEABase`].
#' * `list`: a named list of gene sets represented as character vectors of gene IDs.
#' This format is much simpler and cannot store the metadata required for automatic
#' mapping of gene IDs.
#'
#' @param deduplUse A character vector of length 1 specifying one of several
#' methods to handle duplicated gene set names.
#' Duplicated gene set names are explicitly forbidden by the
#' [GMT file format specification](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
#' but can nevertheless be encountered in the wild.
#' The available choices are:
#' * `first` (the default): drops all gene sets whose names are [`duplicated`]
#' according to the base R function and retains only the first occurence of a
#' gene set name.
#' * `drop`:  removes *all* gene sets that have a duplicated name, including its
#' first occurrence.
#' * `union`: replaces gene sets with duplicated names by a single gene set
#' containing the union of all their gene IDs.
#' * `smallest`: drops gene sets with duplicated names and retains only the
#' smallest of them, i.e. the one with the fewest gene IDs.  If there are
#' several smallest gene sets, the first will be selected.
#' * `largest`: drops gene sets with duplicated names and retains only the
#' largest of them, i.e. the one with the most gene IDs.  If there are
#' several largest gene sets, the first will be selected.
#'
#' @param ... Further arguments passed on to `readLines()`
#' 
#' @return The gene sets imported from the GMT file, with duplicate gene sets
#' resolved according to argument `deduplUse` and in the format determined by
#' argument `valueType`.
#' 
#' @seealso [`readLines`], [`GeneSetCollection`], [`getGmt`]
#'
#' @aliases readGMT
#' @name readGMT
#' @rdname readGMT
#' @export
#' 
readGMT <- function (con,
                     sep = "\t",
                     geneIdType = "auto",
                     collectionType = NullCollection(), 
                     valueType = c("GeneSetCollection", "list"),
                     deduplUse = c("first", "drop", "union", "smallest", "largest"),
                     ...) {
    valueType <- match.arg(valueType)
    
    ## from GSEABase::getGmt()
    lines <- strsplit(readLines(con, ...), sep)
    if (any(sapply(lines, length) < 2)) {
        txt <- paste("all records in the GMT file must have >= 2 fields", 
                     "\n  first invalid line:  %s\n", collapse = "")
        .stopf(txt, lines[sapply(lines, length) < 2][[1]])
    }
    dups <- new.env(parent = emptyenv())
    lines <- lapply(lines, function(elt, dups) {
        if (any(d <- duplicated(elt[-(1:2)]))) {
            dups[[elt[[1]]]] <- unique(elt[-(1:2)][d])
            elt <- c(elt[1:2], unique(elt[-(1:2)]))
        }
        elt
    }, dups)
    if (length(dups)) 
        .warningf("%d record(s) contain duplicate ids: %s", length(dups), 
                  paste(selectSome(sort(ls(dups))), collapse = ", "))

    ## our small addition to tolerate duplicate gene set names
    lines <- deduplicateGmtLines(lines, deduplUse)

    if(geneIdType == "auto") {
        geneIdType <- guessGeneIdType(lapply(lines, tail, -2))
    }

    ## on second thoughts, another small addition: let the user choose the return type
    if(valueType == "GeneSetCollection") {
        ## from GSEABase::getGmt()
        template <- GeneSet(geneIdType = geneIdType, collectionType = collectionType)
        return(GeneSetCollection(lapply(lines, function(line) {
            initialize(template, geneIds = unlist(line[-(1:2)]), 
                       setName = line[[1]], shortDescription = line[[2]], 
                       setIdentifier = .uniqueIdentifier())
        })))
    } else if(valueType == "list") {
        gs <- lapply(lines, tail, -2)
        names(gs) <- sapply(lines, head, 1)

        ## even more thoughts, now that we make use of gene ID metadata in lists
        if(!is.null(geneIdType)) {
            gsvaAnnotation(gs) <- geneIdType
        }
        
        return(gs)
    }
}


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
setMethod("mapGeneSetsToAnno", signature(geneSets="list", anno="NULL"),
          function(geneSets, anno, verbose=FALSE) {
              return(geneSets)
          })

setMethod("mapGeneSetsToAnno", signature(geneSets="list", anno="character"),
          function(geneSets, anno, verbose=FALSE) {
              gsc <- geneIdsToGeneSetCollection(geneIdsList=geneSets)
              return(mapGeneSetsToAnno(gsc, anno))
          })

setMethod("mapGeneSetsToAnno",
          signature(geneSets="list", anno="GeneIdentifierType"),
          function(geneSets, anno, verbose=FALSE) {
              gsc <- geneIdsToGeneSetCollection(geneIdsList=geneSets)
              return(mapGeneSetsToAnno(gsc, anno))
          })

setMethod("mapGeneSetsToAnno",
          signature(geneSets="GeneSetCollection", anno="NULL"),
          function(geneSets, anno, verbose=FALSE) {
              return(geneSets)
          })

#' @importFrom cli cli_alert_info cli_alert_warning
setMethod("mapGeneSetsToAnno",
          signature(geneSets="GeneSetCollection", anno="character"),
          function(geneSets, anno, verbose=FALSE) {
              if(.isAnnoPkgValid(anno)) {
                  if(!.isAnnoPkgInstalled(anno)) {
                      msg <- "Please install the annotation package %s"
                      stop(sprintf(msg, anno))
                  }

                  if (verbose)
                      cli_alert_info("Mapping identifiers")

                  mappedGeneSets <- mapIdentifiers(geneSets,
                                                   AnnoOrEntrezIdentifier(anno))
                  rval <- geneIds(mappedGeneSets)

              } else {
                  if (verbose) {
                      msg <- paste("No annotation metadata available in the",
                                   "input expression data object")
                      cli_alert_warning(msg)
                      msg <- paste("Attempting to directly match identifiers",
                                   "in expression data to gene sets")
                      cli_alert_warning(msg)
                  }

                  rval <- geneIds(geneSets)
              }

              return(rval)
          })

#' @importFrom cli cli_alert_info cli_alert_warning
setMethod("mapGeneSetsToAnno",
          signature(geneSets="GeneSetCollection",
                    anno="GeneIdentifierType"),
          function(geneSets, anno, verbose=FALSE) {
              annoDb <- annotation(anno)

              if(.isAnnoPkgValid(annoDb)) {
                  if(!.isAnnoPkgInstalled(annoDb)) {
                      msg <- "Please install the annotation package %s"
                      stop(sprintf(msg, annoDb))
                  }

                  if (verbose)
                      cli_alert_info("Mapping identifiers")

                  mappedGeneSets <- mapIdentifiers(geneSets, anno)
                  rval <- geneIds(mappedGeneSets)

              } else {
                  if (verbose) {
                      msg <- paste("No annotation metadata available in the",
                                   "input expression data object")
                      cli_alert_warning(msg)
                      msg <- paste("Attempting to directly match identifiers",
                                   "in expression data to gene sets")
                      cli_alert_warning(msg)
                  }

                  rval <- geneIds(geneSets)
              }

              return(rval)
          })
