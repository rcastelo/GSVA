
#' @title Gene Set Variation Analysis
#' 
#' @description Run Gene Set Variation Analysis on a matrix with gene sets in a list.
#' 
#' @param expr A parameter object of one of the following classes:
#' * A [`gsvaParam`] object built using the constructor function [`gsvaParam`].
#'   This object will trigger `gsva()` to use the GSVA algorithm by
#'   Hänzelmann et al. (2013).
#' * A [`plageParam`] object built using the constructor function [`plageParam`].
#'   This object will trigger `gsva()` to use the PLAGE algorithm by
#'   Tomfohr et al. (2005).
#' * A [`zscoreParam`] object built using the constructor function [`zscoreParam`]
#'   This object will trigger `gsva()` to use the combined z-score algorithm by
#'   Lee et al. (2008).
#' * A [`ssgseaParam`] object built using the constructor function [`ssgseaParam`]
#'   This object will trigger `gsva()` to use the ssGSEA algorithm by
#'   Barbie et al. (2009).
#'
#' @param gset.idx.list Dummy parameter, only present for backward compatibility,
#' do not use it. It will be removed once the deprecated version of 'gsva()'
#' is defunct.
#' 
#' @param verbose Gives information about each calculation step. Default: `FALSE`.
#' 
#' @param BPPARAM An object of class [`BiocParallelParam`] specifiying parameters
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

#' @aliases gsva,plageParam,missing-method
#' @rdname gsva
#' @exportMethod
setMethod("gsva", signature(expr="plageParam", gset.idx.list="missing"),
          function(expr, gset.idx.list,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr # for backward compatibility with the old API only
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              
              ## filter genes according to various criteria,
              ## e.g., constant expression
              dataMatrix <- .filterFeatures_newAPI(dataMatrix)

              anno <- annotation(exprData, get_annotation(param))
              geneSets <- get_geneSets(param)
              geneSets <- mapGeneSetsToAnno(geneSets, anno)
              
              ## map to the actual features for which expression data is available
              mappedGeneSets <- .mapGeneSetsToFeatures(geneSets, rownames(dataMatrix))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mappedGeneSets <- filterGeneSets(mappedGeneSets,
                                               min.sz=get_minSize(param),
                                               max.sz=get_maxSize(param))

              rval <- .gsva_newAPI(expr = dataMatrix,
                                   gset.idx.list = mappedGeneSets,
                                   param = param,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, exprData)
              
              return(rval)
          })


#' @aliases gsva,zscoreParam,missing-method
#' @rdname gsva
#' @exportMethod
setMethod("gsva", signature(expr="zscoreParam", gset.idx.list="missing"),
          function(expr, gset.idx.list,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr # for backward compatibility with the old API only
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              
              ## filter genes according to various criteria,
              ## e.g., constant expression
              dataMatrix <- .filterFeatures_newAPI(dataMatrix)

              anno <- annotation(exprData, get_annotation(param))
              geneSets <- get_geneSets(param)
              geneSets <- mapGeneSetsToAnno(geneSets, anno)
              
              ## map to the actual features for which expression data is available
              mappedGeneSets <- .mapGeneSetsToFeatures(geneSets, rownames(dataMatrix))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mappedGeneSets <- filterGeneSets(mappedGeneSets,
                                               min.sz=get_minSize(param),
                                               max.sz=get_maxSize(param))

              rval <- .gsva_newAPI(expr = dataMatrix,
                                   gset.idx.list = mappedGeneSets,
                                   param = param,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, exprData)
              
              return(rval)
          })


#' @aliases gsva,ssgseaParam,missing-method
#' @rdname gsva
#' @exportMethod
setMethod("gsva", signature(expr="ssgseaParam", gset.idx.list="missing"),
          function(expr, gset.idx.list,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr # for backward compatibility with the old API only
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              
              ## filter genes according to various criteria,
              ## e.g., constant expression
              dataMatrix <- .filterFeatures_newAPI(dataMatrix,
                                                   dropConstantRows=FALSE)

              anno <- annotation(exprData, get_annotation(param))
              geneSets <- get_geneSets(param)
              geneSets <- mapGeneSetsToAnno(geneSets, anno)
              
              ## map to the actual features for which expression data is available
              mappedGeneSets <- .mapGeneSetsToFeatures(geneSets, rownames(dataMatrix))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mappedGeneSets <- filterGeneSets(mappedGeneSets,
                                               min.sz=get_minSize(param),
                                               max.sz=get_maxSize(param))

              rval <- .gsva_newAPI(expr = dataMatrix,
                                   gset.idx.list = mappedGeneSets,
                                   param = param,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, exprData)
              
              return(rval)
          })


#' @aliases gsva,gsvaParam,missing-method
#' @rdname gsva
#' @exportMethod
setMethod("gsva", signature(expr="gsvaParam", gset.idx.list="missing"),
          function(expr, gset.idx.list,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr # for backward compatibility with the old API only
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              
              ## filter genes according to various criteria,
              ## e.g., constant expression
              dataMatrix <- .filterFeatures_newAPI(dataMatrix)

              anno <- annotation(exprData, get_annotation(param))
              geneSets <- get_geneSets(param)
              geneSets <- mapGeneSetsToAnno(geneSets, anno)
              
              ## map to the actual features for which expression data is available
              mappedGeneSets <- .mapGeneSetsToFeatures(geneSets, rownames(dataMatrix))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mappedGeneSets <- filterGeneSets(mappedGeneSets,
                                               min.sz=get_minSize(param),
                                               max.sz=get_maxSize(param))

              rval <- .gsva_newAPI(expr = dataMatrix,
                                   gset.idx.list = mappedGeneSets,
                                   param = param,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, exprData)
              
              return(rval)
          })



### ----- GSVA utility functions, slightly adapted for new API -----
.filterFeatures_newAPI <- function(expr, dropConstantRows = TRUE) {
    ## filter out genes with constant expression values
    ## DelayedMatrixStats::rowSds() works for both base and 
    ## DelayedArray matrices
    sdGenes <- DelayedMatrixStats::rowSds(expr)
    ## the following fixes this bug, see issues
    ## https://github.com/rcastelo/GSVA/issues/54
    ## https://github.com/HenrikBengtsson/matrixStats/issues/204
    sdGenes[sdGenes < 1e-10] <- 0
    if(any(sdGenes == 0) || any(is.na(sdGenes))) {
        ## CHECK: unless we actually drop, we do all this only for this warning, right?
        warning(sum(sdGenes == 0 | is.na(sdGenes)),
                " rows with constant values throughout the samples.")
        if(dropConstantRows) {
            warning("Rows with constant values are discarded.")
            expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
        }
    } 

    if(nrow(expr) < 2)
        stop("Less than two genes in the input assay object\n")
    
    ## CHECK: is this the right place to check this?
    if(is.null(rownames(expr)))
        stop("The input assay object doesn't have rownames\n")
    
    return(expr)
}


.gsva_newAPI <- function(expr, gset.idx.list, param,
                         rnaseq = FALSE,
                         verbose=TRUE,
                         BPPARAM=SerialParam(progressbar=verbose)) {
    
    if(is(expr, "DelayedArray")){
        warning("Using 'DelayedArray' objects as input is still in an experimental stage.")

#         return(.gsvaDelayedArray(expr, gset.idx.list, method, kcdf, rnaseq, abs.ranking,
#                                  parallel.sz, mx.diff, tau, kernel, ssgsea.norm, verbose, BPPARAM))
    }

    if (length(gset.idx.list) == 0)
        stop("The gene set list is empty! Filter may be too stringent.")

    if (any(lengths(gset.idx.list) == 1))
        warning("Some gene sets have size one. Consider setting 'minSize > 1'.")

    parallel.sz <- if (inherits(BPPARAM, "SerialParam")) 1L else bpnworkers(BPPARAM)
    
    if (!inherits(BPPARAM, "SerialParam") && verbose)
        cat(sprintf("Setting parallel calculations through a %s back-end\nwith workers=%d and tasks=100.\n",
                    class(BPPARAM), bpnworkers(BPPARAM)))

    if (inherits(param, "ssgseaParam")) {
        if(verbose)
            cat("Estimating ssGSEA scores for", length(gset.idx.list),"gene sets.\n")

        return(ssgsea(expr, gset.idx.list, alpha=get_alpha(param), parallel.sz=parallel.sz,
                      normalization=do_normalize(param), verbose=verbose, BPPARAM=BPPARAM))
    }

    ## CHECK: AFAICS rnaseq is either FALSE because kcdf=="Gaussian" by default (and kcdf shouldn't be considered
    ## for methods other than 'gsva')or it is missing and consequently FALSE by the default of this function; and
    ## the only way rnaseq could be TRUE would be to set kcdf = "Poisson", regardless of the method?!?
    if (inherits(param, "zscoreParam")) {
        if (rnaseq)
            stop("rnaseq=TRUE does not work with method='zscore'.")

        if(verbose)
            cat("Estimating combined z-scores for", length(gset.idx.list), "gene sets.\n")

        return(zscore(expr, gset.idx.list, parallel.sz, verbose, BPPARAM=BPPARAM))
    }

    if (inherits(param, "plageParam")) {
        if (rnaseq)
            stop("rnaseq=TRUE does not work with method='plage'.")

        if(verbose)
            cat("Estimating PLAGE scores for", length(gset.idx.list),"gene sets.\n")

        return(plage(expr, gset.idx.list, parallel.sz, verbose, BPPARAM=BPPARAM))
    }

    if(verbose)
        cat("Estimating GSVA scores for", length(gset.idx.list),"gene sets.\n")
    
    n.samples <- ncol(expr)
    n.genes <- nrow(expr)
    n.gset <- length(gset.idx.list)

    ## CHECK: isn't this the only place this should be used?!?
    if (get_kcdf(param) == "Gaussian") {
        rnaseq <- FALSE
        kernel <- TRUE
    } else if (get_kcdf(param) == "Poisson") {
        rnaseq <- TRUE
        kernel <- TRUE
    } else
        kernel <- FALSE
    
    es.obs <- matrix(NaN, n.gset, n.samples, dimnames=list(names(gset.idx.list),colnames(expr)))
    colnames(es.obs) <- colnames(expr)
    rownames(es.obs) <- names(gset.idx.list)
    
    es.obs <- compute.geneset.es(expr, gset.idx.list, 1:n.samples,
                                 rnaseq=rnaseq, abs.ranking=get_absRanking(param),
                                 parallel.sz=parallel.sz,
                                 mx.diff=get_maxDiff(param), tau=get_tau(param), kernel=kernel,
                                 verbose=verbose, BPPARAM=BPPARAM)
    
    colnames(es.obs) <- colnames(expr)
    rownames(es.obs) <- names(gset.idx.list)

    es.obs
}


### -----  methods for data pre-/post-processing -----

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

              return(as.matrix(assays(container)[[assay]]))
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

              return(as.matrix(assays(container)[[assay]]))
          })


## wrapData: put the resulting data into the original data container type
setMethod("wrapData", signature("matrix", "matrix"),
          function(dataMatrix, container) {
              return(dataMatrix)
          })

setMethod("wrapData", signature("matrix", "ExpressionSet"),
          function(dataMatrix, container) {
              rval <- new("ExpressionSet", exprs=dataMatrix,
                          phenoData=phenoData(container),
                          experimentData=experimentData(container),
                          annotation="")

              return(rval)
          })

setMethod("wrapData", signature("matrix", "SummarizedExperiment"),
          function(dataMatrix, container) {
              rval <- SummarizedExperiment(assays=SimpleList(es=dataMatrix),
                                           colData=colData(container),
                                           metadata=metadata(container))
              metadata(rval)$annotation <- NULL

              return(rval)
          })

setMethod("wrapData", signature("matrix", "SingleCellExperiment"),
          function(dataMatrix, container) {
              rval <- SingleCellExperiment(assays=SimpleList(es=dataMatrix),
                                           colData=colData(container),
                                           metadata=metadata(container))
              metadata(rval)$annotation <- NULL
              
              return(rval)
          })

## annotation: from BiocGenerics, some more methods for us
setMethod("annotation", signature("SingleCellExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })

setMethod("annotation", signature("SummarizedExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })

setMethod("annotation", signature("GsvaExprData"),
          function(object, default = NULL) {
              return(default)
          })

## annotation package checks
isAnnoPkgValid <- function(ap) {
    return((!is.null(ap)) &&
           (length(ap) == 1) &&
           (is.character(ap)) &&
           (!is.na(ap)) &&
           (nchar(ap)> 0))
}

isAnnoPkgInstalled <- function(ap) {
    ap <- c(ap, paste0(ap, ".db"))
    return(any(ap %in% rownames(installed.packages())))
}

## mapGeneSetsToAnno: translate feature IDs used in gene sets to specified annotation type (if any, and if possible)
setMethod("mapGeneSetsToAnno", signature("list"),
          function(geneSets, anno) {
              return(geneSets)
          })

setMethod("mapGeneSetsToAnno", signature("GeneSetCollection"),
          function(geneSets, anno) {
              if(isAnnoPkgValid(anno)) {
                  if(!isAnnoPkgInstalled(anno))
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


### ----- generate test input data -----
generateTestInputData <- function(p = 10, nGS = 3,
                                  n = 30, nGrp1 = round(n / 2), nGrp2 = n - nGrp1) {
    ## we need to be able to reproduce this
    set.seed(2023-08-18)
    tid <- list()

    geneSets <- split(paste0("g", seq.int(p)),
                      cut(seq.int(p), breaks = nGS, labels = paste0("GeneSet", seq.int(nGS))))
    tid[["geneSets"]] <- geneSets

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p) , paste0("s", 1:n)))
    ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
    y[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] <- y[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] + 2
    tid[["continuousData"]] <- y

    ## same structure with count data
    d <- matrix(rpois(n*p, 42), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p) , paste0("s", 1:n)))
    d[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] <- d[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] + 23
    tid[["countData"]] <- d

    return(tid)
}

