
#' @title Gene Set Variation Analysis
#' 
#' @description Run Gene Set Variation Analysis on a matrix with gene sets in a list.
#' 
#' @param expr A parameter object determining the analysis method to be performed
#'   as well as containing any method-specific parameters.
#' @param minSize Minimum size of the resulting gene sets.
#' @param maxSize Maximum size of the resulting gene sets.
#' @param verbose Gives information about each calculation step. Default: `FALSE`.
#' @param BPPARAM An object of class [`BiocParallelParam`] specifiying parameters
#'   related to the parallel execution of some of the tasks and calculations
#'   within this function.
#' 
#' @return A gene-set by sample matrix (of `matrix` or [`dgCMatrix-class`] type, 
#'   depending on the input) of GSVA enrichment scores.
#' 
#' @seealso [`plageParam`], [`zscoreParam`], [`ssgseaParam`], [`gsvaParam`]
#' @aliases gsva,plageParam,missing-method
#' @rdname gsva
#' 
#' @exportMethod
setMethod("gsva", signature(expr="plageParam", gset.idx.list="missing"),
          function(expr, gset.idx.list,
                   minSize=1,
                   maxSize=Inf,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, annotation)
              
              ## filter genes according to various criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(exprData)
                      } else {
                          Biobase::annotation(exprData, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     minSize=max(1, minSize),
                                                     maxSize=maxSize)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   param = param,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, exprData)
              
              return(rval)
          })


#' @title Gene Set Variation Analysis
#' 
#' @description Run combined z-score analysis on a matrix with gene sets in a list.
#' 
#' @aliases gsva,zscoreParam,missing-method
#' @rdname gsva
#' 
#' @exportMethod
setMethod("gsva", signature(expr="zscoreParam", gset.idx.list="missing"),
          function(expr, gset.idx.list,
                   minSize=1,
                   maxSize=Inf,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, annotation)
              
              ## filter genes according to verious criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(exprData)
                      } else {
                          Biobase::annotation(exprData, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     minSize=max(1, minSize),
                                                     maxSize=maxSize)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   param = param,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, exprData)

              return(rval)
          })


#' @title Gene Set Variation Analysis
#' 
#' @description Run ssGSEA on a matrix with gene sets in a list.
#' 
#' @aliases gsva,ssgseaParam,missing-method
#' @rdname gsva
#' 
#' @exportMethod
setMethod("gsva", signature(expr="ssgseaParam", gset.idx.list="missing"),
          function(expr, gset.idx.list, param,
                   minSize=1,
                   maxSize=Inf,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, annotation)
              
              ## filter genes according to verious criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix, dropConstantRows = FALSE)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(exprData)
                      } else {
                          Biobase::annotation(exprData, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     minSize=max(1, minSize),
                                                     maxSize=maxSize)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   param = param,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, exprData)
              
              return(rval)
          })


#' @title Gene Set Variation Analysis
#' 
#' @description Run GSVA on a matrix with gene sets in a list.
#' 
#' @aliases gsva,gsvaParam,missing-method
#' @rdname gsva
#' 
#' @exportMethod
setMethod("gsva", signature(expr="gsvaParam", gset.idx.list="missing"),
          function(expr, gset.idx.list,
                   minSize=1,
                   maxSize=Inf,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              param <- expr
              
              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, annotation)
              
              ## filter genes according to verious criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(exprData)
                      } else {
                          Biobase::annotation(exprData, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     minSize=max(1, minSize),
                                                     maxSize=maxSize)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
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
          function(container, unused) {
              return(container)
          })

setMethod("unwrapData", signature("dgCMatrix"),
          function(container, unused) {
              return(container)
          })

setMethod("unwrapData", signature("ExpressionSet"),
          function(container, unused) {
              return(exprs(container))
          })

setMethod("unwrapData", signature("SummarizedExperiment"),
          function(container, assay) {
              if (length(assays(container)) == 0L)
                  stop("The input SummarizedExperiment object has no assay data.")

              if (missing(assay)) {
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

              if (missing(assay)) {
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

