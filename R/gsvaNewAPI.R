
###
### gsva() S4 methods for running all analysis methods
###

##' @title Gene Set Variation Analysis
##' @description Run PLAGE on a matrix with gene sets in a list.
##' @describeIn gsvaNewAPI Run PLAGE on a matrix with gene sets in a list.
##' @param expr Gene expression data which can be given either as a
##'   [`SummarizedExperiment`], [`SingleCellExperiment`]
##'   [`ExpressionSet`] object, or as a matrix of expression
##'   values where rows correspond to genes and columns correspond to samples.
##'   This matrix can be also in a sparse format, as a [`dgCMatrix-class`], or
##'   as an on-disk backend representation, such as [`HDF5Array`] .
##' @param gset.idx.list Gene sets provided either as a `list` object or as a
##'   [`GeneSetCollection`] object.
##' @param param A parameter object determining the analysis method to be performed
##'   as well as containing any method-specific parameters.
##' @param annotation In the case of calling `gsva()` on a
##'   [`SummarizedExperiment`] or [`SingleCellExperiment`] object,
##'   the `annotation` argument can be used to select the assay
##'   containing the molecular data we want as input to the `gsva()`
##'   function, otherwise the first assay is selected.
##'   In the case of calling `gsva()` with expression data in
##'   a `matrix` and gene sets as a [`GeneSetCollection`]
##'   object, the `annotation` argument can be used to supply
##'   the name of the Bioconductor package that contains
##'   annotations for the class of gene identifiers occurring in
##'   the row names of the expression data matrix.
##'   In the case of calling `gsva()` on a
##'   [`ExpressionSet`] object, the `annotation` argument
##'   is ignored. See details information below.
##' @param min.sz Minimum size of the resulting gene sets.
##' @param max.sz Maximum size of the resulting gene sets.
##' @param parallel.sz Number of threads of execution to use when doing the calculations in parallel.
##'   The argument BPPARAM allows one to set the parallel back-end and fine
##'   tune its configuration.
##' @param verbose Gives information about each calculation step. Default: `FALSE`.
##' @param BPPARAM An object of class [`BiocParallelParam`] specifiying parameters
##'   related to the parallel execution of some of the tasks and calculations within this function.
##' @return A gene-set by sample matrix (of `matrix` or [`dgCMatrix-class`] type, 
##'   depending on the input) of GSVA enrichment scores.
##' @seealso [`plageParam`], [`zscoreParam`], [`ssgseaParam`], [`gsvaParam`]
setMethod("gsva", signature(expr="missing", gset.idx.list="missing", param = "plageParam"),
          function(expr, gset.idx.list, param,
                   annotation, 
                   min.sz=1,
                   max.sz=Inf,
                   parallel.sz=1L, 
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              message("¡Hola PLAGE!")
              dataSet <- get_dataSet(param)
              dataMatrix <- unwrapData(dataSet, annotation)
              
              ## filter genes according to various criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(dataSet)
                      } else {
                          Biobase::annotation(dataSet, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     min.sz=max(1, min.sz),
                                                     max.sz=max.sz)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   param = param,
                                   parallel.sz = parallel.sz,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, dataSet)
              
              return(rval)
          })


##' @title Gene Set Variation Analysis
##' @description Run z-score analysis on a matrix with gene sets in a list.
##' @describeIn gsvaNewAPI Run z-score analysis on a matrix with gene sets in a list.
setMethod("gsva", signature(expr="missing", gset.idx.list="missing", param = "zscoreParam"),
          function(expr, gset.idx.list, param,
                   annotation, 
                   min.sz=1,
                   max.sz=Inf,
                   parallel.sz=1L, 
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              message("¡Hola z-Score!")
              dataSet <- get_dataSet(param)
              dataMatrix <- unwrapData(dataSet, annotation)
              
              ## filter genes according to verious criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(dataSet)
                      } else {
                          Biobase::annotation(dataSet, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     min.sz=max(1, min.sz),
                                                     max.sz=max.sz)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   param = param,
                                   parallel.sz = parallel.sz,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, dataSet)

              return(rval)
          })


##' @title Gene Set Variation Analysis
##' @description Run ssGSEA on a matrix with gene sets in a list.
##' @describeIn gsvaNewAPI Run ssGSEA on a matrix with gene sets in a list.
setMethod("gsva", signature(expr="missing", gset.idx.list="missing", param = "ssgseaParam"),
          function(expr, gset.idx.list, param,
                   annotation, 
                   min.sz=1,
                   max.sz=Inf,
                   parallel.sz=1L, 
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              message("¡Hola ssGSEA!")
              dataSet <- get_dataSet(param)
              dataMatrix <- unwrapData(dataSet, annotation)
              
              ## filter genes according to verious criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix, dropConstantRows = FALSE)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(dataSet)
                      } else {
                          Biobase::annotation(dataSet, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     min.sz=max(1, min.sz),
                                                     max.sz=max.sz)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   param = param,
                                   parallel.sz = parallel.sz,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, dataSet)
              
              return(rval)
          })


##' @title Gene Set Variation Analysis
##' @description Run GSVA on a matrix with gene sets in a list.
##' @describeIn gsvaNewAPI Run GSVA on a matrix with gene sets in a list.
setMethod("gsva", signature(expr="missing", gset.idx.list="missing", param = "gsvaParam"),
          function(expr, gset.idx.list, param,
                   annotation, 
                   min.sz=1,
                   max.sz=Inf,
                   parallel.sz=1L, 
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              message("¡Hola GSVA!")
              dataSet <- get_dataSet(param)
              dataMatrix <- unwrapData(dataSet, annotation)
              
              ## filter genes according to verious criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(dataMatrix)

              anno <- if(missing(annotation)) {
                          Biobase::annotation(dataSet)
                      } else {
                          Biobase::annotation(dataSet, annotation)
                      }
              
              gset.idx.list <- get_geneSets(param)
              gset.idx.list <- mapGeneSetsToAnno(gset.idx.list, anno)
              
              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     min.sz=max(1, min.sz),
                                                     max.sz=max.sz)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   param = param,
                                   parallel.sz = parallel.sz,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)
              rval <- wrapData(rval, dataSet)
              
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
                         parallel.sz=1L,
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
        warning("Some gene sets have size one. Consider setting 'min.sz > 1'.")

    parallel.sz <- as.integer(parallel.sz)
    if (parallel.sz < 1L)
        parallel.sz <- 1L
    
    ## because we keep the argument 'parallel.sz' for backwards compatibility
    ## we need to harmonize it with the contents of BPPARAM
    if (parallel.sz > 1L && inherits(BPPARAM, "SerialParam")) {
        BPPARAM=MulticoreParam(progressbar=verbose, workers=parallel.sz, tasks=100)
    } else if (parallel.sz == 1L && !inherits(BPPARAM, "SerialParam")) {
        parallel.sz <- bpnworkers(BPPARAM)
    } else if (parallel.sz > 1L && !inherits(BPPARAM, "SerialParam")) {
        bpworkers(BPPARAM) <- parallel.sz
    }

    if (!inherits(BPPARAM, "SerialParam") && verbose)
        cat(sprintf("Setting parallel calculations through a %s back-end\nwith workers=%d and tasks=100.\n",
                    class(BPPARAM), parallel.sz))

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
                                 rnaseq=rnaseq, abs.ranking=get_abs.ranking(param),
                                 parallel.sz=parallel.sz,
                                 mx.diff=get_mx.diff(param), tau=get_tau(param), kernel=kernel,
                                 verbose=verbose, BPPARAM=BPPARAM)
    
    colnames(es.obs) <- colnames(expr)
    rownames(es.obs) <- names(gset.idx.list)

    es.obs
}


### ----- generics and methods for data pre-/post-processing

## unwrapData: extract a data matrix from a container object
setGeneric("unwrapData",
           function(container, ...) standardGeneric("unwrapData"))

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
setGeneric("wrapData",
           function(dataMatrix, container) standardGeneric("wrapData"))

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

setMethod("annotation", signature("GsvaDataSet"),
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
setGeneric("mapGeneSetsToAnno",
           function(geneSets, ...) standardGeneric("mapGeneSetsToAnno"))

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

