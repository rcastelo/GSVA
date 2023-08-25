
###
### 2023-08-21  axel: gsva() S4 methods for running all analysis methods
###


##
## PLAGE with a matrix of data and a list of sets
##
##' Estimates GSVA enrichment scores.
##'
##' GSVA assesses the relative enrichment of gene sets across samples using
##' a non-parametric approach. Conceptually, GSVA transforms a p-gene by n-sample
##' gene expression matrix into a g-geneset by n-sample pathway enrichment matrix.
##' This facilitates many forms of statistical analysis in the 'space' of pathways
##' rather than genes, providing a higher level of interpretability.
##' 
##' By default, `gsva()` will try to match the identifiers in `expr` to
##' the identifiers in `gset.idx.list` just as they are, unless the
##' `annotation` argument is set.
##' 
##' The `gsva()` function first maps the identifiers in the gene sets in
##' `gset.idx.list` to the identifiers in the input expression data `expr`.
##' When the input gene sets in `gset.idx.list` is provided as a `list`
##' object, `gsva()` will try to match the identifiers in `expr` directly
##' to the identifiers in `gset.idx.list` just as they are. Because unmatching
##' identifiers will be discarded in both, `expr` and `gset.idx.list`,
##' `gsva()` may prompt an error if no identifiers can be matched as in the case
##' of different types of identifiers (e.g., gene symbols vs Entrez identitifers).
##' 
##' However, then the input gene sets in `gset.idx.list` is provided as a
##' `GeneSetCollection` object, `gsva()` will try to automatically convert
##' those identifiers to the type of identifier in the input expression data `expr`.
##' Such an automatic conversion, however, will only occur in three scenarios: 1. when
##' `expr` is an `ExpressionSet` object with an appropriately set
##' `annotation` slot; 2. when `expr` is a `SummarizedExperiment` or a
##' `SingleCellExperiment` object with an appropriately set `annotation` slot
##' in the metadata of `expr`; 3. when `expr` is a `matrix` or a 
##' `dgCMatrix` and the `annotation` argument of the `gsva()` function
##' is set to the name of the annotation package that provides
##' the relationships between the type of identifiers in `expr` and `gset.idx.list`.
##' 
##' The collection of gene sets resulting from the previous identifier matching,
##' can be further filtered to require a minimun and/or maximum size by using the
##' arguments `min.sz` and `max.sz`.
##' 
##' If you use GSVA in your research, please cite also the corresponding method as
##' described in the `method` parameter.
##' 
##' @title Gene Set Variation Analysis
##' @param expr Gene expression data which can be given either as a
##'   `SummarizedExperiment`, `SingleCellExperiment`
##'   `ExpressionSet` object, or as a matrix of expression
##'   values where rows correspond to genes and columns correspond to samples.
##'   This matrix can be also in a sparse format, as a `dgCMatrix`, or
##'   as an on-disk backend representation, such as `HDF5Array` .
##' @param gset.idx.list Gene sets provided either as a `list` object or as a
##'   `GeneSetCollection` object.
##' @param param A parameter object determining the analysis method to be performed
##'   as well as containing any method-specific parameters.
##' @param annotation In the case of calling `gsva()` on a
##'   `SummarizedExperiment` or `SingleCellExperiment` object,
##'   the `annotation` argument can be used to select the assay
##'   containing the molecular data we want as input to the `gsva()`
##'   function, otherwise the first assay is selected.
##'   In the case of calling `gsva()` with expression data in
##'   a `matrix` and gene sets as a `GeneSetCollection`
##'   object, the `annotation` argument can be used to supply
##'   the name of the Bioconductor package that contains
##'   annotations for the class of gene identifiers occurring in
##'   the row names of the expression data matrix.
##'   In the case of calling `gsva()` on a
##'   `ExpressionSet` object, the `annotation` argument
##'   is ignored. See details information below.
##' @param min.sz Minimum size of the resulting gene sets.
##' @param max.sz Maximum size of the resulting gene sets.
##' @param parallel.sz Number of threads of execution to use when doing the calculations in parallel.
##'   The argument BPPARAM allows one to set the parallel back-end and fine
##'   tune its configuration.
##' @param verbose Gives information about each calculation step. Default: `FALSE`.
##' @param BPPARAM An object of class \linkS4class{BiocParallelParam} specifiying parameters
##'   related to the parallel execution of some of the tasks and calculations within this function.
##' @return A gene-set by sample matrix (of `matrix` or `dgCMatrix` type, 
##'   depending on the input) of GSVA enrichment scores.
##' 
##' @rdname gsvaNewAPI
##' @export 
setMethod("gsva", signature(expr="matrix", gset.idx.list="list", param = "plageParam"),
          function(expr, gset.idx.list, param,
                   annotation, 
                   min.sz=1,
                   max.sz=Inf,
                   parallel.sz=1L, 
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              ## filter genes according to verious criteria,
              ## e.g., constant expression
              expr <- .filterFeatures_newAPI(expr)

              ## map to the actual features for which expression data is available
              mapped.gset.idx.list <- .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
              
              ## remove gene sets from the analysis for which no features are available
              ## and meet the minimum and maximum gene-set size specified by the user
              mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                                     min.sz=max(1, min.sz),
                                                     max.sz=max.sz)

              rval <- .gsva_newAPI(expr = expr,
                                   gset.idx.list = mapped.gset.idx.list,
                                   method = "plage",
                                   rnaseq = FALSE, 
                                   parallel.sz = parallel.sz,
                                   verbose = verbose,
                                   BPPARAM = BPPARAM)

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


.gsva_newAPI <- function(expr, gset.idx.list, method,
                         kcdf=c("Gaussian", "Poisson", "none"),
                         rnaseq=FALSE,
                         parallel.sz=1L,
                         kernel=TRUE,
                         verbose=TRUE,
                         BPPARAM=SerialParam(progressbar=verbose)) {
    
    if(is(expr, "DelayedArray")){
        warning("Using 'DelayedArray' objects as input is still in an experimental stage.")

        return(.gsvaDelayedArray(expr, gset.idx.list, method, kcdf, rnaseq, abs.ranking,
                                 parallel.sz, mx.diff, tau, kernel, ssgsea.norm, verbose, BPPARAM))
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

    if (method == "ssgsea") {
        if(verbose)
            cat("Estimating ssGSEA scores for", length(gset.idx.list),"gene sets.\n")

        return(ssgsea(expr, gset.idx.list, alpha=tau, parallel.sz=parallel.sz,
                      normalization=ssgsea.norm, verbose=verbose, BPPARAM=BPPARAM))
    }

    if (method == "zscore") {
        if (rnaseq)
            stop("rnaseq=TRUE does not work with method='zscore'.")

        if(verbose)
            cat("Estimating combined z-scores for", length(gset.idx.list), "gene sets.\n")

        return(zscore(expr, gset.idx.list, parallel.sz, verbose, BPPARAM=BPPARAM))
    }

    if (method == "plage") {
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
    
    es.obs <- matrix(NaN, n.gset, n.samples, dimnames=list(names(gset.idx.list),colnames(expr)))
    colnames(es.obs) <- colnames(expr)
    rownames(es.obs) <- names(gset.idx.list)
    
    es.obs <- compute.geneset.es(expr, gset.idx.list, 1:n.samples,
                                 rnaseq=rnaseq, abs.ranking=abs.ranking,
                                 parallel.sz=parallel.sz,
                                 mx.diff=mx.diff, tau=tau, kernel=kernel,
                                 verbose=verbose, BPPARAM=BPPARAM)
    
    colnames(es.obs) <- colnames(expr)
    rownames(es.obs) <- names(gset.idx.list)

    es.obs
}



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
    tid[["mxFloat"]] <- y

    ## same structure with count data
    d <- matrix(rpois(n*p, 42), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p) , paste0("s", 1:n)))
    d[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] <- d[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] + 23
    tid[["mxCount"]] <- d

    return(tid)
}

