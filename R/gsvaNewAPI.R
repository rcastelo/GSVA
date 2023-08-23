
###
### 2023-08-21  axel: gsva() S4 methods for running all analysis methods
###


###
### PLAGE with a matrix of data and a list of sets
###
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

