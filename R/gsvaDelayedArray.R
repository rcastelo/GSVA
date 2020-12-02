.gsvaDelayedArray <- function(expr, gset.idx.list,
                  method=c("gsva", "ssgsea", "zscore", "plage"),
                  kcdf=c("Gaussian", "Poisson", "none"),
                  rnaseq=FALSE,
                  abs.ranking=FALSE,
                  parallel.sz=1L,
                  mx.diff=TRUE,
                  tau=1,
                  kernel=TRUE,
                  ssgsea.norm=TRUE,
                  verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {
  
  if (length(gset.idx.list) == 0)
    stop("The gene set list is empty! Filter may be too stringent.")
  
  if (any(lengths(gset.idx.list) == 1))
    warning("Some gene sets have size one. Consider setting 'min.sz > 1'.")
  
  parallel.sz <- as.integer(parallel.sz)
  if (parallel.sz < 1L)
    parallel.sz <- 1L
  
  ## because we keep the argument 'parallel.sz' for backwards compatibility
  ## we need to harmonize it with the contents of BPPARAM
  if (parallel.sz > 1L && class(BPPARAM) == "SerialParam") {
    BPPARAM=MulticoreParam(progressbar=verbose, workers=parallel.sz, tasks=100)
  } else if (parallel.sz == 1L && class(BPPARAM) != "SerialParam") {
    parallel.sz <- bpnworkers(BPPARAM)
  } else if (parallel.sz > 1L && class(BPPARAM) != "SerialParam") {
    bpworkers(BPPARAM) <- parallel.sz
  }
  
  if (class(BPPARAM) != "SerialParam" && verbose)
    cat(sprintf("Setting parallel calculations through a %s back-end\nwith workers=%d and tasks=100.\n",
                class(BPPARAM), parallel.sz))

  if (method == "zscore") {
    if (rnaseq)
      stop("rnaseq=TRUE does not work with method='zscore'.")
    
    if(verbose)
      cat("Estimating combined z-scores for", length(gset.idx.list), "gene sets.\n")
    
    return(zscoreDelayed(expr, gset.idx.list, parallel.sz, verbose, BPPARAM=BPPARAM))
  }
  
  if (method == "plage") {
    if (rnaseq)
      stop("rnaseq=TRUE does not work with method='plage'.")
    
    if(verbose)
      cat("Estimating PLAGE scores for", length(gset.idx.list),"gene sets.\n")
    
    return(plageDelayed(expr, gset.idx.list, parallel.sz, verbose, BPPARAM=BPPARAM))
  }
  
}

h5BackendRealization <- function(gSetIdx, FUN, Z) {
  
  FUN <- match.fun(FUN)
  
  # step 1: create realization sink
  sink <- HDF5Array::HDF5RealizationSink(dim = c(1L, ncol(Z)))
  # step 2: create grid over sink
  sink_grid <- DelayedArray::rowAutoGrid(sink, nrow = 1)
  # step 3: create block using FUN and write it on sink
  block <- FUN(gSetIdx, Z)
  block <- matrix(block, 1, length(block))
  sink <- DelayedArray::write_block(sink, sink_grid[[1L]], block)
  # step 4: close the sink as an hdf5Array
  DelayedArray::close(sink)
  res <- as(sink, "DelayedArray")
  
  res
  
}

rightsingularsvdvectorgset <- function(gSetIdx, Z) {
  s <- svd(Z[gSetIdx, ])
  s$v[, 1]
}

plageDelayed <- function(X, geneSets, parallel.sz, verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {

  Z <- t(DelayedArray::scale(t(X)))
  
  es <- bplapply(geneSets, h5BackendRealization, rightsingularsvdvectorgset,
                 Z, BPPARAM=BPPARAM)
  
  es <- do.call(rbind, es)
  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)
  
  es <- as(es, "HDF5Array")
  
  es
}

combinezDelayed <- function(gSetIdx, Z){
  DelayedMatrixStats::colSums2(Z[gSetIdx,]) / sqrt(length(gSetIdx))
}

zscoreDelayed <- function(X, geneSets, parallel.sz, verbose=TRUE,
                          BPPARAM=SerialParam(progressbar=verbose)){
  
  Z <- t(DelayedArray::scale(t(X)))
  
  es <- bplapply(geneSets, h5BackendRealization, combinezDelayed,
                 Z, BPPARAM = BPPARAM)
  
  es <- do.call(rbind, es)
  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)
  
  es <- as(es, "HDF5Array")
  
  es
}