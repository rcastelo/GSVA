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
  
  if (method == "ssgsea") {
    if(verbose)
      cat("Estimating ssGSEA scores for", length(gset.idx.list),"gene sets.\n")
    
    return(ssgseaDelayed(expr, gset.idx.list, alpha=tau, parallel.sz=parallel.sz,
                  normalization=ssgsea.norm, verbose=verbose, BPPARAM=BPPARAM))
  }

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

#### rank function for hdf5 files using sink and grid methods
rankHDF5 <- function(X){
  sink <- HDF5RealizationSink(dim(X))
  sink_grid <- colAutoGrid(sink, ncol = 1)
  X_grid <- colAutoGrid(X, ncol=1)
  
  
  FUN <- function(sink_grid, sink){
    bid <- currentBlockId()
    X_viewport <- X_grid[[bid]]
    block <- read_block(X, X_viewport)
    block <- as.integer(rank(block))
    block <- matrix(block, length(block), 1)
    write_block(sink, sink_grid, block)
  }
  
  sink <- viewportReduce(FUN, sink_grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

## slightly modified .fastRndWalk() for porpuse of only 
## receiving a vector and not a matrix column for sums
.fastRndWalk2 <- function(gSetIdx, geneRanking, ra_block) {
  n <- length(geneRanking)
  k <- length(gSetIdx)
  idxs <- sort.int(fastmatch::fmatch(gSetIdx, geneRanking))
  stepCDFinGeneSet2 <-
    sum(ra_block[geneRanking[idxs]] * (n - idxs + 1)) /
    sum((ra_block[geneRanking[idxs]]))    
  stepCDFoutGeneSet2 <- (n * (n + 1) / 2 - sum(n - idxs + 1)) / (n - k)
  walkStat <- stepCDFinGeneSet2 - stepCDFoutGeneSet2
  walkStat
}

delayedGeneRanking <- function(r_block, ra_block, geneSets, BPPARAM){
  geneRanking <- order(r_block, decreasing=TRUE)
  res <- bplapply(geneSets, .fastRndWalk2, geneRanking, ra_block, BPPARAM = BPPARAM)
  return(unlist(res))
}

ssgseaDelayed <- function(X, geneSets, alpha=0.25, parallel.sz,
                   normalization=TRUE, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {

  n <- ncol(X)
  
  R <- rankHDF5(X)
  
  Ra <- abs(R)^alpha
  
  # step1: creating a grid for both DelayedArray objects
  # that will be iterated
  R_grid <- colAutoGrid(R, ncol=1)
  Ra_grid <- colAutoGrid(Ra, ncol=1)
  
  # step2: create a sink and a grid for it
  sink <- HDF5RealizationSink(dim = c(1L, ncol(R)))
  sink_grid <- colAutoGrid(sink, ncol=1)
  
  # step3: function that will read blocks of R and
  # Ra grids, apply transformation and write into
  # sink object
  FUN <- function(sink_grid, sink, geneSets){
    bid <- currentBlockId()
    r_block <- read_block(R, R_grid[[bid]])
    ra_block <- read_block(Ra, Ra_grid[[bid]])
    res_block <- delayedGeneRanking(r_block, ra_block, geneSets, BPPARAM)
    res_block <- matrix(res_block, 1, length(res_block))
    write_block(sink, sink_grid, res_block)
  }
  
  sink <- viewportReduce(FUN, sink_grid, sink, geneSets)
  
  # step4: close sink as a DelayedArray object
  close(sink)
  es <- as(sink, "DelayedArray")
  
  if (normalization) {
    ## normalize enrichment scores by using the entire data set, as indicated
    ## by Barbie et al., 2009, online methods, pg. 2
    sink <- HDF5RealizationSink(dim(es))
    sink_grid <- colAutoGrid(sink, ncol=1)
    es_grid <- colAutoGrid(es, ncol=1)
    for(bid in seq_along(sink_grid)){
      block <- read_block(es, es_grid[[bid]])
      block <- block / (range(es)[2] - range(es)[1])
      write_block(sink, sink_grid[[bid]], block)
    }
    close(sink)
    es <- as(sink, "DelayedArray")
  }
  
  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)
  
  es <- as(es, "HDF5Array")
  es
}

