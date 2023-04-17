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

  stop("Not yet implemented for 'method=\"gsva\"'")
}

h5BackendRealization <- function(gSetIdx, FUN, Z, bpp) {
  
  FUN <- match.fun(FUN)
  
  # step 1: create realization sink
  sink <- HDF5RealizationSink(dim = c(1L, ncol(Z)))
  # step 2: create grid over sink
  sink_grid <- rowAutoGrid(sink, nrow = 1)
  # step 3: create block using FUN and write it on sink
  block <- FUN(gSetIdx, Z, bpp)
  block <- matrix(block, 1, length(block))
  sink <- DelayedArray::write_block(sink, sink_grid[[1L]], block)
  # step 4: close the sink as an hdf5Array
  DelayedArray::close(sink)
  res <- as(sink, "DelayedArray")
  
  res
  
}

rightsvdDelayed <- function(gSetIdx, Z, bpp) {
  s <-runRandomSVD(Z[gSetIdx,], k=1, nu=0, BPPARAM = bpp)
  s$v[, 1]
}

plageDelayed <- function(X, geneSets, parallel.sz, verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {

  Z <- t(DelayedArray::scale(t(X)))
  
  es <- bplapply(geneSets, h5BackendRealization, rightsvdDelayed,
                 Z, bpp=BPPARAM, BPPARAM=BPPARAM)
  
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
  
  es <- bplapply(geneSets, function(gSetIdx, Z){
    x <- combinezDelayed(gSetIdx, Z)
    x <- matrix(x, 1, length(x))
    x <- writeHDF5Array(x)
  }, Z, BPPARAM = BPPARAM)
  es <- do.call(DelayedArray::rbind, es)
  es <- as(es, "HDF5Array")
  es
}

#### rank function for hdf5 files using sink and grid methods
rankHDF5 <- function(X){
  sink <- HDF5RealizationSink(dim(X))
  grid <- defaultAutoGrid(sink, block.shape="first-dim-grows-first")
  
  colRanks_byBlock <- function(grid, sink){
    block <- read_block(X, grid)
    block <- t(colRanks(block, ties.method = "average"))
    mode(block) <- "integer"
    write_block(sink, grid, block)
  }
  
  sink <- gridReduce(colRanks_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

## slightly modified .fastRndWalk() for porpoise of only 
## receiving a vector and not a matrix column for sums
.fastRndWalk2 <- function(gSetIdx, geneRanking, ra_block) {
  n <- length(geneRanking)
  k <- length(gSetIdx)
  idxs <- sort.int(match(gSetIdx, geneRanking))
  stepCDFinGeneSet2 <-
    sum(ra_block[geneRanking[idxs]] * (n - idxs + 1)) /
    sum((ra_block[geneRanking[idxs]]))    
  stepCDFoutGeneSet2 <- (n * (n + 1) / 2 - sum(n - idxs + 1)) / (n - k)
  walkStat <- stepCDFinGeneSet2 - stepCDFoutGeneSet2
  walkStat
}


ssgseaDelayed <- function(X, geneSets, alpha=0.25, parallel.sz,
                   normalization=TRUE, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {
  
  n <- ncol(X)
  
  R <- rankHDF5(X)
  
  Ra <- abs(R)^alpha
  
  es <- bplapply(as.list(1:n), function(j) {
    geneRanking <- order(R[, j], decreasing=TRUE)
    colRa <- Ra[,j]
    sink <- HDF5RealizationSink(c(length(names(geneSets)), 1L))
    sink_grid <- colAutoGrid(sink, ncol=1)
    es_sample <- lapply(geneSets, .fastRndWalk2, geneRanking, colRa)
    sink <- DelayedArray::write_block(sink, sink_grid[[1L]], do.call("rbind", es_sample))
    DelayedArray::close(sink)
    res <- as(sink, "DelayedArray")
    res
  }, BPPARAM=BPPARAM)
  
  es <- do.call(DelayedArray::cbind, es)
  es
  
  if (normalization) {
    ## normalize enrichment scores by using the entire data set, as indicated
    ## by Barbie et al., 2009, online methods, pg. 2
    sink <- HDF5RealizationSink(dim(es))
    sink_grid <- defaultAutoGrid(sink, block.shape="first-dim-grows-first")
    es_grid <- defaultAutoGrid(es, block.shape="first-dim-grows-first")
    fin <- range(es)[2]
    ini <- range(es)[1]
    for(bid in seq_along(sink_grid)){
      block <- read_block(es, es_grid[[bid]])
      block <- apply(block, 2, function(x) x / ( fin - ini))
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

