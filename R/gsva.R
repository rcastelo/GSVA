##
## function: gsva
## purpose: main function of the package which estimates activity
##          scores for each given gene-set

#' @title Gene Set Variation Analysis
#' @description This is the old manual page of the defunct version
#' of the function `gsva()`.
#'
#' @name gsva-defunct
## #' @usage gsva(expr, gset.idx.list, ...)
#' @seealso [`GSVA-pkg-defunct`]
#' @keywords internal
NULL

#' @rdname GSVA-pkg-defunct
#' @section Details:
#' Instead of `gsva(expr=., gset.idx.list=., method=., ...)`, use a
#' method-specific parameter object,
#' see [`plageParam`] [`zscoreParam`] [`ssgseaParam`] [`gsvaParam`],
#' followed by a call to the new `gsva()` function, see [`gsva`].
#' @export
setMethod("gsva", signature(param="missing"), function(param, ...) {
    .Defunct(new = "gsva()", package="GSVA",
             msg="Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').")
})

#' @rdname GSVA-pkg-defunct
#' @export
setMethod("gsva", signature(param="matrix"), function(param, ...) {
    .Defunct(new = "gsva()", package="GSVA",
             msg="Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').")
})

#' @rdname GSVA-pkg-defunct
#' @export
setMethod("gsva", signature(param="ExpressionSet"), function(param, ...) {
    .Defunct(new = "gsva()", package="GSVA",
             msg="Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').")
})

#' @rdname GSVA-pkg-defunct
#' @export
setMethod("gsva", signature(param="SummarizedExperiment"), function(param, ...) {
    .Defunct(new = "gsva()", package="GSVA",
             msg="Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').")
})

#' @rdname GSVA-pkg-defunct
#' @export
setMethod("gsva", signature(param="dgCMatrix"), function(param, ...) {
    .Defunct(new = "gsva()", package="GSVA",
             msg="Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').")
})

#' @rdname GSVA-pkg-defunct
#' @export
setMethod("gsva", signature(param="SingleCellExperiment"), function(param, ...) {
    .Defunct(new = "gsva()", package="GSVA",
             msg="Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').")
})


compute.gene.density <- function(expr, sample.idxs, rnaseq=FALSE, kernel=TRUE){
    n.test.samples <- ncol(expr)
    n.genes <- nrow(expr)
    n.density.samples <- length(sample.idxs)
    
    gene.density <- NA
    if (kernel) {
        A = .Call("matrix_density_R",
                  as.double(t(expr[ ,sample.idxs, drop=FALSE])),
                  as.double(t(expr)),
                  n.density.samples,
                  n.test.samples,
                  n.genes,
                  as.integer(rnaseq))
	
        gene.density <- t(matrix(A, n.test.samples, n.genes))
    } else {
        gene.density <- t(apply(expr, 1, function(x, sample.idxs) {
            f <- ecdf(x[sample.idxs])
            f(x)
        }, sample.idxs))
        gene.density <- log(gene.density / (1-gene.density))
    }

    return(gene.density)	
}

compute.geneset.es <- function(expr, gset.idx.list, sample.idxs, rnaseq=FALSE,
                               abs.ranking, parallel.sz=1L, 
                               mx.diff=TRUE, tau=1, kernel=TRUE,
                               verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose)) {
    num_genes <- nrow(expr)
    if (verbose) {
        if (kernel) {
            if (rnaseq)
                cat("Estimating ECDFs with Poisson kernels\n")
            else
                cat("Estimating ECDFs with Gaussian kernels\n")
        } else
            cat("Estimating ECDFs directly\n")
    }

    ## open parallelism only if ECDFs have to be estimated for
    ## more than 100 genes on more than 100 samples
    if (parallel.sz > 1 && length(sample.idxs) > 100 && nrow(expr) > 100) {
        if (verbose)
            cat(sprintf("Estimating ECDFs in parallel on %d cores\n", as.integer(parallel.sz)))
        iter <- function(Y, n_chunks=BiocParallel::multicoreWorkers()) {
            idx <- splitIndices(nrow(Y), min(nrow(Y), n_chunks))
            i <- 0L
            function() {
                if (i == length(idx))
                    return(NULL)
                i <<- i + 1L
                Y[idx[[i]], , drop=FALSE]
            }
        }
        gene.density <- bpiterate(iter(expr, 100),
                                  compute.gene.density,
                                  sample.idxs=sample.idxs,
                                  rnaseq=rnaseq, kernel=kernel,
                                  REDUCE=rbind, reduce.in.order=TRUE,
                                  BPPARAM=BPPARAM)
    } else 
        gene.density <- compute.gene.density(expr, sample.idxs, rnaseq, kernel)
    
    compute_rank_score <- function(sort_idx_vec){
        tmp <- rep(0, num_genes)
        tmp[sort_idx_vec] <- abs(seq(from=num_genes,to=1) - num_genes/2)
        return (tmp)
    }
    
    rank.scores <- rep(0, num_genes)
    sort.sgn.idxs <- apply(gene.density, 2, order, decreasing=TRUE) # n.genes * n.samples
    
    rank.scores <- apply(sort.sgn.idxs, 2, compute_rank_score)

    m <- bplapply(gset.idx.list, ks_test_m,
                  gene.density=rank.scores,
                  sort.idxs=sort.sgn.idxs,
                  mx.diff=mx.diff, abs.ranking=abs.ranking,
                  tau=tau, verbose=verbose,
                  BPPARAM=BPPARAM)
    m <- do.call("rbind", m)
    colnames(m) <- colnames(expr)

    return (m)
}


ks_test_m <- function(gset_idxs, gene.density, sort.idxs, mx.diff=TRUE,
                      abs.ranking=FALSE, tau=1, verbose=TRUE){
	
	n.genes <- nrow(gene.density)
	n.samples <- ncol(gene.density)
	n.geneset <- length(gset_idxs)

	geneset.sample.es <- .Call("ks_matrix_R",
			                       as.double(gene.density),
			                       as.integer(sort.idxs),
			                       n.genes,
			                       as.integer(gset_idxs),
			                       n.geneset,
			                       as.double(tau),
			                       n.samples,
			                       as.integer(mx.diff),
                             as.integer(abs.ranking))

	return(geneset.sample.es)
}


## ks-test in R code - testing only
ks_test_Rcode <- function(gene.density, gset_idxs, tau=1, make.plot=FALSE){
	
	n.genes = length(gene.density)
	n.gset = length(gset_idxs)
	
	sum.gset <- sum(abs(gene.density[gset_idxs])^tau)
	
	dec = 1 / (n.genes - n.gset)
	
	sort.idxs <- order(gene.density,decreasing=T)
	offsets <- sort(match(gset_idxs, sort.idxs))
	
	last.idx = 0
	values <- rep(NaN, length(gset_idxs))
	current = 0
	for(i in seq_along(offsets)){
		current = current + abs(gene.density[sort.idxs[offsets[i]]])^tau / sum.gset - dec * (offsets[i]-last.idx-1)
		
		values[i] = current
		last.idx = offsets[i]
	}
	check_zero = current - dec * (n.genes-last.idx)
	#if(check_zero > 10^-15){ 
	#	stop(paste=c("Expected zero sum for ks:", check_zero))
	#}
	if(make.plot){ plot(offsets, values,type="l") } 
	
	max.idx = order(abs(values),decreasing=T)[1]
	mx.value <- values[max.idx]
	
	return (mx.value)
}


.rndWalk <- function(gSetIdx, geneRanking, j, R, alpha) {
  indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  stepCDFinGeneSet <- cumsum((abs(R[geneRanking, j])^alpha * 
                      indicatorFunInsideGeneSet)) /
                      sum((abs(R[geneRanking, j])^alpha *
                      indicatorFunInsideGeneSet))
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet) /
                       sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

  sum(walkStat) 
}

## optimized version of the function .rndWalk by Alexey Sergushichev
## https://github.com/rcastelo/GSVA/pull/15
## based on his paper https://doi.org/10.1101/060012
.fastRndWalk <- function(gSetIdx, geneRanking, j, Ra) {
    n <- length(geneRanking)
    k <- length(gSetIdx)
    idxs <- sort.int(match(gSetIdx, geneRanking))
    
    stepCDFinGeneSet2 <- 
        sum(Ra[geneRanking[idxs], j] * (n - idxs + 1)) /
        sum((Ra[geneRanking[idxs], j]))    
    
    
    stepCDFoutGeneSet2 <- (n * (n + 1) / 2 - sum(n - idxs + 1)) / (n - k)
    
    walkStat <- stepCDFinGeneSet2 - stepCDFoutGeneSet2

    walkStat
}


ssgsea <- function(X, geneSets, alpha=0.25,
                   normalization=TRUE, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {

  n <- ncol(X)
  
  print("Calculating ranks...")
  
  R <- t(sparseMatrixStats::colRanks(X, ties.method="average"))
  mode(R) <- "integer"

  print("Calculating absolute values from ranks...")
  
  Ra <- abs(R)^alpha
  
  es <- bplapply(as.list(1:n), function(j) {
    geneRanking <- order(R[, j], decreasing=TRUE)
    es_sample <- lapply(geneSets, .fastRndWalk, geneRanking, j, Ra)
    
    unlist(es_sample)
  }, BPPARAM=BPPARAM)
  
  es <- do.call("cbind", es)
  
  
  if (normalization) {
    print("Normalizing...")
    ## normalize enrichment scores by using the entire data set, as indicated
    ## by Barbie et al., 2009, online methods, pg. 2
    es <- es[, 1:n, drop=FALSE] / (range(es)[2] - range(es)[1])
  }
  
  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)
  
  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)
  
  es
}

combinez <- function(gSetIdx, j, Z) sum(Z[gSetIdx, j]) / sqrt(length(gSetIdx))

zscore <- function(X, geneSets, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {
    if(is(X, "dgCMatrix")){
        .sparseScaleMessage()
        Z <- t(.sparseColumnApplyAndReplace(t(X), FUN=scale))
    } else {
        Z <- t(scale(t(X)))
    }

    es <- bplapply(geneSets, function(gSetIdx) {
        colSums(Z[gSetIdx, ]) / sqrt(length(gSetIdx))
    }, BPPARAM=BPPARAM)
    
    es <- do.call("rbind", es)
    es
}

rightsingularsvdvectorgset <- function(gSetIdx, Z) {
  if(is(Z, "dgCMatrix")){
    s <- BiocSingular::runExactSVD(Z[gSetIdx, ])
  } else {
    s <- svd(Z[gSetIdx, ])
  }
  s$v[, 1]
}

plage <- function(X, geneSets, verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {
    if(is(X, "dgCMatrix")){
        .sparseScaleMessage()
        Z <- t(.sparseColumnApplyAndReplace(t(X), FUN=scale))

        es <- bplapply(geneSets, rightsingularsvdvectorgset, Z,
                       BPPARAM=BPPARAM)
        
        es <- do.call(rbind, es)
    } else {
        
        Z <- t(scale(t(X)))
        
        es <- bplapply(geneSets, rightsingularsvdvectorgset, Z,
                       BPPARAM=BPPARAM)
        
        es <- do.call(rbind, es)
        
        ## why these extra steps for dense matrices?
        ## svd() removes all dimnames, while BiocSingular::runExactSVD()
        ## retains them.  colnames must be set; however, I don't see a reason
        ## to set rownames nor recreating the matrix if only one gene set.
        ## hmmm...
        if (length(geneSets) == 1)
            es <- matrix(es, nrow=1)
        
        rownames(es) <- names(geneSets)
        colnames(es) <- colnames(X)
    }
    
    es
}


#' @title Filter gene sets
#' 
#' @description Filters gene sets through a given minimum and maximum set size.
#' 
#' This function filters the input gene sets according to a given minimum and
#' maximum set size.
#' 
#' @aliases filterGeneSets
#'
#' @name filterGeneSets
#' 
#' @rdname filterGeneSets
#'
#' @param gSets Gene sets given either as a `list` or a
#' `GeneSetCollection` object.
#' 
#' @param minSize Minimum size.
#' 
#' @param maxSize Maximum size.
#' 
#' @return A collection of gene sets that meet the given minimum and maximum
#' set size.
#' 
#' @author J. Guinney
#' 
#' @seealso [`computeGeneSetsOverlap`]
#' 
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' 
#' @keywords Gene set
#' 
#' @examples
#' geneSets <- list(set1=as.character(1:4), set2=as.character(4:10))
#' filterGeneSets(geneSets, minSize=5)
NULL

#' @aliases filterGeneSets,list-method
#' @rdname filterGeneSets
#' @exportMethod filterGeneSets
setMethod("filterGeneSets", signature(gSets="list"),
          function(gSets, minSize=1, maxSize=Inf) {
	gSetsLen <- lengths(gSets)
	return (gSets[gSetsLen >= minSize & gSetsLen <= maxSize])	
})

#' @aliases filterGeneSets,GeneSetCollection-method
#' @rdname filterGeneSets
#' @exportMethod filterGeneSets
setMethod("filterGeneSets", signature(gSets="GeneSetCollection"),
          function(gSets, minSize=1, maxSize=Inf) {
  filterGeneSets(geneIds(gSets), minSize, maxSize)
})


#' @title Compute gene-sets overlap
#' 
#' @description Calculates the overlap among every pair of gene-sets given as
#' input.
#' 
#' This function calculates the overlap between every pair of gene sets of the
#' input argument `gSets`. Before this calculation takes place, the gene
#' sets in `gSets` are firstly filtered to discard genes that do not match
#' to the identifiers in `uniqGenes`. Secondly, they are further filtered
#' to meet the minimum and/or maximum size specified with the arguments
#' `minSize` and `maxSize`. The overlap between two gene sets is
#' calculated as the number of common genes between the two gene sets divided
#' by the smallest size of the two gene sets.
#' 
#' @aliases computeGeneSetsOverlap
#'
#' @name computeGeneSetsOverlap
#'
#' @rdname computeGeneSetsOverlap
#' 
#' @param gSets Gene sets given either as a `list` or a
#' `GeneSetCollection` object.
#' 
#' @param uniqGenes Vector of unique genes to be considered when calculating
#' the overlaps.
#' 
#' @param minSize Minimum size.
#' 
#' @param maxSize Maximum size.
#' 
#' @return A gene-set by gene-set matrix of the overlap among every pair of
#' gene sets.
#' 
#' @author J. Guinney
#' 
#' @seealso [`filterGeneSets`]
#' 
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' 
#' @keywords Gene set
#' 
#' @examples
#' geneSets <- list(set1=as.character(1:4), set2=as.character(4:10))
#' computeGeneSetsOverlap(geneSets, unique(unlist(geneSets)))
NULL

#' @aliases computeGeneSetsOverlap,list,character-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="character"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
  totalGenes <- length(uniqGenes)

  ## map to the actual features for which expression data is available
  gSets <- .mapGeneSetsToFeatures(gSets, uniqGenes)

  lenGsets <- lengths(gSets)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,list,ExpressionSet-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
  uniqGenes <- featureNames(uniqGenes)
  totalGenes <- length(uniqGenes)

  ## map to the actual features for which expression data is available
  gSets <- .mapGeneSetsToFeatures(gSets, uniqGenes)

  lenGsets <- lengths(gSets)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,GeneSetCollection,character-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="character"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,GeneSetCollection,ExpressionSet-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
  ## map gene identifiers of the gene sets to the features in the chip
  ## Biobase::annotation() is necessary to disambiguate from the
  ## 'annotation' argument
  gSets <- mapIdentifiers(gSets, AnnoOrEntrezIdentifier(Biobase::annotation(uniqGenes)))
  
  uniqGenes <- featureNames(uniqGenes)

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

.computeGeneSetsOverlap <- function(gSetsMembershipMatrix, minSize=1, maxSize=Inf) {
  ## gSetsMembershipMatrix should be a (genes x gene-sets) incidence matrix

  lenGsets <- colSums(gSetsMembershipMatrix)

  szFilterMask <- lenGsets >= max(1, minSize) & lenGsets <= maxSize
  if (!any(szFilterMask))
    stop("No gene set meets the minimum and maximum size filter\n")

  gSetsMembershipMatrix <- gSetsMembershipMatrix[, szFilterMask]
  lenGsets <- lenGsets[szFilterMask]

  totalGsets <- ncol(gSetsMembershipMatrix)

  M <- t(gSetsMembershipMatrix) %*% gSetsMembershipMatrix

  M1 <- matrix(lenGsets, nrow=totalGsets, ncol=totalGsets,
               dimnames=list(colnames(gSetsMembershipMatrix), colnames(gSetsMembershipMatrix)))
  M2 <- t(M1)
  M.min <- matrix(0, nrow=totalGsets, ncol=totalGsets)
  M.min[M1 < M2] <- M1[M1 < M2]
  M.min[M2 <= M1] <- M2[M2 <= M1]
  overlapMatrix <- M / M.min

  return (overlapMatrix)
}

## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.isPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
  (name %in% loadedNamespaces())
}

##
## ARE THESE FUNCTIONS STILL NECESSARY ?????
##

##a <- replicate(1000, compute.null.enrichment(10000,50,make.plot=F))

compute.null.enrichment <- function(n.genes, n.geneset, make.plot=FALSE){
	ranks <- (n.genes/2) - rev(1:n.genes)
	#null.gset.idxs <- seq(1, n.genes, by=round(n.genes / n.geneset))
	null.gset.idxs <- sample(n.genes, n.geneset)
	null.es <- ks_test_Rcode(ranks, null.gset.idxs,make.plot=make.plot)
	return (null.es)
}


load.gmt.data <- function(gmt.file.path){
	tmp <- readLines(gmt.file.path)
	gsets <- list()
	for(i in 1:length(tmp)){
		t <- strsplit(tmp[i],'\t')[[1]]
		gsets[[t[1]]] <- t[3:length(t)]
	}
	return (gsets)
}

compute.gset.overlap.score <- function(gset.idxs){
	n <- length(gset.idxs)
	mx.idx <- max(unlist(gset.idxs, use.names=F))
	l <- c(sapply(gset.idxs, length))
	
	gset.M <- matrix(0, nrow=mx.idx, ncol=n)
	for(i in 1:n){
		gset.M[gset.idxs[[i]],i] = 1
	}
	M <- t(gset.M) %*% gset.M
	
	M1 <- matrix(l, nrow=n, ncol=n)
	M2 <- t(M1)
	M.min <- matrix(0, nrow=n, ncol=n)
	M.min[M1 < M2] <- M1[M1 < M2]
	M.min[M2 <= M1] <- M2[M2 <= M1]
	M.score <- M / M.min
	return (M.score)
}
