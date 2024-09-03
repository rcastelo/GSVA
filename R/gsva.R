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


compute.gene.cdf <- function(expr, sample.idxs, rnaseq=FALSE, kernel=TRUE,
                             sparse=FALSE, verbose=TRUE) {
    n.test.samples <- ncol(expr)
    n.genes <- nrow(expr)
    n.density.samples <- length(sample.idxs)
    
    gene.cdf <- NA
    if (kernel) {
        if (is(expr, "dgCMatrix")) {
            if (sparse)
                gene.cdf <- .kcdfvals_sparse_to_sparse(expr[, sample.idxs, drop=FALSE],
                                                       !rnaseq, verbose)
            else
                gene.cdf <- .kcdfvals_sparse_to_dense(expr[, sample.idxs, drop=FALSE],
                                                      !rnaseq, verbose)
        } else if (is.matrix(expr)) {
            A = .Call("matrix_density_R",
                      as.double(t(expr[ ,sample.idxs, drop=FALSE])),
                      as.double(t(expr)),
                      n.density.samples,
                      n.test.samples,
                      n.genes,
                      as.integer(rnaseq),
                      verbose)
            gene.cdf <- t(matrix(A, n.test.samples, n.genes))
        } else
            stop(sprintf("Matrix class %s cannot be handled yet.", class(expr)))
    } else {
        if (is(expr, "dgCMatrix")) {
            if (sparse)
                gene.cdf <- .ecdfvals_sparse_to_sparse(expr[, sample.idxs, drop=FALSE],
                                                       verbose)
            else
                gene.cdf <- .ecdfvals_sparse_to_dense(expr[, sample.idxs, drop=FALSE],
                                                      verbose)
        } else if (is.matrix(expr))
            gene.cdf <- .ecdfvals_dense_to_dense(expr[, sample.idxs, drop=FALSE],
                                                 verbose)
        else
            stop(sprintf("Matrix class %s cannot be handled yet.", class(expr)))
    }

    return(gene.cdf)	
}

zorder_rankstat <- function(z, p) {
  ## calculation of the ranks by expression-level statistic
  zord <- apply(z, 2, order, decreasing=TRUE)

  ## calculation of the rank-order statistic
  zrs <- apply(zord, 2, function(x, p)
               do.call("[<-", list(rep(0, p), x, abs(p:1-p/2))), p)

  list(Zorder=zord, ZrankStat=zrs)
}

.gsvaRndWalk <- function(gSetIdx, geneRanking, rankStat) {
    n <- length(geneRanking)
    k <- length(gSetIdx)

    stepCDFinGeneSet <- integer(n)
    stepCDFinGeneSet[gSetIdx] <- rankStat[geneRanking[gSetIdx]]
    stepCDFinGeneSet <- cumsum(stepCDFinGeneSet)
    stepCDFinGeneSet <- stepCDFinGeneSet / stepCDFinGeneSet[n]

    stepCDFoutGeneSet <- rep(1L, n)
    stepCDFoutGeneSet[gSetIdx] <- 0L
    stepCDFoutGeneSet <- cumsum(stepCDFoutGeneSet)
    stepCDFoutGeneSet <- stepCDFoutGeneSet / stepCDFoutGeneSet[n]

    walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

    walkStat
}

#' @importFrom Matrix nnzero
.sufficient_ssize <- function(expr, kcdf.min.ssize) {
  ## in the sparse case stored in a 'dgCMatrix', by now,
  ## use the average nonzero values per row
  if (is(expr, "dgCMatrix"))
    return((nnzero(expr) / nrow(expr)) >= kcdf.min.ssize)

  ## in every other case, including the dense case, by now,
  ## just look at the number of columns
  return(ncol(expr) >= kcdf.min.ssize)
}

#' @importFrom parallel splitIndices
#' @importFrom IRanges IntegerList match
#' @importFrom cli cli_alert_info cli_progress_bar
#' @importFrom cli cli_progress_update cli_progress_done
compute.geneset.es <- function(expr, gset.idx.list, sample.idxs, kcdf,
                               kcdf.min.ssize, abs.ranking, parallel.sz=1L,
                               mx.diff=TRUE, tau=1, sparse=FALSE,
                               verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose)) {
    num_genes <- nrow(expr)

    kernel <- rnaseq <- FALSE
    if (kcdf == "auto") {
        if (verbose)
            cli_alert_info("kcdf='auto' (default)")
        if (!.sufficient_ssize(expr, kcdf.min.ssize)) {
            kernel <- TRUE
            if (is(expr, "dgCMatrix")) { ## dgCMatrix does not store integers
                                         ## so we check them with x == floor(x)
                sam <- sample(expr@x, size=min(1000, length(expr@x)),
                              replace=FALSE)
                rnaseq <- all(sam == floor(sam))
            } else if (is.integer(expr[1, 1]))
                rnaseq <- TRUE
        }
    } else {
        if (kcdf == "Gaussian") {
            kernel <- TRUE
            rnaseq <- FALSE
        } else if (kcdf == "Poisson") {
            kernel <- TRUE
            rnaseq <- TRUE
        } else
            kernel <- FALSE
    }
    if (verbose) {
        if (kernel) {
            if (rnaseq)
                cli_alert_info("Row-wise ECDF estimation with Poisson kernels")
            else
                cli_alert_info("Row-wise ECDF estimation with Gaussian kernels")
        } else
            cli_alert_info("Direct row-wise ECDFs estimation")
    }

    ## open parallelism only if ECDFs have to be estimated for
    ## more than 100 genes on more than 100 samples
    if (parallel.sz > 1 && length(sample.idxs) > 100 && nrow(expr) > 100) {
        iter <- function(Y, idpb, n_chunks=BiocParallel::multicoreWorkers()) {
            idx <- splitIndices(nrow(Y), min(nrow(Y), n_chunks))
            i <- 0L
            function() {
                if (i == length(idx))
                    return(NULL)
                i <<- i + 1L
                if (!is.null(idpb))
                    cli_progress_update(id=idpb, set=i)
                Y[idx[[i]], , drop=FALSE]
            }
        }
        if (verbose) {
            msg <- sprintf("Estimating ECDFs with %d cores",
                           as.integer(parallel.sz))
            idpb <- cli_progress_bar(msg, total=100)
        }
        gene.density <- bpiterate(iter(expr, idpb, 100),
                                  compute.gene.cdf,
                                  sample.idxs=sample.idxs,
                                  rnaseq=rnaseq, kernel=kernel,
                                  sparse=sparse, verbose=FALSE,
                                  REDUCE=rbind, reduce.in.order=TRUE,
                                  BPPARAM=BPPARAM)
        if (verbose)
            cli_progress_done(idpb)
    } else
        gene.density <- compute.gene.cdf(expr, sample.idxs, rnaseq, kernel,
                                         sparse, verbose)
    
    gset.idx.list <- IntegerList(gset.idx.list)
    n <- ncol(expr)
    es <- NULL
    if (n > 10 && bpnworkers(BPPARAM) > 1) {
        es <- bplapply(as.list(1:n), function(j, Z) {
            gene_ord_rnkstat <- list()
            if (is(Z, "dgCMatrix"))
                gene_ord_rnkstat <- .order_rankstat_sparse_to_sparse(Z, j)
            else
                gene_ord_rnkstat <- .order_rankstat(Z[, j])
            geneRanking <- gene_ord_rnkstat[[1]]
            rankStat <- gene_ord_rnkstat[[2]]

            geneSetsRankIdx <- match(gset.idx.list, geneRanking)
            .gsva_score_genesets(as.list(geneSetsRankIdx), geneRanking, rankStat,
                                 mx.diff, abs.ranking, tau)
        }, Z=gene.density, BPPARAM=BPPARAM)
    } else {
        idpb <- NULL
        if (verbose)
            idpb <- cli_progress_bar("Calculating GSVA scores", total=n)
        es <- lapply(as.list(1:n), function(j, Z) {
            gene_ord_rnkstat <- list()
            if (is(Z, "dgCMatrix"))
                gene_ord_rnkstat <- .order_rankstat_sparse_to_sparse(Z, j)
            else
                gene_ord_rnkstat <- .order_rankstat(Z[, j])
            geneRanking <- gene_ord_rnkstat[[1]]
            rankStat <- gene_ord_rnkstat[[2]]

            geneSetsRankIdx <- match(gset.idx.list, geneRanking)
            sco <- .gsva_score_genesets(as.list(geneSetsRankIdx), geneRanking, rankStat,
                                        mx.diff, abs.ranking, tau)
            if (verbose)
                cli_progress_update(id=idpb)
            sco
        }, Z=gene.density)
        if (verbose)
            cli_progress_done(idpb)
    }
    es <- do.call("cbind", es)

    return(es)
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
## with further optimizations with Bob Policastro discussed in
## https://github.com/rcastelo/GSVA/issues/71
.fastRndWalk <- function(gSetIdx, geneRanking, j, Ra) {
    n <- length(geneRanking)
    k <- length(gSetIdx)
    
    stepCDFinGeneSet <- 
        sum(Ra[geneRanking[gSetIdx], j] * (n - gSetIdx + 1)) /
        sum((Ra[geneRanking[gSetIdx], j]))    
    
    stepCDFoutGeneSet <- (n * (n + 1) / 2 - sum(n - gSetIdx + 1)) / (n - k)
    
    walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

    walkStat
}

.fastRndWalkNArm <- function(gSetIdx, geneRanking, j, Ra,
                             anyna=FALSE, na_use="everything", minSize=1,
                             wna_env=new.env()) {
    n <- length(geneRanking)
    if (anyna && na_use == "na.rm")
        gSetIdx <- na.omit(gSetIdx)
    k <- length(gSetIdx)
    
    walkStat <- NA_real_
    if (k >= minSize) {
      stepCDFinGeneSet <- 
          sum(Ra[geneRanking[gSetIdx], j] * (n - gSetIdx + 1)) /
          sum((Ra[geneRanking[gSetIdx], j]))    
    
      stepCDFoutGeneSet <- (n * (n + 1) / 2 - sum(n - gSetIdx + 1)) / (n - k)
    
      walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
    } else if (!get("w", envir=wna_env)) ## warn only once. it can only happen
      assign("w", TRUE, envir=wna_env)   ## with anyna=TRUE and na_use="na.rm"

    walkStat
}

#' @importFrom IRanges IntegerList match
ssgsea <- function(X, geneSets, alpha=0.25,
                   normalization=TRUE,
                   check_na=FALSE,
                   any_na=FALSE,
                   na_use=c("everything", "all.obs", "na.rm"),
                   minSize=1, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {
  na_use <- match.arg(na_use)

  n <- ncol(X)
  
  if (verbose)
    print("Calculating ranks...")
  
  R <- t(sparseMatrixStats::colRanks(X, ties.method="average"))
  mode(R) <- "integer"

  if (verbose)
    print("Calculating absolute values from ranks...")
  
  Ra <- abs(R)^alpha
  
  wna_env <- new.env()
  assign("w", FALSE, envir=wna_env)
  geneSets <- IntegerList(geneSets)
  es <- bplapply(as.list(1:n), function(j) {
    if (any_na && na_use == "na.rm") {
      geneRanking <- order(R[, j], decreasing=TRUE, na.last=NA)
      geneSetsRankIdx <- match(geneSets, geneRanking)
      es_sample <- sapply(geneSetsRankIdx, .fastRndWalkNArm,
                          geneRanking, j, Ra, any_na, na_use,
                          minSize, wna_env, USE.NAMES=FALSE)
    } else {
      geneRanking <- order(R[, j], decreasing=TRUE)
      geneSetsRankIdx <- match(geneSets, geneRanking)
      es_sample <- sapply(geneSetsRankIdx, .fastRndWalk,
                          geneRanking, j, Ra)
    }
    es_sample
  }, BPPARAM=BPPARAM)
  es <- do.call("cbind", es)
  
  if (any_na && na_use =="na.rm")
    if (get("w", envir=wna_env))
      warning(sprintf(paste("one or more gene sets had less than %d genes after",
                            "after removing NA values, the corresponding",
                            "enrichment scores have been set to NA."), minSize))
  
  if (normalization) {
    if (verbose)
      cat("Normalizing...\n")
    ## normalize enrichment scores by using the entire data set, as indicated
    ## by Barbie et al., 2009, online methods, pg. 2

    ## consider calculating the range on the fly to avoid having to do this
    ## on a large matrix
    rng <- NULL
    if (any_na)
      rng <- range(es, na.rm=TRUE) ## discard always NA values to calculate the
                                   ## normalization factor to enable either the
                                   ## propogation of NA values or to avoid them
    else
      rng <- range(es) ## na.rm increases execution time and memory consumption

    if (any(is.na(rng) | !is.finite(rng)))
      stop(paste("Cannot calculate normalizing factor for the enrichment scores in",
                 "ssGSEA, likely due to NA values in the input expression data."))
    es <- es[, 1:n, drop=FALSE] / (rng[2] - rng[1])
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
