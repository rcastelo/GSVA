##
## function: gsva
## purpose: main function of the package which estimates activity
##          scores for each given gene-set

setGeneric("gsva", function(expr, gset.idx.list, ...) standardGeneric("gsva"))

setMethod("gsva", signature(expr="ExpressionSet", gset.idx.list="list"),
          function(expr, gset.idx.list,
  abs.ranking=FALSE,
  min.sz=1,
  max.sz=Inf,
  no.bootstraps=0, 
  bootstrap.percent = .632, 
  parallel.sz=0, 
  parallel.type="SOCK",
  verbose=TRUE,
  mx.diff=TRUE)
{
  ## map to the actual features for which expression data is available
  mapped.gset.idx.list <- lapply(gset.idx.list,
                                 function(x, y) na.omit(match(x, y)),
                                 featureNames(expr))

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                         min.sz=max(1, min.sz),
                                         max.sz=max.sz)

  eSco <- GSVA:::.gsva(Biobase::exprs(expr), mapped.gset.idx.list, abs.ranking,
                       no.bootstraps, bootstrap.percent, parallel.sz, parallel.type,
                       verbose, mx.diff)
  eScoEset <- expr
  eScoEset <- Biobase::`exprs<-`(eScoEset, eSco$es.obs)
  eScoEset <- Biobase::`annotation<-`(eScoEset, "")

	return(list(es.obs=eScoEset,
				      bootstrap=eSco$bootstrap,
              p.vals.sign=eSco$p.vals.sign))
})

setMethod("gsva", signature(expr="ExpressionSet", gset.idx.list="GeneSetCollection"),
          function(expr, gset.idx.list,
  abs.ranking=FALSE,
  min.sz=1,
  max.sz=Inf,
  no.bootstraps=0, 
  bootstrap.percent = .632, 
  parallel.sz=0, 
  parallel.type="SOCK",
  verbose=TRUE,
  mx.diff=TRUE)
{
  if (verbose)
    cat("Mapping identifiers between gene sets and feature names\n")

  ## map gene identifiers of the gene sets to the features in the chip
  mapped.gset.idx.list <- mapIdentifiers(gset.idx.list,
                                         AnnotationIdentifier(annotation(expr)))
  
  ## map to the actual features for which expression data is available
  tmp <- lapply(geneIds(mapped.gset.idx.list),
                                 function(x, y) na.omit(match(x, y)),
                                 featureNames(expr))
  names(tmp) <- names(mapped.gset.idx.list)
  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  mapped.gset.idx.list <- filterGeneSets(tmp,
                                         min.sz=max(1, min.sz),
                                         max.sz=max.sz)

  eSco <- GSVA:::.gsva(Biobase::exprs(expr), mapped.gset.idx.list, abs.ranking,
                       no.bootstraps, bootstrap.percent, parallel.sz, parallel.type,
                       verbose, mx.diff)
  eScoEset <- expr
  eScoEset <- Biobase::`exprs<-`(eScoEset, eSco$es.obs)
  eScoEset <- Biobase::`annotation<-`(eScoEset, "")

	return(list(es.obs=eScoEset,
				      bootstrap=eSco$bootstrap,
              p.vals.sign=eSco$p.vals.sign))
})

setMethod("gsva", signature(expr="matrix", gset.idx.list="list"),
          function(expr, gset.idx.list,
  abs.ranking=FALSE,
  min.sz=1,
  max.sz=Inf,
  no.bootstraps=0, 
  bootstrap.percent = .632, 
  parallel.sz=0, 
  parallel.type="SOCK",
  verbose=TRUE,
  mx.diff=TRUE)
{
  mapped.gset.idx.list <- lapply(gset.idx.list,
                                 function(x ,y) na.omit(match(x, y)),
                                 rownames(expr))

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  mapped.gset.idx.list <- filterGeneSets(mapped.gset.idx.list,
                                         min.sz=max(1, min.sz),
                                         max.sz=max.sz)

  GSVA:::.gsva(expr, mapped.gset.idx.list, abs.ranking, no.bootstraps,
               bootstrap.percent, parallel.sz, parallel.type,
               verbose, mx.diff)
})

.gsva <- function(expr, gset.idx.list,
  abs.ranking=FALSE,
  no.bootstraps=0, 
  bootstrap.percent = .632, 
  parallel.sz=0, 
  parallel.type="SOCK",
  verbose=TRUE,
  mx.diff=TRUE)
{
	
	if(length(gset.idx.list) == 0){
		stop("The gene set list is empty!  Filter may be too stringent.")
	}
	
	if(verbose)
		cat("Testing", length(gset.idx.list),"gene sets.\n")
	
	if(parallel.sz > 0 && no.bootstraps > 0){
		if((no.bootstraps %% parallel.sz) != 0){
			stop("'parrallel.sz' must be an integer divisor of 'no.bootsraps'" )
		}
	}
	n.samples <- ncol(expr)
	n.genes <- nrow(expr)
	n.gset <- length(gset.idx.list)
	
	es.obs <- matrix(NaN, n.gset, n.samples, dimnames=list(names(gset.idx.list),colnames(expr)))
	colnames(es.obs) <- colnames(expr)
	rownames(es.obs) <- names(gset.idx.list)
	
	
	if (verbose)
    cat("Computing observed enrichment scores\n")
	es.obs <- compute.geneset.es(expr, gset.idx.list, 1:n.samples,
                               abs.ranking,parallel.sz,
                               parallel.type,verbose=verbose, mx.diff=mx.diff)
	
	# es.bootstraps -> n.gset by n.samples by n.resamples
	es.bootstraps=NULL
	p.vals.wilcoxon=NULL
	p.vals.sign=NULL
	
	if(no.bootstraps > 0){
		if(verbose) cat("Computing bootstrap enrichment scores\n")
		bootstrap.nsamples <- floor(bootstrap.percent * n.samples)
		
		p.vals.sign <- matrix(NaN, n.gset, n.samples,dimnames=list(names(gset.idx.list),colnames(expr)))
		
		es.bootstraps <- array(NaN, c(n.gset, n.samples, no.bootstraps))
		if(parallel.sz > 0){
			
		  if(!GSVA:::.isPackageLoaded("snow")) {
			  stop("Please load the 'snow' library")
		  }
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      clSetupRNG <- get("clusterSetupRNG", mode="function")
      clEvalQ <- get("clusterEvalQ", mode="function")
      clExport <- get("clusterExport", mode="function")
      stopCl <- get("stopCluster", mode="function")
			
			cl <- makeCl(parallel.sz, type = parallel.type) 
			.GlobalEnv[["expr"]] <- expr
			.GlobalEnv[["bootstrap.nsamples"]] <- bootstrap.nsamples
			.GlobalEnv[["n.samples"]] <- n.samples
			.GlobalEnv[["gset.idx.list"]] <- gset.idx.list
			clExport(cl,"expr")
			clExport(cl,"bootstrap.nsamples")
			clExport(cl, "n.samples")
			clExport(cl, "gset.idx.list")
			clEvalQ(cl, library(GSVA))
			
			clSetupRNG(cl)
			
			if(verbose) cat("Parallel bootstrap...\n")
			## parallelized bootstrap
			n.cycles <- floor(no.bootstraps / parallel.sz)
			for(i in 1:n.cycles){
				if(verbose) cat("bootstrap cycle ", i, "\n")
				r <- clEvalQ(cl, compute.geneset.es(expr, gset.idx.list, 
								sample(n.samples, bootstrap.nsamples, replace=T),
								abs.ranking))
				for(j in 1:length(r)){
					es.bootstraps[,,(parallel.sz * (i-1) + j)] <- r[[j]]
				}	
			}
			stopCl(cl)
		}else{
			if(verbose) cat("Sequential bootstrap...\n")
			for(i in 1:no.bootstraps){
				es.bootstraps[,,i] <- compute.geneset.es(expr, gset.idx.list,
						sample(n.samples, bootstrap.nsamples, replace=T),
						abs.ranking)
			}
		}
	
		
		for(i in 1:n.gset){
			
			for(j in 1:n.samples){
				# non-parametric test if median of empirical dist is 0 
				if(es.obs[i,j] > 0){
					p.vals.sign[i,j] <- (1 + sum(es.bootstraps[i,j,] < 0)) / (1 + no.bootstraps)
				}else{
					p.vals.sign[i,j] <- (1 + sum(es.bootstraps[i,j,] > 0)) / (1 + no.bootstraps)
				}
			}
		}
	}
	
	colnames(es.obs) <- colnames(expr)
	rownames(es.obs) <- names(gset.idx.list)
	return(list(es.obs=es.obs,
				      bootstrap=list(es.bootstraps=es.bootstraps,
              p.vals.sign=p.vals.sign)))
}


compute.gene.density <- function(expr, sample.idxs){
	n.test.samples <- ncol(expr)
	n.genes <- nrow(expr)
	n.density.samples <- length(sample.idxs)
	
	A = .C("assess_matrix_density_R",
			as.double(t(expr[,sample.idxs])),
			as.double(t(expr)),
			R = double(n.test.samples * n.genes),
			n.density.samples,
			n.test.samples,
			n.genes)$R
	
	gene.density <- t(matrix(A, n.test.samples, n.genes))
	return (gene.density)	
}

compute.geneset.es <- function(expr, gset.idx.list, sample.idxs, abs.ranking,
                               parallel.sz=0, parallel.type="SOCK",
                               verbose=FALSE, mx.diff){
	num_genes <- nrow(expr)
	if(verbose) cat("Computing gene densities\n")
	gene.density <- compute.gene.density(expr, sample.idxs)
	
	compute_rank_score <- function(sort_idx_vec){
		tmp <- rep(0, num_genes)
		tmp[sort_idx_vec] <- abs(seq(from=num_genes,to=1) - num_genes/2)
		return (tmp)
	}
	
	rank.scores <- rep(0, num_genes)
	if(abs.ranking){
		sort.sgn.idxs <- apply(abs(gene.density), 2, order, decreasing=TRUE) # n.genes * n.samples	
	}else{
		sort.sgn.idxs <- apply(gene.density, 2, order, decreasing=TRUE) # n.genes * n.samples
	}
	
	rank.scores <- apply(sort.sgn.idxs, 2, compute_rank_score)
	
	haveMulticore <- GSVA:::.isPackageLoaded("multicore")
	haveSnow <- GSVA:::.isPackageLoaded("snow")
	
	if(parallel.sz > 0 || haveMulticore) {
		if(!haveMulticore && !haveSnow) {
			stop("In order to run calculations in parallel either the 'snow', or the 'multicore' library, should be loaded first")
		}

    if (!haveMulticore) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      parSapp <- get("parSapply", mode="function")
      clEvalQ <- get("clusterEvalQ", mode="function")
      stopCl <- get("stopCluster", mode="function")

      if (verbose)
        cat("Allocating cluster\n")
		  cl <- makeCl(parallel.sz, type = parallel.type) 
		  clEvalQ(cl, library(GSVA))
		  if (verbose) {
		    cat("Evaluating parallel ks-tests...\n")
	      if(mx.diff) {
          cat("Taking diff of max KS.\n")
        } else{
          cat("Evaluting max KS.\n")
        }
      }
	
		  m <- t(parSapp(cl, gset.idx.list, ks_test_m, 
						  gene.density=rank.scores, 
						  sort.idxs=sort.sgn.idxs,
						  mx.diff=mx.diff, verbose=FALSE))
		  if(verbose)
        cat("Cleaning up\n")
		  stopCl(cl)

    } else {             ## use multicore

      mclapp <- get('mclapply', envir=getNamespace('multicore'))
      detCor <- get('detectCores', envir=getNamespace('multicore'))
      nCores <- detCor()
      options(cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(cores=parallel.sz)

      pb <- NULL
      if (verbose){
        cat("Using multicore with", getOption("cores"), "cores\n")
        assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
        assign("nGeneSets", ceiling(length(gset.idx.list) / getOption("cores")), envir=globalenv())
        assign("iGeneSet", 0, envir=globalenv())
      }

      m <- mclapp(gset.idx.list, GSVA:::ks_test_m,
                  gene.density=rank.scores,
                  sort.idxs=sort.sgn.idxs,
                  mx.diff=mx.diff, verbose=verbose)
      m <- do.call("rbind", m)
      colnames(m) <- colnames(expr)

      if (verbose) {
        close(get("progressBar", envir=globalenv()))
      }
    }

	} else {
		if (verbose) {
      cat("Evaluating ks-tests\n")
	    if (mx.diff) {
        cat("Taking diff of max KS.\n")
      } else{
        cat("Evaluting max KS.\n")
      }
    }
    pb <- NULL
    if (verbose){
      assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
      assign("nGeneSets", length(gset.idx.list), envir=globalenv())
      assign("iGeneSet", 0, envir=globalenv())
    }

		m <- t(sapply(gset.idx.list, ks_test_m, rank.scores, sort.sgn.idxs,
                  mx.diff=mx.diff, verbose=verbose))

    if (verbose) {
      close(get("progressBar", envir=globalenv()))
    }
	}
	return (m)
}


ks_test_m <- function(gset_idxs, gene.density, sort.idxs, tau.factor=1, mx.diff=TRUE, verbose=TRUE){
	
	n.genes <- nrow(gene.density)
	n.samples <- ncol(gene.density)
	n.geneset <- length(gset_idxs)

	geneset.sample.es = .C("ks_matrix_R",
			as.double(gene.density),
			R = double(n.samples),
			as.integer(sort.idxs),
			n.genes,
			as.integer(gset_idxs),
			n.geneset,
			as.double(tau.factor),
			n.samples,
			as.integer(mx.diff))$R

  if (verbose) {
    assign("iGeneSet", get("iGeneSet", envir=globalenv()) + 1, envir=globalenv())
    setTxtProgressBar(get("progressBar", envir=globalenv()),
                      get("iGeneSet", envir=globalenv()) / get("nGeneSets", envir=globalenv()))
  }
	
	return (geneset.sample.es)
}


## ks-test in R code - testing only
ks_test_Rcode <- function(gene.density, gset_idxs, tau.factor=1, make.plot=FALSE){
	
	n.genes = length(gene.density)
	n.gset = length(gset_idxs)
	
	sum.gset <- sum(abs(gene.density[gset_idxs])^tau.factor)
	
	dec = 1 / (n.genes - n.gset)
	
	sort.idxs <- order(gene.density,decreasing=T)
	offsets <- sort(match(gset_idxs, sort.idxs))
	
	last.idx = 0
	values <- rep(NaN, length(gset_idxs))
	current = 0
	for(i in seq_along(offsets)){
		current = current + abs(gene.density[sort.idxs[offsets[i]]])^tau.factor / sum.gset - dec * (offsets[i]-last.idx-1)
		
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

setGeneric("filterGeneSets", function(gSets, ...) standardGeneric("filterGeneSets"))

setMethod("filterGeneSets", signature(gSets="list"),
          function(gSets, min.sz=1, max.sz=Inf) {
	gSetsLen <- sapply(gSets,length)
	return (gSets[gSetsLen >= min.sz & gSetsLen <= max.sz])	
})

setMethod("filterGeneSets", signature(gSets="GeneSetCollection"),
          function(gSets, min.sz=1, max.sz=Inf) {
	gSetsLen <- sapply(geneIds(gSets),length)
	return (gSets[gSetsLen >= min.sz & gSetsLen <= max.sz])	
})



setGeneric("computeGeneSetsOverlap", function(gSets, uniqGenes=unique(unlist(gSets, use.names=FALSE)), ...) standardGeneric("computeGeneSetsOverlap"))

setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="character"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  totalGenes <- length(uniqGenes)

  ## map to the features requested
  gSets <- lapply(gSets, function(x, y) as.vector(na.omit(match(x, y))), uniqGenes)

  lenGsets <- sapply(gSets, length)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  GSVA:::.computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  uniqGenes <- Biobase::featureNames(uniqGenes)
  totalGenes <- length(uniqGenes)

  ## map to the actual features for which expression data is available
  gSets <- lapply(gSets, function(x, y) as.vector(na.omit(match(x, y))), uniqGenes)

  lenGsets <- sapply(gSets, length)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  GSVA:::.computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="character"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  GSVA:::.computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  ## map gene identifiers of the gene sets to the features in the chip
  gSets <- mapIdentifiers(gSets, AnnotationIdentifier(annotation(uniqGenes)))
  
  uniqGenes <- Biobase::featureNames(uniqGenes)

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  GSVA:::.computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})

.computeGeneSetsOverlap <- function(gSetsMembershipMatrix, min.sz=1, max.sz=Inf) {
  ## gSetsMembershipMatrix should be a (genes x gene-sets) incidence matrix

  lenGsets <- colSums(gSetsMembershipMatrix)

  szFilterMask <- lenGsets >= max(1, min.sz) & lenGsets <= max.sz
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
