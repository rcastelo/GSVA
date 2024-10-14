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


compute.gene.cdf <- function(expr, sample.idxs, Gaussk=TRUE, kernel=TRUE,
                             sparse=FALSE, verbose=TRUE) {
    n.test.samples <- ncol(expr)
    n.genes <- nrow(expr)
    n.density.samples <- length(sample.idxs)
    
    gene.cdf <- NA
    if (kernel) {
        if (is(expr, "dgCMatrix")) {
            if (sparse)
                gene.cdf <- .kcdfvals_sparse_to_sparse(expr[, sample.idxs, drop=FALSE],
                                                       Gaussk, verbose)
            else
                gene.cdf <- .kcdfvals_sparse_to_dense(expr[, sample.idxs, drop=FALSE],
                                                      Gaussk, verbose)
        } else if (is.matrix(expr)) {
            A = .Call("matrix_density_R",
                      as.double(t(expr[ ,sample.idxs, drop=FALSE])),
                      as.double(t(expr)),
                      n.density.samples,
                      n.test.samples,
                      n.genes,
                      as.integer(Gaussk),
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

.parse_kcdf_param <- function(expr, kcdf, kcdf.min.ssize, sparse, verbose) {
    kernel <- FALSE
    Gaussk <- TRUE  ## default (TRUE) is a Gaussian kernel, Poisson otherwise (FALSE)
    if (kcdf == "auto") {
        if (verbose)
            cli_alert_info("kcdf='auto' (default)")
        if (!.sufficient_ssize(expr, kcdf.min.ssize)) {
            kernel <- TRUE
            if (is(expr, "dgCMatrix")) { ## dgCMatrix does not store integers
                                         ## so we check them with x == floor(x)
                sam <- sample(expr@x, size=min(1000, length(expr@x)),
                              replace=FALSE)
                Gaussk <- any((sam < 0) | (sam != floor(sam)))
            } else if (is.integer(expr[1, 1]))
                Gaussk <- FALSE
        }
    } else {
        if (kcdf == "Gaussian") {
            kernel <- TRUE
            Gaussk <- TRUE
        } else if (kcdf == "Poisson") {
            kernel <- TRUE
            Gaussk <- FALSE
        } else
            kernel <- FALSE
    }

    if (verbose) {
        if (is(expr, "dgCMatrix") && sparse)
            cli_alert_info("GSVA sparse algorithm")
        else
            cli_alert_info("GSVA dense (classical) algorithm")
        if (kernel) {
            if (Gaussk)
                cli_alert_info("Row-wise ECDF estimation with Gaussian kernels")
            else
                cli_alert_info("Row-wise ECDF estimation with Poisson kernels")
        } else
            cli_alert_info("Direct row-wise ECDFs estimation")
    }

    list(kernel=kernel, Gaussk=Gaussk)
}

#' @importFrom parallel splitIndices
#' @importFrom IRanges IntegerList match
#' @importFrom BiocParallel bpnworkers
#' @importFrom cli cli_alert_info cli_progress_bar
#' @importFrom cli cli_progress_update cli_progress_done
compute.geneset.es <- function(expr, gset.idx.list, sample.idxs, kcdf,
                               kcdf.min.ssize, abs.ranking, parallel.sz=1L,
                               mx.diff=TRUE, tau=1, sparse=FALSE,
                               verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose)) {

    kcdfparam <- .parse_kcdf_param(expr, kcdf, kcdf.min.ssize, sparse, verbose)
    kernel <- kcdfparam$kernel
    Gaussk <- kcdfparam$Gaussk

    ## open parallelism only if ECDFs have to be estimated for
    ## more than 100 genes on more than 100 samples
    if (parallel.sz > 1 && length(sample.idxs) > 100 && nrow(expr) > 100) {
        iter <- function(Y, idpb, n_chunks=bpnworkers(BPPARAM)) {
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
                                  Gaussk=Gaussk, kernel=kernel,
                                  sparse=sparse, verbose=FALSE,
                                  REDUCE=rbind, reduce.in.order=TRUE,
                                  BPPARAM=BPPARAM)
        if (verbose)
            cli_progress_done(idpb)
    } else
        gene.density <- compute.gene.cdf(expr, sample.idxs, Gaussk, kernel,
                                         sparse, verbose)
    
    gset.idx.list <- IntegerList(gset.idx.list)
    n <- ncol(expr)
    es <- NULL
    if (n > 10 && bpnworkers(BPPARAM) > 1) {
        if (verbose) {
            msg <- sprintf("Calculating GSVA scores with %d cores",
                           as.integer(parallel.sz))
            cli_alert_info(msg)
        }
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

##
## implement calculations separating ranks from random walks
##

## BEGIN exported methods (to be moved to 'gsvaNewAPI.R')

#' @title GSVA ranks and scores
#'
#' @description Calculate GSVA scores in two steps: (1) calculate GSVA
#' ranks; and (2) calculate GSVA scores using the previously calculated
#' ranks.
#'
#' @param param A [`gsvaParam`] object built using the constructor function
#' [`gsvaParam`].
#'
#' @param verbose Gives information about each calculation step. Default: `TRUE`.
#'
#' @param BPPARAM An object of class [`BiocParallelParam`] specifying parameters
#'   related to the parallel execution of some of the tasks and calculations
#'   within this function.
#'
#' @return In the case of the `gsvaRanks()` method, a matrix of GSVA rank
#' values per column.
#'
#' @seealso [`gsvaParam`], [`gsva`]
#'
#' @aliases gsvaRanks,gsvaParam-method
#' @name gsvaRanks
#' @rdname gsvaRanks
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' [DOI](https://doi.org/10.1186/1471-2105-14-7)
#'
#' @examples
#' library(GSVA)
#'
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#' nGrp1 <- 15 ## number of samples in group 1
#' nGrp2 <- n - nGrp1 ## number of samples in group 2
#'
#' ## consider three disjoint gene sets
#' geneSets <- list(gset1=paste0("g", 1:3),
#'                  gset2=paste0("g", 4:6),
#'                  gset3=paste0("g", 7:10))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate GSVA ranks
#' gsva_ranks <- gsvaRanks(gsvapar)
#' gsva_ranks
#' ## calculate GSVA scores
#' gsva_es <- gsvaScores(gsvapar, gsva_ranks)
#' gsva_es
#'
#' ## calculate now GSVA scores in a single step
#' gsva_es1 <- gsva(gsvapar)
#'
#' ## both approaches give the same result with the same input gene sets
#' all.equal(gsva_es1, gsva_es)
#'
#' ## however, results will be (obviously) different with different gene sets
#' geneSets2 <- list(gset1=paste0("g", 3:6),
#'                   gset2=paste0("g", c(1, 2, 7, 8)))
#'
#' ## note that there is no need to calculate the GSVA ranks again
#' gsvaScores(gsvapar, gsva_ranks, geneSets2)
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @exportMethod gsvaRanks
setMethod("gsvaRanks", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              filteredDataMatrix <- .filterGenes(dataMatrix, removeConstant=TRUE,
                                                 removeNzConstant=TRUE)

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if (verbose)
                  cli_alert_info(sprintf("Calculating GSVA ranks"))

              psz <- if(inherits(BPPARAM, "SerialParam")) 1L else bpnworkers(BPPARAM)

              gsva_rnk <- .compute_gsva_ranks(expr=filteredDataMatrix,
                                              kcdf=get_kcdf(param),
                                              kcdf.min.ssize=get_kcdfNoneMinSampleSize(param),
                                              sparse=get_sparse(param),
                                              verbose=verbose,
                                              BPPARAM=BPPARAM)

              rownames(gsva_rnk) <- rownames(filteredDataMatrix)
              colnames(gsva_rnk) <- colnames(filteredDataMatrix)

              if (verbose)
                  cli_alert_success("Calculations finished")

              return(gsva_rnk)
          })

.check_geneSets_minSize_maxSize_tau <- function(geneSets, minSize, maxSize, tau) {
  if (all(!is.na(geneSets))) {
      if (!is.list(geneSets) && !is(geneSets, "GeneSetCollection"))
          cli_abort(c("x"="'geneSets' must be either a list or a 'GeneSetCollection' object"))

      if (length(geneSets) == 0)
          cli_abort(c("x"="'geneSets' has length 0"))
  }

  if (length(minSize) != 1)
      cli_abort(c("x"="'minSize' must be of length 1"))
  if (length(maxSize) != 1)
      cli_abort(c("x"="'maxSize' must be of length 1"))

  if ((is.na(minSize) && !is.na(maxSize)) || ## here assuming length 'minSize' and 'maxSize'
      (!is.na(minSize) && is.na(maxSize)))   ## is 1, otherwise 'is.na()' would return > 1 value
      cli_abort(c("x"="'minSize' and 'maxSize' should be either both NA or both non-NA"))

  if (!is.na(minSize) && !is.na(maxSize)) {
      if (!is.integer(minSize) && !is.numeric(minSize))
          cli_abort(c("x"="'minSize' must be a positive integer value"))
      if (!is.integer(maxSize) && !is.numeric(maxSize))
          cli_abort(c("x"="'maxSize' must be a positive integer value"))
      if (minSize < 1)
          cli_abort(c("x"="'minSize' must be a positive integer value"))
      if (maxSize < 1)
          cli_abort(c("x"="'maxSize' must be a positive integer value"))
      if (maxSize < minSize)
          cli_abort(c("x"="'maxSize' must be at least 'minSize' or greater"))
  }

  if (length(tau) != 1)
      cli_abort(c("x"="'tau' must be of length 1"))
  if (!is.na(tau)) {
    if (!is.integer(tau) && !is.numeric(tau))
          cli_abort(c("x"="'tau' must be a numeric value"))
  }
}

.check_maxDiff_absRanking_sparse <- function(maxDiff, absRanking, sparse) {
  if (length(maxDiff) != 1)
      cli_abort(c("x"="'maxDiff' must be of length 1"))

  if (!is.na(maxDiff)) {
    if (!is.logical(maxDiff))
          cli_abort(c("x"="'maxDiff' must be a logical value"))
  }

  if (length(absRanking) != 1)
      cli_abort(c("x"="'absRanking' must be of length 1"))

  if (!is.na(absRanking)) {
    if (!is.logical(absRanking))
          cli_abort(c("x"="'absRanking' must be a logical value"))
  }

  if (length(sparse) != 1)
      cli_abort(c("x"="'sparse' must be of length 1"))

  if (!is.na(sparse)) {
    if (!is.logical(sparse))
          cli_abort(c("x"="'sparse' must be a logical value"))
  }
}


#' @param ranks A matrix-like object storing GSVA ranks calculated with the
#' method [`gsvaRanks`].
#'
#' @param geneSets A collection of gene sets. Must be one of the classes
#' supported by [`GsvaGeneSets-class`]. For a list of these classes, see its
#' help page using `help(GsvaGeneSets)`. By default, this parameter is set to
#' the `NA` missing value, which means that GSVA scores will be calculated
#' using the gene sets specified in the `param` argument. If this parameter is
#' set to a non-missing value corresponding to an object of the classes
#' supported by [`GsvaGeneSets-class`], then GSVA scores will be calculated
#' using the gene sets in this argument, instead of the ones specified in the
#' `param` argument.
#'
#' @param minSize Numeric vector of length 1.  Minimum size of the resulting gene
#' sets after gene identifier mapping. Its default value is `NA`, indicating that
#' this minimum value will be taken from the input `param` argument, otherwise,
#' non-`NA` values override those from the input `param` argument.
#'
#' @param maxSize Numeric vector of length 1.  Minimum size of the resulting gene
#' sets after gene identifier mapping. Its default value is `NA`, indicating that
#' this minimum value will be taken from the input `param` argument, otherwise,
#' non-`NA` values override those from the input `param` argument.
#'
#' @param tau Numeric vector of length 1.  The exponent defining the weight of
#' the tail in the random walk performed by the `GSVA` (Hänzelmann et al.,
#' 2013) method.  The default value is 1 as described in the paper.
#'
#' @param maxDiff Logical vector of length 1 which offers two approaches to
#' calculate the enrichment statistic (ES) from the KS random walk statistic.
#' * `FALSE`: ES is calculated as the maximum distance of the random walk
#' from 0. This approach produces a distribution of enrichment scores that is
#' bimodal, but it can give large enrichment scores to gene sets whose genes
#' are not concordantly activated in one direction only.
#' * `TRUE` (the default): ES is calculated as the magnitude difference between
#' the largest positive and negative random walk deviations. This default value
#' gives larger enrichment scores to gene sets whose genes are concordantly
#' activated in one direction only.
#'
#' @param absRanking Logical vector of length 1 used only when `maxDiff=TRUE`.
#' When `absRanking=FALSE` (default) a modified Kuiper statistic is used to
#' calculate enrichment scores, taking the magnitude difference between the
#' largest positive and negative random walk deviations. When
#' `absRanking=TRUE` the original Kuiper statistic that sums the largest
#' positive and negative random walk deviations is used.
#'
#' @param sparse Logical vector of length 1 used only when the input expression
#' data in `exprData` is stored in a sparse matrix (e.g., a `dgCMatrix` or a
#' `SingleCellExperiment` object storing the expression data in a `dgCMatrix`).
#' In such a case, when `sparse=TRUE` (default), a sparse version of the GSVA
#' algorithm will be applied. Otherwise, when `sparse=FALSE`, the classical
#' (dense) version of the GSVA algorithm will be used.
#'
#' @return In the case of the `gsvaScores()` method, a gene-set by sample matrix
#' of GSVA enrichment scores stored in a ocntainer object of the same type as
#' the input expression data container in the `param` argument.
#'
#' @aliases gsvaScores,gsvaParam,GsvaExprData-method
#' @name gsvaScores
#' @rdname gsvaRanks
#'
#' @importFrom cli cli_alert_info cli_abort cli_alert_success
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @exportMethod gsvaScores
setMethod("gsvaScores", signature(param="gsvaParam", ranks="GsvaExprData"),
          function(param, ranks, geneSets=NA, minSize=NA, maxSize=NA,
                   tau=NA, maxDiff=NA, absRanking=NA, sparse=NA,
                   verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              .check_geneSets_minSize_maxSize_tau(geneSets, minSize, maxSize, tau)

              .check_maxDiff_absRanking_sparse(maxDiff, absRanking, sparse)

              tau <- ifelse(is.na(tau), get_tau(param), tau)
              maxDiff <- ifelse(is.na(maxDiff), get_maxDiff(param), maxDiff)
              absRanking <- ifelse(is.na(absRanking), get_absRanking(param),
                                   absRanking)
              sparse <- ifelse(is.na(sparse), get_sparse(param), sparse)

              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              filteredDataMatrix <- .filterGenes(dataMatrix,
                                                 removeConstant=TRUE,
                                                 removeNzConstant=TRUE)

              if (!identical(rownames(filteredDataMatrix),
                             rownames(unwrapData(ranks)))) {
                  msg <- paste("Rownames in ranks don't match those from the",
                               "input expression data in 'param'")
                  cli_abort(c("x"=msg))
              }

              filteredMappedGeneSets <- .filterAndMapGeneSets(param, geneSets,
                                                              minSize, maxSize,
                                                              filteredDataMatrix,
                                                              verbose)

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if (verbose)
                  cli_alert_info(sprintf("Calculating GSVA scores"))

              gsva_es <- .compute_gsva_scores(R=unwrapData(ranks),
                                              geneSetsIdx=filteredMappedGeneSets,
                                              tau=tau, maxDiff=maxDiff,
                                              absRanking=absRanking,
                                              sparse=sparse, verbose=verbose,
                                              BPPARAM=BPPARAM)

              rownames(gsva_es) <- names(filteredMappedGeneSets)
              colnames(gsva_es) <- colnames(filteredDataMatrix)

              gs <- .geneSetsIndices2Names(indices=filteredMappedGeneSets,
                                           names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), gsva_es, gs)

              if (verbose)
                  cli_alert_success("Calculations finished")

              return(rval)
          })


## END exported methods (to be moved to 'gsvaNewAPI.R')


#' @importFrom cli cli_progress_update
#' @importFrom parallel splitIndices
.row_iter <- function(X, idpb, n_chunks) {
    idx <- splitIndices(nrow(X), min(nrow(X), n_chunks))
    i <- 0L
    function() {
        if (i == length(idx))
            return(NULL)
        i <<- i + 1L
        if (!is.null(idpb))
            cli_progress_update(id=idpb, set=i)
        X[idx[[i]], , drop=FALSE]
    }
}

#' @importFrom cli cli_progress_update
#' @importFrom parallel splitIndices
.col_iter <- function(X, idpb, n_chunks) {
    idx <- splitIndices(ncol(X), min(ncol(X), n_chunks))
    i <- 0L
    function() {
        if (i == length(idx))
            return(NULL)
        i <<- i + 1L
        if (!is.null(idpb))
            cli_progress_update(id=idpb, set=i)
        X[, idx[[i]], drop=FALSE]
    }
}

#' @importFrom IRanges IntegerList match
#' @importFrom BiocParallel bpnworkers
#' @importFrom cli cli_alert_info cli_progress_bar
#' @importFrom cli cli_progress_done
#' @importFrom sparseMatrixStats colRanks
.compute_gsva_ranks <- function(expr, kcdf, kcdf.min.ssize,
                                sparse=FALSE, verbose=TRUE,
                                BPPARAM=SerialParam(progressbar=verbose)) {

    kcdfparam <- .parse_kcdf_param(expr, kcdf, kcdf.min.ssize, sparse, verbose)
    kernel <- kcdfparam$kernel
    Gaussk <- kcdfparam$Gaussk
    Z <- idpb <- NULL

    ## open parallelism only if ECDFs have to be estimated for
    ## more than 100 genes on more than 100 samples
    if (bpnworkers(BPPARAM) > 1 && nrow(expr) > 100 && ncol(expr) > 100) {
        n_chunks <- 10 ## 10 chunks of (nrow(expr) / 10) rows
        if (verbose) {
            msg <- sprintf("Estimating row ECDFs with %d cores",
                           as.integer(bpnworkers(BPPARAM)))
            idpb <- cli_progress_bar(msg, total=n_chunks)
        }
        Z <- bpiterate(.row_iter(expr, idpb, n_chunks),
                       compute.gene.cdf,
                       sample.idxs=seq.int(ncol(expr)),
                       Gaussk=Gaussk, kernel=kernel,
                       sparse=sparse, verbose=FALSE,
                       REDUCE=rbind, reduce.in.order=TRUE,
                       BPPARAM=BPPARAM)
        if (verbose)
            cli_progress_done(idpb)
    } else
        Z <- compute.gene.cdf(expr, seq.int(ncol(expr)), Gaussk, kernel,
                              sparse, verbose)

    R <- NULL
    ## here 'ties.method="last"' allows one to obtain the result
    ## from 'order()' based on ranks
    if (is(Z, "dgCMatrix")) { ## assumes expression values are positive
        if (verbose)
            cli_alert_info("Calculating GSVA column ranks")
        R <- .sparseColumnApplyAndReplace(Z, rank, ties.method="last")
    } else {
        ## open parallelism only if ranks have to be calculated for
        ## more than 10000 genes on more than 1000 samples
        if (bpnworkers(BPPARAM) > 1 && nrow(Z) > 10000 && ncol(Z) > 1000) {
            n_chunks <- 100 ## 100 chunks of (ncol(expr) / 100) columns
            if (verbose)
                idpb <- cli_progress_bar("Calculating GSVA column ranks")
            R <- bpiterate(.col_iter(Z, idpb, n_chunks),
                           colRanks, ties.method="last", preserveShape=TRUE,
                           REDUCE=cbind, reduce.in.order=TRUE,
                           BPPARAM=BPPARAM)
            if (verbose)
                cli_progress_done(idpb) 
        } else {
            if (verbose)
                cli_alert_info("Calculating GSVA column ranks")
            R <- colRanks(Z, ties.method="last", preserveShape=TRUE)
        }
    }

    return(R)
}

.gsvaRndWalk2 <- function(gSetIdx, decOrderStat, symRnkStat, tau) {
    n <- length(decOrderStat)
    k <- length(gSetIdx)
    gSetRnk <- decOrderStat[gSetIdx]

    stepCDFinGeneSet <- integer(n)
    if (tau == 1)
      stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]
    else {
      stepCDFinGeneSet <- numeric(n)
      stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]^tau
    }
    stepCDFinGeneSet <- cumsum(stepCDFinGeneSet)
    stepCDFinGeneSet <- stepCDFinGeneSet / stepCDFinGeneSet[n]

    stepCDFoutGeneSet <- rep(1L, n)
    stepCDFoutGeneSet[gSetRnk] <- 0L
    stepCDFoutGeneSet <- cumsum(stepCDFoutGeneSet)
    stepCDFoutGeneSet <- stepCDFoutGeneSet / stepCDFoutGeneSet[n]

    walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

    walkStat
}

## convert ranks into decreasing order statistics and symmetric rank statistics
.ranks2stats <- function(r, sparse) {
    mask <- r == 0
    p <- length(r)
    r_dense <- r

    if (any(mask)) {                  ## sparse ranks into dense ranks
        nzs <- sum(mask)
        r_dense[!mask] <- r_dense[!mask] + nzs ## shift ranks of nonzero values
        r_dense[mask] <- seq.int(nzs)          ## zeros get increasing ranks
    }

    dos <- p - r_dense + 1            ## dense ranks into decreasing order stats
    srs <- numeric(p)

    if (any(mask) && sparse) {
        r[!mask] <- r[!mask] + 1      ## shift ranks of nonzero values by one
        r[mask] <- 1                  ## all zeros get the same first rank
        srs <- abs(max(r)/2 - r)
    } else
        srs <- abs(p/2 - r_dense)

    list(dos=dos, srs=srs)
}

#' @importFrom BiocParallel bpnworkers
.compute_gsva_scores <- function(R, geneSetsIdx, tau=1, maxDiff=TRUE,
                                 absRanking=FALSE, sparse=TRUE, verbose=TRUE,
                                 BPPARAM=SerialParam(progressbar=verbose)) {
    n <- ncol(R)
    es <- NULL
    if (!is(R, "dgCMatrix"))
        sparse <- FALSE

    if (n > 10 && bpnworkers(BPPARAM) > 1) {
        if (verbose) {
            msg <- sprintf("Calculating GSVA scores with %d cores",
                           as.integer(bpnworkers(BPPARAM)))
            idpb <- cli_progress_bar(msg, total=n)
        }
        es <- bplapply(as.list(1:n), function(j, R) {
            rnkstats <- .ranks2stats(R[, j], sparse)
            md <- lapply(geneSetsIdx, function(gSetIdx, decOrdStat, symRnkStat) {
              walkStat <- .gsvaRndWalk2(gSetIdx, decOrdStat, symRnkStat, tau)
              maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
              maxDev
            }, rnkstats$dos, rnkstats$srs)
            if (maxDiff && absRanking)
                md[, 2] <- -1 * md[, 2]
            sco <- rowSums(md)
            if (!maxDiff)
                sco <- md[cbind(1:length(sco), ifelse(sco > 0, 1, 2))]
            sco
        }, R=R, BPPARAM=BPPARAM)
        es
    } else {
        idpb <- NULL
        if (verbose)
          idpb <- cli_progress_bar("Calculating GSVA scores", total=n)
        es <- lapply(as.list(1:n), function(j, R, idpb) {
            rnkstats <- .ranks2stats(R[, j], sparse)
            md <- lapply(geneSetsIdx, function(gSetIdx, decOrdStat, symRnkStat) {
              walkStat <- .gsvaRndWalk2(gSetIdx, decOrdStat, symRnkStat, tau)
              maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
              maxDev
            }, rnkstats$dos, rnkstats$srs)
            md <- do.call("rbind", md)
            if (maxDiff && absRanking)
                md[, 2] <- -1 * md[, 2]
            sco <- rowSums(md)
            if (!maxDiff)
                sco <- md[cbind(1:length(sco), ifelse(sco > 0, 1, 2))]
            if (verbose)
                cli_progress_update(id=idpb)
            sco
        }, R=R, idpb=idpb)
        if (verbose)
            cli_progress_done(id=idpb)
    }
    es <- do.call("cbind", es)

    return(es)
}

##
## functions interfacing C code
##

.ecdfvals_sparse_to_sparse <- function(X, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_sparse_to_sparse_R", X, Xrsp, verbose)
}

.ecdfvals_sparse_to_dense <- function(X, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_sparse_to_dense_R", X, Xrsp, verbose)
}

.ecdfvals_dense_to_dense <- function(X, verbose) {
  stopifnot(is.matrix(X)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_dense_to_dense_R", X, verbose)
}

.kcdfvals_sparse_to_sparse <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_sparse_to_sparse_R", X, Xrsp, Gaussk, verbose)
}

.kcdfvals_sparse_to_dense <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_sparse_to_dense_R", X, Xrsp, Gaussk, verbose)
}

.order_rankstat <- function(x) {
  stopifnot(is.numeric(x)) ## QC
  .Call("order_rankstat_R", x)
}

.gsva_rnd_walk <- function(gsetIdx, geneRanking, rankStat) {
  stopifnot(is.integer(gsetIdx)) ## QC
  stopifnot(is.integer(geneRanking)) ## QC
  ## stopifnot(is.integer(rankStat)) ## QC
  stopifnot(is.numeric(rankStat)) ## QC
  .Call("gsva_rnd_walk_R", gsetIdx, geneRanking, rankStat)
}

.gsva_score_genesets <- function(geneSetsRankIdx, geneRanking, rankStat,
                                 maxDiff, absRnk, tau) {
  stopifnot(is.list(geneSetsRankIdx)) ## QC
  stopifnot(length(geneSetsRankIdx) > 0) ## QC
  stopifnot(is.integer(geneSetsRankIdx[[1]])) ## QC
  stopifnot(is.integer(geneRanking)) ## QC
  ## stopifnot(is.integer(rankStat)) ## QC
  stopifnot(is.numeric(rankStat)) ## QC
  stopifnot(is.logical(maxDiff)) ## QC
  stopifnot(is.logical(absRnk)) ## QC
  stopifnot(is.numeric(tau)) ## QC
  .Call("gsva_score_genesets_R", geneSetsRankIdx, geneRanking, rankStat,
        maxDiff, absRnk, tau)
}

.order_rankstat_sparse_to_dense <- function(X, j) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  .Call("order_rankstat_sparse_to_dense_R", X, j)
}

.order_rankstat_sparse_to_sparse <- function(X, j) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  .Call("order_rankstat_sparse_to_sparse_R", X, j)
}



##
## we may not need anymore these two functions below
##

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
