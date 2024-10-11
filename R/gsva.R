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

#' @title GSVA ranks
#'
#' @description Calculates GSVA ranks. This functions allows one to
#' split the calculations of the GSVA method in two steps: (1) calculate
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
#' @return A matrix of GSVA rank values per column.
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
#' geneSets <- list(set1=paste("g", 1:3, sep=""),
#'                  set2=paste("g", 4:6, sep=""),
#'                  set3=paste("g", 7:10, sep=""))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate GSVA ranks
#' gsva_ranks <- gsvaRanks(gsvapar)
#' gsva_ranks
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

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       verbose)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if (verbose)
                  cli_alert_info(sprintf("Calculating GSVA ranks"))

              gsvaRanks <- .compute_gsva_ranks(expr=filteredDataMatrix,
                                               kcdf=get_kcdf(param),
                                               kcdf.min.ssize=get_kcdfNoneMinSampleSize(param),
                                               parallel.sz=if (inherits(BPPARAM, "SerialParam"))
                                                               1L else bpnworkers(BPPARAM),
                                               sparse=get_sparse(param),
                                               verbose=verbose,
                                               BPPARAM=BPPARAM)

              rownames(gsvaRanks) <- rownames(filteredDataMatrix)
              colnames(gsvaRanks) <- colnames(filteredDataMatrix)

              return(gsvaRanks)
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
.compute_gsva_ranks <- function(expr, kcdf, kcdf.min.ssize, parallel.sz=1L,
                                sparse=FALSE, verbose=TRUE,
                                BPPARAM=SerialParam(progressbar=verbose)) {

    kcdfparam <- .parse_kcdf_param(expr, kcdf, kcdf.min.ssize, sparse, verbose)
    kernel <- kcdfparam$kernel
    Gaussk <- kcdfparam$Gaussk
    Z <- idpb <- NULL

    ## open parallelism only if ECDFs have to be estimated for
    ## more than 100 genes on more than 100 samples
    if (parallel.sz > 1 && nrow(expr) > 100 && ncol(expr) > 100 &&
        bpnworkers(BPPARAM) > 1) {
        n_chunks <- 10 ## 10 chunks of (nrow(expr) / 10) rows
        if (verbose) {
            msg <- sprintf("Estimating row ECDFs with %d cores",
                           as.integer(parallel.sz))
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
        if (parallel.sz > 1 && nrow(Z) > 10000 && ncol(Z) > 1000) {
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

## assumes all ranks are different
.gsvaRndWalk2 <- function(gSetIdx, R_j) {
    stopifnot(all(!duplicated(R_j))) ## QC
    n <- length(R_j)
    k <- length(gSetIdx)
    rnkStat <- abs(n/2-R_j)             ## symmetric rank statistic
    gSetRnk <- n-R_j[gSetIdx]+1         ## decreasing order statistic

    stepCDFinGeneSet <- integer(n)
    stepCDFinGeneSet[gSetRnk] <- rnkStat[gSetIdx]
    stepCDFinGeneSet <- cumsum(stepCDFinGeneSet)
    stepCDFinGeneSet <- stepCDFinGeneSet / stepCDFinGeneSet[n]

    stepCDFoutGeneSet <- rep(1L, n)
    stepCDFoutGeneSet[gSetRnk] <- 0L
    stepCDFoutGeneSet <- cumsum(stepCDFoutGeneSet)
    stepCDFoutGeneSet <- stepCDFoutGeneSet / stepCDFoutGeneSet[n]

    walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

    walkStat
}

.ranks_sparse_to_dense <- function(r) {
    mask <- r == 0
    if (any(mask)) {
        nzeros <- sum(mask)
        r[!mask] <- r[!mask] + nzeros ## shift ranks of nonzero values
        r[mask] <- 1:nzeros           ## zeros get increasing ranks
    }
    r
}

.compute_gsva_scores <- function(R, gset.idx.list, abs.ranking, parallel.sz=1L,
                                 mx.diff=TRUE, tau=1, sparse=FALSE, verbose=TRUE,
                                 BPPARAM=SerialParam(progressbar=verbose)) {
    stopifnot(tau==1) ## MOMENTARILY!!!
    n <- ncol(R)
    es <- NULL
    if (n > 10 && bpnworkers(BPPARAM) > 1) {
        if (verbose) {
            msg <- sprintf("Calculating GSVA scores with %d cores",
                           as.integer(parallel.sz))
            idpb <- cli_progress_bar(msg, total=n)
        }
        es <- bplapply(as.list(1:n), function(j, R) {
            md <- lapply(gset.idx.list, function(gSetIdx, R, j) {
              R_j <- .ranks_sparse_to_dense(R[, j])
              walkStat <- .gsvaRndWalk2(gSetIdx, R_j)
              maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
              maxDev
            }, R, j)
            if (mx.diff && abs.ranking)
                md[, 2] <- -1 * md[, 2]
            sco <- rowSums(md)
            if (!mx.diff)
                sco <- md[cbind(1:length(sco), ifelse(sco > 0, 1, 2))]
            sco
        }, R=R, BPPARAM=BPPARAM)
        es
    } else {
        idpb <- NULL
        if (verbose)
          idpb <- cli_progress_bar("Calculating GSVA scores", total=n)
        es <- lapply(as.list(1:n), function(j, R, idpb) {
            md <- lapply(gset.idx.list, function(gSetIdx, R, j) {
              R_j <- .ranks_sparse_to_dense(R[, j])
              walkStat <- .gsvaRndWalk2(gSetIdx, R_j)
              maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
              maxDev
            }, R, j)
            md <- do.call("rbind", md)
            if (mx.diff && abs.ranking)
                md[, 2] <- -1 * md[, 2]
            sco <- rowSums(md)
            if (!mx.diff)
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
