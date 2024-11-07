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
                             sparse=FALSE, any_na=FALSE,
                             na_use=c("everything", "all.obs", "na.rm"),
                             verbose=TRUE) {
    na_use <- match.arg(na_use)
    n.test.samples <- ncol(expr)
    n.genes <- nrow(expr)
    n.density.samples <- length(sample.idxs)

    if (any_na && na_use == "all.obs") {
        msg <- paste("missing values present in the input expression data and",
                     "'use=\"all.obs\".")
        cli_abort(c("x"=msg))
    }
    
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
                      any_na,
                      as.integer(factor(na_use,
                                        levels=c("everything", "all.obs",
                                                 "na.rm"))),
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
        } else if (is.matrix(expr)) {
            if (any_na)
                gene.cdf <- .ecdfvals_dense_to_dense_nas(expr[, sample.idxs, drop=FALSE],
                                                         verbose)
            else
                gene.cdf <- .ecdfvals_dense_to_dense(expr[, sample.idxs, drop=FALSE],
                                                 verbose)
        } else
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

## here gSetIdx contains the positions in the decreasing gene ranking
## and rankStat contains the rank statistic value in the original
## gene order of the data
.gsvaRndWalk_rankingpos_notau <- function(gSetIdx, geneRanking, rankStat) {
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
                               kcdf.min.ssize, abs.ranking,
                               mx.diff=TRUE, tau=1, sparse=FALSE,
                               verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose)) {

    kcdfparam <- .parse_kcdf_param(expr, kcdf, kcdf.min.ssize, sparse, verbose)
    kernel <- kcdfparam$kernel
    Gaussk <- kcdfparam$Gaussk

    ## open parallelism only if ECDFs have to be estimated for
    ## more than 100 genes on more than 100 samples
    if (bpnworkers(BPPARAM) > 1 && length(sample.idxs) > 100 &&
        nrow(expr) > 100) {
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
                           as.integer(bpnworkers(BPPARAM)))
            idpb <- cli_progress_bar(msg, total=100)
        }
        gene.density <- bpiterate(iter(expr, idpb, 100),
                                  compute.gene.cdf,
                                  sample.idxs=sample.idxs,
                                  Gaussk=Gaussk, kernel=kernel,
                                  sparse=sparse, any_na=FALSE,
                                  na_use="everything", verbose=FALSE,
                                  REDUCE=rbind, reduce.in.order=TRUE,
                                  BPPARAM=BPPARAM)
        if (verbose)
            cli_progress_done(idpb)
    } else
        gene.density <- compute.gene.cdf(expr, sample.idxs, Gaussk, kernel,
                                         sparse, any_na=FALSE,
                                         na_use="everything", verbose)
    
    gset.idx.list <- IntegerList(gset.idx.list)
    n <- ncol(expr)
    es <- NULL
    if (n > 10 && bpnworkers(BPPARAM) > 1) {
        if (verbose) {
            msg <- sprintf("Calculating GSVA scores with %d cores",
                           as.integer(bpnworkers(BPPARAM)))
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
#' @param param A [`gsvaParam-class`] object built using the constructor
#' function [`gsvaParam`].
#'
#' @param verbose Gives information about each calculation step. Default: `TRUE`.
#'
#' @param BPPARAM An object of class [`BiocParallelParam`] specifying parameters
#'   related to the parallel execution of some of the tasks and calculations
#'   within this function.
#'
#' @return In the case of the `gsvaRanks()` method, an object of class
#' [`gsvaRanksParam-class`].
#'
#' @seealso [`gsvaParam-class`], [`gsvaRanksParam-class`], [`gsva`]
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
#' ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
#' y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2
#'
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate GSVA ranks
#' gsvarankspar <- gsvaRanks(gsvapar)
#' gsvarankspar
#' ## calculate GSVA scores
#' gsva_es <- gsvaScores(gsvarankspar)
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
#' geneSets(gsvarankspar) <- geneSets2
#' gsvaScores(gsvarankspar)
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom BiocParallel bpnworkers
#' @exportMethod gsvaRanks
setMethod("gsvaRanks", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose && gsva_global$show_start_and_end_messages) {
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))
              }

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

              kcdfminssize <-get_kcdfNoneMinSampleSize(param)
              gsvarnks <- .compute_gsva_ranks(expr=filteredDataMatrix,
                                              kcdf=get_kcdf(param),
                                              kcdf.min.ssize=kcdfminssize,
                                              sparse=get_sparse(param),
                                              any_na=anyNA(param),
                                              na_use=get_NAuse(param),
                                              verbose=verbose,
                                              BPPARAM=BPPARAM)

              rownames(gsvarnks) <- rownames(filteredDataMatrix)
              colnames(gsvarnks) <- colnames(filteredDataMatrix)

              rnkcontainer <- wrapData(get_exprData(param), gsvarnks)
              rval <- new("gsvaRanksParam",
                          exprData=rnkcontainer, geneSets=get_geneSets(param),
                          assay="gsvaranks", annotation=get_annotation(param),
                          minSize=get_minSize(param), maxSize=get_maxSize(param),
                          kcdf=get_kcdf(param),
                          kcdfNoneMinSampleSize=get_kcdfNoneMinSampleSize(param),
                          tau=get_tau(param), maxDiff=get_maxDiff(param),
                          absRanking=get_absRanking(param),
                          sparse=get_sparse(param), checkNA=get_checkNA(param),
                          didCheckNA=get_didCheckNA(param), anyNA=anyNA(param),
                          use=get_NAuse(param))

              if (verbose && gsva_global$show_start_and_end_messages)
                  cli_alert_success("Calculations finished")

              return(rval)
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

.check_maxDiff_absRanking <- function(maxDiff, absRanking) {
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
}


#' @param param A parameter object of the [`gsvaRanksParam-class`] class.
#'
#' @return In the case of the `gsvaScores()` method, a gene-set by sample matrix
#' of GSVA enrichment scores stored in a container object of the same type as
#' the input ranks data container. If
#' the input was a base matrix or a [`dgCMatrix-class`] object, then the output will
#' be a base matrix object with the gene sets employed in the calculations
#' stored in an attribute called `geneSets`. If the input was an
#' [`ExpressionSet`] object, then the output will be also an [`ExpressionSet`]
#' object with the gene sets employed in the calculations stored in an
#' attributed called `geneSets`. If the input was an object of one of the
#' classes described in [`GsvaExprData`], such as a [`SingleCellExperiment`],
#' then the output will be of the same class, where enrichment scores will be
#' stored in an assay called `es` and the gene sets employed in the
#' calculations will be stored in the `rowData` slot of the object under the
#' column name `gs`.
#'
#' @aliases gsvaScores,gsvaRanksParam-method
#' @name gsvaScores
#' @rdname gsvaRanks
#'
#' @importFrom cli cli_alert_info cli_abort cli_alert_success
#' @importFrom BiocParallel bpnworkers
#' @exportMethod gsvaScores
setMethod("gsvaScores", signature(param="gsvaRanksParam"),
          function(param, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose && gsva_global$show_start_and_end_messages) {
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))
              }

              ## assuming rows in the rank data have been already filtered
              exprData <- get_exprData(param)
              filteredDataMatrix <- unwrapData(exprData, get_assay(param))
              filteredMappedGeneSets <- .filterAndMapGeneSets(param=param,
                                                              filteredDataMatrix=filteredDataMatrix,
                                                              verbose=verbose)

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if (verbose)
                  cli_alert_info(sprintf("Calculating GSVA scores"))

              gsva_es <- .compute_gsva_scores(R=filteredDataMatrix,
                                              geneSetsIdx=filteredMappedGeneSets,
                                              tau=get_tau(param),
                                              maxDiff=get_maxDiff(param),
                                              absRanking=get_absRanking(param),
                                              sparse=get_sparse(param),
                                              any_na=anyNA(param),
                                              na_use=get_NAuse(param),
                                              minSize=get_minSize(param),
                                              verbose=verbose, BPPARAM=BPPARAM)

              rownames(gsva_es) <- names(filteredMappedGeneSets)
              colnames(gsva_es) <- colnames(filteredDataMatrix)

              gs <- .geneSetsIndices2Names(indices=filteredMappedGeneSets,
                                           names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), gsva_es, gs)

              if (verbose && gsva_global$show_start_and_end_messages)
                  cli_alert_success("Calculations finished")

              return(rval)
          })

#' @title GSVA enrichment data and visualization
#'
#' @description Extract and plot enrichment data from GSVA scores.
#'
#' @param param A [`gsvaRanksParam-class`] object obtained with the method
#' [`gsvaRanks`].
#'
#' @param column The column for which we want to retrieve the enrichment data.
#' This parameter is only available in the `gsvaEnrichment()` method.
#'
#' @param geneSet Either a positive integer number between 1 and the number of
#' available gene sets in `param`, or a character string with the name of
#' one of the gene sets available in `param`.
#'
#' @param plot A character string indicating whether an enrichment plot should
#' be produced using either base R graphics (`plot="base"`) or the ggplot2 package
#' (`plot="ggplot"`), or not (`plot="no"`). In the latter case, the enrichment
#' data will be returned. By default `plot="auto"`, which implies that if this
#' method is called from an interactive session, a plot using base R graphics
#' will be produced and, otherwise, the enrichment data is returned.
#'
#' @param ... Further arguments passed to the `plot()` function when the
#' previous parameter `plot="base"`.
#'
#' @return When `plot="no"`, this method returns the enrichment data. When
#' `plot="ggplot"`, this method returns a `ggplot` object. When `plot="base"`
#' no value is returned.
#'
#' @aliases gsvaEnrichment,gsvaRanksParam-method
#' @name gsvaEnrichment
#' @rdname gsvaEnrichment
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
#' ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
#' y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2
#'
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate GSVA ranks
#' gsvarankspar <- gsvaRanks(gsvapar)
#' gsvarankspar
#'
#' ## by default the enrichment data for the first column and the first
#' ## gene set are retrieved
#' gsvaEnrichment(gsvarankspar)
#'
#' @importFrom cli cli_alert_info cli_abort cli_alert_danger
#' @exportMethod gsvaScores
setMethod("gsvaEnrichment", signature(param="gsvaRanksParam"),
          function(param, column=1, geneSet=1,
                   plot=c("auto", "base", "ggplot", "no"), ...)
          {
              plot <- match.arg(plot)

              geneSets <- get_geneSets(param)
              if (length(geneSet) > 1) {
                  msg <- paste("Please provide only the name or position of a",
                               "single gene set.")
                  cli_abort(c("x"=msg))
              }
              if (is.character(geneSet)) {
                  if (!geneSet %in% names(geneSets)) {
                      msg <- paste("Gene set %s is missing from the input",
                                   "parameter object")
                      cli_abort(c("x"=sprintf(msg, geneSet)))
                  }
              } else if (is.integer(geneSet) || is.numeric(geneSet)) {
                  if (geneSet < 1 || geneSet > length(geneSets)) {
                       msg <- paste("When 'geneSet' is numeric, it should be a",
                                    "number between 1 and the number of gene",
                                    "sets (%d).")
                       cli_abort(c("x"=sprintf(msg, length(geneSets))))
                  }
              } else {
                  msg <- paste("'geneSet' should be either numeric or",
                               "character.")
                  cli_abort(c("x"=msg))
              }

              tau <- get_tau(param)
              maxDiff <- get_maxDiff(param)
              absRanking <- get_absRanking(param)
              sparse <- get_sparse(param)
              any_na <- anyNA(param)
              na_use <- get_NAuse(param)
              minsize <- get_minSize(param)

              exprData <- get_exprData(param)
              filteredDataMatrix <- unwrapData(exprData, get_assay(param))

              ## no need for verbosity when mapping a single gene set
              filteredMappedGeneSets <- .filterAndMapGeneSets(param, wgset=geneSet,
                                                              filteredDataMatrix,
                                                              verbose=FALSE)

              geneSetIdx <- filteredMappedGeneSets[[1]]
              edata <- .gsva_enrichment_data(R=filteredDataMatrix,
                                             column=column,
                                             geneSetIdx=geneSetIdx,
                                             maxDiff=maxDiff,
                                             absRanking=absRanking,
                                             tau=tau,
                                             sparse=sparse,
                                             any_na=any_na,
                                             na_use=na_use,
                                             minSize=minsize)

              if (plot == "no" || (plot == "auto" && !interactive()))
                  return(edata)

              if (plot == "auto" || plot == "base")
                  .plot_enrichment_base(edata, ...) 
              else { ## plot == "ggplot"
                  instpkgs <- installed.packages(noCache=TRUE)[, "Package"]
                  if (!"ggplot2" %in% instpkgs)
                      cli_alert_danger("Please install the ggplot2 package")
                  else
                      .plot_enrichment_ggplot(edata)
              }
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
                                sparse, any_na, na_use, verbose,
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
                       sparse=sparse, any_na=any_na,
                       na_use=na_use, verbose=FALSE,
                       REDUCE=rbind, reduce.in.order=TRUE,
                       BPPARAM=BPPARAM)
        if (verbose)
            cli_progress_done(idpb)
    } else
        Z <- compute.gene.cdf(expr, seq.int(ncol(expr)), Gaussk, kernel,
                              sparse, any_na=any_na, na_use=na_use, verbose)

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

## here gSetIdx, decOrderStat and symRnkStat contain the positions with respect
## to the original order of genes in the data
.gsvaRndWalk <- function(gSetIdx, decOrderStat, symRnkStat, tau) {
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
    stepCDFoutGeneSet <- rep(1L, n)
    stepCDFoutGeneSet[gSetRnk] <- 0L
    stepCDFoutGeneSet <- cumsum(stepCDFoutGeneSet)

    walkStat <- rep(NA_real_, n)
    if (stepCDFinGeneSet[n] > 0 && stepCDFoutGeneSet[n] > 0) {
        stepCDFinGeneSet <- stepCDFinGeneSet / stepCDFinGeneSet[n]
        stepCDFoutGeneSet <- stepCDFoutGeneSet / stepCDFoutGeneSet[n]

        walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
    }

    walkStat
}

.gsva_score_genesets_Rimp <- function(geneSetsIdx, decOrdStat, symRnkStat,
                                      maxDiff, absRanking, tau, any_na, na_use,
                                      minSize) {
   md <- lapply(geneSetsIdx, function(gSetIdx, decOrdStat, symRnkStat) {
             maxDev <- c(NA_real_, NA_real_)
             if (any_na) {
                 walkStat <- .gsvaRndWalk_nas(gSetIdx, decOrdStat, symRnkStat,
                                              tau, na_use, minSize)
                 if (any(!is.na(walkStat))) {
                     if (na_use == "na.rm")
                         maxDev <- c(max(c(0, max(walkStat, na.rm=TRUE))),
                                     min(c(0, min(walkStat, na.rm=TRUE))))
                     else
                         maxDev <- c(max(c(0, max(walkStat))),
                                     min(c(0, min(walkStat))))
                 }
             } else {
                 walkStat <- .gsvaRndWalk(gSetIdx, decOrdStat, symRnkStat, tau)
                 maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
             }
             maxDev
         }, decOrdStat, symRnkStat)
   md <- do.call("rbind", md)
   if (maxDiff && absRanking)
       md[, 2] <- -1 * md[, 2]
   sco <- rowSums(md)
   if (!maxDiff) {
       if (any_na) {
         mask <- is.na(sco)
         sco[!mask] <- md[cbind(1:sum(!mask), ifelse(sco[!mask] > 0, 1, 2))]
       } else
         sco <- md[cbind(1:length(sco), ifelse(sco > 0, 1, 2))]
   }
   sco
}

## here gSetIdx, decOrderStat and symRnkStat contain the positions with respect
## to the original order of genes in the data
.gsvaRndWalk_nas <- function(gSetIdx, decOrderStat, symRnkStat, tau, na_use,
                             minSize=1L, wna_env) {
    n <- length(decOrderStat)
    gSetRnk <- decOrderStat[gSetIdx]

    if (any(is.na(gSetRnk))) {
        if (na_use == "everything")
            return(rep(NA_real_, n))
        else if (na_use == "all.obs")
            cli_abort(c("x"="Input GSVA ranks have NA values."))
        else if (na_use == "na.rm") {
            gSetIdx <- gSetIdx[!is.na(gSetRnk)]
            gSetRnk <- gSetRnk[!is.na(gSetRnk)]
        }
    }
    k <- length(gSetIdx)

    walkStat <- rep(NA_real_, n)
    if (k >= minSize) {

        stepCDFinGeneSet <- integer(n)
        if (tau == 1)
          stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]
        else {
          stepCDFinGeneSet <- numeric(n)
          stepCDFinGeneSet[gSetRnk] <- symRnkStat[gSetIdx]^tau
        }

        stepCDFinGeneSet <- cumsum(stepCDFinGeneSet)
        stepCDFoutGeneSet <- rep(1L, n)
        stepCDFoutGeneSet[gSetRnk] <- 0L
        stepCDFoutGeneSet <- cumsum(stepCDFoutGeneSet)

        if (stepCDFinGeneSet[n] > 0 && stepCDFinGeneSet[n] > 0) {
            stepCDFinGeneSet <- stepCDFinGeneSet / stepCDFinGeneSet[n]
            stepCDFoutGeneSet <- stepCDFoutGeneSet / stepCDFoutGeneSet[n]

            walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
        }
    } else if (!get("w", envir=wna_env)) ## warn only once. it can only happen
        assign("w", TRUE, envir=wna_env) ## with na_use="na.rm"

    walkStat
}

## convert ranks into decreasing order statistics and symmetric rank statistics
.ranks2stats <- function(r, sparse) {
    mask <- r == 0
    p <- length(r)
    r_dense <- as.integer(r)          ## assume ranks are integer

    if (any(mask)) {                  ## sparse ranks into dense ranks
        nzs <- sum(mask)
        r_dense[!mask] <- r_dense[!mask] + nzs ## shift ranks of nonzero values
        r_dense[mask] <- seq.int(nzs)          ## zeros get increasing ranks
    }

    dos <- p - r_dense + 1L           ## dense ranks into decreasing order stats
    srs <- numeric(p)

    if (any(mask) && sparse) {
        r[!mask] <- r[!mask] + 1      ## shift ranks of nonzero values by one
        r[mask] <- 1                  ## all zeros get the same first rank
        srs <- abs(max(r)/2 - r)
    } else
        srs <- abs(p/2 - r_dense)

    list(dos=dos, srs=srs)
}

## convert ranks into decreasing order statistics and symmetric rank statistics
## skipping NA values
.ranks2stats_nas <- function(r, sparse) {
    na_mask <- is.na(r)

    if (all(na_mask))
        return(list(dos=rep(NA, length(r)), srs=rep(NA, length(r))))

    n_nas <- sum(na_mask)
    mask <- !na_mask & r == 0
    p <- length(r)
    r_dense <- as.integer(r)          ## assume ranks are integer

    if (any(mask)) {                  ## sparse ranks into dense ranks
        nzs <- sum(mask)
        r_dense[!mask] <- r_dense[!mask] + nzs ## shift ranks of nonzero values
        r_dense[mask] <- seq.int(nzs)          ## zeros get increasing ranks
    }

    dos <- p - n_nas - r_dense + 1L   ## dense ranks into decreasing order stats
    srs <- numeric(p)

    if (any(mask) && sparse) {
        r[!mask] <- r[!mask] + 1L     ## shift ranks of nonzero values by one
        r[mask] <- 1L                 ## all zeros get the same first rank
        srs <- abs(max(r, na.rm=TRUE)/2 - r)
    } else
        srs <- abs((p - n_nas)/2 - r_dense)

    list(dos=dos, srs=srs)
}


#' @importFrom cli cli_alert_warning
#' @importFrom BiocParallel bpnworkers
.compute_gsva_scores <- function(R, geneSetsIdx, tau, maxDiff, absRanking,
                                 sparse, any_na, na_use, minSize, verbose,
                                 BPPARAM=SerialParam(progressbar=verbose)) {
    n <- ncol(R)
    es <- NULL
    if (!is(R, "dgCMatrix"))
        sparse <- FALSE
    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    if (n > 10 && bpnworkers(BPPARAM) > 1) {
        if (verbose) {
            msg <- sprintf("Calculating GSVA scores with %d cores",
                           as.integer(bpnworkers(BPPARAM)))
            idpb <- cli_alert_info(msg)
        }
        es <- bplapply(as.list(1:n), function(j, R) {
            rnkstats <- NULL
            if (any_na)
                rnkstats <- .ranks2stats_nas(R[, j], sparse)
            else
                rnkstats <- .ranks2stats(R[, j], sparse)
            sco <- .gsva_score_genesets(geneSetsIdx, decOrdStat=rnkstats$dos,
                                        symRnkStat=rnkstats$srs, maxDiff,
                                        absRanking, tau, any_na, na_use, minSize,
                                        wna_env)
            sco
        }, R=R, BPPARAM=BPPARAM)
        es
    } else {
        idpb <- NULL
        if (verbose)
          idpb <- cli_progress_bar("Calculating GSVA scores", total=n)
        es <- lapply(as.list(1:n), function(j, R, idpb) {
            rnkstats <- NULL
            if (any_na)
                rnkstats <- .ranks2stats_nas(R[, j], sparse)
            else
                rnkstats <- .ranks2stats(R[, j], sparse)
            sco <- .gsva_score_genesets(geneSetsIdx, decOrdStat=rnkstats$dos,
                                        symRnkStat=rnkstats$srs, maxDiff,
                                        absRanking, tau, any_na, na_use, minSize,
                                        wna_env)
            if (verbose)
                cli_progress_update(id=idpb)
            sco
        }, R=R, idpb=idpb)
        if (verbose)
            cli_progress_done(id=idpb)
    }
    es <- do.call("cbind", es)

    if (any_na && na_use == "na.rm")
        if (get("w", envir=wna_env)) {
            msg <- sprintf(paste("NA enrichment scores in gene sets with less",
                                 "than %d genes after removing missing values"),
                           minSize)
            cli_alert_warning(msg)
        }

    return(es)
}

.gsva_enrichment_data <- function(R, column, geneSetIdx, maxDiff,
                                  absRanking, tau, sparse, any_na,
                                  na_use, minSize) {
    n <- ncol(R)
    es <- NULL
    if (!is(R, "dgCMatrix"))
        sparse <- FALSE
    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    if (any_na) {
        rnkstats <- .ranks2stats_nas(R[, column], sparse)
        walkStat <- .gsvaRndWalk_nas(geneSetIdx, rnkstats$dos, rnkstats$srs,
                                     tau, na_use, minSize, wna_env=wna_env)
        maxDev <- whMaxDev <- c(NA, NA)
        if (any(!is.na(walkStat))) {
          if (na_use == "na.rm")
              maxDev <- c(max(c(0, max(walkStat, na.rm=TRUE))),
                          min(c(0, min(walkStat, na.rm=TRUE))))
          else
              maxDev <- c(max(c(0, max(walkStat))),
                          min(c(0, min(walkStat))))
        }
        if (length(which.max(walkStat)) > 0)
            whMaxDev[1] <- which.max(walkStat)
        if (length(which.min(walkStat)) > 0)
            whMaxDev[2] <- which.min(walkStat)
    } else {
        rnkstats <- .ranks2stats(R[, column], sparse)
        walkStat <- .gsvaRndWalk(geneSetIdx, rnkstats$dos, rnkstats$srs, tau)
        maxDev <- c(max(c(0, max(walkStat))), min(c(0, min(walkStat))))
        whMaxDev <- c(which.max(walkStat), which.min(walkStat))
        whMaxDev[maxDev == 0] <- NA
    }
    
    if (any_na && na_use == "na.rm")
        if (get("w", envir=wna_env)) {
            msg <- sprintf(paste("Gene set has fewer than %d genes after",
                                 "removing missing values, no enrichment data",
                                 "available"),
                           minSize)
            cli_alert_warning(msg)
            return(list())
        }

    if (maxDiff && absRanking)
        maxDev[2] <- -1 * maxDev[2]
    sco <- sum(maxDev)
    if (!maxDiff) {
        if (any_na) {
            if (!is.na(sco)) {
                sco <- maxDev[1]
                if (abs(maxDev[2]) > maxDev[1])
                    sco <- maxDev[2]
            }
        }
    }

    edat <- data.frame(rank=seq.int(nrow(R)),
                       stat=walkStat)
    rownames(edat)[na.omit(rnkstats$dos)] <- rownames(R)[!is.na(rnkstats$dos)]

    geneSetIdx <- geneSetIdx[!is.na(rnkstats$dos[geneSetIdx])]
    gsetrnk <- rnkstats$dos[geneSetIdx]
    lepos <- leneg <- NA
    if (!is.na(whMaxDev[1]))
        lepos <- geneSetIdx[gsetrnk <= whMaxDev[1]]
    if (!is.na(whMaxDev[2])) {
        if (!is.na(whMaxDev[1]) && whMaxDev[2] < whMaxDev[1]) {
            mask <- gsetrnk >= whMaxDev[2] & gsetrnk <= whMaxDev[1]
            lepos <- leneg <- geneSetIdx[mask]
        } else
            leneg <- geneSetIdx[gsetrnk >= whMaxDev[2]]
    }
    if (all(!is.na(lepos)))
        lepos <- rownames(R)[lepos]
    if (all(!is.na(leneg)))
        leneg <- rownames(R)[leneg]

    res <- list(stats=edat,
                gsetrnk=gsetrnk,
                maxPos=maxDev[1],
                whichMaxPos=whMaxDev[1],
                maxNeg=maxDev[2],
                whichMaxNeg=whMaxDev[2],
                leadingEdgePos=lepos,
                leadingEdgeNeg=leneg,
                score=sco,
                tau=tau,
                maxDiff=maxDiff,
                absRanking=absRanking,
                sparse=sparse)

    return(res)
}

#' @importFrom graphics abline grid lines segments
.plot_enrichment_base <- function(edata, ...) {
    ylim <- range(edata$stats$stat)
    hgsetticks <- (ylim[2] - ylim[1]) * 0.1
    plot(edata$stats, type="l", lwd=2, las=1, panel.first=grid(),
         xlab="Gene Ranking", ylab="Random Walk Statistic", col="green", ...)
    abline(h=0, lwd=2, lty=2, col="grey")
    lines(edata$stats, lwd=2, col="green")
    segments(edata$gsetrnk, -hgsetticks/2, edata$gsetrnk, hgsetticks/2, lwd=2)
    if (!is.na(edata$whichMaxPos) &&
        (edata$maxDiff || edata$maxPos >= abs(edata$maxNeg)))
        segments(edata$whichMaxPos, 0, edata$whichMaxPos, edata$maxPos,
                 lwd=2, lty=2, col="darkred")
    if (!is.na(edata$whichMaxNeg) &&
        (edata$maxDiff || edata$maxPos < abs(edata$maxNeg)))
        segments(edata$whichMaxNeg, 0, edata$whichMaxNeg, edata$maxNeg,
                 lwd=2, lty=2, col="darkred")
}

#' @importFrom cli cli_abort
#' @importFrom utils globalVariables
.plot_enrichment_ggplot <- function(edata, ...) {
    if (!.isPackageLoaded("ggplot2")) {
        loaded <- suppressPackageStartupMessages(requireNamespace("ggplot2"))
        if (!loaded)
            cli_abort(c("x"="ggplot2 could not be loaded"))
    }

    ylim <- range(edata$stats$stat)
    hgsetticks <- (ylim[2] - ylim[1]) * 0.1
    gsetticks <- data.frame(gsetrnk=edata$gsetrnk)
    ## from https://stackoverflow.com/a/39877048
    fintticks <- function(x) unique(floor(pretty(seq(min(x),
                                    (max(x) + 1) * 1.1))))
    ggplot2::ggplot(data=edata$stats) +
        ggplot2::scale_x_continuous(breaks=fintticks) +
        ggplot2::geom_line(ggplot2::aes_string(x="rank", y="stat"), color="green") +
        ggplot2::geom_segment(data=gsetticks,
                     mapping=ggplot2::aes_string(x="gsetrnk", y=-hgsetticks/2,
                                 xend="gsetrnk", yend=hgsetticks/2),
                     linewidth=1) +
        ggplot2::geom_hline(yintercept=0, colour="grey", linetype="dashed") +
        { if (!is.na(edata$whichMaxPos) &&
              (edata$maxDiff || edata$maxPos >= abs(edata$maxNeg)))
              ggplot2::geom_segment(data=data.frame(whichMaxPos=edata$whichMaxPos,
                                           maxPos=edata$maxPos),
                           mapping=ggplot2::aes_string(x="whichMaxPos", y=0,
                                       xend="whichMaxPos", yend="maxPos"),
                           colour="darkred", linetype="dashed") } +
        { if (!is.na(edata$whichMaxPos) &&
              (edata$maxDiff || edata$maxPos < abs(edata$maxNeg)))
              ggplot2::geom_segment(data=data.frame(whichMaxNeg=edata$whichMaxNeg,
                                           maxNeg=edata$maxNeg),
                           mapping=ggplot2::aes_string(x="whichMaxNeg", y=0,
                                       xend="whichMaxNeg", yend="maxNeg"),
                           colour="darkred", linetype="dashed") } +
        ggplot2::theme(panel.background=ggplot2::element_blank(),
              panel.grid.major=ggplot2::element_line(colour="grey", linetype="dotted"),
              panel.grid.minor=ggplot2::element_line(colour=NA),
              axis.text=ggplot2::element_text(size=12),
              axis.title=ggplot2::element_text(size=14),
              panel.border=ggplot2::element_rect(colour="black", fill=NA)) +
        ggplot2::labs(x="Gene Ranking", y="Random Walk Statistic")
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

.ecdfvals_dense_to_dense_nas <- function(X, verbose) {
  stopifnot(is.matrix(X)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_dense_to_dense_nas_R", X, verbose)
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

.gsva_rnd_walk <- function(gsetIdx, decOrdStat, symRnkStat) {
  stopifnot(is.integer(gsetIdx)) ## QC
  stopifnot(is.integer(decOrdStat)) ## QC
  stopifnot(is.numeric(symRnkStat)) ## QC
  .Call("gsva_rnd_walk_R", gsetIdx, decOrdStat, symRnkStat)
}

#' @importFrom cli cli_abort
.gsva_score_genesets <- function(geneSetsIdx, decOrdStat, symRnkStat,
                                 maxDiff, absRanking, tau, any_na, na_use,
                                 minSize, wna_env) {
  minSize <- as.integer(minSize)
  stopifnot(is.list(geneSetsIdx)) ## QC
  stopifnot(length(geneSetsIdx) > 0) ## QC
  stopifnot(is.integer(geneSetsIdx[[1]])) ## QC
  stopifnot(is.integer(decOrdStat)) ## QC
  stopifnot(is.numeric(symRnkStat)) ## QC
  stopifnot(is.logical(maxDiff)) ## QC
  stopifnot(is.logical(absRanking)) ## QC
  stopifnot(is.numeric(tau)) ## QC
  stopifnot(is.logical(any_na)) ## QC
  stopifnot(is.character(na_use)) ## QC
  stopifnot(is.integer(minSize)) ## QC
  na_use <- as.integer(factor(na_use, levels=c("everything", "all.obs",
                                               "na.rm")))
  sco <- .Call("gsva_score_genesets_R", geneSetsIdx, decOrdStat, symRnkStat,
               maxDiff, absRanking, tau, any_na, na_use, minSize)
  if (any_na) {
    if (na_use == 2 && !is.null(attr(sco, "class")))
        cli_abort(c("x"="Input GSVA ranks have NA values."))

    if (na_use == 3 && !is.null(attr(sco, "class")))
        assign("w", TRUE, envir=wna_env)
  }

  sco
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
