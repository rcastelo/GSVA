#' @importFrom S4Arrays is_sparse
#' @importFrom DelayedArray seed
#' @importFrom cli cli_abort
compute.gene.cdf <- function(expr, Gaussk=TRUE, kernel=TRUE,
                             sparse=FALSE, any_na=FALSE,
                             na_use=c("everything", "all.obs", "na.rm"),
                             grid=NULL, verbose=TRUE, BPPARAM=NULL) {
    na_use <- match.arg(na_use)
    n.test.samples <- ncol(expr)
    n.genes <- nrow(expr)
    n.density.samples <- ncol(expr)

    if (any_na && na_use == "all.obs") {
        msg <- paste("missing values present in the input expression data and",
                     "'use=\"all.obs\".")
        cli_abort(c("x"=msg))
    }
    
    gene.cdf <- NA
    if (kernel) { ## kernel ECDF estimation
        if (is(expr, "dgCMatrix")) {
            if (sparse)
                gene.cdf <- .kcdfvals_sparse_to_sparse(expr, Gaussk, verbose)
            else
                gene.cdf <- .kcdfvals_sparse_to_dense(expr, Gaussk, verbose)
        } else if (is(expr, "SVT_SparseArray")) {
            if (sparse)
                gene.cdf <- .kcdfvals_svt_to_svt(expr, Gaussk, verbose)
            else
                gene.cdf <- .kcdfvals_svt_to_dense(expr, Gaussk, verbose)
        } else if (is(expr, "DelayedMatrix")) {
            if (sparse)
                gene.cdf <- .kcdfvals_sparseh5_to_sparseh5(expr, Gaussk=Gaussk,
                                                           grid=grid,
                                                           verbose)
            else {
                if (is_sparse(expr)) ## input HDF5 may be sparse or not
                    gene.cdf <- .kcdfvals_sparseh5_to_denseh5(expr, Gaussk=Gaussk,
                                                              grid=grid,
                                                              verbose)
                else
                    gene.cdf <- .kcdfvals_denseh5_to_denseh5(expr,
                                                             grid=grid,
                                                             verbose)
            }
        } else if (is.matrix(expr)) {
            A = .Call("matrix_density_R",
                      as.double(t(expr)),
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
            cli_abort(c("x"=sprintf("Matrix class %s cannot be handled yet.",
                                    class(expr))))
    } else { ## direct ECDF estimation
        if (is(expr, "dgCMatrix")) {
            if (sparse)
                gene.cdf <- .ecdfvals_sparse_to_sparse(expr, verbose)
            else
                gene.cdf <- .ecdfvals_sparse_to_dense(expr, verbose)
        } else if (is(expr, "SVT_SparseArray")) {
            if (sparse)
                gene.cdf <- .ecdfvals_svt_to_svt(expr, verbose)
            else
                gene.cdf <- .ecdfvals_svt_to_dense(expr, verbose)
        } else if (is(expr, "DelayedMatrix")) {
            if (sparse)
                gene.cdf <- .ecdfvals_sparseh5_to_sparseh5(expr, grid=grid,
                                                           verbose=verbose)
            else {
                if (is_sparse(expr)) ## input HDF5 may be sparse or not
                    gene.cdf <- .ecdfvals_sparseh5_to_denseh5(expr, grid=grid,
                                                              verbose)
                else
                    gene.cdf <- .ecdfvals_denseh5_to_denseh5(expr, grid=grid,
                                                             verbose)
            }
        } else if (is.matrix(expr)) {
            if (any_na)
                gene.cdf <- .ecdfvals_dense_to_dense_nas(expr, verbose)
            else
                gene.cdf <- .ecdfvals_dense_to_dense(expr, verbose)
        } else
            cli_abort(c("x"=sprintf("Input container class %s cannot be handled yet.",
                                    class(expr))))
    }

    return(gene.cdf)	
}

## here 'ties.method="last"' allows one to obtain the result
## from 'order()' based on ranks
## pending how to propagate verbosity if necessary
compute.col.ranks <- function(Z, ties.method="last", verbose=TRUE) {
    R <- NULL

    if (is(Z, "dgCMatrix")) { ## assumes expression values are positive
        R <- .sparseColumnApplyAndReplace(Z, rank, ties.method=ties.method)
    } else if (is(Z, "SVT_SparseArray")) {
        R <- .colRanks_SVT_SparseArray(Z, ties.method=ties.method)
    } else if (is(Z, "DelayedMatrix")) {
        R <- .colRanksHDF5(Z, ties.method=ties.method)
    } else {
        R <- colRanks(Z, ties.method="last", preserveShape=TRUE)
    }

    return(R)
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
  ## in the sparse case stored in a 'dgCMatrix' or a 'SVT_SparseArray',
  ## by now, use the average nonzero values per row
  if (is(expr, "dgCMatrix") || is(expr, "SVT_SparseArray"))
    return((nnzero(expr) / nrow(expr)) >= kcdf.min.ssize)

  ## in every other case, including the dense case, by now,
  ## just look at the number of columns
  return(ncol(expr) >= kcdf.min.ssize)
}

#' @importFrom S4Arrays is_sparse
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
        is_sparse_matrix <- is(expr, "dgCMatrix") ||
                            is(expr, "SVT_SparseArray") ||
                            (is(expr, "DelayedMatrix") && is_sparse(expr))
        if (is_sparse_matrix && sparse)
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
#' @param BPPARAM An object of class `BiocParallelParam` specifying parameters
#' related to the parallel execution of some of the tasks and calculations
#' within this function.
#'
#' @param maxmem A vector of length 1 either specifying a number in bytes, or
#' a character string with either the word `auto` (default), or a number
#' followed by a suffix indicating kilobytes (K), megabytes (M), gigabytes (G)
#' or terabytes (T), which GSVA will use to attempt bounding the maximum amount
#' of main memory used across all threads of execution to that given quantity.
#' By default `maxmem="auto"`, indicating that the maximum memory will be the
#' 90% of the total main memory, as calculated by [`Sys.meminfo()`][memuse::Sys.meminfo].
#' To avoid setting any bound on the maximum memory, please use `maxmem=Inf`.
#' Note that the amount of main memory used in an R session or script may depend
#' on other commands and packages used in that same session or script.
#'
#' @return In the case of the `gsvaRanks()` method, an object of class
#' [`gsvaRanksParam-class`].
#'
#' @seealso [`gsvaParam-class`], [`gsvaRanksParam-class`], [`gsva`],
#' [`BiocParallelParam`][BiocParallel::BiocParallelParam-class],
#' [`dgCMatrix`][Matrix::dgCMatrix-class],
#' \code{\link[Biobase]{ExpressionSet}},
### we are using the plain Rd above because
###  #' [`ExpressionSet`][Biobase::ExpressionSet-class],
### results in the following R CMD check NOTE:
### Non-topic package-anchored link(s) in Rd file 'gsvaRanks.Rd':
###  ‘[Biobase:class.ExpressionSet]{ExpressionSet}’
#' [`SingleCellExperiment`][SingleCellExperiment::SingleCellExperiment-class]
#'
#' @aliases gsvaRanks,gsvaParam-method
#' @name gsvaRanks
#' @rdname gsvaRanks
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' \doi{10.1186/1471-2105-14-7}
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
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto")
          {
              if (verbose && gsva_global$show_start_and_end_messages) {
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))
              }

              exprData <- get_exprData(param)
              dataMatrix <- unwrapData(exprData, get_assay(param))
              maxmem <- .check_maxmem(maxmem, verbose)
              ondisk <- .check_ondisk(param, maxmem, verbose)

              if (is(dataMatrix, "DelayedMatrix") && ondisk == "no") {
                  if (verbose)
                      cli_alert_info("Loading input expression data into main memory")
                  if (is_sparse(dataMatrix)) {
                      if (nzcount(param) < .Machine$integer.max)
                          dataMatrix <- as(dataMatrix, "dgCMatrix")
                      else
                          dataMatrix <- as(dataMatrix, "SVT_SparseArray")
                  } else
                      dataMatrix <- as.matrix(dataMatrix)
              }

              filtDataMatrix <- dataMatrix
              if (get_filterRows(param))
                  filtDataMatrix <- .filterGenes(dataMatrix, anyNA(param),
                                                 removeConstant=TRUE,
                                                 removeNzConstant=TRUE,
                                                 verbose, BPPARAM=BPPARAM,
                                                 maxmem=maxmem)
              else if (verbose)
                  cli_alert_warning("Skipping filtering of constant rows (filterRows=FALSE)")
              
              if (verbose)
                  cli_alert_info(sprintf("Calculating GSVA ranks"))

              kcdfminssize <-get_kcdfNoneMinSampleSize(param)
              gsvarnks <- .compute_gsva_ranks(expr=filtDataMatrix,
                                              kcdf=get_kcdf(param),
                                              kcdf.min.ssize=kcdfminssize,
                                              sparse=get_sparse(param),
                                              any_na=anyNA(param),
                                              na_use=get_NAuse(param),
                                              verbose=verbose,
                                              BPPARAM=BPPARAM,
                                              maxmem=maxmem)

              rownames(gsvarnks) <- rownames(filtDataMatrix)
              colnames(gsvarnks) <- colnames(filtDataMatrix)

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
                          use=get_NAuse(param), filterRows=get_filterRows(param),
                          nzcount=nzcount(param), ondisk=get_ondisk(param))

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
#' the input was a base matrix or a `dgCMatrix` object, then the output will
#' be a base matrix object with the gene sets employed in the calculations
#' stored in an attribute called `geneSets`. If the input was an
#' `ExpressionSet` object, then the output will be also an `ExpressionSet`
#' object with the gene sets employed in the calculations stored in an
#' attributed called `geneSets`. If the input was an object of one of the
#' classes described in [`GsvaExprData`], such as a `SingleCellExperiment`,
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
                   BPPARAM=SerialParam(progressbar=verbose),
                   maxmem="auto") {
              if (verbose && gsva_global$show_start_and_end_messages) {
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))
              }

              ## assuming rows in the rank data have been already filtered
              exprData <- get_exprData(param)
              filtDataMatrix <- unwrapData(exprData, get_assay(param))

              filtMappedGeneSets <- .filterAndMapGeneSets(param=param,
                                           filteredDataMatrix=filtDataMatrix,
                                           verbose=verbose)

              sparse <- get_sparse(param)
              if (sparse && !is_sparse(filtDataMatrix))
                  sparse <- FALSE

              if (verbose) {
                if (sparse)
                    cli_alert_info("GSVA sparse algorithm")
                  else
                    cli_alert_info("GSVA dense (classical) algorithm")
              }

              maxmem <- .check_maxmem(maxmem, verbose)
              ondisk <- .check_ondisk(param, maxmem, verbose)

              if (is(filtDataMatrix, "DelayedMatrix") && ondisk == "no") {
                  if (verbose)
                      cli_alert_info("Loading input expression data into main memory")
                  if (is_sparse(filtDataMatrix)) {
                      if (nzcount(param) < .Machine$integer.max)
                          filtDataMatrix <- as(filtDataMatrix, "dgCMatrix")
                      else
                          filtDataMatrix <- as(filtDataMatrix, "SVT_SparseArray")
                  } else
                      filtDataMatrix <- as.matrix(filtDataMatrix)
              }

              if (bpnworkers(BPPARAM) > 1 && nrow(filtDataMatrix) > 100 &&
                  ncol(filtDataMatrix) > 100) {
                  if (verbose) {
                      msg <- sprintf("Calculating GSVA scores with %d cores",
                                     as.integer(bpnworkers(BPPARAM)))
                      cli_alert_info(msg)
                  }
              } else {
                  if (verbose)
                      cli_alert_info("Calculating GSVA scores")
                  BPPARAM <- NULL
              }

              ondisk <- FALSE
              esreqmem <- as.numeric(length(filtMappedGeneSets)) *
                          as.numeric(ncol(filtDataMatrix)) * 8 ## 8 bytes per double
              if (esreqmem > maxmem) {
                cli_alert_warning("The resulting matrix of enrichment scores will not fit")
                cli_alert_warning("in the given maximum main memory size, the returned")
                cli_alert_warning("object will use an on-disk data structure")
                ondisk <- TRUE
              }

              gsva_es <- .processMatrixCols(filtDataMatrix,
                                            FUN=.compute_gsva_scores,
                                            geneSetsIdx=filtMappedGeneSets,
                                            tau=get_tau(param),
                                            maxDiff=get_maxDiff(param),
                                            absRanking=get_absRanking(param),
                                            sparse=sparse, any_na=anyNA(param),
                                            na_use=get_NAuse(param),
                                            minSize=get_minSize(param),
                                            ondisk=ondisk, verbose=verbose,
                                            minparrows=100, minparcols=100,
                                            BPPARAM=BPPARAM,
                                            maxmem=ceiling(maxmem/100)) ## use
                                            ## of memory increases here about
                                            ## 10-fold over block size memory

              rownames(gsva_es) <- names(filtMappedGeneSets)
              colnames(gsva_es) <- colnames(filtDataMatrix)

              gs <- .geneSetsIndices2Names(indices=filtMappedGeneSets,
                                           names=rownames(filtDataMatrix))
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
#' \doi{10.1186/1471-2105-14-7}
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
              filtDataMatrix <- unwrapData(exprData, get_assay(param))

              ## no need for verbosity when mapping a single gene set
              filtMappedGeneSets <- .filterAndMapGeneSets(param,
                                           wgset=geneSet,
                                           filteredDataMatrix=filtDataMatrix,
                                           verbose=FALSE)

              geneSetIdx <- filtMappedGeneSets[[1]]
              edata <- .gsva_enrichment_data(R=filtDataMatrix,
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

#' @importFrom cli cli_progress_update
#' @importFrom parallel splitIndices
.col_iter_idx <- function(X, idpb, n_chunks) {
    idx <- splitIndices(ncol(X), min(ncol(X), n_chunks))
    i <- 0L
    function() {
        if (i == length(idx))
            return(NULL)
        i <<- i + 1L
        if (!is.null(idpb))
            cli_progress_update(id=idpb, set=i)
        idx[[i]]
    }
}

#' @importFrom IRanges IntegerList match
#' @importFrom BiocParallel bpnworkers
#' @importFrom cli cli_alert_info cli_progress_bar
#' @importFrom cli cli_progress_done cli_abort
#' @importFrom BiocGenerics rbind cbind
#' @importFrom sparseMatrixStats colRanks
.compute_gsva_ranks <- function(expr, kcdf, kcdf.min.ssize,
                                sparse, any_na, na_use, verbose,
                                BPPARAM=NULL, maxmem=Inf) {

    kcdfparam <- .parse_kcdf_param(expr, kcdf, kcdf.min.ssize, sparse, verbose)
    kernel <- kcdfparam$kernel
    Gaussk <- kcdfparam$Gaussk

    if (bpnworkers(BPPARAM) > 1 && nrow(expr) > 100 && ncol(expr) > 100) {
        if (verbose) {
            msg <- sprintf("Calculating row ECDFs with %d cores",
                           as.integer(bpnworkers(BPPARAM)))
            cli_alert_info(msg)
        }
    } else {
        if (verbose)
            cli_alert_info("Calculating row ECDFs")
        BPPARAM <- NULL
    }

    Z <- .processMatrixRows(expr, FUN=compute.gene.cdf, Gaussk=Gaussk,
                            kernel=kernel, sparse=sparse, any_na=any_na,
                            na_use=na_use, verbose=verbose, minparrows=100,
                            minparcols=100, BPPARAM=BPPARAM, maxmem=maxmem)

    if (!is.null(BPPARAM) && verbose) {
        msg <- sprintf("Calculating column ranks with %d cores",
                       as.integer(bpnworkers(BPPARAM)))
        cli_alert_info(msg)
    } else if (verbose)
        cli_alert_info("Calculating column ranks")
 
    ## here 'ties.method="last"' allows one to obtain the result
    ## from 'order()' based on ranks
    R <- .processMatrixCols(Z, FUN=compute.col.ranks, ties.method="last",
                            verbose=verbose, minparrows=100, minparcols=100,
                            BPPARAM=BPPARAM, maxmem=maxmem)

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
## r is a matrix of features x samples/cells

#' @importFrom MatrixGenerics colMaxs
#' @importFrom BiocGenerics which
.ranks2stats_block <- function(r, sparse) {
    stopifnot(length(dim(r)) == 2) ## QC
    mask <- unname(as.matrix(r)) == 0L
    p <- nrow(r)
    r_dense <- as.matrix(unname(r))           ## convert to dense
    mode(r_dense) <- "integer"                ## assume ranks are integer
    wh <- NULL

    if (any(mask)) {                          ## sparse ranks into dense ranks
        nzs <- colSums(mask)
        mode(nzs) <- "integer"
        nzsmat <- matrix(nzs, nrow=nrow(r), ncol=ncol(r), byrow=TRUE)
        wh <- which(mask, arr.ind=TRUE)
        nzsmat[wh] <- 0L
        r_dense <- r_dense + nzsmat           ## shift ranks of nonzero values
        r_dense[wh] <- unlist(sapply(nzs,     ## zeros get increasing ranks
                                     seq.int))
    }

    dos <- p - r_dense + 1L                   ## dense ranks into
                                              ## decreasing order stats
    srs <- NULL
    if (any(mask) && sparse) {
        r[wh] <- 1L                   ## all zeros get the same first rank
        wh <- which(!mask, arr.ind=TRUE)
        r[wh] <- r[wh] + 1L       ## shift ranks of nonzero values by one
        maxrmat <- matrix(colMaxs(r)/2, nrow=nrow(r), ncol=ncol(r), byrow=TRUE)
        srs <- as.matrix(abs(maxrmat - r))
    } else
        srs <- abs(p/2 - r_dense)

    list(dos=dos, srs=srs)
}

## convert ranks into decreasing order statistics and symmetric rank statistics
## skipping NA values, r is a matrix of features x samples/cells
.ranks2stats_nas_block <- function(r, sparse) {
    stopifnot(length(dim(r)) == 2) ## QC
    mask <- unname(as.matrix(r)) == 0L
    na_mask <- is.na(mask)

    if (all(na_mask))
        return(list(dos=matrix(NA_integer_, nrow(r), ncol(r)),
                    srs=matrix(NA_real_, nrow(r), ncol(r))))

    n_nas <- colSums(na_mask)
    mode(n_nas) <- "integer"
    mask <- !na_mask & mask
    p <- nrow(r)
    r_dense <- as.matrix(r)
    mode(r_dense) <- "integer"          ## assume ranks are integer

    if (any(mask)) {                    ## sparse ranks into dense ranks
        nzs <- colSums(mask)
        mode(nzs) <- "integer"
        nzsmat <- matrix(nzs, nrow=nrow(r), ncol=ncol(r), byrow=TRUE)
        nzsmat[mask] <- 0L
        r_dense <- r_dense + nzsmat            ## shift ranks of nonzero values
        r_dense[!mask] <- r_dense[!mask] + nzs ## shift ranks of nonzero values
        r_dense[mask] <- unlist(sapply(nzs,    ## zeros get increasing ranks
                                       seq.int))
    }

    n_nasmat <-  matrix(n_nas, nrow=nrow(r), ncol=ncol(r), byrow=TRUE)
    dos <- p - n_nasmat - r_dense + 1L   ## dense ranks into decreasing order stats

    srs <- NULL
    if (any(mask) && sparse) {
        r[!mask] <- r[!mask] + 1L     ## shift ranks of nonzero values by one
        r[mask] <- 1L                 ## all zeros get the same first rank
        maxrmat <- matrix(colMaxs(r, na.rm=TRUE)/2, nrow=nrow(r), ncol=ncol(r),
                          byrow=TRUE)
        srs <- abs(maxrmat - r)
    } else
        srs <- abs((p - n_nasmat)/2 - r_dense)

    list(dos=dos, srs=srs)
}


## convert ranks into decreasing order statistics and symmetric rank statistics
## skipping NA values
.ranks2stats_nas <- function(r, sparse) {
    na_mask <- is.na(r)

    if (all(na_mask))
        return(list(dos=rep(NA, length(r)), srs=rep(NA, length(r))))

    n_nas <- sum(na_mask)
    mode(n_nas) <- "integer"
    mask <- !na_mask & r == 0
    p <- length(r)
    r_dense <- as.integer(r)          ## assume ranks are integer

    if (any(mask)) {                  ## sparse ranks into dense ranks
        nzs <- sum(mask)
        mode(nzs) <- "integer"
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


#' @importFrom cli cli_alert_info cli_alert_warning
#' @importFrom BiocParallel bpnworkers
#' @importFrom S4Arrays is_sparse refdim
.compute_gsva_scores <- function(R, geneSetsIdx, tau, maxDiff, absRanking,
                                 sparse, any_na, na_use, minSize, ondisk,
                                 verbose) {
    p <- nrow(R)
    n <- ncol(R)
    es <- NULL
    if (sparse && !is_sparse(R))
        sparse <- FALSE

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)
    es <- NULL
    if (is(R, "DelayedMatrix") || ondisk) {
        sink <- HDF5RealizationSink(c(length(geneSetsIdx), ncol(R)),
                                    as.sparse=FALSE) ## GSVA scores are dense
        grid <- DummyArrayGrid(dim(R))
        grid_es <- DummyArrayGrid(dim(sink))

        if (length(grid) != length(grid_es) ||
            refdim(grid)[2] != refdim(grid_es)[2] ||
            dim(grid)[2] != dim(grid_es)[2]) {
            cli_abort(c("x"="Grid column blocks for ranks should match grid column blocks for enrichment scores"))
        }

        ## avp - ArrayViewport for reaching the (possibly sparse) rank matrix
        ## avp_es - ArrayViewport for writing the enrichment dense scores matrix
        colScores_byBlock <- function(avp, avp_es, sink) {
            block <- read_block(R, avp)
            rnkstats <- NULL
            if (any_na)
                rnkstats <- .ranks2stats_nas_block(block, sparse)
            else
                rnkstats <- .ranks2stats_block(block, sparse)
            block <- .gsva_score_genesets(NULL, geneSetsIdx,
                                          decOrdStat=rnkstats$dos,
                                          symRnkStat=rnkstats$srs,
                                          maxDiff, absRanking, tau,
                                          any_na, na_use, minSize,
                                          wna_env, verbose=verbose)
            write_block(sink, avp_es, block)
        }

        nblock <- length(grid)
        for (bid in seq_len(nblock))
            sink <- colScores_byBlock(grid[[bid]], grid_es[[bid]], sink)
        close(sink)
        es <- as(sink, "DelayedArray")

    } else {
        rnkstats <- NULL
         if (any_na)
             rnkstats <- .ranks2stats_nas_block(R, sparse)
         else
             rnkstats <- .ranks2stats_block(R, sparse)

         es <- .gsva_score_genesets(NULL, geneSetsIdx,
                                    decOrdStat=rnkstats$dos,
                                    symRnkStat=rnkstats$srs,
                                    maxDiff, absRanking, tau, any_na,
                                    na_use, minSize, wna_env, verbose=verbose)
    }

    if (any_na && na_use == "na.rm")
        if (get("w", envir=wna_env)) {
            msg <- sprintf(paste("NA enrichment scores in gene sets with less",
                                 "than %d genes after removing missing values"),
                           minSize)
            cli_alert_warning(msg)
        }

    return(es)
}

#' @importFrom S4Arrays is_sparse
.gsva_enrichment_data <- function(R, column, geneSetIdx, maxDiff,
                                  absRanking, tau, sparse, any_na,
                                  na_use, minSize) {
    n <- ncol(R)
    es <- NULL
    if (!is_sparse(R))
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
    .data <- get(".data")
    ggplot2::ggplot(data=edata$stats) +
        ggplot2::scale_x_continuous(breaks=fintticks) +
        ggplot2::geom_line(ggplot2::aes(x=.data$rank, y=.data$stat), color="green") +
        ggplot2::geom_segment(data=gsetticks,
                     mapping=ggplot2::aes(x=.data$gsetrnk, y=-hgsetticks/2,
                                 xend=.data$gsetrnk, yend=hgsetticks/2),
                     linewidth=1) +
        ggplot2::geom_hline(yintercept=0, colour="grey", linetype="dashed") +
        { if (!is.na(edata$whichMaxPos) &&
              (edata$maxDiff || edata$maxPos >= abs(edata$maxNeg)))
              ggplot2::geom_segment(data=data.frame(whichMaxPos=edata$whichMaxPos,
                                           maxPos=edata$maxPos),
                           mapping=ggplot2::aes(x=.data$whichMaxPos, y=0,
                                       xend=.data$whichMaxPos, yend=.data$maxPos),
                           colour="darkred", linetype="dashed") } +
        { if (!is.na(edata$whichMaxPos) &&
              (edata$maxDiff || edata$maxPos < abs(edata$maxNeg)))
              ggplot2::geom_segment(data=data.frame(whichMaxNeg=edata$whichMaxNeg,
                                           maxNeg=edata$maxNeg),
                           mapping=ggplot2::aes(x=.data$whichMaxNeg, y=0,
                                       xend=.data$whichMaxNeg, yend=.data$maxNeg),
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

.fetch_row_nzvals <- function(X, i, whimin1=NULL) {
  stopifnot(is(X, "SVT_SparseArray")) ## QC
  stopifnot(is.numeric(i)) ## QC
  if (!is.null(whimin1)) {
      stopifnot(is.numeric(whimin1)) ## QC
      whimin1 <- as.integer(whimin1)
      stopifnot(length(whimin1) == ncol(X)) ## QC
  }
  .Call("fetch_row_nzvals_R", X, as.integer(i), whimin1)
}

.ecdfvals_svt_to_dense <- function(X, verbose) {
  stopifnot(is(X, "SVT_SparseArray")) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_svt_to_dense_R", X, verbose)
}

.ecdfvals_svt_to_sparse <- function(X, verbose) {
  stopifnot(is(X, "SVT_SparseArray")) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_svt_to_sparse_R", X, verbose)
}

.ecdfvals_svt_to_svt <- function(X, verbose) {
  stopifnot(is(X, "SVT_SparseArray")) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_svt_to_svt_R", X, verbose)
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom DelayedArray seed gridReduce
.ecdfvals_sparseh5_to_sparseh5 <- function(X, grid=NULL, verbose=FALSE) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=TRUE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowEcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .ecdfvals_svt_to_svt(block, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowEcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom S4Arrays DummyArrayGrid
#' @importFrom DelayedArray seed gridReduce
.ecdfvals_sparseh5_to_denseh5 <- function(X, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowEcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .ecdfvals_svt_to_dense(block, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowEcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce
.ecdfvals_denseh5_to_denseh5 <- function(X, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowEcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .ecdfvals_dense_to_dense(block, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowEcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

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

.kcdfvals_svt_to_dense <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "SVT_SparseArray")) ## QC
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_svt_to_dense_R", X, Gaussk, verbose)
}

.kcdfvals_svt_to_svt <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "SVT_SparseArray")) ## QC
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_svt_to_svt_R", X, Gaussk, verbose)
}

.kcdfvals_sparse_to_sparse <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_sparse_to_sparse_R", X, Xrsp, Gaussk, verbose)
}

.kcdfvals_sparse_to_dense <- function(X, Gaussk, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(Gaussk)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  .Call("kcdfvals_sparse_to_dense_R", X, Xrsp, Gaussk, verbose)
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce
.kcdfvals_sparseh5_to_sparseh5 <- function(X, Gaussk, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=TRUE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowKcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .kcdfvals_svt_to_svt(block, Gaussk=Gaussk, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowKcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce
.kcdfvals_sparseh5_to_denseh5 <- function(X, Gaussk, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowKcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- .kcdfvals_svt_to_dense(block, Gaussk=Gaussk, verbose=verbose)
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowKcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
}

#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray seed rowAutoGrid blockReduce
.kcdfvals_denseh5_to_denseh5 <- function(X, Gaussk, grid=NULL, verbose) {
  stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

  sink <- HDF5RealizationSink(dim(X), as.sparse=FALSE)
  if (is.null(grid))
      grid <- DummyArrayGrid(dim(X))

  rowKcdf_byBlock <- function(grid, sink) {
    block <- read_block(X, grid)
    block <- t(matrix(.Call("matrix_density_R",
                            as.double(t(X)),
                            as.double(t(X)),
                            ncol(X),
                            ncol(X),
                            nrow(X),
                            as.integer(Gaussk),
                            FALSE, 1L,
                            verbose), ncol(X), nrow(X)))
    write_block(sink, grid, block)
  }
  sink <- gridReduce(rowKcdf_byBlock, grid, sink)
  close(sink)
  res <- as(sink, "DelayedArray")
  res
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
.gsva_score_genesets <- function(colIdx, geneSetsIdx, decOrdStat, symRnkStat,
                                 maxDiff, absRanking, tau, any_na, na_use,
                                 minSize, wna_env, verbose) {
  minSize <- as.integer(minSize)
  stopifnot(is.null(colIdx) || is.integer(colIdx)) ## QC
  stopifnot(is.list(geneSetsIdx)) ## QC
  stopifnot(length(geneSetsIdx) > 0) ## QC
  stopifnot(is.integer(geneSetsIdx[[1]])) ## QC
  stopifnot(is.integer(decOrdStat)) ## QC
  stopifnot(is.numeric(symRnkStat)) ## QC
  stopifnot(all(dim(decOrdStat) == dim(symRnkStat))) ## QC
  stopifnot(is.logical(maxDiff)) ## QC
  stopifnot(is.logical(absRanking)) ## QC
  stopifnot(is.numeric(tau)) ## QC but it still might be an integer!!
  stopifnot(is.logical(any_na)) ## QC
  stopifnot(is.character(na_use)) ## QC
  stopifnot(is.integer(minSize)) ## QC
  stopifnot(is.logical(verbose)) ## QC
  na_use <- as.integer(factor(na_use, levels=c("everything", "all.obs",
                                               "na.rm")))
  if (is.null(colIdx)) {
    if (is.null(dim(decOrdStat)))
      decOrdStat <- matrix(decOrdStat, ncol=1)
    if (is.null(dim(symRnkStat)))
      symRnkStat <- matrix(symRnkStat, ncol=1)
  } else {
    decOrdStat <- decOrdStat[, colIdx, drop=FALSE]
    symRnkStat <- symRnkStat[, colIdx, drop=FALSE]
  }
  sco <- .Call("gsva_score_genesets_R", geneSetsIdx, decOrdStat, symRnkStat,
               maxDiff, absRanking, as.double(tau), any_na, na_use, minSize,
               verbose)
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

## calculate ranks using on an SVT_SparseArray object
.colRanks_SVT_SparseArray <- function(X, ties.method="last") {
    R <- X
    rnks <- lapply(lapply(X@SVT, "[[", 1), rank, ties.method=ties.method)
    R@type <- "integer" ## rank() w/ ties.method="last" returns integer
    R@SVT <- mapply(list, rnks, lapply(X@SVT, "[[", 2), SIMPLIFY=FALSE)
    R
}

## calculate ranks using an HDF5 backend

#' @importFrom MatrixGenerics colRanks
#' @importFrom BiocParallel SerialParam
.colRanksHDF5 <- function(X, grid=NULL, ties.method="last") {
    stopifnot(is(X, "DelayedMatrix") || is(X, "HDF5Matrix")) ## QC

    sink <- HDF5RealizationSink(dim(X), as.sparse=is_sparse(X))
    if (is.null(grid))
        grid <- DummyArrayGrid(dim(X))

    colRanks_byBlock <- function(grid, sink) {
        block <- read_block(X, grid)
        if (is(block, "SVT_SparseArray")) {
            block <- .colRanks_SVT_SparseArray(block, ties.method=ties.method)
        } else {
            block <- colRanks(block, ties.method=ties.method,
                              preserveShape=TRUE)
            if (ties.method == "last")
                mode(block) <- "integer"
        }
        write_block(sink, grid, block)
    }

    sink <- gridReduce(colRanks_byBlock, grid, sink)
    close(sink)
    res <- as(sink, "DelayedArray")
    res
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
