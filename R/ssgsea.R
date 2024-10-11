##
## methods for the ssGSEA method from Barbie et al. (2009)
##

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
#' @importFrom BiocParallel bpnworkers
#' @importFrom sparseMatrixStats colRanks
#' @importFrom cli cli_alert_info cli_alert_warning
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
ssgsea <- function(X, geneSets, alpha=0.25,
                   normalization=TRUE,
                   check_na=FALSE,
                   any_na=FALSE,
                   na_use=c("everything", "all.obs", "na.rm"),
                   minSize=1, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {
  na_use <- match.arg(na_use)

  if (verbose)
    cli_alert_info("Calculating ranks")
  
  R <- t(colRanks(X, ties.method="average"))
  mode(R) <- "integer"

  if (verbose)
    cli_alert_info("Calculating rank weights")
  
  Ra <- abs(R)^alpha
  
  wna_env <- new.env()
  assign("w", FALSE, envir=wna_env)
  geneSets <- IntegerList(geneSets)
  n <- ncol(X)
  es <- NULL
  if (n > 10 && bpnworkers(BPPARAM) > 1) {
  
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
  } else {
    idpb <- NULL
    if (verbose)
      idpb <- cli_progress_bar("Calculating ssGSEA scores", total=n)
    es <- lapply(as.list(1:n), function(j) {
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
      if (verbose)
        cli_progress_update(id=idpb)
      es_sample
    })
    if (verbose)
      cli_progress_done(idpb)
  }
  es <- do.call("cbind", es)
  
  if (any_na && na_use =="na.rm")
    if (get("w", envir=wna_env)) {
      msg <- sprintf(paste("NA enrichment scores in gene sets with less than",
                           "%d genes after removing missing values", minSize))
      cli_alert_warning(msg)
    }
  
  if (normalization) {
    if (verbose)
      cli_alert_info("Normalizing ssGSEA scores")
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
