##
## methods for the PLAGE method from Tomfohr et al. (2005)
##

#' @importFrom cli cli_progress_update
rightsingularsvdvectorgset <- function(gSetIdx, Z, verbose, idpb) {
  if(is(Z, "dgCMatrix")){
    s <- BiocSingular::runExactSVD(Z[gSetIdx, ])
  } else {
    s <- svd(Z[gSetIdx, ])
  }
  if (verbose)
      cli_progress_update(id=idpb)
  s$v[, 1]
}

#' @importFrom cli cli_alert_info
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom BiocParallel bpnworkers
plage <- function(X, geneSets, verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {
    Z <- NULL
    if (is(X, "dgCMatrix")){
        if (verbose)
            cli_alert_info("Centering and scaling non-zero values")

        Z <- t(.sparseColumnApplyAndReplace(t(X), FUN=scale))
    } else if (is.matrix(X)) {
        if (verbose)
            cli_alert_info("Centering and scaling values")

        Z <- t(scale(t(X)))
    } else
        stop(sprintf("Matrix class %s cannot be handled yet.", class(X)))

    if (bpnworkers(BPPARAM) > 1) {
        if (verbose)
            cli_progress_bar("Calculating PLAGE scores")
        es <- bplapply(geneSets, rightsingularsvdvectorgset, Z,
                       verbose=FALSE, idpb=character(0), BPPARAM=BPPARAM)
    } else {
        idpb <- NULL
        if (verbose)
            idpb <- cli_progress_bar("Calculating PLAGE scores",
                                     total=length(geneSets))
        es <- lapply(geneSets, rightsingularsvdvectorgset, Z,
                     verbose, idpb)
        if (verbose)
            cli_progress_done(idpb)
    }
        
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
    
    es
}
