##
## methods for the z-score method from Lee et al. (2008)
##

combinez <- function(gSetIdx, j, Z) sum(Z[gSetIdx, j]) / sqrt(length(gSetIdx))

#' @importFrom cli cli_alert_info
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom BiocParallel bpnworkers
#' @importFrom Matrix colSums
zscore <- function(X, geneSets, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {
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

    if (bpnworkers(BPPARAM) > 1)
        es <- bplapply(geneSets, function(gSetIdx) {
                           colSums(Z[gSetIdx, , drop=FALSE]) / sqrt(length(gSetIdx))
                       }, BPPARAM=BPPARAM)
    else {
        idpb <- NULL
        if (verbose)
            idpb <- cli_progress_bar("Calculating Z-scores",
                                     total=length(geneSets))
        es <- lapply(geneSets, function(gSetIdx, verbose, idpb) {
                           if (verbose)
                               cli_progress_update(idpb)
                           colSums(Z[gSetIdx, , drop=FALSE]) / sqrt(length(gSetIdx))
                     }, verbose, idpb)
        if (verbose)
            cli_progress_done(idpb)
    }
    
    es <- do.call(rbind, es)

    es
}
