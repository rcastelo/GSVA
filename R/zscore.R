##
## methods for the z-score method from Lee et al. (2008)
##

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
