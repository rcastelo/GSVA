##
## methods for the PLAGE method from Tomfohr et al. (2005)
##

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
