test_sparseMethods <- function(){
    message("Running unit tests for sparse methods")
    
    set.seed(123)
    
    m <- matrix(runif(100), 10, 10)
    colnames(m) <- paste0("cell_", 1:10)
    rownames(m) <- paste0("gene_", 1:10)
    
    gene.sets <- list("my_list1"= paste0("gene_", 1:2),
                      "my_list2"= paste0("gene_", 3:4))
    
    suppressPackageStartupMessages(library(Matrix))
    M <- Matrix(m, sparse=TRUE)

    suppressPackageStartupMessages(library(cli)) ## for cli_fmt()

    out <- cli_fmt(mg <- gsva(gsvaParam(m, gene.sets), verbose=TRUE))
    out <- cli_fmt(Mg <- gsva(gsvaParam(M, gene.sets, sparse=FALSE, verbose=TRUE), verbose=TRUE))
    checkEqualsNumeric(mg, Mg)
    
    out <- cli_fmt(mp <- gsva(plageParam(m, gene.sets), verbose=TRUE))
    out <- cli_fmt(Mp <- gsva(plageParam(M, gene.sets), verbose=TRUE))
    checkEqualsNumeric(mp, Mp)
    
    out <- cli_fmt(mz <- gsva(zscoreParam(m, gene.sets), verbose=TRUE))
    out <- cli_fmt(Mz <- gsva(zscoreParam(M, gene.sets), verbose=TRUE))
    checkEqualsNumeric(mz, Mz)
    
    out <- cli_fmt(ms <- gsva(ssgseaParam(m, gene.sets), verbose=TRUE))
    out <- cli_fmt(Ms <- gsva(ssgseaParam(M, gene.sets), verbose=TRUE))
    checkEqualsNumeric(ms, Ms)
}

test_ecdfvals <- function() {
    message("Running unit tests for ECDF values calculations")

    ecdfvals_dense <- function(X) t(apply(X, 1, function(rx) ecdf(rx)(rx)))
    ecdfvals_sparse_to_sparse <- function(X) {
        stopifnot(is(X, "dgCMatrix")) ## QC
        for (i in 1:nrow(X)) {
            rx <- X[i, , drop=FALSE]
            vals <- sort(unique(rx@x))
            mt <- match(rx@x, vals)
            tab <- tabulate(mt, nbins=length(vals))
            ecdfvals <- cumsum(tab) / nnzero(rx)
            X[i, which(diff(rx@p) > 0)] <- ecdfvals[mt]
        }
        X
    }

    suppressPackageStartupMessages({
        library(Matrix)
	library(HDF5Array)
	library(SparseArray)
    })

    n <- 100
    p <- 100
    z <- numeric(p * n)
    nnz <- ceiling(0.05 * p * n) ## 5% nonzero values
    z[sample(1:(p*n), size=nnz, replace=FALSE)] <- rnorm(nnz)
    zz <- matrix(z, nrow=p, ncol=n)
    zzs <- Matrix(zz, sparse=TRUE)
    res_R_dense <- ecdfvals_dense(zz)
    res_C_dense <- GSVA:::compute.gene.cdf(zz, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_dense)

    res_C_sparse_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_sparse_to_dense)

    res_R_sparse_to_sparse <- ecdfvals_sparse_to_sparse(zzs)
    res_C_sparse_to_sparse <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparse_to_sparse)

    zzs <- SparseArray(zzs)
    res_C_svt_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_svt_to_dense)
    res_C_svt_to_svt <- GSVA:::.ecdfvals_svt_to_svt(zzs, FALSE)
    checkEqualsNumeric(SparseArray(res_C_sparse_to_sparse), res_C_svt_to_svt)

    zz <- as(zz, "HDF5Array")
    res_C_denseh5_to_denseh5 <- GSVA:::compute.gene.cdf(zz, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_denseh5_to_denseh5)
    zzs <- as(zzs, "HDF5Array")
    res_C_sparseh5_to_denseh5 <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_sparseh5_to_denseh5)
    res_C_sparseh5_to_sparseh5 <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparseh5_to_sparseh5)
}

test_kcdfvals <- function() {
    message("Running unit tests for KCDF values calculations")

    kcdfegaussianvals_dense_to_dense <- function(x) {
        x <- as.matrix(x)
        t(apply(x, 1, function(rx, n) {
           bw <- sd(rx) / 4
	   if (is.na(bw) || bw == 0)
               bw = 0.001
           kecdf <- rowSums(pnorm(outer(rx, rx, FUN="-")/bw)/n)
	   -log((1-kecdf)/kecdf)
        }, ncol(x)))
    }

    kcdfegaussianvals_sparse_to_dense <- function(x) {
        x <- as.matrix(x)
        t(apply(x, 1, function(rx, n) {
           bw <- sd(rx) / 4
	   if (is.na(bw) || bw == 0)
               bw = 0.001
           rowSums(pnorm(outer(rx, rx, FUN="-")/bw)/n)
        }, ncol(x)))
    }

    kcdfgaussianvals_sparse_to_sparse <- function(X) {
        stopifnot(is(X, "dgCMatrix")) ## QC
        for (i in 1:nrow(X)) {
            rx <- X[i, , drop=FALSE]
            bw <- sd(rx@x) / 4
	    if (is.na(bw) || bw == 0)
                bw = 0.001
            kcdfvals <- rowSums(pnorm(outer(rx@x, rx@x, FUN="-")/bw)/nnzero(rx))
            X[i, which(diff(rx@p) > 0)] <- kcdfvals
        }
        X
    }

    suppressPackageStartupMessages({
        library(Matrix)
	library(HDF5Array)
	library(SparseArray)
    })

    n <- 100
    p <- 100
    z <- numeric(p * n)
    nnz <- ceiling(0.05 * p * n) ## 5% nonzero values
    z[sample(1:(p*n), size=nnz, replace=FALSE)] <- rnorm(nnz)
    zz <- matrix(z, nrow=p, ncol=n)
    zzs <- Matrix(zz, sparse=TRUE)

    res_R_sparse_to_dense <- kcdfegaussianvals_sparse_to_dense(zzs)
    res_C_sparse_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=FALSE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_sparse_to_dense, tolerance=0.00001)

    res_R_sparse_to_sparse <- kcdfgaussianvals_sparse_to_sparse(zzs)
    res_C_sparse_to_sparse <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparse_to_sparse, tolerance=0.0001)

    zzs <- SparseArray(zzs)
    res_C_svt_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=FALSE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_svt_to_dense, tolerance=0.00001)
    res_C_svt_to_svt <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=TRUE)
    checkEqualsNumeric(SparseArray(res_C_sparse_to_sparse), res_C_svt_to_svt)

    res_R_dense_to_dense <- kcdfegaussianvals_dense_to_dense(zz)
    zz <- as(zz, "HDF5Array")
    res_C_denseh5_to_denseh5 <- GSVA:::compute.gene.cdf(zz, Gaussk=TRUE, kernel=TRUE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense_to_dense, res_C_denseh5_to_denseh5, tolerance=0.0001)
    zzs <- as(zzs, "HDF5Array")
    res_C_sparseh5_to_denseh5 <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=FALSE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_sparseh5_to_denseh5, tolerance=0.00001)
    res_C_sparseh5_to_sparseh5 <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparseh5_to_sparseh5, tolerance=0.0001)
}
