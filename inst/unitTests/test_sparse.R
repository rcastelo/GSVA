test_sparseMethods <- function(){
    message("Running unit tests for sparse methods")
    
    suppressPackageStartupMessages({
        library(Matrix)
        library(cli) ## for cli_fmt()
    })

    p <- 50 ## number of genes
    n <- 100 ## number of samples/cells
    m <- 5 ## number of gene sets

    ## create random gene sets of size between 3 and 10
    set.seed(123)
    gssizes <- sample(3:10, size=m, replace=TRUE)
    gsets <- lapply(gssizes, function(x) sample(paste0("g", 1:p), size=x, replace=FALSE))
    names(gsets) <- paste0("gs", 1:m)


    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    s <- ceiling(0.15 * p * n)
    sam <- sample(1:(p * n), size=s, replace=FALSE)
    x <- numeric(p * n)
    x[sam] <- runif(s)
    m <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    M <- Matrix(m, sparse=TRUE)
    
    checkException(mg <- gsva(gsvaParam(x, gsets), verbose=TRUE))
    out <- cli_fmt(mg <- gsva(gsvaParam(m, gsets), verbose=TRUE))
    out <- cli_fmt(Mg <- gsva(gsvaParam(M, gsets, sparse=FALSE), verbose=TRUE))
    checkEqualsNumeric(mg, Mg)

    out <- cli_fmt(mp <- gsva(plageParam(m, gsets), verbose=TRUE))
    out <- cli_fmt(Mp <- gsva(plageParam(M, gsets), verbose=TRUE))
    checkEqualsNumeric(mp, Mp)
    
    out <- cli_fmt(mz <- gsva(zscoreParam(m, gsets), verbose=TRUE))
    out <- cli_fmt(Mz <- gsva(zscoreParam(M, gsets), verbose=TRUE))
    checkEqualsNumeric(mz, Mz)
    
    out <- cli_fmt(ms <- gsva(ssgseaParam(m, gsets), verbose=TRUE))
    out <- cli_fmt(Ms <- gsva(ssgseaParam(M, gsets), verbose=TRUE))
    checkEqualsNumeric(ms, Ms)

    M[1, 1] <- NA
    checkException(Mg <- gsva(gsvaParam(M, gsets), verbose=TRUE))
    
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

    ## on-disk input is processed by blocks of rows, read into main memory, and
    ## its output written on disk, see GSVA:::.ondisk_blocks()
    zz <- as(zz, "HDF5Array")
    res_C_denseh5_to_denseh5 <- GSVA:::.processMatrixRows(zz, GSVA:::compute.gene.cdf, Gaussk=FALSE, kernel=FALSE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_dense, res_C_denseh5_to_denseh5)
    zzs <- as(zzs, "HDF5Array")
    res_C_sparseh5_to_denseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=FALSE, kernel=FALSE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_dense, res_C_sparseh5_to_denseh5)
    res_C_sparseh5_to_sparseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=FALSE, kernel=FALSE, sparse=TRUE, verbose=FALSE, sinkout=TRUE)
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
    ## on-disk input is processed by blocks of rows, read into main memory, and
    ## its output written on disk, see GSVA:::.ondisk_blocks()
    zz <- as(zz, "HDF5Array")
    res_C_denseh5_to_denseh5 <- GSVA:::.processMatrixRows(zz, GSVA:::compute.gene.cdf, Gaussk=TRUE, kernel=TRUE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_dense_to_dense, res_C_denseh5_to_denseh5, tolerance=0.0001)
    zzs <- as(zzs, "HDF5Array")
    res_C_sparseh5_to_denseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=TRUE, kernel=TRUE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_sparseh5_to_denseh5, tolerance=0.00001)
    res_C_sparseh5_to_sparseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=TRUE, kernel=TRUE, sparse=TRUE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparseh5_to_sparseh5, tolerance=0.0001)
}

test_svt_lacunar <- function() {
    message("Running unit tests for ECDF and KCDF values on lacunar SVT leaves")

    suppressPackageStartupMessages({
        library(Matrix)
        library(SparseArray)
    })

    ## sparse count data, where the first column has only ones as nonzero
    ## values, which SparseArray stores in a 'lacunar' leaf without values,
    ## and the fourth column has only zeros
    m <- matrix(0L, nrow=6, ncol=5)
    m[c(1, 3, 5), 1] <- 1L
    m[c(1, 2, 4), 2] <- c(2L, 3L, 1L)
    m[c(1, 3, 5, 6), 3] <- c(1L, 4L, 2L, 7L)
    m[c(2, 3, 4, 6), 5] <- c(5L, 1L, 2L, 3L)

    for (type in c("integer", "double")) {
        storage.mode(m) <- type
        svt <- as(m, "SVT_SparseMatrix")
        ## check that the first column is stored in a lacunar leaf
        checkTrue(is.null(svt@SVT[[1]][[1]]))
        dgc <- as(m, "dgCMatrix")

        res_svt <- GSVA:::.ecdfvals_svt_to_svt(svt, FALSE)
        checkEqualsNumeric(as.matrix(GSVA:::.ecdfvals_sparse_to_sparse(dgc,
                                                                       FALSE)),
                           as.matrix(res_svt))

        for (Gaussk in c(TRUE, FALSE)) {
            res_svt <- GSVA:::.kcdfvals_svt_to_svt(svt, Gaussk, FALSE)
            res_dgc <- GSVA:::.kcdfvals_sparse_to_sparse(dgc, Gaussk, FALSE)
            checkEqualsNumeric(as.matrix(res_dgc), as.matrix(res_svt))
        }
    }

    ## a matrix with only zeros, whose SVT is NULL
    svt <- as(matrix(0L, nrow=3, ncol=2), "SVT_SparseMatrix")
    checkEqualsNumeric(matrix(0, nrow=3, ncol=2),
                       as.matrix(GSVA:::.ecdfvals_svt_to_svt(svt, FALSE)))
    checkEqualsNumeric(matrix(0, nrow=3, ncol=2),
                       as.matrix(GSVA:::.kcdfvals_svt_to_svt(svt, TRUE, FALSE)))
}
