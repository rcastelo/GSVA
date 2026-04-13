test_sparseMethods <- function(){
    message("Running unit tests for sparse methods.")
    
    set.seed(123)
    
    m <- matrix(runif(100), 10, 10)
    colnames(m) <- paste0("cell_", 1:10)
    rownames(m) <- paste0("gene_", 1:10)
    
    gene.sets <- list("my_list1"= paste0("gene_", 1:2),
                      "my_list2"= paste0("gene_", 3:4))
    
    M <- as(as(as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix")

    mg <- gsva(gsvaParam(m, gene.sets), verbose=FALSE)
    Mg <- gsva(gsvaParam(M, gene.sets, sparse=FALSE), verbose=FALSE)
    checkEqualsNumeric(mg, Mg)
    
    mp <- gsva(plageParam(m, gene.sets), verbose=FALSE)
    Mp <- gsva(plageParam(M, gene.sets), verbose=FALSE)
    checkEqualsNumeric(mp, Mp)
    
    mz <- gsva(zscoreParam(m, gene.sets), verbose=FALSE)
    Mz <- gsva(zscoreParam(M, gene.sets), verbose=FALSE)
    checkEqualsNumeric(mz, Mz)
    
    ms <- gsva(ssgseaParam(m, gene.sets), verbose=FALSE)
    Ms <- gsva(ssgseaParam(M, gene.sets), verbose=FALSE)
    checkEqualsNumeric(ms, Ms)
}

test_sparse_ecdfvals <- function() {
    message("Running unit tests for sparse ECDF values calculations.")

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

    suppressPackageStartupMessages(library(Matrix))
    suppressPackageStartupMessages(library(SparseArray))
    n <- 100
    p <- 100
    z <- numeric(p * n)
    nnz <- ceiling(0.05 * p * n) ## 5% nonzero values
    z[sample(1:(p*n), size=nnz, replace=FALSE)] <- rnorm(nnz)
    zz <- matrix(z, nrow=p, ncol=n)
    zzs <- Matrix(zz, sparse=TRUE)
    res_R_dense <- ecdfvals_dense(zz)
    res_C_dense_to_dense <- GSVA:::.ecdfvals_dense_to_dense(zz, FALSE)
    checkEqualsNumeric(res_R_dense, res_C_dense_to_dense)

    res_C_sparse_to_dense <- GSVA:::.ecdfvals_sparse_to_dense(zzs, FALSE)
    checkEqualsNumeric(res_R_dense, res_C_sparse_to_dense)

    res_R_sparse_to_sparse <- ecdfvals_sparse_to_sparse(zzs)
    res_C_sparse_to_sparse <- GSVA:::.ecdfvals_sparse_to_sparse(zzs, FALSE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparse_to_sparse)

    zzs <- SparseArray(zzs)
    res_C_svt_to_dense <- GSVA:::.ecdfvals_svt_to_dense(zzs, FALSE)
    checkEqualsNumeric(res_R_dense, res_C_svt_to_dense)
    res_C_svt_to_sparse <- GSVA:::.ecdfvals_svt_to_sparse(zzs, FALSE)
    checkEqualsNumeric(res_C_sparse_to_sparse, res_C_svt_to_sparse)
    res_C_svt_to_svt <- GSVA:::.ecdfvals_svt_to_svt(zzs, FALSE)
    checkEqualsNumeric(SparseArray(res_C_sparse_to_sparse), res_C_svt_to_svt)
}

test_sparse_kcdfvals <- function() {
    message("Running unit tests for sparse KCDF values calculations.")

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

    suppressPackageStartupMessages(library(Matrix))
    suppressPackageStartupMessages(library(SparseArray))
    n <- 100
    p <- 100
    z <- numeric(p * n)
    nnz <- ceiling(0.05 * p * n) ## 5% nonzero values
    z[sample(1:(p*n), size=nnz, replace=FALSE)] <- rnorm(nnz)
    zz <- matrix(z, nrow=p, ncol=n)
    zzs <- Matrix(zz, sparse=TRUE)

    res_R_sparse_to_dense <- kcdfegaussianvals_sparse_to_dense(zzs)
    res_C_sparse_to_dense <- GSVA:::.kcdfvals_sparse_to_dense(zzs, TRUE, FALSE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_sparse_to_dense, tolerance=0.001)

    res_R_sparse_to_sparse <- kcdfgaussianvals_sparse_to_sparse(zzs)
    res_C_sparse_to_sparse <- GSVA:::.kcdfvals_sparse_to_sparse(zzs, TRUE, FALSE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparse_to_sparse, tolerance=0.001)

    zzs <- SparseArray(zzs)
    res_C_svt_to_dense <- GSVA:::.kcdfvals_svt_to_dense(zzs, TRUE, FALSE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_svt_to_dense, tolerance=0.001)
    res_C_svt_to_svt <- GSVA:::.kcdfvals_svt_to_svt(zzs, TRUE, FALSE)
    checkEqualsNumeric(SparseArray(res_C_sparse_to_sparse), res_C_svt_to_svt)
}
