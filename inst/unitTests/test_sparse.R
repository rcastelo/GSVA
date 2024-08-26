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

text_sparse_ecdfvals <- function() {
    message("Running unit tests for sparse ECDF values calculations.")

    ecdfvals_dense <- function(X) t(apply(X, 1, function(rx) ecdf(rx)(rx)))
    ecdfvals_sparse_to_sparse <- function(X) {
        for (i in 1:nrow(X)) {
            rx <- X[i, , drop=FALSE]
            vals <- unique(sort(rx@x))
            mt <- match(rx@x, vals)
            tab <- tabulate(mt, nbins=length(vals))
            ecdfvals <- cumsum(tab) / nnzero(rx)
            X[i, rep.int(1:n, diff(rx@p))] <- ecdfvals[mt]
        }
        X
    }

    n <- 100
    p <- 100
    z <- numeric(p * n)
    nnz <- ceiling(0.05 * p * n) ## 5% nonzero values
    z[sample(1:(p*n), size=nnz, replace=FALSE)] <- rnorm(nnz)
    zz <- matrix(z, nrow=p, ncol=n)
    zzs <- Matrix(zz, sparse=TRUE)
    res_R_dense <- ecdfvals_dense(zz)
    res_C_dense_to_dense <- GSVA:::.ecdfvals_dense_to_dense(zz)
    checkEqualsNumeric(res_R_dense, res_C_dense_to_dense)

    res_C_sparse_to_dense <- GSVA:::.ecdfvals_sparse_to_dense(zzs)
    checkEqualsNumeric(res_R_dense, res_C_sparse_to_dense)

    res_R_sparse_to_sparse <- ecdfvals_sparse_to_sparse(zzs)
    res_C_sparse_to_sparse <- GSVA:::.ecdfvals_sparse_to_sparse(zzs)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparse_to_sparse)
}
