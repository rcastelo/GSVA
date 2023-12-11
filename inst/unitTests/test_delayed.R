test_delayed <- function() {
    DEACTIVATED("HDF5Array is not yet supported by the new API.")

    message("Running unit tests for DelayedMatrix.")
    
    ## set.seed(1001)
    m <- matrix(runif(10000), 1000, 100)
    rownames(m) <- paste0("gene_", 1:1000)
    colnames(m) <- paste0("cell_", 1:100)
    
    gene.sets <- list("my_list1"= paste0("gene_", 1:100),
                      "my_list2"= paste0("gene_", 401:500))
    
    h5 <- as(m, "HDF5Array")
    
    x1 <- gsva(m, gene.sets, method="plage", verbose=FALSE)
    x2 <- gsva(h5, gene.sets, method="plage", verbose=FALSE)
    
    ## SVD sign is arbitrary and implementation-dependent
    checkEqualsNumeric(abs(as.vector(t(x1))), abs(as.vector(t(x2))))
    
    x1 <- gsva(m, gene.sets, method="zscore", verbose=FALSE)
    x2 <- gsva(h5, gene.sets, method="zscore", verbose=FALSE)
    
    checkEqualsNumeric(as.vector(t(x1)), as.vector(t(x2)))
    
    x1 <- gsva(m, gene.sets, method="ssgsea", verbose=FALSE)
    x2 <- gsva(h5, gene.sets, method="ssgsea", verbose=FALSE)
    
    checkEqualsNumeric(as.vector(t(x1)), as.vector(t(x2)))
}
