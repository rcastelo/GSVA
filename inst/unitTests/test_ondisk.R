test_ondisk <- function() {

    message("Running unit tests for ondisk input.")

    suppressPackageStartupMessages({
        library(Matrix)
        library(HDF5Array)
        library(SparseArray)
    })

    p <- 50 ## number of genes
    n <- 100 ## number of samples
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
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    M <- Matrix(y, sparse=TRUE)
    H5 <- as(M, "HDF5Matrix")

    ## estimate GSVA enrichment scores with and without HDF5 input and check that they are identical
    es_noh5 <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE)
    es_h5 <- gsva(gsvaParam(H5, gsets, verbose=FALSE), verbose=FALSE)
    checkIdentical(es_noh5, es_h5)

    ## estimate GSVA enrichment scores with HDF5 input and output and check that they are identical
    es_h5ondisk <- gsva(gsvaParam(H5, gsets, ondisk="yes", verbose=FALSE), verbose=FALSE)
    es_h5ondiskmat <- as.matrix(es_h5ondisk)
    checkEqualsNumeric(es_noh5, es_h5ondiskmat)

    ## estimate ssGSEA enrichment scores with and without HDF5 input and check that they are identical
    es_noh5 <- gsva(ssgseaParam(M, gsets, verbose=FALSE), verbose=FALSE)
    es_h5 <- gsva(ssgseaParam(H5, gsets, verbose=FALSE), verbose=FALSE)
    checkIdentical(es_noh5, es_h5)

    ## estimate ssGSEA enrichment scores with HDF5 input and output and check that they are identical
    es_h5ondisk <- gsva(ssgseaParam(H5, gsets, ondisk="yes", verbose=FALSE),
			verbose=FALSE)
    es_h5ondiskmat <- as.matrix(es_h5ondisk)
    checkEqualsNumeric(es_noh5, es_h5ondiskmat)

    ## estimate Z-scores enrichment scores with and without HDF5 input and check that they are identical
    es_noh5 <- gsva(zscoreParam(M, gsets, verbose=FALSE), verbose=TRUE)
    es_h5 <- gsva(zscoreParam(H5, gsets, verbose=FALSE), verbose=FALSE)
    ## not identical also due to the rowSds() vs sd() differences
    checkEqualsNumeric(es_noh5, es_h5)

    ## estimate Z-scores enrichment scores with HDF5 input and output and check that they are identical
    es_h5ondisk <- gsva(zscoreParam(H5, gsets, ondisk="yes", verbose=FALSE),
			verbose=FALSE)
    es_h5ondiskmat <- as.matrix(es_h5ondisk)
    ## not identical also due to the rowSds() vs sd() differences
    checkEqualsNumeric(es_noh5, es_h5ondiskmat)

    ## test the block processing of a small toy HDF5 input and output by
    ## setting a small block size and maximum available memory
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(1024)
    es_noh5 <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE)
    es_chunks <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=TRUE, maxmem="25K")
    checkEqualsNumeric(es_noh5, es_chunks)
    setAutoBlockSize(oldautoblocksize)
}
