test_parallel <- function() {
    message("Running unit tests for parallel execution.")

    library(BiocParallel)
    library(Matrix)

    p <- 150 ## number of genes
    n <- 150 ## number of samples
    m <- 10 ## number of gene sets

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

    ## estimate GSVA enrichment scores with and without parallel execution
    ## and check that they are identical
    es_serial <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE)
    es_parallel <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE,
			BPPARAM=MulticoreParam(workers=2))
    checkIdentical(es_serial, es_parallel)
}
