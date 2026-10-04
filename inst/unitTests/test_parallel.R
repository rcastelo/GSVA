test_parallel <- function() {
    message("Running unit tests for parallel execution.")

    suppressPackageStartupMessages({
        library(BiocParallel)
        library(Matrix)
    })

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

    M[1, 2] <- NA
    checkException(gsva(gsvaParam(M, gsets, kcdf="Gaussian", verbose=FALSE),
			verbose=FALSE, BPPARAM=MulticoreParam(workers=2)))
}

test_parallel_bpparam_unchanged <- function() {

    message("Running unit tests for parallel execution keeping BPPARAM unchanged")

    suppressPackageStartupMessages({
        library(BiocParallel)
        library(SpatialExperiment)
    })

    set.seed(123)
    p <- 200
    n <- 150
    y <- matrix(rnorm(p * n), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    gsets <- list(gs1=paste0("g", 1:20), gs2=paste0("g", 21:60),
                  gs3=paste0("g", 61:120))

    ## with verbose=TRUE, the parallel execution sets progressbar=TRUE in
    ## BPPARAM, which is a reference class object, and this setting is
    ## restored when the calculations finish, also when spatCor() loops over
    ## several samples
    xy <- cbind(x=runif(n), y=runif(n))
    rownames(xy) <- colnames(y)
    spe <- SpatialExperiment(assays=list(es=y[1:3, ]),
                             sample_id=rep(c("A", "B"), each=n/2),
                             spatialCoords=xy)
    bpparam <- if (.Platform$OS.type == "unix")
                   MulticoreParam(workers=2, progressbar=FALSE)
               else
                   SnowParam(workers=2, progressbar=FALSE)
    suppressMessages(spatCor(spe, assay="es", verbose=TRUE, BPPARAM=bpparam))
    checkTrue(!bpprogressbar(bpparam))

    ## and in gsvaMap(), where BTPARAM gets the value of the verbose argument
    btparam <- BatchtoolsParam(workers=1, progressbar=FALSE,
                               resources=list(ncpus=1, memory="1G"))
    suppressMessages(gsvaMap(gsvaRowNorm, gsvaParam(y, gsets), verbose=TRUE,
                             BTPARAM=btparam))
    checkTrue(!bpprogressbar(btparam))
}
