test_gsvaRowNorm <- function() {
    message("Running unit tests for GSVA row normalization")

    p <- 10 ## number of genes
    n <- 30 ## number of samples
    nGrp1 <- 15 ## number of samples in group 1
    nGrp2 <- n - nGrp1 ## number of samples in group 2

    ## consider three disjoint gene sets
    gsets <- list(set1=paste("g", 1:3, sep=""),
                  set2=paste("g", 4:6, sep=""),
                  set3=paste("g", 7:10, sep=""))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    y <- abs(y) ## make all values positive to avoid issues with log-transformations

    ## build GSVA parameter object, specifying row normalization by ECDF
    gsvapar <- gsvaParam(y, gsets, rowNorm="ecdf", verbose=FALSE)

    ## calculate row-normalized expression values with ECDF
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## the current default of dense matrices and ECDF, calculates the
    ## log of the odds of the normalized ECDF values, so we need to
    ## transform those back to ECDF values first before we check that
    ## they are in the range [0, 1]
    z <- exp(t(gsvarownorm)) / (1 + exp(t(gsvarownorm)))
    checkTrue(all(z >= 0 & z <= 1))

    ## build GSVA parameter object, specifying row normalization by CLR
    gsvapar <- gsvaParam(y, gsets, rowNorm="clr", verbose=FALSE)

    ## calculate row-normalized expression values with CLR
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## calculate geometric mean of the normalized rows
    gsvarownorm_gmean <- exp(rowMeans(log(gsvarownorm)))

    ## check that the geometric mean of the normalized rows is 1
    checkEqualsNumeric(gsvarownorm_gmean, rep(1, p))

    ## repeat now the tests with sparse matrices
    suppressPackageStartupMessages(library(Matrix))
    set.seed(123)
    s <- ceiling(0.15 * p * n)
    sam <- sample(1:(p * n), size=s, replace=FALSE)
    x <- numeric(p * n)
    x[sam] <- runif(s)
    m <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    M <- Matrix(m, sparse=TRUE)
    
    ## build GSVA parameter object, specifying row normalization by ECDF
    gsvapar <- gsvaParam(M, gsets, rowNorm="ecdf", verbose=FALSE)

    ## calculate row-normalized expression values with ECDF
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## the sparse regime directly returns the ECDF values
    checkTrue(all(gsvarownorm >= 0 & gsvarownorm <= 1))

    ## build GSVA parameter object, specifying row normalization by CLR
    gsvapar <- gsvaParam(M, gsets, rowNorm="clr", verbose=FALSE)

    ## calculate row-normalized expression values with CLR
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## calculate geometric mean of nonzero values in the normalized rows
    gsvarownorm <- as.matrix(gsvarownorm) ## assignment below not supported for sparse matrices
    gsvarownorm[gsvarownorm == 0] <- NA
    gsvarownorm_gmean <- exp(rowMeans(log(gsvarownorm), na.rm=TRUE))

    ## check that the geometric mean of the normalized rows is 1
    checkEqualsNumeric(gsvarownorm_gmean, rep(1, p))
}
