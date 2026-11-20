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

    ## the direct ECDF ('kcdf="none"') on dense integer values, e.g., counts,
    ## with and without missing values, gives the same values as on doubles
    set.seed(123)
    yi <- matrix(rpois(n*p, 5), nrow=p, ncol=n, dimnames=dimnames(y))
    yd <- yi
    storage.mode(yd) <- "double"
    for (use in c("everything", "na.rm")) {
        if (use == "na.rm") {
            yi[sample(length(yi), 10)] <- NA
            yd[is.na(yi)] <- NA
        }
        gsvarownorm_int <- suppressWarnings(
            gsvaRowNorm(gsvaParam(yi, gsets, rowNorm="ecdf", kcdf="none",
                                  use=use), verbose=FALSE))
        gsvarownorm_dbl <- suppressWarnings(
            gsvaRowNorm(gsvaParam(yd, gsets, rowNorm="ecdf", kcdf="none",
                                  use=use), verbose=FALSE))
        checkEqualsNumeric(as.vector(gsvarownorm_int),
                           as.vector(gsvarownorm_dbl))
    }
}

## rowNorm="auto" selects CLR for sparse input processed with sparse=TRUE whose
## nonzero values are positive, including integer counts, and ECDF otherwise,
## informing of the selection, and stores the selected method
test_rowNormAuto <- function() {
    message("Running unit tests for rowNorm='auto'")

    suppressPackageStartupMessages({
        library(Matrix)
        library(SparseArray)
        library(HDF5Array)
        library(cli)
    })

    set.seed(123)
    p <- 40; n <- 30
    cnt <- matrix(rpois(p * n, 0.8), p, n,
                  dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    gsets <- list(set1=paste0("g", 1:5), set2=paste0("g", 6:12))
    rownorm <- function(x, ...)
        GSVA:::.get_rowNorm(gsvaParam(x, gsets, verbose=FALSE, ...))

    ## dense input
    checkIdentical(rownorm(cnt + 1), "ecdf")
    out <- cli_fmt(gsvaParam(cnt + 1, gsets))
    checkTrue(any(grepl("rowNorm='auto' selects 'ecdf'", out)))

    ## sparse input in main memory with positive nonzero values, also integer
    ## counts, unless processed as dense
    S <- as(log1p(cnt), "dgCMatrix")
    checkIdentical(rownorm(S), "clr")
    checkIdentical(rownorm(as(cnt, "dgCMatrix")), "clr")
    checkIdentical(rownorm(as(S, "SVT_SparseMatrix")), "clr")
    checkIdentical(rownorm(S, sparse=FALSE), "ecdf")
    out <- cli_fmt(gsvaParam(S, gsets))
    checkTrue(any(grepl("rowNorm='auto' selects 'clr'", out)))

    ## a nonpositive nonzero value selects ECDF, while the method given
    ## explicitly is kept
    N <- S
    N[2, 3] <- -1
    checkIdentical(rownorm(N), "ecdf")
    checkIdentical(rownorm(S, rowNorm="ecdf"), "ecdf")
    checkIdentical(rownorm(S, rowNorm="none"), "none")

    ## the smallest nonzero value, from all of them or from a sample of the
    ## data on disk, whose number of nonzero values its seeds store
    nzc <- GSVA:::.estimate_nzcount(S, NA, FALSE, nzmin=TRUE)
    checkEqualsNumeric(attr(nzc, "nzmin"), min(S@x))
    checkTrue(attr(nzc, "nzminexact"))
    H <- SummarizedExperiment::SummarizedExperiment(
             assays=list(counts=writeTENxMatrix(S)))
    nzc <- GSVA:::.estimate_nzcount(H, "counts", FALSE, nzmin=TRUE)
    checkEqualsNumeric(as.vector(nzc), nnzero(S))
    checkTrue(attr(nzc, "nzmin") > 0)
    checkTrue(!attr(nzc, "nzminexact"))
    checkIdentical(rownorm(H), "clr")

    ## CLR on nonpositive nonzero values gives an error suggesting ECDF
    err <- tryCatch(gsvaRowNorm(gsvaParam(N, gsets, rowNorm="clr",
                                          verbose=FALSE), verbose=FALSE),
                    error=conditionMessage)
    checkTrue(is.character(err) && grepl("rowNorm=\"ecdf\"", err, fixed=TRUE))
}
