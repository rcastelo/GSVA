test_gsvaRanks <- function() {
    message("Running unit tests for GSVA ranks.")

    p <- 10 ## number of genes
    n <- 30 ## number of samples
    nGrp1 <- 15 ## number of samples in group 1
    nGrp2 <- n - nGrp1 ## number of samples in group 2

    ## consider three disjoint gene sets
    geneSets <- list(set1=paste("g", 1:3, sep=""),
                     set2=paste("g", 4:6, sep=""),
                     set3=paste("g", 7:10, sep=""))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(y, geneSets)

    ## calculate GSVA scores in one step
    gsva_es1 <- gsva(gsvapar, verbose=FALSE)

    ## calculate GSVA scores in two steps
    ## first calculate GSVA ranks
    gsvarankspar <- gsvaRanks(gsvapar, verbose=FALSE)

    ## second calculate GSVA scores using GSVA ranks
    gsva_es2 <- gsvaScores(gsvarankspar, verbose=FALSE)

    ## both approaches to calculate GSVA scores must give
    ## the same result with the same input gene sets
    checkEqualsNumeric(gsva_es1, gsva_es2)
    checkTrue(all.equal(gsva_es1, gsva_es2))
}
