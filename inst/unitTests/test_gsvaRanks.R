test_gsvaRanks <- function() {
    message("Running unit tests for GSVA ranks.")

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

    ## build GSVA parameter object
    gsvapar <- gsvaParam(y, gsets)

    ## calculate GSVA scores in one step
    gsva_es1 <- gsva(gsvapar, verbose=FALSE)

    ## calculate GSVA scores in three steps
    ## first calculate row-normalized expression values
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## second calculate GSVA column ranks
    gsvaranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)

    ## third calculate GSVA scores from column ranks
    gsva_es2 <- gsvaColScores(gsvaranks, verbose=FALSE)

    ## both approaches to calculate GSVA scores must give
    ## the same result with the same input gene sets
    checkEqualsNumeric(gsva_es1, gsva_es2)

    ## check that gsvaEnrichment() works
    gsvaenrich <- gsvaEnrichment(gsvaranks, plot="no")
    checkEqualsNumeric(gsva_es1[1, 1], gsvaenrich$score)
    gsvaenrich2 <- gsvaEnrichment(gsvaranks,
                                  geneSet=c("g1", "g4", "g7"),
                                  plot="no")
    checkTrue(!is.na(gsvaenrich2$score))

    ## test the ggplotting from gsvaEnrichment()
    ggp <- gsvaEnrichment(gsvaranks, plot="ggplot")
    checkTrue(is(ggp, "ggplot"))
    checkTrue(identical(gsvaenrich$stats, ggp@data))

    ## calculate again row-normalized expression values,
    ## but this time in chunks
    gsvarownorm1 <- gsvaRowNorm(gsvapar, first=1, last=5, verbose=FALSE)
    gsvarownorm2 <- gsvaRowNorm(gsvapar, first=6, last=10, verbose=FALSE)
    checkEqualsNumeric(gsvarownorm, rbind(gsvarownorm1, gsvarownorm2))
    checkException(gsvaRowNorm(gsvapar, first=10, last=6, verbose=FALSE))
    checkException(gsvaRowNorm(gsvapar, first=11, last=20, verbose=FALSE))

    ## calculate again GSVA column ranks, but this time in chunks
    gsvaranks1 <- gsvaColRanks(gsvarownorm, first=1, last=10, verbose=FALSE)
    gsvaranks2 <- gsvaColRanks(gsvarownorm, first=11, last=20, verbose=FALSE)
    gsvaranks3 <- gsvaColRanks(gsvarownorm, first=21, last=30, verbose=FALSE)
    checkEqualsNumeric(gsvaranks, cbind(gsvaranks1, gsvaranks2, gsvaranks3))

    ## calculate again GSVA scores from column ranks, but this time in chunks
    gsva_es_c1 <- gsvaColScores(gsvaranks, first=1, last=10, verbose=FALSE)
    gsva_es_c2 <- gsvaColScores(gsvaranks, first=11, last=20, verbose=FALSE)
    gsva_es_c3 <- gsvaColScores(gsvaranks, first=21, last=30, verbose=FALSE)
    checkEqualsNumeric(gsva_es1, cbind(gsva_es_c1, gsva_es_c2, gsva_es_c3))
}
