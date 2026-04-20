test_ranksserialization <- function() {
    p <- 10 ## number of genes
    n <- 30 ## number of samples
    nGrp1 <- 15 ## number of samples in group 1
    nGrp2 <- n - nGrp1 ## number of samples in group 2

    ## consider three disjoint gene sets
    geneSets <- list(gset1=paste0("g", 1:3),
                     gset2=paste0("g", 4:6),
                     gset3=paste0("g", 7:10))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(y, geneSets)

    ## calculate GSVA ranks
    gsvarankspar <- gsvaRanks(gsvapar)

    ## calculate GSVA scores
    es <- gsvaScores(gsvarankspar)
 
    ## save the GSVA ranks to disk
    dir <- tempfile()     
    saveHDF5GSVAranks(gsvarankspar, dir)               
                                               
    ## load the GSVA ranks from disk               
    loaded_gsvarankspar <- loadHDF5GSVAranks(dir)      
                                               
    ## check that the loaded ranks provide the
    ## same scores as the original ranks
    loaded_es <- gsvaScores(loaded_gsvarankspar)

    checkTrue(identical(es, loaded_es))
}
