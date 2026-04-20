test_ranksserialization <- function() {

    message("Running unit tests for ranks serialization")

    suppressPackageStartupMessages(library(GSEABase))
    suppressPackageStartupMessages(library(SummarizedExperiment))

    p <- 10 ## number of genes
    n <- 30 ## number of samples
    nGrp1 <- 15 ## number of samples in group 1
    nGrp2 <- n - nGrp1 ## number of samples in group 2

    ## consider three disjoint gene sets
    gsets <- list(gset1=paste0("g", 1:3),
                  gset2=paste0("g", 4:6),
                  gset3=paste0("g", 7:10))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    se <- SummarizedExperiment(assays=list(counts=y))
    gsvaAnnotation(se) <- SymbolIdentifier("org.Hs.eg.db")

    ## build GSVA parameter object
    gsvapar <- gsvaParam(se, gsets)

    ## calculate GSVA ranks
    gsvarankspar <- gsvaRanks(gsvapar)

    ## calculate GSVA scores
    es <- gsvaScores(gsvarankspar)
 
    ## save the GSVA ranks to disk
    dir <- tempfile()     
    checkException(saveHDF5GSVAranks(gsvapar, dir))
    saveHDF5GSVAranks(gsvarankspar, dir)
                                               
    ## load the GSVA ranks from disk               
    loaded_gsvarankspar <- loadHDF5GSVAranks(dir)      
                                               
    ## check that the loaded ranks provide the
    ## same scores as the original ranks
    loaded_es <- gsvaScores(loaded_gsvarankspar)

    checkTrue(identical(es, loaded_es))

    ## tweak the input object to test the error handling of the loading function
    assay(gsvarankspar@exprData, "gsvaranks") <- NULL
    checkException(saveHDF5GSVAranks(gsvarankspar, dir))
}
