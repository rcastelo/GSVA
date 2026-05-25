test_ranksserialization <- function() {

    message("Running unit tests for ranks serialization")

    suppressPackageStartupMessages(library(Matrix))
    suppressPackageStartupMessages(library(GSEABase))
    suppressPackageStartupMessages(library(SummarizedExperiment))

    p <- 10 ## number of genes
    n <- 30 ## number of samples

    ## consider three disjoint gene sets
    gsets <- list(gset1=paste0("g", 1:3),
                  gset2=paste0("g", 4:6),
                  gset3=paste0("g", 7:10))

    ## build a random sparse count matrix with 85% sparsity
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.15)) ## 85% sparsity
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- matrix(cnt, nrow=p, ncol=n,
                  dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    cnt <- Matrix(cnt, sparse=TRUE)

    se <- SummarizedExperiment(assays=list(counts=cnt))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(se, gsets, verbose=FALSE)

    ## calculate row-normalized expression values
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## calculate GSVA column ranks
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)

    ## calculate GSVA scores
    es <- gsvaColScores(gsvacolranks, verbose=FALSE)
 
    ## save the GSVA ranks to disk
    dir <- tempfile()     
    checkException(saveHDF5GSVAranks(gsvapar, dir))
    saveHDF5GSVAranks(gsvacolranks, dir)
                                               
    ## load the GSVA ranks from disk               
    loaded_gsvacolranks <- loadHDF5GSVAranks(dir)      
                                               
    ## check that the loaded ranks provide the
    ## same scores as the original ranks
    loaded_es <- gsvaColScores(loaded_gsvacolranks, verbose=FALSE)

    checkTrue(identical(es, loaded_es))

    ## check it again saving and loading the ranks stored in a
    ## non-'SummarizedExperiment' object
    gsvapar <- gsvaParam(cnt, gsets, verbose=FALSE)
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    saveHDF5GSVAranks(gsvacolranks, dir, replace=TRUE)
    loaded_gsvacolranks <- loadHDF5GSVAranks(dir)      
    loaded_es <- gsvaColScores(loaded_gsvacolranks, verbose=FALSE)
    checkEqualsNumeric(assay(es), loaded_es)
}
