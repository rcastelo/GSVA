test_ranksserialization <- function() {

    message("Running unit tests for ranks serialization")

    suppressPackageStartupMessages({
        library(Matrix)
        library(GSEABase)
        library(SummarizedExperiment)
    })

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
 
    ## check out that 'saveHDF5GSVA()' throws an error when the input
    ## object is not one of the classes of the union of classes defined
    ## by 'GsvaExprData'
    rnormdir <- tempfile()
    checkException(saveHDF5GSVA(gsvapar, rnormdir))

    ## save the GSVA row-normalized expression values to disk
    savedrnormdir <- saveHDF5GSVA(gsvarownorm, rnormdir)
    checkTrue(file.exists(savedrnormdir) && identical(savedrnormdir, rnormdir))
                                               
    ## save the GSVA rank values to disk
    ranksdir <- tempfile()
    savedranksdir <- saveHDF5GSVA(gsvacolranks, ranksdir)
    checkTrue(file.exists(savedranksdir) && identical(savedranksdir, ranksdir))
                                               
    ## load the GSVA row-normalized values from disk               
    loaded_gsvarownorm <- loadHDF5GSVA(savedrnormdir)      
                                               
    ## load the GSVA ranks from disk               
    loaded_gsvacolranks <- loadHDF5GSVA(savedranksdir)      
                                               
    ## check that the loaded row-normalized values provide the
    ## same ranks as the ones calculated from the original values
    gsvacolranks_from_loaded_gsvarownorm <- gsvaColRanks(loaded_gsvarownorm, verbose=FALSE)
    checkEqualsNumeric(assay(gsvacolranks, "gsvaranks"),
		       assay(gsvacolranks_from_loaded_gsvarownorm, "gsvaranks"))

    ## check that the loaded ranks provide the
    ## same scores as the original ranks
    es_from_loaded_gsvaranks <- gsvaColScores(loaded_gsvacolranks, verbose=FALSE)
    checkTrue(identical(es, es_from_loaded_gsvaranks))

    ## check it again saving and loading the ranks stored in a
    ## non-'SummarizedExperiment' object
    gsvapar <- gsvaParam(cnt, gsets, verbose=FALSE)
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    savedrnormdir <- saveHDF5GSVA(gsvarownorm, rnormdir, replace=TRUE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    savedranksdir <- saveHDF5GSVA(gsvacolranks, ranksdir, replace=TRUE)

    loaded_gsvarownorm <- loadHDF5GSVA(savedrnormdir)
    gsvacolranks_from_loaded_gsvarownorm <- gsvaColRanks(loaded_gsvarownorm, verbose=FALSE)
    checkEqualsNumeric(gsvacolranks, gsvacolranks_from_loaded_gsvarownorm)

    loaded_gsvacolranks <- loadHDF5GSVA(savedranksdir)
    es_from_loaded_gsvaranks <- gsvaColScores(loaded_gsvacolranks, verbose=FALSE)
    checkEqualsNumeric(assay(es), es_from_loaded_gsvaranks)
}
