test_genesets <- function() {
    message("Running unit tests for input gene sets.")

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
    
    ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
    y[gsets$set1, (nGrp1+1):n] <- y[gsets$set1, (nGrp1+1):n] + 2

    ## estimate GSVA enrichment scores with gene sets input as a list
    gsvapar <- gsvaParam(y, gsets, verbose=FALSE)
    es.mat <- gsva(gsvapar, verbose=FALSE)

    ## check that gene sets are correctly propagated to the parameter object
    checkIdentical(geneSets(gsvapar), gsets)
    ## check that gene sets are correctly propagated to the output enrichment score matrix
    checkIdentical(geneSets(es.mat), gsets)

    ## convert input gene sets into a GeneSetCollection object
    gsc <- geneIdsToGeneSetCollection(gsets)

    ## estimate GSVA enrichment scores with gene sets input as a GeneSetCollection object
    es.mat2 <- gsva(gsvaParam(y, gsc), verbose=FALSE)

    checkTrue(identical(es.mat, es.mat2))
}

test_geneSetDeDuplication <- function() {
    message("Running unit tests for deduplicating gene sets.")

    gsets <- list(gs1=LETTERS[1:3], gs2=LETTERS[4:6], gs2=LETTERS[5:8])

    checkTrue(length(deduplicateGeneSets(gsets)) == 2)
    checkTrue(length(deduplicateGeneSets(gsets, deduplUse="drop")) == 1)
    checkTrue(all(lengths(deduplicateGeneSets(gsets, deduplUse="union")) == c(3, 5)))
    checkTrue(all(lengths(deduplicateGeneSets(gsets, deduplUse="smallest")) == c(3, 3)))
    checkTrue(all(lengths(deduplicateGeneSets(gsets, deduplUse="largest")) == c(3, 4)))
}

test_readGMT <- function() {
    message("Running unit tests for reading GMT files with readGMT()")

    library(GSVAdata)
    fname <- system.file("extdata", "c2.subsetdups.v7.5.symbols.gmt.gz",
                        package="GSVAdata")
    warn <- FALSE
    suppressWarnings(c2.dupgenesets <- readGMT(fname, deduplUse="union", valueType="list"))
    checkTrue(!any(duplicated(names(c2.dupgenesets))))
}
