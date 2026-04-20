test_genesets <- function() {
    message("Running unit tests for input gene sets")

    suppressPackageStartupMessages(library(GSEABase))
    suppressPackageStartupMessages(library(SummarizedExperiment))

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
    gsvaAnnotation(y) <- SymbolIdentifier("org.Hs.eg.db")
    se <- SummarizedExperiment(assays=list(logcounts=y))
    gsvaAnnotation(se) <- SymbolIdentifier("org.Hs.eg.db")
    
    ## check error is thrown when trying to access missing gene sets
    checkException(geneSets(y))

    ## estimate GSVA enrichment scores with gene sets input as a list
    gsvapar <- gsvaParam(y, gsets, verbose=FALSE)
    es.mat <- gsva(gsvapar, verbose=FALSE)

    ## check that gene sets are correctly propagated to the parameter object
    checkIdentical(geneSets(gsvapar), gsets)
    ## check that gene set sizes are correctly propagated to the parameter object
    checkIdentical(geneSetSizes(gsvapar), lengths(gsets))
    ## check that gene sets are correctly propagated to the output enrichment
    ## score matrix
    checkIdentical(geneSets(es.mat), gsets)
    ## check that gene set sizes are correctly propagated to the output
    ## enrichment score matrix
    checkIdentical(geneSetSizes(es.mat), geneSetSizes(gsvapar))

    ## convert input gene sets into a GeneSetCollection object
    gsc <- geneIdsToGeneSetCollection(gsets)

    ## check gene identifier metadata
    checkTrue(identical(gsvaAnnotation(gsc), SymbolIdentifier()))

    ## check filtering method
    checkTrue(identical(filterGeneSets(gsc), geneIds(gsc)))

    ## estimate GSVA enrichment scores with gene sets input as a GeneSetCollection object
    es.mat2 <- gsva(gsvaParam(y, gsc), verbose=FALSE)
    checkTrue(identical(es.mat, es.mat2))

    ## check that when input expression data has no rownames and gene sets
    ## are made out of indexes to the rows, the results do not change
    gsets <- lapply(gsets, function(x) as.numeric(sub("g", "", x)))
    rownames(y) <- NULL
    gsvapar <- gsvaParam(y, gsets, verbose=FALSE)
    es.mat3 <- gsva(gsvapar, verbose=FALSE)
    attr(es.mat, "geneSets") <- attr(es.mat3, "geneSets") <- NULL
    checkTrue(identical(es.mat, es.mat3))
    gsets$gset3 <- c(gsets$gset3, 11)
    gsvapar <- gsvaParam(y, gsets, verbose=FALSE)
    es.mat4 <- gsva(gsvapar, verbose=FALSE)
    attr(es.mat4, "geneSets") <- NULL
    checkTrue(identical(es.mat, es.mat4))
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

    checkException(readGMT(c("nonexistent_file1.gmt", "nonexistent_file2.gmt")))
    checkException(readGMT("nonexistent_file.gmt"))

    suppressPackageStartupMessages(library(GSVAdata))
    fname <- system.file("extdata", "c2.subsetdups.v7.5.symbols.gmt.gz",
                        package="GSVAdata")
    suppressWarnings(c2.dupgenesets <- readGMT(fname, deduplUse="union",
                                               valueType="list"))
    checkTrue(!any(duplicated(names(c2.dupgenesets))))
    suppressWarnings(c2.dupgenesets <- readGMT(fname, deduplUse="smallest",
                                               valueType="GeneSetCollection"))
    checkTrue(!any(duplicated(names(c2.dupgenesets))))
    suppressWarnings(c2.dupgenesets <- readGMT(fname, deduplUse="largest",
                                               valueType="list"))
    suppressWarnings(c2.dupgenesets <- readGMT(fname, deduplUse="drop",
                                               valueType="list"))
    checkTrue(!any(duplicated(names(c2.dupgenesets))))

    gsets <- c2.dupgenesets[1:2]
    fname <- tempfile()
    con <- file(fname, "w")
    writeLines(c(names(gsets)[1],
		 paste(names(gsets)[2], "desc2", paste(gsets[[2]], collapse="\t"), sep="\t")), con)
    close(con)
    checkException(gsets.read <- readGMT(fname, deduplUse="drop", valueType="list"))

    gsets[[1]][2] <- gsets[[1]][1]
    con <- file(fname, "w")
    writeLines(c(paste(names(gsets)[1], "desc1", paste(paste0("ENSG", gsets[[1]]), collapse="\t"), sep="\t"),
		 paste(names(gsets)[2], "desc2", paste(paste0("ENSG", gsets[[2]]), collapse="\t"), sep="\t")), con)
    close(con)
    library(cli)
    gsets.read <- readGMT(fname, deduplUse="drop", valueType="list")
    checkTrue(!any(duplicated(gsets.read[[1]])))
}
