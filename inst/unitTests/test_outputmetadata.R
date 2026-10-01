## metadata recorded in GSVA output objects: 'gsvaVersion', the version of GSVA
## that produced them, and 'ranksNrow', the number of rows on which column
## ranks were calculated, which is checked to prevent calculating GSVA scores
## from ranks after removing rows
getoutmd <- function(x, name) {
    if (is(x, "SummarizedExperiment"))
        metadata(x)[[name]]
    else
        attr(x, name, exact=TRUE)
}

test_outputmetadata <- function() {

    message("Running unit tests for metadata in GSVA output objects")

    suppressPackageStartupMessages({
        library(Matrix)
        library(Biobase)
        library(SummarizedExperiment)
        library(SingleCellExperiment)
    })

    p <- 40 ## number of genes
    n <- 60 ## number of samples
    gsets <- list(gset1=paste0("g", 1:10),
                  gset2=paste0("g", 11:25),
                  gset3=paste0("g", 26:40))
    set.seed(123)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.3))
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- Matrix(matrix(cnt, nrow=p, ncol=n, dimnames=dimnames(y)),
                  sparse=TRUE)

    gsvaversion <- as.character(packageVersion("GSVA"))

    getvals <- function(x, a) {
        if (is(x, "SummarizedExperiment"))
            x <- assay(x, a)
        else if (is(x, "ExpressionSet"))
            x <- exprs(x)
        as.matrix(x)
    }

    inputs <- list(y, SummarizedExperiment(assays=list(exprs=y)),
                   SingleCellExperiment(assays=list(counts=cnt)),
                   ExpressionSet(y))
    for (input in inputs) {
        gsvapar <- gsvaParam(input, gsets, verbose=FALSE)
        gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
        gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
        es <- gsvaColScores(gsvacolranks, verbose=FALSE)
        es2 <- gsva(gsvapar, verbose=FALSE)

        ## all outputs record the version of GSVA, and only the column
        ## ranks record the number of rows on which they were calculated
        for (x in list(gsvarownorm, gsvacolranks, es, es2))
            checkIdentical(gsvaversion, getoutmd(x, "gsvaVersion"))
        checkIdentical(unname(nrow(gsvacolranks)),
                       getoutmd(gsvacolranks, "ranksNrow"))
        checkTrue(is.null(getoutmd(gsvarownorm, "ranksNrow")))
        checkTrue(is.null(getoutmd(es, "ranksNrow")))
        checkTrue(is.null(getoutmd(es2, "ranksNrow")))

        if (is(input, "SummarizedExperiment")) {
            ## removing rows after calculating the ranks gives an error,
            ## instead of calculating GSVA scores from invalid ranks
            subranks <- gsvacolranks[1:20, ]
            msg <- tryCatch(gsvaColScores(subranks, verbose=FALSE),
                            error=conditionMessage)
            checkTrue(is.character(msg) && grepl("calculated on", msg))
            checkException(gsvaEnrichment(subranks, column=1, geneSet=1,
                                          plot="no"), silent=TRUE)

            ## removing columns does not affect the ranks
            es3 <- gsvaColScores(gsvacolranks[, 1:10], verbose=FALSE)
            checkEqualsNumeric(getvals(es, "es")[, 1:10], getvals(es3, "es"))

            ## ranks without the number of rows, as produced by previous
            ## versions of GSVA, are not checked
            oldranks <- gsvacolranks
            metadata(oldranks)$ranksNrow <- NULL
            es4 <- gsvaColScores(oldranks, verbose=FALSE)
            checkEqualsNumeric(getvals(es, "es"), getvals(es4, "es"))
        }

        ## the metadata is kept when saving and loading the output
        if (!is(input, "ExpressionSet")) {
            rankspath <- saveHDF5GSVA(gsvacolranks, tempfile())
            loadedranks <- loadHDF5GSVA(rankspath)
            checkIdentical(gsvaversion, getoutmd(loadedranks, "gsvaVersion"))
            checkIdentical(unname(nrow(gsvacolranks)), getoutmd(loadedranks,
                                                        "ranksNrow"))
            unlink(rankspath, recursive=TRUE)
        }
    }
}

test_outputmetadatamapreduce <- function() {

    message("Running unit tests for metadata in GSVA map-reduce output")

    suppressPackageStartupMessages({
        library(DelayedArray)
        library(SummarizedExperiment)
        library(BiocParallel)
    })

    p <- 40 ## number of genes
    n <- 90 ## number of samples
    gsets <- list(gset1=paste0("g", 1:10),
                  gset2=paste0("g", 11:25),
                  gset3=paste0("g", 26:40))
    set.seed(123)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))

    gsvaversion <- as.character(packageVersion("GSVA"))
    ## gsvaMap() is run in this R process with one worker, which is much
    ## faster than starting batchtools jobs, and small blocks, together with a
    ## finite maximum memory, split the input into several chunks
    btpar <- BatchtoolsParam(workers=1, resources=list(ncpus=1, memory="1K"))
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(p * 8 * 30) ## blocks of 30 columns
    on.exit(setAutoBlockSize(oldautoblocksize), add=TRUE)

    for (input in list(y, SummarizedExperiment(assays=list(exprs=y)))) {
        gsvapar <- gsvaParam(input, gsets, verbose=FALSE)
        gsvarownorm <- gsvaReduce(gsvaMap(gsvaRowNorm, gsvapar, verbose=FALSE,
                                          BTPARAM=btpar), verbose=FALSE)
        checkIdentical(gsvaversion, getoutmd(gsvarownorm, "gsvaVersion"))
        checkTrue(is.null(getoutmd(gsvarownorm, "ranksNrow")))

        gsvacolranks <- gsvaReduce(gsvaMap(gsvaColRanks, gsvarownorm,
                                           verbose=FALSE, BTPARAM=btpar),
                                   verbose=FALSE)
        checkIdentical(gsvaversion, getoutmd(gsvacolranks, "gsvaVersion"))
        checkIdentical(unname(nrow(gsvacolranks)),
                       getoutmd(gsvacolranks, "ranksNrow"))

        es <- gsvaReduce(gsvaMap(gsvaColScores, gsvacolranks, verbose=FALSE,
                                 BTPARAM=btpar), verbose=FALSE)
        checkIdentical(gsvaversion, getoutmd(es, "gsvaVersion"))
        checkTrue(is.null(getoutmd(es, "ranksNrow")))
    }
}
