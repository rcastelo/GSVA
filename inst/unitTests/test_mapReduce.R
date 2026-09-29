test_mapReduce <- function() {

    message("Running unit tests for map reduce")

    suppressPackageStartupMessages({
        library(Matrix)
        library(GSEABase)
        library(BiocParallel)
        library(SingleCellExperiment)
        library(scrapper)
    })

    p <- 10 ## number of genes
    n <- 30 ## number of samples

    ## consider three disjoint gene sets
    gsets <- list(gset1=paste0("g", 1:3),
                  gset2=paste0("g", 4:6),
                  gset3=paste0("g", 7:10))

    ## build a random sparse count matrix with 85% sparsity
    set.seed(123)
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.15)) ## 85% sparsity
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- matrix(cnt, nrow=p, ncol=n,
                  dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    cnt <- Matrix(cnt, sparse=TRUE)

    sce <- SingleCellExperiment(assays=list(counts=cnt))

    ## process it as if it were single-cell RNA-seq data
    sce <- quickRnaQc.se(sce, subsets=list(mito=rep(FALSE, nrow(sce))))
    sce <- sce[, sce$keep]
    sce <- normalizeRnaCounts.se(sce, size.factors=sce$sum)

    ## build GSVA parameter object
    gsvapar <- gsvaParam(sce, gsets, verbose=FALSE)

    ## calculate row-normalized expression values without map-reduce
    gsvarnorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## calculate row-normalized expression values with map-reduce
    gsvarnorm2 <- gsvaReduce(gsvaMap(gsvaRowNorm, gsvapar, verbose=FALSE),
                             verbose=FALSE)

    ## check that both approaches yield same row-normalized expression values
    checkEqualsNumeric(assay(gsvarnorm, "gsvarnorm"),
                       assay(gsvarnorm2, "gsvarnorm"))

    ## calculate column rank values without map-reduce
    gsvaranks <- gsvaColRanks(gsvarnorm, verbose=FALSE)

    ## calculate column rank values with map-reduce
    ## gsvaMap() saves its results in the working directory of the registry,
    ## which by default is the current working directory
    wd <- tempfile("gsvamapwd")
    dir.create(wd)
    btregargs <- batchtoolsRegistryargs(work.dir=wd)
    btpar <- BatchtoolsParam(workers=1, registryargs=btregargs) ## just to force unit testing internal MAP_FUN_WRAPPER()
    ## two workers to have more than one chunk
    btpar2 <- BatchtoolsParam(workers=2, registryargs=btregargs)
    gsvamapranks <- gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE, BTPARAM=btpar)
    gsvaredranks <- gsvaReduce(gsvamapranks, verbose=FALSE)

    ## check that both approaches yield the same column rank values
    checkEqualsNumeric(assay(gsvaranks, "gsvaranks"),
                       assay(gsvaredranks, "gsvaranks"))

    ## check that gsvaReduce() correctly reduces chunks provided in
    ## a non-sequential order
    p <- sample(seq_len(length(gsvamapranks))) ## permutation
    gsvamapranks2 <- gsvamapranks[p]
    attributes(gsvamapranks2) <- attributes(gsvamapranks)
    gsvaredranks <- gsvaReduce(gsvamapranks2, verbose=FALSE)

    ## check that both approaches yield the same column rank values
    checkEqualsNumeric(assay(gsvaranks, "gsvaranks"),
                       assay(gsvaredranks, "gsvaranks"))

    ## calculate column rank values with map-reduce returning paths to results
    gsvamapranksfls <- gsvaMap(gsvaColRanks, gsvarnorm, output="HDF5",
                               verbose=FALSE, BTPARAM=btpar)
    gsvaredranksfls <- gsvaReduce(gsvamapranksfls, verbose=FALSE)

    ## check that this approach also yields the same column rank values
    checkEqualsNumeric(assay(gsvaranks, "gsvaranks"),
                       assay(gsvaredranksfls, "gsvaranks"))

    ## calculate column GSVA scores without map-reduce
    gsvaes <- gsvaColScores(gsvaranks, verbose=FALSE)

    ## calculate column GSVA scores with map-reduce
    gsvaesmapred <- gsvaReduce(gsvaMap(gsvaColScores, gsvaranks, verbose=FALSE),
                               verbose=FALSE)

    ## check that both approaches yield the same column GSVA scores
    checkEqualsNumeric(assay(gsvaes, "es"),
                       assay(gsvaesmapred, "es"))

    ## calculate column GSVA scores with map-reduce on mapped ranks
    gsvaesmaprnkred <- gsvaReduce(gsvaMap(gsvaColScores, gsvamapranks,
                                          verbose=FALSE, BTPARAM=btpar),
                                  verbose=FALSE)

    ## check that we obtain the same column GSVA scores as before
    checkEqualsNumeric(assay(gsvaes, "es"),
                       assay(gsvaesmaprnkred, "es"))

    ## calculate column GSVA scores with map-reduce on mapped ranks stored in
    ## temporary files
    gsvaesmaprnkflsred <- gsvaReduce(gsvaMap(gsvaColScores, gsvamapranksfls,
                                             verbose=FALSE, BTPARAM=btpar),
                                     verbose=FALSE)

    ## check that we obtain the same column GSVA scores as before
    checkEqualsNumeric(assay(gsvaes, "es"),
                       assay(gsvaesmaprnkflsred, "es"))

    ## calculate column GSVA scores with map-reduce on mapped ranks stored in
    ## temporary files, returning paths to results
    gsvaesmaprnkflsredfls <- gsvaReduce(gsvaMap(gsvaColScores, gsvamapranksfls,
                                                output="HDF5", verbose=FALSE,
                                                BTPARAM=btpar2),
                                        verbose=FALSE)

    ## check that we obtain the same column GSVA scores as before
    checkEqualsNumeric(assay(gsvaes, "es"),
                       assay(gsvaesmaprnkflsredfls, "es"))

    ## check with input expression data stored in a matrix
    expr <- as(logcounts(sce), "matrix")
    gsvapar <- gsvaParam(expr, gsets, verbose=FALSE)
    gsvaes <- gsva(gsvapar, verbose=FALSE)
    ## set output="HDF5" once to test stripping of attributes and wrapping into an SE for saving
    gsvarnorm <- gsvaReduce(gsvaMap(gsvaRowNorm, gsvapar, output="HDF5", verbose=FALSE, BTPARAM=btpar),
                            verbose=FALSE)
    gsvaranks <- gsvaReduce(gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE, BTPARAM=btpar),
                            verbose=FALSE)
    ## returning paths also when mapping by columns an input that is not a
    ## SummarizedExperiment, with two workers to have more than one chunk
    gsvaranksfls <- gsvaReduce(gsvaMap(gsvaColRanks, gsvarnorm, output="HDF5",
                                       verbose=FALSE, BTPARAM=btpar2),
                               verbose=FALSE)
    checkEqualsNumeric(gsvaranks, gsvaranksfls)
    gsvaes2 <- gsvaReduce(gsvaMap(gsvaColScores, gsvaranks, verbose=FALSE, BTPARAM=btpar),
                          verbose=FALSE)
    checkEqualsNumeric(gsvaes, gsvaes2)

    unlink(wd, recursive=TRUE)
}

test_mapReduceParquet <- function() {

    if (!requireNamespace("arrow", quietly=TRUE)) {
        message("Skipping unit tests for map reduce with Parquet files (no 'arrow')")
        return(invisible(TRUE))
    }

    message("Running unit tests for map reduce with Parquet files")

    suppressPackageStartupMessages({
        library(Matrix)
        library(SummarizedExperiment)
        library(BiocParallel)
    })

    p <- 40 ## number of genes
    n <- 90 ## number of samples
    gsets <- list(gset1=paste0("g", 1:10),
                  gset2=paste0("g", 11:25),
                  gset3=paste0("g", 26:40))

    set.seed(123)
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.15)) ## 85% sparsity
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- Matrix(matrix(cnt, nrow=p, ncol=n,
                         dimnames=list(paste0("g", 1:p), paste0("s", 1:n))),
                  sparse=TRUE)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n, dimnames=dimnames(cnt))

    getvals <- function(x, a)
        as.matrix(if (is(x, "SummarizedExperiment")) assay(x, a) else x)

    ## gsvaMap() saves its results in the working directory of the registry
    wd <- tempfile("gsvamapwd")
    dir.create(wd)
    btpar <- BatchtoolsParam(workers=2,
                             registryargs=batchtoolsRegistryargs(work.dir=wd))
    isparquet <- function(paths)
        all(vapply(paths, function(x) grepl("\\.parquet$", x) &&
                                      GSVA:::.is_parquet_file(x), logical(1)))

    checkException(gsvaMap(gsvaRowNorm, gsvaParam(y, gsets, verbose=FALSE),
                           output="csv", verbose=FALSE), silent=TRUE)

    for (input in list(SummarizedExperiment(assays=list(counts=cnt)), y)) {
        gsvapar <- gsvaParam(input, gsets, verbose=FALSE)
        gsvarnorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
        gsvaranks <- gsvaColRanks(gsvarnorm, verbose=FALSE)
        gsvaes <- gsvaColScores(gsvaranks, verbose=FALSE)

        ## row-normalized values mapped by rows
        rnormpaths <- gsvaMap(gsvaRowNorm, gsvapar, output="Parquet",
                              verbose=FALSE, BTPARAM=btpar)
        checkTrue(length(rnormpaths) > 1 && isparquet(rnormpaths))
        gsvarnorm2 <- gsvaReduce(rnormpaths, verbose=FALSE)
        checkEqualsNumeric(getvals(gsvarnorm, "gsvarnorm"),
                           getvals(gsvarnorm2, "gsvarnorm"))

        ## column ranks mapped by columns, given in a non-sequential order
        rankspaths <- gsvaMap(gsvaColRanks, gsvarnorm, output="Parquet",
                              verbose=FALSE, BTPARAM=btpar)
        checkTrue(length(rankspaths) > 1 && isparquet(rankspaths))
        rankspaths2 <- rev(rankspaths)
        attributes(rankspaths2) <- attributes(rankspaths)
        gsvaranks2 <- gsvaReduce(rankspaths2, verbose=FALSE)
        checkEqualsNumeric(getvals(gsvaranks, "gsvaranks"),
                           getvals(gsvaranks2, "gsvaranks"))
        rnks2 <- if (is(gsvaranks2, "SummarizedExperiment"))
                     assay(gsvaranks2, "gsvaranks") else gsvaranks2
        checkTrue(GSVA:::.is_parquet_backed(rnks2))

        ## scores from reduced ranks read from the Parquet files, and from
        ## the Parquet files given as input to gsvaMap()
        gsvaes2 <- gsvaColScores(gsvaranks2, verbose=FALSE)
        checkEqualsNumeric(getvals(gsvaes, "es"), getvals(gsvaes2, "es"))
        gsvaes3 <- gsvaReduce(gsvaMap(gsvaColScores, rankspaths,
                                      verbose=FALSE), verbose=FALSE)
        checkEqualsNumeric(getvals(gsvaes, "es"), getvals(gsvaes3, "es"))

        unlink(c(unlist(rnormpaths), unlist(rankspaths)))
    }

    unlink(wd, recursive=TRUE)
}

test_mapReduceScoresOutput <- function() {

    message("Running unit tests for map reduce saving GSVA scores")

    suppressPackageStartupMessages({
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

    getvals <- function(x, a)
        as.matrix(if (is(x, "SummarizedExperiment")) assay(x, a) else x)

    formats <- "HDF5"
    if (requireNamespace("arrow", quietly=TRUE))
        formats <- c(formats, "Parquet")

    ## gsvaMap() saves its results in the working directory of the registry
    wd <- tempfile("gsvamapwd")
    dir.create(wd)
    btpar <- BatchtoolsParam(workers=2,
                             registryargs=batchtoolsRegistryargs(work.dir=wd))

    for (input in list(y, SummarizedExperiment(assays=list(exprs=y)))) {
        gsvaranks <- gsvaColRanks(gsvaRowNorm(gsvaParam(input, gsets,
                                                        verbose=FALSE),
                                              verbose=FALSE),
                                  verbose=FALSE)
        gsvaes <- gsvaColScores(gsvaranks, verbose=FALSE)

        for (fmt in formats) {
            espaths <- gsvaMap(gsvaColScores, gsvaranks, output=fmt,
                               verbose=FALSE, BTPARAM=btpar)
            checkTrue(length(espaths) > 1)
            checkTrue(all(vapply(espaths, is.character, logical(1))))
            gsvaes2 <- gsvaReduce(espaths, verbose=FALSE)
            checkEqualsNumeric(getvals(gsvaes, "es"), getvals(gsvaes2, "es"))
            if (!is(input, "SummarizedExperiment"))
                checkIdentical(attr(gsvaes, "geneSets"),
                               attr(gsvaes2, "geneSets"))
            unlink(unlist(espaths), recursive=TRUE)
        }
    }

    unlink(wd, recursive=TRUE)
}
