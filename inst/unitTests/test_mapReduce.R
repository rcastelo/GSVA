## batchtools waits 5 seconds before removing a registry, and checks the state
## of the jobs every 5 seconds or more, which makes each call to gsvaMap()
## using batchtools last about 10 seconds, even with the small data in these
## tests. Setting the option BIOCPARALLEL_BATCHTOOLS_REMOVE_REGISTRY_WAIT=0
## and giving the registry the configuration file below shorten those waits,
## without changing how the jobs are run
fastRegistryargs <- function(...) {
    conf <- tempfile("batchtools.conf", fileext=".R")
    writeLines("sleep <- 0.1", conf)
    batchtoolsRegistryargs(conf.file=conf, ...)
}

test_mapReduce <- function() {

    message("Running unit tests for map reduce")

    oldopt <- options(BIOCPARALLEL_BATCHTOOLS_REMOVE_REGISTRY_WAIT=0)
    on.exit(options(oldopt), add=TRUE)

    suppressPackageStartupMessages({
        library(DelayedArray)
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

    ## calculate row-normalized expression values with map-reduce, using
    ## batchtools with two workers, to test the parallel route
    gsvarnorm2 <- gsvaReduce(gsvaMap(gsvaRowNorm, gsvapar, verbose=FALSE,
                                     BTPARAM=BatchtoolsParam(workers=2,
                                         registryargs=fastRegistryargs())),
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
    btregargs <- fastRegistryargs(work.dir=wd)
    ## except for two calls using batchtools with two workers, to test the
    ## parallel route, gsvaMap() is run in this R process with one worker,
    ## which is much faster than starting batchtools jobs, and small blocks,
    ## together with a finite maximum memory, split the input into several
    ## chunks, as two workers would do
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(nrow(gsvarnorm) * 8 * 5) ## blocks of 5 columns
    on.exit(setAutoBlockSize(oldautoblocksize), add=TRUE)
    btpar <- BatchtoolsParam(workers=1, resources=list(ncpus=1, memory="1K"),
                             registryargs=btregargs)
    btpar2 <- BatchtoolsParam(workers=2, registryargs=btregargs)
    gsvamapranks <- gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE, BTPARAM=btpar)
    checkTrue(length(gsvamapranks) > 1)
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
    gsvaesmapred <- gsvaReduce(gsvaMap(gsvaColScores, gsvaranks, verbose=FALSE,
                                       BTPARAM=btpar),
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
    ## temporary files, returning paths to results, using batchtools with two
    ## workers, which load and save files, to test the parallel route
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
    ## returning paths also when mapping by columns, in several chunks, an
    ## input that is not a SummarizedExperiment
    gsvamapranksfls2 <- gsvaMap(gsvaColRanks, gsvarnorm, output="HDF5",
                                verbose=FALSE, BTPARAM=btpar)
    checkTrue(length(gsvamapranksfls2) > 1)
    gsvaranksfls <- gsvaReduce(gsvamapranksfls2, verbose=FALSE)
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
        library(DelayedArray)
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
    ## gsvaMap() is run in this R process with one worker, which is much
    ## faster than starting batchtools jobs, and small blocks, together with a
    ## finite maximum memory, split the input into several chunks
    btpar <- BatchtoolsParam(workers=1, resources=list(ncpus=1, memory="1K"),
                             registryargs=batchtoolsRegistryargs(work.dir=wd))
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(p * 8 * 30) ## blocks of 30 columns
    on.exit(setAutoBlockSize(oldautoblocksize), add=TRUE)
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
                                      verbose=FALSE, BTPARAM=btpar),
                              verbose=FALSE)
        checkEqualsNumeric(getvals(gsvaes, "es"), getvals(gsvaes3, "es"))

        unlink(c(unlist(rnormpaths), unlist(rankspaths)))
    }

    unlink(wd, recursive=TRUE)
}

test_mapReduceScoresOutput <- function() {

    message("Running unit tests for map reduce saving GSVA scores")

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

    getvals <- function(x, a)
        as.matrix(if (is(x, "SummarizedExperiment")) assay(x, a) else x)

    formats <- "HDF5"
    if (requireNamespace("arrow", quietly=TRUE))
        formats <- c(formats, "Parquet")

    ## gsvaMap() saves its results in the working directory of the registry
    wd <- tempfile("gsvamapwd")
    dir.create(wd)
    ## gsvaMap() is run in this R process with one worker, which is much
    ## faster than starting batchtools jobs, and small blocks, together with a
    ## finite maximum memory, split the input into several chunks
    btpar <- BatchtoolsParam(workers=1, resources=list(ncpus=1, memory="1K"),
                             registryargs=batchtoolsRegistryargs(work.dir=wd))
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(p * 8 * 30) ## blocks of 30 columns
    on.exit(setAutoBlockSize(oldautoblocksize), add=TRUE)

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

test_mapReduceRedo <- function() {

    message("Running unit tests for map reduce resubmitting failed chunks")

    oldopt <- options(BIOCPARALLEL_BATCHTOOLS_REMOVE_REGISTRY_WAIT=0)
    on.exit(options(oldopt), add=TRUE)

    suppressPackageStartupMessages({
        library(DelayedArray)
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

    ## catch the warning given by gsvaMap() when some chunks fail
    mapwarn <- function(expr) {
        w <- NULL
        res <- withCallingHandlers(expr, warning=function(cond) {
            w <<- cond
            invokeRestart("muffleWarning")
        })
        list(res=res, warning=w)
    }

    ## gsvaMap() saves its results in the working directory of the registry,
    ## where a second call with the same input and output format would find
    ## the results of the first one, so each call starting a new calculation
    ## gets a new working directory. gsvaMap() is run in this R process with
    ## one worker, and small blocks, together with a finite maximum memory,
    ## split the input into several chunks, except for the calls using
    ## batchtools with two workers, to test the parallel route, where a large
    ## maximum memory splits the input into two chunks. With
    ## stop.on.error=FALSE, chunks following a failed one in the same job are
    ## also computed, instead of being reported as failed
    wds <- character(0)
    on.exit(unlink(wds, recursive=TRUE), add=TRUE)
    newbtpar <- function(workers=1, memory="1K") {
        wd <- tempfile("gsvamapwd")
        dir.create(wd)
        wds <<- c(wds, wd)
        regargs <- fastRegistryargs(work.dir=wd,
                                    file.dir=file.path(wd, "registry"))
        BatchtoolsParam(workers=workers, registryargs=regargs,
                        resources=list(ncpus=1, memory=memory),
                        stop.on.error=FALSE)
    }
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(p * 8 * 10) ## blocks of 10 columns
    on.exit(setAutoBlockSize(oldautoblocksize), add=TRUE)

    ## output is saved to disk by default only with workload managers
    checkIdentical(GSVA:::.default_map_output(list(cluster="multicore")),
                   "object")
    checkIdentical(GSVA:::.default_map_output(list(cluster="slurm")), "HDF5")

    gsvapar <- gsvaParam(y, gsets, verbose=FALSE)
    gsvarnorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    gsvaranks <- gsvaColRanks(gsvarnorm, verbose=FALSE)
    gsvaes <- gsvaColScores(gsvaranks, verbose=FALSE)
    btparranks <- newbtpar()
    rankspaths <- gsvaMap(gsvaColRanks, gsvarnorm, output="HDF5",
                          verbose=FALSE, BTPARAM=btparranks)
    checkTrue(all(file.exists(unlist(rankspaths))))
    checkTrue(!any(grepl("partial$", list.files(wds[1]))))

    ## a second call with the same input and output format in the same
    ## directory refuses to overwrite the results of the first one
    checkException(gsvaMap(gsvaColRanks, gsvarnorm, output="HDF5",
                           verbose=FALSE, BTPARAM=btparranks), silent=TRUE)

    ## simulate the failure of some chunks with in-memory results, and check
    ## that only those chunks are computed again, with the chunk boundaries of
    ## the previous call, even if the new BTPARAM argument, with a maximum
    ## memory that would put the input in a single chunk, splits it otherwise
    btpar <- newbtpar()
    mapout <- gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE, BTPARAM=btpar)
    nchunks <- length(mapout)
    checkTrue(nchunks > 2)
    whfail <- c(2L, nchunks)
    partial <- mapout
    partial[whfail] <- list(simpleError("simulated failure"))
    checkException(gsvaReduce(partial, verbose=FALSE), silent=TRUE)
    redone <- gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE,
                      BTPARAM=newbtpar(memory="10G"), MAPREDO=partial)
    checkIdentical(length(redone), nchunks)
    checkIdentical(redone[-whfail], mapout[-whfail])
    checkEqualsNumeric(gsvaranks, gsvaReduce(redone, verbose=FALSE))

    ## simulate the failure of some chunks in each of the three steps, with
    ## results saved to disk, by deleting their files, as if their jobs had
    ## been killed. Only those chunks are computed again, saving them in the
    ## directory of the previous call, given either the output of the previous
    ## call, or the path to its directory or manifest, as if the R session of
    ## the previous call had ended
    steps <- list(list(FUN=gsvaRowNorm, input=gsvapar, ref=gsvarnorm),
                  list(FUN=gsvaColRanks, input=gsvarnorm, ref=gsvaranks),
                  list(FUN=gsvaColScores, input=gsvaranks, ref=gsvaes),
                  list(FUN=gsvaColScores, input=rankspaths, ref=gsvaes))
    for (st in steps) {
        btpar <- newbtpar()
        wd <- wds[length(wds)]
        mapout <- gsvaMap(st$FUN, st$input, output="HDF5", verbose=FALSE,
                          BTPARAM=btpar)
        paths <- unlist(mapout)
        nchunks <- length(mapout)
        checkTrue(nchunks > 2)
        checkTrue(all(dirname(paths) == normalizePath(wd)))
        manifest <- list.files(wd, pattern="_manifest\\.rds$",
                               full.names=TRUE)
        checkIdentical(length(manifest), 1L)
        whfail <- c(2L, nchunks)
        mtimes <- file.mtime(paths[-whfail])
        for (redo in list("output", "dir", "manifest")) {
            unlink(paths[whfail], recursive=TRUE)
            partial <- mapout
            partial[whfail] <- list(simpleError("simulated failure"))
            mapredo <- switch(redo, output=partial, dir=wd, manifest=manifest)
            redone <- gsvaMap(st$FUN, st$input, verbose=FALSE,
                              BTPARAM=newbtpar(memory="10G"), MAPREDO=mapredo)
            checkIdentical(redone, mapout)
            checkTrue(all(file.exists(paths)))
        }
        checkIdentical(file.mtime(paths[-whfail]), mtimes)
        checkEqualsNumeric(st$ref, gsvaReduce(redone, verbose=FALSE))
        ## nothing to resubmit
        checkIdentical(gsvaMap(st$FUN, st$input, verbose=FALSE,
                               BTPARAM=btpar, MAPREDO=redone), redone)
        checkIdentical(gsvaMap(st$FUN, st$input, verbose=FALSE,
                               BTPARAM=btpar, MAPREDO=wd), redone)
    }

    ## a partial output cannot be the input of the next step
    partial <- rankspaths
    partial[[2]] <- simpleError("simulated failure")
    checkException(gsvaMap(gsvaColScores, partial, verbose=FALSE,
                           BTPARAM=btpar), silent=TRUE)

    ## MAPREDO must come from a call with the same FUN, inputData and output
    checkException(gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE,
                           BTPARAM=btpar, MAPREDO=list(1, 2)), silent=TRUE)
    checkException(gsvaMap(gsvaColScores, gsvaranks, verbose=FALSE,
                           BTPARAM=btpar, MAPREDO=partial), silent=TRUE)
    gsvarnorm2 <- gsvarnorm
    colnames(gsvarnorm2)[1] <- "other"
    checkException(gsvaMap(gsvaColRanks, gsvarnorm2, verbose=FALSE,
                           BTPARAM=btpar, MAPREDO=partial), silent=TRUE)
    checkException(gsvaMap(gsvaColRanks, gsvarnorm, output="object",
                           verbose=FALSE, BTPARAM=btpar, MAPREDO=partial),
                   silent=TRUE)
    ## and a path must have a manifest matching FUN and inputData
    checkException(gsvaMap(gsvaColRanks, gsvarnorm2, verbose=FALSE,
                           BTPARAM=btpar, MAPREDO=wds[1]), silent=TRUE)
    checkException(gsvaMap(gsvaColScores, gsvaranks, verbose=FALSE,
                           BTPARAM=btpar, MAPREDO=wds[1]), silent=TRUE)
    checkException(gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE,
                           BTPARAM=btpar, MAPREDO=tempfile()), silent=TRUE)
    checkException(gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE,
                           BTPARAM=btpar, MAPREDO=c(wds[1], wds[2])),
                   silent=TRUE)
    ## while the right ones resubmit the failed chunk
    unlink(rankspaths[[2]], recursive=TRUE)
    redone <- gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE, BTPARAM=btpar,
                      MAPREDO=partial)
    checkIdentical(redone, rankspaths)

    ## the same input saved in two formats in the same directory requires
    ## setting the output format to choose between their manifests
    if (requireNamespace("arrow", quietly=TRUE)) {
        pqpaths <- gsvaMap(gsvaColRanks, gsvarnorm, output="Parquet",
                           verbose=FALSE, BTPARAM=btparranks)
        checkTrue(all(grepl("\\.parquet$", unlist(pqpaths))))
        checkException(gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE,
                               BTPARAM=btpar, MAPREDO=wds[1]), silent=TRUE)
        unlink(pqpaths[[1]])
        redone <- gsvaMap(gsvaColRanks, gsvarnorm, output="Parquet",
                          verbose=FALSE, BTPARAM=btpar, MAPREDO=wds[1])
        checkIdentical(redone, pqpaths)
        checkIdentical(gsvaMap(gsvaColRanks, gsvarnorm, output="HDF5",
                               verbose=FALSE, BTPARAM=btpar, MAPREDO=wds[1]),
                       rankspaths)
    }

    ## real failure of a chunk whose input cannot be read, running in this R
    ## process and through batchtools, which is fixed before resubmitting it
    tmppath <- paste0(rankspaths[[2]], "_moved")
    for (workers in 1:2) {
        bp <- newbtpar(workers=workers)
        file.rename(rankspaths[[2]], tmppath)
        mw <- mapwarn(gsvaMap(gsvaColScores, rankspaths, verbose=FALSE,
                              BTPARAM=bp))
        checkTrue(!is.null(mw$warning))
        checkIdentical(which(GSVA:::.map_failed(mw$res)), 2L)
        file.rename(tmppath, rankspaths[[2]])
        redone <- gsvaMap(gsvaColScores, rankspaths, verbose=FALSE,
                          BTPARAM=bp, MAPREDO=mw$res)
        checkTrue(!any(GSVA:::.map_failed(redone)))
        checkEqualsNumeric(gsvaes, gsvaReduce(redone, verbose=FALSE))
    }

    ## jobs saving the same chunk use different temporary names, also when
    ## they run in forked processes, where tempfile() gives names that only
    ## differ in the process id
    if (.Platform$OS.type == "unix") {
        fname <- file.path(tempdir(), "chunk_1_10")
        tmpnames <- unlist(parallel::mclapply(1:4, function(i)
                                                  GSVA:::.unique_tmpname(fname),
                                              mc.cores=4))
        tmpnames <- c(tmpnames, GSVA:::.unique_tmpname(fname))
        checkTrue(!anyDuplicated(tmpnames))
        checkTrue(all(dirname(tmpnames) == dirname(fname)))
        checkTrue(all(startsWith(basename(tmpnames), "chunk_1_10.")))
        checkTrue(all(grepl("\\.partial$", tmpnames)))
    }

    ## another job, such as one left running by a call to gsvaMap() whose R
    ## session ended, saves the same chunk while this job is saving it,
    ## simulated by creating the result of that other job when this job has
    ## saved its own under a temporary name. An HDF5 directory saved by the
    ## other job is kept, while a Parquet file saved by the other job is
    ## replaced, and the temporary output of this job is not left on disk
    node <- gsub(".", "\\.", Sys.info()[["nodename"]], fixed=TRUE)
    final <- function(tmpname)
        sub(paste0("\\.", node, "\\.[0-9a-f]+\\.partial$"), "", tmpname)
    formats <- list(list(output="HDF5", fun="saveHDF5GSVA",
                         other=quote({
                             fname <- final(dir)
                             dir.create(fname)
                             file.copy(list.files(dir, full.names=TRUE),
                                       fname, recursive=TRUE)
                             file.create(file.path(fname, "otherjob"))
                         })))
    if (requireNamespace("arrow", quietly=TRUE))
        formats <- c(formats,
                     list(list(output="Parquet", fun="saveParquetGSVA",
                               other=quote(writeLines("otherjob",
                                                      final(file))))))
    for (fmt in formats) {
        bp <- newbtpar()
        wd <- wds[length(wds)]
        ## count the calls to the tracer, to check that it was run
        ncalls <- new.env()
        ncalls$n <- 0L
        tracer <- substitute({
            final <- FINAL
            OTHER
            assign("n", NCALLS$n + 1L, envir=NCALLS)
        }, list(FINAL=final, OTHER=fmt$other, NCALLS=ncalls))
        suppressMessages(trace(fmt$fun, where=asNamespace("GSVA"),
                               print=FALSE, exit=tracer))
        mapout <- tryCatch(gsvaMap(gsvaColRanks, gsvarnorm, output=fmt$output,
                                   verbose=FALSE, BTPARAM=bp),
                           finally=suppressMessages(untrace(fmt$fun,
                                                    where=asNamespace("GSVA"))))
        checkTrue(!any(GSVA:::.map_failed(mapout)))
        checkIdentical(ncalls$n, length(mapout))
        checkTrue(!any(grepl("partial$", list.files(wd))))
        if (fmt$output == "HDF5")
            checkTrue(all(file.exists(file.path(unlist(mapout), "otherjob"))))
        checkEqualsNumeric(gsvaranks, gsvaReduce(mapout, verbose=FALSE))
    }

    ## a job killed by the workload manager, simulated by a batchtools job
    ## killing its own process, in a call to gsvaMap() whose registry
    ## directory was left by a previous call whose R session ended
    if (.Platform$OS.type == "unix") {
        bp <- newbtpar(workers=2, memory="10G")
        wd <- wds[length(wds)]
        regdir <- bp$registryargs$file.dir
        dir.create(regdir)
        suppressMessages(trace("MAP_FUN_WRAPPER", where=asNamespace("GSVA"),
                               print=FALSE,
                               tracer=quote(if (is(X, "gsvaMapChunk") &&
                                                IRanges::start(X$chunk) > 1)
                                                tools::pskill(Sys.getpid(),
                                                              tools::SIGKILL))))
        mw <- tryCatch(mapwarn(gsvaMap(gsvaColRanks, gsvarnorm, output="HDF5",
                                       verbose=FALSE, BTPARAM=bp)),
                       finally=suppressMessages(untrace("MAP_FUN_WRAPPER",
                                                where=asNamespace("GSVA"))))
        checkTrue(!is.null(mw$warning))
        checkIdentical(which(GSVA:::.map_failed(mw$res)), 2L)
        checkTrue(file.exists(mw$res[[1]]))
        ## the registry directory given in BTPARAM is restored
        checkIdentical(bp$registryargs$file.dir, regdir)
        redone <- gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE,
                          BTPARAM=newbtpar(), MAPREDO=wd)
        checkTrue(!any(GSVA:::.map_failed(redone)))
        checkEqualsNumeric(gsvaranks, gsvaReduce(redone, verbose=FALSE))
    }
}
