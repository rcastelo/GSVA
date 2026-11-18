source(system.file("unitTests", "test_helpers.R", package = "GSVA"))

test_parallel <- function() {
    message("Running unit tests for parallel execution.")

    suppressPackageStartupMessages({
        library(BiocParallel)
        library(Matrix)
    })

    p <- 150 ## number of genes
    n <- 150 ## number of samples
    m <- 10 ## number of gene sets

    ## create random gene sets of size between 3 and 10
    set.seed(123)
    gssizes <- sample(3:10, size=m, replace=TRUE)
    gsets <- lapply(gssizes, function(x) sample(paste0("g", 1:p), size=x, replace=FALSE))
    names(gsets) <- paste0("gs", 1:m)

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    s <- ceiling(0.15 * p * n)
    sam <- sample(1:(p * n), size=s, replace=FALSE)
    x <- numeric(p * n)
    x[sam] <- runif(s)
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    M <- Matrix(y, sparse=TRUE)

    ## estimate GSVA enrichment scores with and without parallel execution
    ## and check that they are identical
    es_serial <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE, device="cpu")
    es_parallel <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE,
                        BPPARAM=MulticoreParam(workers=2), device="cpu")
    checkIdentical(es_serial, es_parallel)

    ## CPU vs GPU. the GPU scores the gene sets in this R process, without
    ## the workers of a MulticoreParam, which are forked from it after it has
    ## initialized CUDA, and calculate the row normalization and the ranks
    .run_with_gpu({
        es_serial_gpu <- gsva(gsvaParam(M, gsets, verbose=FALSE),
                              verbose=FALSE, device="gpu")
        .check_gpu_equiv(es_serial, es_serial_gpu)
        es_parallel_gpu <- gsva(gsvaParam(M, gsets, verbose=FALSE),
                                verbose=FALSE,
                                BPPARAM=MulticoreParam(workers=2),
                                device="gpu")
        checkIdentical(es_serial_gpu, es_parallel_gpu)

        ## a process forked from this one cannot use CUDA, which this one
        ## has initialized: device="auto" selects the CPU there, and
        ## device="gpu" gives an error
        if (.Platform$OS.type == "unix") {
            res <- bplapply(1:2, function(i)
                list(state=GSVA:::.gsva_cuda_state(),
                     auto=gsva(gsvaParam(M, gsets, verbose=FALSE),
                               verbose=FALSE),
                     gpu=tryCatch(gsva(gsvaParam(M, gsets, verbose=FALSE),
                                       verbose=FALSE, device="gpu"),
                                  error=function(e) e)),
                BPPARAM=MulticoreParam(workers=2))
            checkIdentical(res[[1]]$state, "init")
            checkIdentical(res[[1]]$auto, es_serial)
            checkTrue(is(res[[1]]$gpu, "error"))
        }

        ## the CPU scores all gene sets of data with missing values, with
        ## the workers of BPPARAM, also with device="gpu"
        yna <- y
        yna[1, 2] <- NA
        es_na_cpu <- gsva(gsvaParam(yna, gsets, verbose=FALSE),
                          verbose=FALSE, device="cpu")
        msgs <- character(0)
        es_na_gpu <- withCallingHandlers(
            gsva(gsvaParam(yna, gsets, verbose=FALSE), verbose=TRUE,
                 BPPARAM=MulticoreParam(workers=2, progressbar=FALSE),
                 device="gpu"),
            message=function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
            })
        checkIdentical(es_na_cpu, es_na_gpu)
        checkTrue(any(grepl("scores on the CPU", msgs)))
        checkTrue(!any(grepl("on the GPU", msgs)))
    })

    M[1, 2] <- NA
    checkException(gsva(gsvaParam(M, gsets, kcdf="Gaussian", verbose=FALSE),
			verbose=FALSE, BPPARAM=MulticoreParam(workers=2), device="cpu"))
}

test_device_selection <- function() {
    message("Running unit tests for the selection of the GPU or the CPU")

    ## data with missing values are scored on the CPU
    checkIdentical(GSVA:::.resolve_device("auto", TRUE), "cpu")
    checkIdentical(GSVA:::.resolve_device("cpu", FALSE), "cpu")

    ## the jobs of gsvaMap() can use the GPU only when they are not processes
    ## forked from this R session; the fields of BatchtoolsParam objects that
    ## these functions use are taken from lists, to avoid needing Slurm
    mc <- list(cluster="multicore", resources=list())
    sl <- list(cluster="slurm", resources=list())
    slgpu <- list(cluster="slurm", resources=list(gres="gpu:1"))
    checkIdentical(GSVA:::.map_device("cpu", mc), "cpu")
    checkIdentical(GSVA:::.map_device("gpu", sl), "gpu")
    checkIdentical(GSVA:::.map_device("gpu", list(cluster="socket")), "gpu")
    checkException(GSVA:::.map_device("gpu", mc), silent=TRUE)

    ## jobs of Slurm using the GPU without requesting it give a warning, and
    ## those requesting it without using it, a message
    .conds <- function(expr) {
        conds <- character(0)
        withCallingHandlers(expr,
            warning=function(w) {
                conds <<- c(conds, "warning")
                invokeRestart("muffleWarning")
            },
            message=function(m) {
                conds <<- c(conds, "message")
                invokeRestart("muffleMessage")
            })
        conds
    }
    checkIdentical(.conds(GSVA:::.check_map_gres(TRUE, sl, TRUE)), "warning")
    checkIdentical(.conds(GSVA:::.check_map_gres(TRUE, slgpu, TRUE)),
                   character(0))
    checkIdentical(.conds(GSVA:::.check_map_gres(FALSE, slgpu, TRUE)),
                   "message")
    checkIdentical(.conds(GSVA:::.check_map_gres(FALSE, slgpu, FALSE)),
                   character(0))
    checkIdentical(.conds(GSVA:::.check_map_gres(FALSE, sl, TRUE)),
                   character(0))
    checkIdentical(.conds(GSVA:::.check_map_gres(TRUE, mc, TRUE)),
                   character(0))

    ## gsvaMap() gives the error before sending any job, and has no "auto"
    ## device
    if (.Platform$OS.type == "unix") {
        set.seed(123)
        y <- matrix(rnorm(200 * 20), nrow=200, ncol=20,
                    dimnames=list(paste0("g", 1:200), paste0("s", 1:20)))
        gsets <- list(gs1=paste0("g", 1:20), gs2=paste0("g", 21:60))
        rnk <- gsvaColRanks(gsvaRowNorm(gsvaParam(y, gsets), verbose=FALSE),
                            verbose=FALSE)
        btmc <- BiocParallel::BatchtoolsParam(workers=1, cluster="multicore",
                                              progressbar=FALSE,
                                              resources=list(ncpus=1,
                                                             memory="2G"))
        checkException(gsvaMap(gsvaColScores, rnk, verbose=FALSE,
                               BTPARAM=btmc, device="gpu"), silent=TRUE)
        checkException(gsvaMap(gsvaColScores, rnk, verbose=FALSE,
                               BTPARAM=btmc, device="auto"), silent=TRUE)
    }
}

test_parallel_bpparam_unchanged <- function() {

    message("Running unit tests for parallel execution keeping BPPARAM unchanged")

    suppressPackageStartupMessages({
        library(BiocParallel)
        library(SpatialExperiment)
    })

    set.seed(123)
    p <- 200
    n <- 150
    y <- matrix(rnorm(p * n), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    gsets <- list(gs1=paste0("g", 1:20), gs2=paste0("g", 21:60),
                  gs3=paste0("g", 61:120))

    ## with verbose=TRUE, the parallel execution sets progressbar=TRUE in
    ## BPPARAM, which is a reference class object, and this setting is
    ## restored when the calculations finish, also when spatCor() loops over
    ## several samples
    xy <- cbind(x=runif(n), y=runif(n))
    rownames(xy) <- colnames(y)
    spe <- SpatialExperiment(assays=list(es=y[1:3, ]),
                             sample_id=rep(c("A", "B"), each=n/2),
                             spatialCoords=xy)
    bpparam <- if (.Platform$OS.type == "unix")
                   MulticoreParam(workers=2, progressbar=FALSE)
               else
                   SnowParam(workers=2, progressbar=FALSE)
    suppressMessages(spatCor(spe, assay="es", verbose=TRUE, BPPARAM=bpparam))
    checkTrue(!bpprogressbar(bpparam))

    ## and in gsvaMap(), where BTPARAM gets the value of the verbose argument
    btparam <- BatchtoolsParam(workers=1, progressbar=FALSE,
                               resources=list(ncpus=1, memory="2G"))
    suppressMessages(gsvaMap(gsvaRowNorm, gsvaParam(y, gsets), verbose=TRUE,
                             BTPARAM=btparam))
    checkTrue(!bpprogressbar(btparam))
}
