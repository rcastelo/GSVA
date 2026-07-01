#' @title MapReduce parallelization for HPC environments
#'
#' @description The functions `gsvaMap()` and `gsvaReduce()` allow one
#' to run GSVA calculations across multiple compute nodes in a high-performance
#' computing (HPC) environment. The `gsvaMap()` function is used to launch the
#' GSVA calculations on each compute node, while the `gsvaReduce()` function is
#' used to combine the results from all nodes into a single object.
#'
#' The specific HPC backend used for parallelization is determined by the
#' `BTPARAM` argument to the `gsvaMap()` function, which should be an object of
#' class [`BatchtoolsParam`][BiocParallel::BatchtoolsParam-class]. The function
#' `gsvaBatchtoolsSlurmParam()` is provided to create a `BatchtoolsParam` object
#' with some sensible defaults for running GSVA calculations on a Slurm cluster.
#' The calculations across independent compute nodes without shared memory are
#' enabled by using on-disk data structures that store enrichment scores and
#' intermediate results in HDF5 files located in a filesystem path specified in
#' the `dir` argument of `gsvaBatchtoolsSlurmParam()`, which defaults to a
#' directory named "GSVAOUTPUT" in the current working directory from where the
#' R session calling `gsvaMap()` was launched. The user must ensure that this
#' path is reachable by all compute nodes in the HPC environment, and must
#' manually delete its contents after the GSVA calculations are finished.
#'
#' @param FUN In `gsvaMap()`, function to map to the data in the `inputData`
#' argument.
#'
#' @param inputData In `gsvaMap()`, input data for performing the calculations
#' in parallel. It should be an object either of class [`gsvaParam`], or one of
#' the classes supported by [`GsvaExprData-class`] as output from
#' [`gsvaRowNorm`] or [`gsvaColRanks`]; for a list of these classes consult
#' `class ? GsvaExprData`. It can also be a list object with the output of
#' `gsvaMap()` itself, to proceed through the three-step pipeline of calculating
#' row-normalized expression values, column ranks, and GSVA scores without
#' having to call to `gsvaReduce()` in between.
#'
#' @param returnPath In `gsvaMap()`, if `TRUE`, the output of the function will
#' be a list of file paths where the resulting objects have been serialized
#' using [`saveHDF5GSVA`], instead returning the list of resulting objects
#' themselves, which is the default behavior (`FALSE`).
#'
#' @param mapOutput In `gsvaReduce()`, the output of `gsvaMap()`.
#'
#' @param verbose Gives information about the progress of the calculations.
#' Default: `TRUE`.
#'
#' @param dir In `gsvaBatchtoolsSlurmParam()`, path to a directory where the
#' output of the GSVA calculations will be saved. Default: "GSVAOUTPUT" in the
#' current working directory.
#'
#' @param partition In `gsvaBatchtoolsSlurmParam()`, name of the Slurm partition
#' to use for the GSVA calculations. No default value, the user must provide a
#' valid partition name for the Slurm cluster.
#'
#' @param walltime In `gsvaBatchtoolsSlurmParam()`, maximum wall time in seconds
#' for the GSVA calculations. Default: 600 seconds (10 minutes).
#'
#' @param nodes In `gsvaBatchtoolsSlurmParam()`, number of independent compute
#' nodes to distribute the GSVA calculations (tasks) across. Default: 1.
#'
#' @param ncpus_per_task In `gsvaBatchtoolsSlurmParam()`, number of CPU cores
#' to use for each independent task executed within a compute node. Default: 1.
#'
#' @param mem In `gsvaBatchtoolsSlurmParam()`, amount of memory to allocate for
#' each independent task executed within a compute node. Default: "10G".
#'
#' @param BTPARAM In `gsvaMap()`, an object of class
#' [`BatchtoolsParam`][BiocParallel::BatchtoolsParam-class] specifying
#' parameters for parallel execution in an HPC environment. By default, it is
#' set to a `BatchtoolsParam` object with 2 workers and a progress bar enabled,
#' and this will start a multicore execution using CPU cores in the compute node
#' where `gsvaMap()` has been called, i.e., by default it will not deploy an HPC
#' environment. For that purpose, users should either create a `BatchtoolsParam`
#' object themselves with appropriate arguments or, if an SLURM HPC environment
#' is available, they may use the wrapper function `gsvaBatchtoolsSlurmParam()`,
#' which is provided to create a `BatchtoolsParam` object with some sensible
#' defaults for running GSVA calculations on a SLURM cluster.
#'
#' @return The `gsvaMap()` function returns either a list of objects with the
#' results of the GSVA calculations for each compute node, or a list of file paths
#' where the results are saved. The `gsvaReduce()` function returns a single
#' object that combines the results from all compute nodes. The
#' `gsvaBatchtoolsSlurmParam()` function returns a `BatchtoolsParam` object with
#' some sensible defaults for running GSVA calculations on a SLURM cluster.
#'
#' @examples
#'
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#' nGrp1 <- 15 ## number of samples in group 1
#' nGrp2 <- n - nGrp1 ## number of samples in group 2
#'
#' ## consider three disjoint gene sets
#' geneSets <- list(gset1=paste0("g", 1:3),
#'                  gset2=paste0("g", 4:6),
#'                  gset3=paste0("g", 7:10))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#'
#' ## calculate row-normalized expression values in parallel across multiple
#' ## compute nodes in a high-performance computing (HPC) environment
#' gsvarnorm <- gsvaReduce(gsvaMap(gsvaRowNorm, gsvapar))
#'
#' ## calculate column GSVA ranks in parallel across multiple
#' ## compute nodes in a high-performance computing (HPC) environment
#' gsvaranks <- gsvaReduce(gsvaMap(gsvaColRanks, gsvarnorm))
#'
#' ## calculate column GSVA scores in parallel across multiple
#' ## compute nodes in a high-performance computing (HPC) environment
#' gsvaes <- gsvaReduce(gsvaMap(gsvaColScores, gsvaranks))
#'
#' ## the example below assumes that a SLURM HPC environment is available
#' ## with a partition named 'short' and that the user has write access to a
#' ## filesystem path called 'GSVAOUTPUT' in the current working directory
#' ## where this script is run and where the GSVA calculations will be saved.
#' ## The user must ensure that this path is reachable by all compute nodes in
#' ## the HPC environment, and must manually delete its contents after the GSVA
#' ## calculations are finished.
#' \dontrun{
#' gsvabtpar <- gsvaBatchtoolsSlurmParam(partition="short")
#' gsvaes <- gsvaReduce(gsvaMap(gsvaColScores, gsvaranks, BTPARAM=gsvabtpar))
#' }
#'
#' @importFrom BiocParallel BatchtoolsParam bpnworkers MulticoreParam SnowParam
#' @importFrom BiocParallel bplapply
#' @importFrom IRanges IRanges start end
#' @rdname map-reduce
#' @export gsvaMap
gsvaMap <- function(FUN, inputData, returnPath=FALSE, verbose=TRUE,
                    BTPARAM=BatchtoolsParam(workers=2, progressbar=verbose)) {

    FUN <- match.fun(FUN)

    .check_FUN_inputData(FUN, inputData)

    BTPARAM <- .check_batchtools_param(BTPARAM, verbose)

    nworkers <- bpnworkers(BTPARAM)
    ncpus <- BTPARAM$resources$ncpus
    maxmem <- .memtext2bytes(BTPARAM$resources$memory)

    gridsizefun <- .colgridsize
    splitinrangesfun <- .splitColsInRanges

    funargs <- list()
    assay <- NA_character_

    if (identical(FUN, gsvaRowNorm)) {

        funargs <- c(funargs, list(param=inputData,
                                   dropExistingAssays=TRUE,
                                   errorOnTooFewRows=FALSE))
        gridsizefun <- .rowgridsize
        splitinrangesfun <- .splitRowsInRanges
        assay <- get_assay(inputData)

    } else if (identical(FUN, gsvaColRanks)) {

        if (!"gsvarnorm" %in% gsvaAssayNames(inputData))
            cli_abort(c("x"=paste("FUN=gsvaColRanks requires inputData with",
                                  "row-normalized expression values.")))

        funargs <- c(funargs, list(rowNormExprData=inputData,
                                   dropExistingAssays=TRUE))
        assay <- "gsvarnorm"

    } else if (identical(FUN, gsvaColScores)) {

        if (!is.list(inputData)) {
            if (!"gsvaranks" %in% gsvaAssayNames(inputData))
                cli_abort(c("x"=paste("FUN=gsvaColScores requires inputData",
                                      "with column rank values.")))

            funargs <- c(funargs, list(rankExprData=inputData))
        } else
            funargs <- c(funargs, list(recompute_nzcount=TRUE))
        assay <- "gsvaranks"

    } else
        cli_abort(c("x"="Internal error, invalid FUN argument."))

    path2save <- ""
    if (returnPath)
        path2save <- path.expand(BTPARAM$registryargs$work.dir)

    X <- inputData
    totalInputDim <- NULL
    if (!is.list(X)) {
        grid <- gridsizefun(unwrapData(get_exprData(inputData), assay),
                            nworkers, maxmem)
        X <- splitinrangesfun(grid)
        totalInputDim <- dim(get_exprData(inputData))
    } else {
        if (is.null(attributes(X)$totalInputDim))
            cli_abort(c("x"=paste("If inputData is a list, it must contain",
                                  "the attribute 'totalInputDim'.")))
        totalInputDim <- attributes(X)$totalInputDim
    }

    res <- do.call("bplapply", args=c(list(X=X, FUN=MAP_FUN_WRAPPER,
                                           WRAPPED_FUN=FUN, path2save=path2save,
                                           ncpus=ncpus, maxmem=maxmem,
                                           BPPARAM=BTPARAM), funargs))
    attributes(res)$totalInputDim <- totalInputDim

    return(res)
}

#' @importFrom cli cli_abort
#' @importFrom BiocGenerics rbind cbind
#' @rdname map-reduce
#' @export gsvaReduce
gsvaReduce <- function(mapOutput, verbose=TRUE) {
    if (!is.list(mapOutput))
        cli_abort(c("x"="argument 'mapOutput' must be a list."))

    cls <- unique(lapply(mapOutput, class))
    if (length(cls) > 1)
        cli_abort(c("x"="All inputs must be of the same class."))

    totalInputDim <- attributes(mapOutput)$totalInputDim
    if (is.null(totalInputDim))
        cli_abort(c("x"=paste("The input list argument in 'mapOutput' must",
                              "contain the attribute 'totalInputDim'.")))

    if (is.character(mapOutput[[1]])) {
        mapOutput <- lapply(mapOutput, function(x) {
            if (!dir.exists(x))
                cli_abort(c("x"="Cannot find {x}."))
            loadHDF5GSVA(x)
        })
    }

    param <- .pull_param(mapOutput[[1]])
    nrmdata <- .pull_nonrestrict_metadata(mapOutput[[1]])
    rmdata <- .pull_restrict_metadata_list(mapOutput)
    ord <- .check_and_order_restrict_metadata(rmdata, totalInputDim)
    mapOutput <- .strip_metadata(mapOutput)

    if (is.null(rmdata[[1]]$whdim))
        cli_abort(c("x"=paste("The input list argument in 'mapOutput' must",
                              "contain the 'restrict' metadata with the",
                              "element 'whdim'.")))

    bfun <- "rbind"
    if (rmdata[[1]]$whdim == 2)
        bfun <- "cbind"
    res <- do.call(bfun, mapOutput[ord])
    res <- .add_metadata(res, param, nrmdata)

    rem <- vapply(rmdata, function(x) x$rem, numeric(1))
    if (dim(res)[rmdata[[1]]$whdim]+sum(rem) != totalInputDim[rmdata[[1]]$whdim])
        cli_abort(c("x"=paste("The combined output object does not match the",
                              "expected dimensions of the input data.")))

    return(res)
}

#' @importFrom BiocParallel BatchtoolsParam batchtoolsRegistryargs
#' @importFrom cli cli_alert_warning cli_abort
#' @rdname map-reduce
#' @export gsvaBatchtoolsSlurmParam
gsvaBatchtoolsSlurmParam <- function(dir="GSVAOUTPUT", partition, walltime=600,
                                     nodes=2, ncpus_per_task=2, mem="10G") { # nocov start

    if (dir.exists(dir))
        cli_alert_warning("The directory {dir} already exists.")
    else
        dir.create(dir, recursive=TRUE)

    dir <- normalizePath(dir, mustWork=TRUE)

    if (missing(partition) || is.null(partition) || !is.character(partition))
        cli_abort(c("x"=paste("You must provide a valid partition name",
                              "for the Slurm cluster.")))

    ## automatically created HDF5 datasets should be stored in a filesystem
    ## path that is reachable by all compute nodes, instead of the default
    ## tempdir() path, which is local to each compute node. this also means
    ## that, at least by now, the user must manually delete those HDF5 files
    ## after the GSVA calculations are finished
    con <- file(file.path(dir, "gsvainit.R"))
    hdf5_dump_dir <- file.path(dir, "HDF5Array_dump")
    writeLines(c("library(HDF5Array)",
                 sprintf("setHDF5DumpDir(\"%s\")", hdf5_dump_dir)),
               con)
    close(con)

    registryargs <- batchtoolsRegistryargs(file.dir=file.path(dir, "registry"),
                                           work.dir=dir, packages="GSVA",
                                           source="gsvainit.R")

    BTPARAM <- BatchtoolsParam(workers=nodes, cluster="slurm",
                               jobname="gsva",
                               resources=list(ncpus=ncpus_per_task,
                                              partition=partition,
                                              walltime=walltime,
                                              memory=mem),
                               registryargs=registryargs,
                               stop.on.error=TRUE, log=TRUE, logdir=dir)
    return(BTPARAM)
} # nocov end


## private functions

    MAP_FUN_WRAPPER <- function(X, WRAPPED_FUN, path2save, ncpus, maxmem, ...) {
        rng <- X
        res <- whdim <- NULL
        parallelbackend <- MulticoreParam(workers=ncpus)
        if (.Platform$OS.type != "unix")
            parallelbackend <- SnowParam(workers=ncpus)

        if (is(X, "IRanges")) {
            res <- WRAPPED_FUN(..., first=start(rng), last=end(rng),
                               verbose=FALSE,
                               BPPARAM=parallelbackend,
                               maxmem=maxmem)
        } else {
            rem <- 0
            if (is(X, "SummarizedExperiment")) {
                if (is.null(metadata(X)$restrict))
                    cli_abort(c("x"=paste("Input object must contain 'restrict'",
                                          "metadata with chunk boundaries.")))
                rng <- IRanges(start=metadata(X)$restrict$first,
                               end=metadata(X)$restrict$last)
                rem <- metadata(X)$restrict$rem
                whdim <- metadata(X)$restrict$whdim
            } else if (is(X, "GsvaExprData")) {
                if (is.null(attributes(X)$restrict))
                    cli_abort(c("x"=paste("Input object must contain 'restrict'",
                                          "metadata with chunk boundaries.")))
                rng <- IRanges(start=attributes(X)$restrict$first,
                               end=attributes(X)$restrict$last)
                rem <- attributes(X)$restrict$rem
                whdim <- attributes(X)$restrict$whdim
            } ## if X is path WRAPPED_FUN loads object and metadata from disk

            res <- WRAPPED_FUN(X, ..., verbose=FALSE,
                               BPPARAM=parallelbackend,
                               maxmem=maxmem)

            if (is(X, "SummarizedExperiment")) {
                metadata(res)$restrict <- list(first=start(rng),
                                               last=end(rng),
                                               rem=rem,
                                               whdim=whdim)
            } else if (is(X, "GsvaExprData")) {
                attributes(res)$restrict <- list(first=start(rng),
                                                 last=end(rng),
                                                 rem=rem,
                                                 whdim=whdim)
            } else { ## X is a path, restrict metadata is in the loaded object
                if (is(res, "SummarizedExperiment")) {
                    rng <- IRanges(start=metadata(res)$restrict$first,
                                   end=metadata(res)$restrict$last)
                } else if (is(res, "GsvaExprData")) {
                    rng <- IRanges(start=attributes(res)$restrict$first,
                                   end=attributes(res)$restrict$last)
                }
            }
        }
        if (nchar(path2save) > 0) {
            fname <- file.path(path2save, sprintf("%s_%d_%d",
                                                  basename(tempfile()),
                                                  start(rng), end(rng)))
            if (dir.exists(fname))
                cli_abort(c("x"=paste("cannot save results to {fname} because",
                                      "it already exists. You probably should",
                                      "delete the contents of {path2save}.")))
            res <- saveHDF5GSVA(res, fname)
        }
        return(res)
    }

#' @importFrom S4Vectors metadata
.pull_nonrestrict_metadata <- function(x) {
    if (is(x, "SummarizedExperiment")) {
        nrmdata <- metadata(x)
        nrmdata$restrict <- NULL
    } else {
        nrmdata <- list()
        if (!is.null(attributes(x)$assay))
            nrmdata$assay <- attributes(x)$assay
        if (!is.null(attributes(x)$geneSets))
            nrmdata$geneSets <- attributes(x)$geneSets
    }

    nrmdata
}

#' @importFrom S4Vectors "metadata<-"
.add_metadata <- function(x, param, mdata) {
    if (is(x, "SummarizedExperiment")) {
        metadata(x) <- mdata
        metadata(x)$gsvaParam <- .gsvaParam_as_list(param)
    } else {
        attributes(x) <- c(attributes(x), mdata)
        attributes(x)$gsvaParam <- .gsvaParam_as_list(param)
    }

    return(x)
}

.pull_restrict_metadata <- function(x) {
    rmdt <- NULL
    if (is(x, "SummarizedExperiment")) {
        rmdt <- metadata(x)$restrict
    } else {
        rmdt <- attributes(x)$restrict
    }

    return(rmdt)
}

#' @importFrom cli cli_abort
.pull_restrict_metadata_list <- function(args) {
    rmdt <- NULL
    if (is(args[[1]], "SummarizedExperiment")) {
        rmdt <- lapply(args, function(x) {
        metadata(x)$restrict
        })
    } else {
	    rmdt <- lapply(args, function(x) {
            attributes(x)$restrict
        })
    }
    if (any(vapply(rmdt, is.null, logical(1))))
        cli_abort(c("x"=paste("Missing 'restrict' metadata in at least",
                              "one of the inputs.")))

    return(rmdt)
}

.strip_metadata <- function(inputargs) {
    outputargs <- NULL
    if (is(inputargs[[1]], "SummarizedExperiment")) {
        outputargs <- lapply(inputargs, function(x) {
            metadata(x) <- list()
	    x
        })
    } else {
        outputargs <- lapply(inputargs, function(x) {
            attributes(x)$assay <- NULL
            attributes(x)$geneSets <- NULL
            attributes(x)$gsvaParam <- NULL
            attributes(x)$restrict <- NULL
            x
        })
    }

    outputargs
}

## check that the input restrict metadata in the input list contains chunk
## boundaries in the elements 'first' and 'last' referring to the same dimension
## stored in 'whdim', which are contiguous, non-overlapping and complete, and
## return a permutation that orders the input restrict metadata in ascending
## order of chunk boundaries defined by 'first' and 'last'.

#' @importFrom cli cli_abort
.check_and_order_restrict_metadata <- function(rmdata, totalInputDim) {
    stopifnot(is.list(rmdata))

    first <- vapply(rmdata, function(x) x$first, numeric(1))
    last <- vapply(rmdata, function(x) x$last, numeric(1))
    whdim <- vapply(rmdata, function(x) x$whdim, numeric(1))
    if (length(unique(whdim)) > 1)
        cli_abort(c("x"=paste("All input restrict metadata must refer to the",
			      "same dimension.")))

    if (whdim[1] != 1 && whdim[1] != 2)
        cli_abort(c("x"=paste("Input restrict metadata must refer to either",
                              "rows (whdim=1) or columns (whdim=2).")))

    if (!is.numeric(first) || !is.numeric(last))
        cli_abort(c("x"=paste("Input restrict metadata must contain numeric",
                              "elements 'first' and 'last'.")))
    if (any(first > last))
        cli_abort(c("x"=paste("Input restrict metadata must contain elements",
                              "'first' and 'last' such that first <= last.")))

    ord <- order(first)
    first <- first[ord]
    last <- last[ord]
    if (first[1] != 1 || last[length(last)] != totalInputDim[whdim[1]] ||
        any(first[-1] != (last+1)[-length(last)]) ||
        sum(last - first + 1) != totalInputDim[whdim[1]]) {
        tot <- totalInputDim[whdim[1]]
        cli_abort(c("x"=paste("Input restrict metadata must contain elements",
                              "'first' and 'last' starting at 1, ending at",
                              "{tot}, contiguous and non-overlapping.")))
    }

    return(ord)
}

#' @importFrom cli cli_abort cli_alert_warning
#' @importFrom BiocParallel bpprogressbar "bpprogressbar<-"
.check_batchtools_param <- function(BTPARAM, verbose) {
    if (!is(BTPARAM, "BatchtoolsParam")) {
        msg <- c("{.arg BTPARAM} must be a {.cls BatchtoolsParam} object.",
                 "x"=paste("You provided an object of class {.cls",
                           "{class(BTPARAM)}}."))
        cli_abort(msg)
    }

    if (bpnworkers(BTPARAM) < 1L) {
        msg <- c("{.arg BTPARAM} must have workers assigned.",
                 "x"=paste("The {.cls BatchtoolsParam} object you provided",
			   "has no workers assigned."))
        cli_abort(msg)
    }

    if (is.null(BTPARAM$resources))
        cli_abort(c("x"="{.arg BTPARAM} must have a resources element."))

    if (is.null(BTPARAM$resources$ncpus)) {
        cli_alert_warning(c("x"=paste("{.arg BTPARAM} has no ncpus element in",
                            "its resources. Assuming ncpus=1")))
        BTPARAM$resources$ncpus <- 1
    } else if (!is.numeric(BTPARAM$resources$ncpus) ||
             BTPARAM$resources$ncpus < 1)
        cli_abort(c("x"=paste("{.arg BTPARAM} must have a positive integer",
                    "ncpus element in its resources.")))

    if (is.null(BTPARAM$resources$memory)) {
        cli_alert_warning(c("x"=paste("{.arg BTPARAM} has no memory element",
                            "in its resources. Assuming memory=1G")))
        BTPARAM$resources$memory <- "1G"
    }

    if (is.null(BTPARAM$registryargs))
        cli_abort(c("x"="{.arg BTPARAM} must have a registryargs element."))
    else {
        if (is.null(BTPARAM$registryargs$work.dir))
            cli_abort(c("x"="{.arg BTPARAM} must have a work.dir element in",
                        "its registryargs element."))
        BTPARAM$registryargs$work.dir <- eval(BTPARAM$registryargs$work.dir)
        if (!dir.exists(BTPARAM$registryargs$work.dir))
            cli_abort(c("x"=paste("{.arg BTPARAM} must have a work.dir element",
                                  "in its registryargs element that points to",
                                  "an existing directory.")))
    }

    if (bpprogressbar(BTPARAM) != verbose)
        bpprogressbar(BTPARAM) <- verbose

    return(BTPARAM)
}

.check_FUN_inputData <- function(FUN, inputData) {

    if (!identical(FUN, gsvaRowNorm) && !identical(FUN, gsvaColRanks) &&
        !identical(FUN, gsvaColScores)) {
        msg <- paste("'FUN' must be one of 'gsvaRowNorm',",
                     "'gsvaColRanks' or 'gsvaColScores'.")
        cli_abort(c("x"=msg))
    }

    if (identical(FUN, gsvaRowNorm) && !is(inputData, "gsvaParam"))
        cli_abort(c("x"=paste("FUN=gsvaRowNorm requires inputData of",
                              "class 'gsvaParam'.")))

    if (identical(FUN, gsvaColRanks) && !is(inputData, "GsvaExprData"))
        cli_abort(c("x"=paste("FUN=gsvaColRanks requires inputData of",
                              "class 'GsvaExprData'.")))

    if (identical(FUN, gsvaColScores) && (!is(inputData, "GsvaExprData") &&
                                          !is.list(inputData)))
        cli_abort(c("x"=paste("FUN=gsvaColScores requires inputData of",
                              "class either 'GsvaExprData' or a 'list'",
                              "output from gsvaMap(gsvaColRanks, ...)")))
}
