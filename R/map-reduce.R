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
#' Each job saves its results under a temporary name ending in `.partial`,
#' which it renames once they are complete, so files or directories with that
#' ending left in that path belong to jobs that did not finish, and can be
#' deleted once no job of those calculations is running.
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
#' @param output In `gsvaMap()`, a character string specifying what each
#' worker returns: either the resulting object (`"object"`, default), or the
#' path to that object after saving it in the working directory of the
#' registry of `BTPARAM` with [`saveHDF5GSVA`] (`"HDF5"`) or with
#' [`saveParquetGSVA`] (`"Parquet"`). In the latter two cases, the output of
#' `gsvaMap()` is a list of paths instead of a list of objects. Note that, by
#' default, the working directory of the registry is the current working
#' directory. Saving in Apache Parquet format requires the package
#' [arrow](https://cran.r-project.org/package=arrow). When this argument is
#' not given, its default value is `"HDF5"` if `BTPARAM` sends the jobs to a
#' workload manager, such as Slurm, and `"object"` otherwise, i.e., when the
#' `cluster` field of `BTPARAM` is `"socket"`, `"multicore"` or
#' `"interactive"`. When `MAPREDO` is given, its default value is the one
#' used to produce `MAPREDO`. The names of the saved files depend only on
#' `FUN`, `inputData`, `output` and the chunk boundaries, and `gsvaMap()`
#' also writes in that directory a manifest file, ending in `_manifest.rds`,
#' which allows the `MAPREDO` argument to find them. For this reason,
#' `gsvaMap()` gives an error when that directory already has results of a
#' previous call with the same `FUN`, `inputData` and `output` arguments,
#' unless `MAPREDO` is given.
#'
#' @param MAPREDO In `gsvaMap()`, the output of a previous call to `gsvaMap()`
#' with the same `FUN` and `inputData` arguments, in which some of the chunks
#' failed. When this argument is given, only the failed chunks are computed
#' again, using the same chunk boundaries as in the previous call, and the
#' result combines them with the chunks that did not fail. This makes it
#' possible to resubmit the failed chunks with a different `BTPARAM` argument,
#' for instance with more memory or a longer wall time. When the previous call
#' saved its results to disk (see the `output` argument), this argument can
#' also be the path to the directory where they were saved, or to the
#' manifest file that `gsvaMap()` writes there, and the chunks whose results
#' are not found in that directory are computed again and saved in it. This
#' allows one to resubmit chunks whose jobs were killed by the workload
#' manager, or that were running when the R session calling `gsvaMap()`
#' ended. If that directory has results of several calls matching `FUN` and
#' `inputData`, saved in different formats, the `output` argument selects
#' one of them. Note that `inputData` is matched by its dimensions, dimension
#' names and GSVA parameters, but not by its values. Default: `NULL`, which
#' computes all chunks.
#'
#' @param mapOutput In `gsvaReduce()`, the output of `gsvaMap()`, which can be
#' a list of objects or a list of paths to GSVA output saved with
#' [`saveHDF5GSVA`] or [`saveParquetGSVA`], such as the one returned by
#' `gsvaMap()` with `output="HDF5"` or `output="Parquet"`.
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
#' where the results are saved. If the calculations in some of the chunks
#' failed, `gsvaMap()` gives a warning and the corresponding elements of this
#' list are condition objects describing the errors. This partial result can
#' be given to the `MAPREDO` argument of a new call to `gsvaMap()` to compute
#' only the failed chunks, but not to `gsvaReduce()`, which gives an error.
#' When results are saved to disk, chunks of jobs that did not deliver any
#' result, such as jobs killed by the workload manager, are also reported as
#' failed, unless no chunk was completed, in which case `gsvaMap()` gives an
#' error. The `gsvaReduce()` function returns a single
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
#' ## to keep this example short, the calculations below are split in chunks
#' ## that run in this R session, by using a single worker. See further below
#' ## how to distribute them across multiple compute nodes in a
#' ## high-performance computing (HPC) environment
#' library(BiocParallel)
#' btpar <- BatchtoolsParam(workers=1, resources=list(ncpus=1, memory="1G"))
#'
#' ## calculate row-normalized expression values
#' gsvarnorm <- gsvaReduce(gsvaMap(gsvaRowNorm, gsvapar, BTPARAM=btpar))
#'
#' ## calculate column GSVA ranks
#' gsvaranks <- gsvaReduce(gsvaMap(gsvaColRanks, gsvarnorm, BTPARAM=btpar))
#'
#' ## calculate column GSVA scores
#' gsvaes <- gsvaReduce(gsvaMap(gsvaColScores, gsvaranks, BTPARAM=btpar))
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
#'
#' ## if some of the chunks fail, gsvaMap() gives a warning and returns a
#' ## partial result, whose failed chunks can be resubmitted, for instance
#' ## with more memory per task, without computing again the other chunks
#' gsvamapes <- gsvaMap(gsvaColScores, gsvaranks, BTPARAM=gsvabtpar)
#' gsvabtpar2 <- gsvaBatchtoolsSlurmParam(partition="short", mem="20G")
#' gsvamapes <- gsvaMap(gsvaColScores, gsvaranks, BTPARAM=gsvabtpar2,
#'                      MAPREDO=gsvamapes)
#' gsvaes <- gsvaReduce(gsvamapes)
#'
#' ## if the R session calling gsvaMap() ended before the jobs finished, the
#' ## chunks without results can be resubmitted from a new R session using
#' ## the directory where gsvaMap() saved its results
#' gsvamapes <- gsvaMap(gsvaColScores, gsvaranks, BTPARAM=gsvabtpar,
#'                      MAPREDO="GSVAOUTPUT")
#' }
#'
#' @importFrom BiocParallel BatchtoolsParam bpnworkers
#' @importFrom IRanges IRanges start end
#' @importFrom cli cli_abort cli_alert_info cli_warn
#' @rdname map-reduce
#' @export gsvaMap
gsvaMap <- function(FUN, inputData, output=c("object", "HDF5", "Parquet"),
                    verbose=TRUE,
                    BTPARAM=BatchtoolsParam(workers=2, progressbar=verbose),
                    MAPREDO=NULL) {

    FUN <- match.fun(FUN)
    .check_FUN_inputData(FUN, inputData)
    funname <- .map_fun_name(FUN)

    BTPARAM <- .check_batchtools_param(BTPARAM, verbose)

    totalInputDim <- NULL
    if (!is.list(inputData))
        totalInputDim <- dim(get_exprData(inputData))
    else {
        if (is.null(attributes(inputData)$totalInputDim))
            cli_abort(c("x"=paste("If inputData is a list, it must contain",
                                  "the attribute 'totalInputDim'.")))
        totalInputDim <- attributes(inputData)$totalInputDim
    }
    fingerprint <- .map_fingerprint(funname, inputData, totalInputDim)

    mapinfo <- NULL
    if (!is.null(MAPREDO)) {
        if (is.character(MAPREDO)) ## path to a manifest or its directory
            MAPREDO <- .read_map_manifest(MAPREDO, funname, fingerprint,
                                          if (missing(output)) NULL
                                          else match.arg(output))
        mapinfo <- .check_MAPREDO(MAPREDO, funname)
    }

    if (missing(output)) {
        output <- .default_map_output(BTPARAM)
        if (!is.null(mapinfo))
            output <- mapinfo$output
    } else
        output <- match.arg(output)
    if (!is.null(mapinfo) && output != mapinfo$output)
        cli_abort(c("x"=paste("{.arg MAPREDO} was produced with",
                              "output=\"{mapinfo$output}\", which differs",
                              "from output=\"{output}\".")))
    if (output == "Parquet") ## fail before sending any job to the workers
        .require_arrow()

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

    X <- inputData
    chunks <- NULL
    if (!is.list(X)) {
        if (is.null(mapinfo)) {
            grid <- gridsizefun(unwrapData(get_exprData(inputData), assay),
                                nworkers, maxmem)
            X <- splitinrangesfun(grid)
            chunks <- data.frame(first=vapply(X, start, integer(1)),
                                 last=vapply(X, end, integer(1)))
        } else { ## the chunks of the previous run, which may have used a
                 ## different number of workers or maximum memory
            chunks <- mapinfo$chunks
            X <- lapply(seq_len(nrow(chunks)), function(i)
                            IRanges(start=chunks$first[i], end=chunks$last[i]))
        }
    } else {
        chunks <- .list_input_chunks(inputData)
        attributes(X) <- NULL
    }

    if (!is.null(mapinfo) && (fingerprint != mapinfo$fingerprint ||
                              length(MAPREDO) != length(X)))
        cli_abort(c("x"=paste("{.arg MAPREDO} was not produced by",
                              "{.fn gsvaMap} with FUN={funname} on the",
                              "same {.arg inputData}.")))

    ## results saved to disk get file names that depend only on the input
    ## and the chunk boundaries, and a manifest that allows MAPREDO to find
    ## them, even after the R session calling gsvaMap() has ended
    outdir <- files <- paths <- NULL
    if (output != "object") {
        if (is.null(mapinfo)) {
            outdir <- normalizePath(path.expand(BTPARAM$registryargs$work.dir),
                                    mustWork=TRUE)
            files <- .map_file_names(funname, fingerprint, output, chunks,
                                     length(X))
        } else {
            outdir <- mapinfo$dir
            files <- mapinfo$files
        }
        paths <- file.path(outdir, files)
    }
    mapinfo <- list(FUN=funname, output=output, fingerprint=fingerprint,
                    chunks=chunks, files=files, dir=outdir)
    if (output != "object")
        .write_map_manifest(mapinfo, totalInputDim,
                            mustNotExist=is.null(MAPREDO))

    idx <- seq_along(X)
    res <- vector("list", length(X))
    if (!is.null(MAPREDO)) {
        res <- MAPREDO
        if (output == "object")
            idx <- which(.map_failed(MAPREDO))
        else { ## the results on disk tell which chunks are missing
            done <- file.exists(paths)
            res[done] <- as.list(paths[done])
            idx <- which(!done)
        }
        if (length(idx) == 0) {
            if (verbose)
                cli_alert_info(paste("No failed chunks in {.arg MAPREDO},",
                                     "nothing to resubmit."))
            attributes(res) <- list(totalInputDim=totalInputDim,
                                    mapInfo=mapinfo)
            return(res)
        }
        if (verbose)
            cli_alert_info(paste("Resubmitting {length(idx)} failed chunk{?s}",
                                 "out of {length(X)}."))
    }

    if (output != "object") ## each worker gets the name of its output file
        X <- lapply(seq_along(X), function(i)
                        structure(list(chunk=X[[i]], fname=paths[i]),
                                  class="gsvaMapChunk"))

    ## a registry left on disk by a call whose R session ended prevents
    ## starting a new one in the same directory
    if (nworkers > 1) {
        regdir <- .free_registry_dir(BTPARAM, verbose)
        if (!is.null(regdir))
            on.exit(BTPARAM$registryargs$file.dir <- regdir, add=TRUE)
    }

    res[idx] <- .map_chunks(X[idx], BTPARAM,
                            c(list(WRAPPED_FUN=FUN, output=output,
                                   ncpus=ncpus, maxmem=maxmem), funargs),
                            paths[idx], outdir)
    attributes(res) <- list(totalInputDim=totalInputDim, mapInfo=mapinfo)

    failed <- .map_failed(res)
    if (any(failed)) {
        nfailed <- sum(failed)
        firsterr <- conditionMessage(res[[which(failed)[1]]])
        redo <- "this result"
        if (output != "object")
            redo <- paste("this result, or the directory", outdir)
        cli_warn(c("!"="{nfailed} out of {length(res)} chunk{?s} failed.",
                   "i"="First error: {firsterr}",
                   "i"=paste("Resubmit only the failed chunks by calling",
                             "{.fn gsvaMap} again with the same {.arg FUN}",
                             "and {.arg inputData}, and {redo} in",
                             "{.arg MAPREDO}.")))
    }

    return(res)
}

#' @importFrom cli cli_abort
#' @importFrom BiocGenerics rbind cbind
#' @rdname map-reduce
#' @export gsvaReduce
gsvaReduce <- function(mapOutput, verbose=TRUE) {
    if (!is.list(mapOutput))
        cli_abort(c("x"="argument 'mapOutput' must be a list."))

    failed <- .map_failed(mapOutput)
    if (any(failed)) {
        nfailed <- sum(failed)
        cli_abort(c("x"=paste("{nfailed} out of {length(mapOutput)}",
                              "chunk{?s} in {.arg mapOutput} failed."),
                    "i"=paste("Resubmit them by calling {.fn gsvaMap} with",
                              "{.arg MAPREDO} set to {.arg mapOutput}.")))
    }

    cls <- unique(lapply(mapOutput, class))
    if (length(cls) > 1)
        cli_abort(c("x"="All inputs must be of the same class."))

    totalInputDim <- attributes(mapOutput)$totalInputDim
    if (is.null(totalInputDim))
        cli_abort(c("x"=paste("The input list argument in 'mapOutput' must",
                              "contain the attribute 'totalInputDim'.")))

    if (is.character(mapOutput[[1]]))
        mapOutput <- lapply(mapOutput, .load_gsva_path, assay="auto",
                            argname="mapOutput", verbose=FALSE)

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
                               stop.on.error=FALSE, log=TRUE, logdir=dir)
    return(BTPARAM)
} # nocov end


## private functions

#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
#' @importFrom IRanges IRanges start end
MAP_FUN_WRAPPER <- function(X, WRAPPED_FUN, output, ncpus, maxmem, ...) {
    fname <- NULL
    if (is(X, "gsvaMapChunk")) { ## output is saved to the file 'fname'
        fname <- X$fname
        X <- X$chunk
        ## completed by another job, such as one left running by a call to
        ## gsvaMap() whose R session ended
        if (file.exists(fname))
            return(fname)
    }
    rng <- X
    res <- whdim <- NULL
    parallelbackend <- SerialParam()
    if (ncpus > 1) {
        parallelbackend <- MulticoreParam(workers=ncpus)
        if (.Platform$OS.type != "unix")
            parallelbackend <- SnowParam(workers=ncpus)
    }

    if (is(X, "IRanges")) { ## input is splitted in chunks
        res <- WRAPPED_FUN(..., first=start(rng), last=end(rng),
                           verbose=FALSE,
                           BPPARAM=parallelbackend,
                           maxmem=maxmem)
    } else {                ## input was already splitted in chunks
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
        } ## if X is a path, restrict metadata is in the loaded object
    }
    if (output != "object") {
        ## save to a temporary name and rename it when complete, so that a
        ## job killed while saving does not leave a result that looks complete.
        ## the temporary name is unique to this job, because a job left
        ## running by a call to gsvaMap() whose R session ended may be saving
        ## the same chunk. renaming fails when another job has already saved
        ## an HDF5 directory with this chunk, in which case only the output of
        ## this job is discarded, while a Parquet file of another job is
        ## replaced by the identical one of this job
        tmpname <- .unique_tmpname(fname)
        if (output == "HDF5")
            saveHDF5GSVA(res, tmpname)
        else
            saveParquetGSVA(res, tmpname)
        if (!suppressWarnings(file.rename(tmpname, fname))) {
            unlink(tmpname, recursive=TRUE)
            if (!file.exists(fname))
                cli_abort(c("x"="cannot save results to {fname}."))
        }
        res <- fname
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
        if (!is.null(attributes(x)$gsvaVersion))
            nrmdata$gsvaVersion <- attributes(x)$gsvaVersion
        if (!is.null(attributes(x)$ranksNrow))
            nrmdata$ranksNrow <- attributes(x)$ranksNrow
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
            attributes(x)$gsvaVersion <- NULL
            attributes(x)$ranksNrow <- NULL
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
                            "in its resources. Using all available memory")))
        BTPARAM$resources$memory <- Inf
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

    if (is.list(inputData) && any(.map_failed(inputData)))
        cli_abort(c("x"=paste("Some chunks in {.arg inputData} failed."),
                    "i"=paste("Resubmit them by calling {.fn gsvaMap} with",
                              "{.arg MAPREDO} set to {.arg inputData}.")))
}

.map_fun_name <- function(FUN) {
    for (funname in c("gsvaRowNorm", "gsvaColRanks", "gsvaColScores"))
        if (identical(FUN, get(funname)))
            return(funname)
    cli_abort(c("x"="Internal error, invalid FUN argument."))
}

## which chunks of a gsvaMap() output failed: those whose element is a
## condition, either a 'bperror' object from bplapply(), or an error caught by
## .map_chunks() when running in this R process
.map_failed <- function(x) {
    vapply(x, function(el) inherits(el, "condition"), logical(1))
}

## map the chunks in 'X' with MAP_FUN_WRAPPER() and the arguments in 'args',
## returning failed chunks as condition objects instead of aborting. When
## results are saved to the files in 'paths', the chunks of jobs that did not
## deliver a result, such as jobs killed by the workload manager, are
## recovered from those files
#' @importFrom BiocParallel bplapply bpnworkers bptry
.map_chunks <- function(X, BTPARAM, args, paths=NULL, outdir=NULL) {
    if (bpnworkers(BTPARAM) > 1) {
        res <- tryCatch(bptry(do.call("bplapply",
                                      args=c(list(X=X, FUN=MAP_FUN_WRAPPER,
                                                  BPPARAM=BTPARAM), args))),
                        error=identity)
        ## bptry() returns, instead of a list of results, a 'bperror' raised
        ## for reasons other than failed chunks, such as a timeout
        if (inherits(res, "condition"))
            res <- .map_recover_from_disk(res, paths, outdir)
    } else ## mainly to be able to unit test this
        res <- lapply(X, function(x)
                          tryCatch(do.call("MAP_FUN_WRAPPER",
                                           args=c(list(X=x), args)),
                                   error=identity))

    res
}

## after the error 'err' raised by bplapply(), build the result of the chunks
## whose output is saved in the files in 'paths', from those files
#' @importFrom cli cli_abort
.map_recover_from_disk <- function(err, paths, outdir) {
    if (is.null(paths))
        stop(err)
    done <- file.exists(paths)
    if (!any(done)) {
        errmsg <- conditionMessage(err)
        cli_abort(c("x"="No chunk was completed: {errmsg}",
                    "i"=paste("Once the problem is fixed, resubmit the",
                              "calculations by calling {.fn gsvaMap} again",
                              "with the same {.arg FUN} and {.arg inputData},",
                              "and MAPREDO=\"{outdir}\".")),
                  parent=err)
    }
    errmsg <- paste("The job of this chunk did not save any result, which",
                    "happens, for instance, when the workload manager kills",
                    "it. The error reported was:", conditionMessage(err))
    res <- lapply(seq_along(paths), function(i) {
        if (done[i])
            paths[i]
        else
            simpleError(errmsg)
    })

    res
}

## chunk boundaries of the elements of a list output from a previous call to
## gsvaMap(), taken from their 'restrict' metadata or, when they are paths,
## from the output of gsvaMap(); NULL when they cannot be known without
## loading the elements from disk
.list_input_chunks <- function(inputData) {
    if (any(vapply(inputData, is.character, logical(1)))) {
        chunks <- attributes(inputData)$mapInfo$chunks
        if (!is.null(chunks) && nrow(chunks) != length(inputData))
            chunks <- NULL
        return(chunks)
    }
    rmdata <- lapply(inputData, .pull_restrict_metadata)
    if (any(vapply(rmdata, is.null, logical(1))))
        return(NULL)

    data.frame(first=as.integer(vapply(rmdata, function(x) x$first,
                                       numeric(1))),
               last=as.integer(vapply(rmdata, function(x) x$last, numeric(1))))
}

## identifier of the results of gsvaMap() saved to disk, given by its input
## and its output format
.map_runid <- function(funname, fingerprint, output) {
    step <- c(gsvaRowNorm="rnorm", gsvaColRanks="ranks",
              gsvaColScores="es")[funname]
    paste(substr(.md5_object(list(fingerprint, output)), 1L, 12L), step,
          sep="_")
}

## names of the files with the results of each chunk, given by the identifier
## of the run and the chunk boundaries or, when these are unknown, the chunk
## number
.map_file_names <- function(funname, fingerprint, output, chunks, nchunks) {
    runid <- .map_runid(funname, fingerprint, output)
    if (!is.null(chunks))
        files <- sprintf("%s_%d_%d", runid, chunks$first, chunks$last)
    else
        files <- sprintf("%s_chunk%d", runid, seq_len(nchunks))
    if (output == "Parquet")
        files <- paste0(files, ".parquet")

    files
}

.map_manifest_path <- function(mapinfo) {
    file.path(mapinfo$dir,
              paste0(.map_runid(mapinfo$FUN, mapinfo$fingerprint,
                                mapinfo$output), "_manifest.rds"))
}

## write the manifest of the results of gsvaMap() saved to disk, which allows
## MAPREDO to find them; the directory is not stored, to allow moving it
#' @importFrom cli cli_abort
#' @importFrom utils packageVersion
.write_map_manifest <- function(mapinfo, totalInputDim, mustNotExist) {
    fname <- .map_manifest_path(mapinfo)
    if (file.exists(fname)) {
        if (mustNotExist) {
            dir <- mapinfo$dir
            runid <- sub("_manifest\\.rds$", "", basename(fname))
            cli_abort(c("x"=paste("The directory {.file {dir}} has results of",
                                  "a previous call to {.fn gsvaMap} with the",
                                  "same {.arg FUN} and {.arg inputData}."),
                        "i"=paste("Set MAPREDO=\"{dir}\" to resubmit only",
                                  "the chunks without results, or delete the",
                                  "files starting with {.val {runid}} to",
                                  "compute them all again.")))
        }
        return(invisible(fname))
    }

    manifest <- mapinfo
    manifest$dir <- NULL
    manifest <- c(list(manifestVersion=1L), manifest,
                  list(totalInputDim=totalInputDim,
                       gsvaVersion=as.character(packageVersion("GSVA"))))
    tmpname <- .unique_tmpname(fname)
    saveRDS(manifest, tmpname)
    file.rename(tmpname, fname)

    invisible(fname)
}

## temporary name in the directory of 'fname', unique to this R process
## across compute nodes, because the name given by tempfile() only includes
## the process id and a random part that is the same in forked processes
.unique_tmpname <- function(fname) {
    tempfile(pattern=paste0(basename(fname), ".", Sys.info()[["nodename"]],
                            "."),
             tmpdir=dirname(fname), fileext=".partial")
}

## read the manifest in 'path', or the one in the directory 'path' that
## matches the given FUN, input fingerprint and, if not NULL, output, and
## build from it and the results saved to disk the output of gsvaMap()
#' @importFrom cli cli_abort
.read_map_manifest <- function(path, funname, fingerprint, output) {
    if (length(path) != 1L || is.na(path))
        cli_abort(c("x"=paste("{.arg MAPREDO} must be either the output of",
                              "a previous call to {.fn gsvaMap}, or a single",
                              "path to a directory or a manifest file with",
                              "results saved by {.fn gsvaMap}.")))
    if (dir.exists(path))
        fnames <- list.files(path, pattern="_manifest\\.rds$",
                             full.names=TRUE)
    else if (file.exists(path))
        fnames <- path
    else
        cli_abort(c("x"="{.arg MAPREDO} {.file {path}} cannot be found."))

    manifests <- lapply(fnames, function(f) tryCatch(readRDS(f),
                                                     error=function(e) NULL))
    match <- vapply(manifests, function(m) {
        is.list(m) && !is.null(m$manifestVersion) &&
            identical(m$FUN, funname) &&
            identical(m$fingerprint, fingerprint) &&
            (is.null(output) || identical(m$output, output))
    }, logical(1))
    if (!any(match))
        cli_abort(c("x"=paste("No manifest of {.fn gsvaMap} results in",
                              "{.file {path}} matches FUN={funname} and",
                              "the given {.arg inputData}",
                              if (!is.null(output)) "and {.arg output}",
                              "")))
    if (sum(match) > 1)
        cli_abort(c("x"=paste("Several manifests of {.fn gsvaMap} results",
                              "in {.file {path}} match FUN={funname} and",
                              "the given {.arg inputData}."),
                    "i"=paste("Set the {.arg output} argument or give the",
                              "path to one of these manifests:",
                              "{.file {basename(fnames[match])}}.")))

    m <- manifests[[which(match)]]
    dir <- dirname(normalizePath(fnames[match]))
    paths <- file.path(dir, m$files)
    res <- lapply(paths, function(p) {
        if (file.exists(p))
            p
        else
            simpleError(sprintf("No result was saved in %s", p))
    })
    attributes(res) <- list(totalInputDim=m$totalInputDim,
                            mapInfo=list(FUN=m$FUN, output=m$output,
                                         fingerprint=m$fingerprint,
                                         chunks=m$chunks, files=m$files,
                                         dir=dir))

    res
}

## when the registry directory of 'BTPARAM' exists, which happens when the R
## session of a previous call to gsvaMap() ended before removing it, set it
## to a new directory, returning the previous one to restore it afterwards
#' @importFrom BiocParallel bpisup
#' @importFrom cli cli_alert_warning
.free_registry_dir <- function(BTPARAM, verbose) {
    regdir <- BTPARAM$registryargs$file.dir
    if (is.null(regdir) || is.na(regdir) || !dir.exists(regdir) ||
        bpisup(BTPARAM) || isTRUE(BTPARAM$saveregistry))
        return(NULL)

    k <- 1L
    while (dir.exists(newdir <- paste0(regdir, "-", k)))
        k <- k + 1L
    BTPARAM$registryargs$file.dir <- newdir
    if (verbose)
        cli_alert_warning(paste("The registry directory {.file {regdir}}",
                                "already exists, probably left by a previous",
                                "call to {.fn gsvaMap} whose R session ended.",
                                "Using {.file {newdir}} instead. Delete",
                                "{.file {regdir}} once its jobs are not",
                                "running anymore."))

    regdir
}

## default value of the 'output' argument of gsvaMap(), which saves results
## to disk when jobs are sent to a workload manager, enabling MAPREDO to
## recover them
.default_map_output <- function(BTPARAM) {
    output <- "object"
    if (!BTPARAM$cluster %in% c("socket", "multicore", "interactive"))
        output <- "HDF5"

    output
}

#' @importFrom cli cli_abort
.check_MAPREDO <- function(MAPREDO, funname) {
    mapinfo <- attributes(MAPREDO)$mapInfo
    if (!is.list(MAPREDO) || is.null(mapinfo))
        cli_abort(c("x"=paste("{.arg MAPREDO} must be the output of a",
                              "previous call to {.fn gsvaMap}.")))
    if (mapinfo$FUN != funname)
        cli_abort(c("x"=paste("{.arg MAPREDO} was produced with",
                              "FUN={mapinfo$FUN}, which differs from",
                              "FUN={funname}.")))

    mapinfo
}

## fingerprint identifying the input of a gsvaMap() call, built from the
## dimensions, dimension names and GSVA parameters of the input data, without
## its values, which can be too large to read, or from the chunk boundaries or
## paths of a list output from a previous call to gsvaMap()
.map_fingerprint <- function(funname, inputData, totalInputDim) {
    if (is.list(inputData)) {
        elems <- lapply(inputData, function(x) {
            if (is.character(x)) ## allow moving the directory with the files
                basename(x)
            else
                .pull_restrict_metadata(x)[c("first", "last", "whdim")]
        })
        fp <- list(funname, totalInputDim, elems)
    } else {
        param <- inputData
        if (!is(inputData, "GsvaMethodParam"))
            param <- .pull_param(inputData)
        exprData <- get_exprData(param)
        fp <- list(funname, class(inputData), dim(exprData),
                   dimnames(exprData), .gsvaParam_as_list(param))
    }

    .md5_object(fp)
}

## MD5 hash of an R object, serialized without its header, which stores the
## version of R that serialized it, and with serialization format version 2,
## which writes ALTREP vectors in the same way as regular ones
#' @importFrom tools md5sum
.md5_object <- function(x) {
    bytes <- serialize(x, connection=NULL, version=2L)
    fname <- tempfile()
    on.exit(unlink(fname))
    writeBin(bytes[-seq_len(14L)], fname)
    unname(md5sum(fname))
}
