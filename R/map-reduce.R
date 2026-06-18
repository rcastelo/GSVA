#' @title MapReduce parallelization for HPC environments
#'
#' @description The functions `gsvaMap()` and `gsvaReduce()` allow one
#' to run GSVA calculations across multiple compute nodes in a high-performance
#' computing (HPC) environment. The `gsvaMap()` function is used to launch the
#' GSVA calculations on each compute node, while the `gsvaReduce()` function is
#' used to combine the results from all nodes into a single object.
#'
#' @param FUN Function to map to the data in the `inputData` argument.
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
#' @param ... In `gsvaReduce()`, the output of `gsvaMap()`.
#'
#' @param verbose Gives information about the progress of the calculations.
#' Default: `TRUE`.
#'
#' @param BTPARAM An object of class
#' [`BatchtoolsParam`][BiocParallel::BatchtoolsParam-class] specifying parameters
#' for parallel execution in an HPC enviroment.
#'
#' @return The `gsvaMap()` function returns either a list of objects with the
#' results of the GSVA calculations for each compute node, or a list of file paths
#' where the results are saved. The `gsvaReduce()` function returns a single
#' object that combines the results from all compute nodes.
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
#' @importFrom BiocParallel BatchtoolsParam bpnworkers MulticoreParam bplapply
#' @rdname map-reduce
#' @export gsvaMap
gsvaMap <- function(FUN, inputData, returnPath=FALSE, verbose=TRUE,
                    BTPARAM=BatchtoolsParam(workers=2, progressbar=verbose)) {

    FUN <- match.fun(FUN)

    if (!identical(FUN, gsvaRowNorm) && !identical(FUN, gsvaColRanks) &&
        !identical(FUN, gsvaColScores)) {
        msg <- paste("'FUN' must be one of 'gsvaRowNorm',",
                     "'gsvaColRanks' or 'gsvaColScores'.")
        cli_abort(c("x"=msg))
    }

    if (!is(inputData, "gsvaParam") &&
        !is(inputData, "GsvaExprData") && !is.list(inputData)) {
        msg <- paste("'inputData' must be an object of class
                     'gsvaParam', 'GsvaExprData' or 'list'.")
        cli_abort(c("x"=msg))
    } else if (is.list(inputData))
        cli_abort(c("x"="Not implemented yet."))

    BTPARAM <- .check_batchtools_param(BTPARAM, verbose)

    nworkers <- bpnworkers(BTPARAM)
    ncpus <- BTPARAM$resources$ncpus
    maxmem <- .memtext2bytes(BTPARAM$resources$memory)

    if (identical(FUN, gsvaRowNorm) && !is(inputData, "gsvaParam"))
        cli_abort(c("x"=paste("FUN=gsvaRowNorm requires inputData of",
                              "class 'gsvaParam'.")))

    if ((identical(FUN, gsvaColRanks) || identical(FUN, gsvaColScores)) &&
        !is(inputData, "GsvaExprData"))
        cli_abort(c("x"=paste("FUN=gsvaColRanks or FUN=gsvaColScores",
                              "requires inputData of class 'GsvaExprData'.")))

    FUN_WRAPPER <- function(rng, WRAPPED_FUN, path2save, ncpus, maxmem, ...) {
       res <- WRAPPED_FUN(..., first=start(rng), last=end(rng), verbose=FALSE,
                          BPPARAM=MulticoreParam(workers=ncpus), maxmem=maxmem)
       if (nchar(path2save) > 0) {
           fname <- file.path(path2save, sprintf("%s_%d_%d",
                                                 basename(tempfile()),
                                                 start(rng), end(rng)))
           if (dir.exists(fname))
               cli_abort(c("x"=paste("cannot save results to {fname} because",
                                     "it already exists.")))
           res <- saveHDF5GSVA(res, fname)
       }
       return(res)
    }

    gridsizefun <- .colgridsize
    splitinrangesfun <- .splitColsInRanges

    funargs <- list()

    if (identical(FUN, gsvaRowNorm)) {

        funargs <- c(funargs, list(param=inputData,
                                   dropExistingAssays=TRUE,
                                   errorOnTooFewRows=FALSE))
        gridsizefun <- .rowgridsize
        splitinrangesfun <- .splitRowsInRanges

    } else if (identical(FUN, gsvaColRanks)) {
        if (!"gsvarnorm" %in% gsvaAssayNames(inputData))
            cli_abort(c("x"=paste("FUN=gsvaColRanks requires inputData with",
                                  "row-normalized expression values.")))

        funargs <- c(funargs, list(rowNormExprData=inputData,
                                   dropExistingAssays=TRUE))

    } else if (identical(FUN, gsvaColScores)) {
        if (!"gsvaranks" %in% gsvaAssayNames(inputData))
            cli_abort(c("x"=paste("FUN=gsvaColScores requires inputData with",
                                  "column rank values.")))

        funargs <- c(funargs, list(rankExprData=inputData))

    } else
        cli_abort(c("x"="Internal error, invalid FUN argument."))

    path2save <- ""
    if (returnPath)
        path2save <- path.expand(BTPARAM$registryargs$work.dir)

    X <- unwrapData(get_exprData(inputData))
    grid <- gridsizefun(X, nworkers, maxmem)
    rcir <- splitinrangesfun(grid)

    do.call("bplapply", args=c(list(X=rcir, FUN=FUN_WRAPPER, WRAPPED_FUN=FUN,
                                    path2save=path2save, ncpus=ncpus,
                                    maxmem=maxmem, BPPARAM=BTPARAM), funargs))

}

#' @importFrom cli cli_abort
#' @importFrom BiocGenerics rbind cbind
#' @rdname map-reduce
#' @export gsvaReduce
gsvaReduce <- function(..., verbose=TRUE) {
    args <- list(...)
    if (length(args) == 1 && is.list(args[[1]]))
        args <- args[[1]]

    cls <- unique(lapply(args, class))
    if (length(cls) > 1)
        cli_abort(c("x"="All inputs must be of the same class."))

    if (is.character(args[[1]])) {
        args <- lapply(args, function(x) {
            if (!dir.exists(x))
                cli_abort(c("x"="Cannot find {x}."))
            loadHDF5GSVA(x)
        })
    }

    param <- .pull_param(args[[1]])
    nrmdata <- .pull_nonrestrict_metadata(args[[1]])
    rmdata <- .pull_restrict_metadata(args)
    ord <- .check_and_order_restrict_metadata(rmdata)
    args <- .strip_metadata(args)
    bfun <- "rbind"
    if (rmdata[[1]]$whdim == 2)
        bfun <- "cbind"
    res <- do.call(bfun, args[ord])
    res <- .add_metadata(res, param, nrmdata)

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

#' @importFrom cli cli_abort
.pull_restrict_metadata <- function(args) {
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
.check_and_order_restrict_metadata <- function(rmdata) {
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
    if (first[1] != 1 || any(first[-1] != (last+1)[-length(last)]))
        cli_abort(c("x"=paste("Input restrict metadata must contain elements",
			      "'first' and 'last' starting at position 1,",
                              "contiguous and non-overlapping.")))

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
