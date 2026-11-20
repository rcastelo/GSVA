## ----- private methods for wrapping/unwrapping assay data in containers -----

## unwrapData: extract a data matrix from a container object
setMethod("unwrapData", signature("matrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("dgCMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("SVT_SparseMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("DelayedMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("ExpressionSet"),
          function(container, assay) {
              return(exprs(container))
          })

#' @importFrom cli cli_abort
setMethod("unwrapData", signature("SummarizedExperiment"),
          function(container, assay) {
              assay <- .check_unwrapping_assay(container,
                                               "SummarizedExperiment", assay)

              return(.unwrap_assay(container, assay))
          })

#' @importFrom cli cli_abort
setMethod("unwrapData", signature("SingleCellExperiment"),
          function(container, assay) {
              assay <- .check_unwrapping_assay(container,
                                               "SingleCellExperiment", assay)

              return(.unwrap_assay(container, assay))
          })

setMethod("unwrapData", signature("SpatialExperiment"),
          function(container, assay) {
              assay <- .check_unwrapping_assay(container,
                                               "SpatialExperiment", assay)

              return(.unwrap_assay(container, assay))
          })

## the assay 'assay' of the SummarizedExperiment object 'container', with
## the dimension names of 'container', as assays(container)[[assay]] returns
## it, but without setting the dimension names of the other assays, since
## setting them on a DelayedArray object validates it by extracting an empty
## array from its seeds, which opens every HDF5 file bound by gsvaReduce(),
## and this can take seconds on the shared filesystems of clusters
.unwrap_assay <- function(container, assay) {
    a <- assays(container, withDimnames=FALSE)[[assay]]
    x_dimnames <- dimnames(container)
    if (is.null(x_dimnames))
        x_dimnames <- list(NULL, NULL)
    a_dimnames <- dimnames(a)
    a_dimnames[1:2] <- x_dimnames
    if (all(vapply(a_dimnames, is.null, logical(1))))
        a_dimnames <- NULL
    if (!identical(dimnames(a), a_dimnames))
        dimnames(a) <- a_dimnames
    a
}

## metadata recorded in every GSVA output object: the version of GSVA that
## produced it, 'gsvaVersion', to handle in the future objects produced by
## different versions of GSVA, e.g., when they are stored on disk, and for
## column ranks, the number of rows on which they were calculated,
## 'ranksNrow', since the ranks are only valid for those rows; for other
## outputs 'ranksNrow' is NULL, which removes it when inherited from the input
#' @importFrom utils packageDescription
.gsva_output_metadata <- function(dataMatrix, assay) {
    list(gsvaVersion=packageDescription("GSVA")[["Version"]],
         ranksNrow=if (assay == "gsvaranks") nrow(dataMatrix) else NULL)
}

.wrapdata_nonSE <- function(dataMatrix, param, assay, first, last, rem, whdim,
                            dropAssays, geneSets) {
    stopifnot(!missing(param))
    stopifnot(!missing(assay))
    stopifnot(!missing(first))
    stopifnot(!missing(last))
    stopifnot(!missing(rem))
    stopifnot(!missing(whdim))
    stopifnot(!missing(dropAssays))
    omd <- .gsva_output_metadata(dataMatrix, assay)
    attr(dataMatrix, "gsvaParam") <- .gsvaParam_as_list(param)
    attr(dataMatrix, "assay") <- assay
    attr(dataMatrix, "gsvaVersion") <- omd$gsvaVersion
    attr(dataMatrix, "ranksNrow") <- omd$ranksNrow
    if (!is.na(first) || !is.na(last))
        attr(dataMatrix, "restrict") <- list(first=first, last=last,
                                             rem=rem, whdim=whdim)
    if (!missing(geneSets))
        attr(dataMatrix, "geneSets") <- geneSets
    return(dataMatrix)
}

## wrapData: put the resulting data and gene sets into the original data
##           container type
## parameters:
##  - container is the original data container object
##  - dataMatrix is a matrix of values to be wrapped into the container
##  - param is the parameter object used for the calculation that produced
##    dataMatrix
##  - assay is the name of the assay to be used in the container object to
##    store dataMatrix
##  - first if not NA, it is the index of the first row/column in the original
##    data container that was used to produce dataMatrix
##  - last if not NA, it is the index of the last row/column in the original
##    data container that was used to produce dataMatrix
##  - rem if not NA, it is the number of rows/columns that were discarded from
##    the original data container, e.g., when filtering out constant rows
##  - whdim if not NA, it is the dimension along which first and last are
##    defined, i.e., 1 for rows and 2 for columns
##  - dropAssays if TRUE and container uses assays, the original assays in the
##    container are dropped and only the new assay is kept. If FALSE, the new
##    assay is added to the existing assays
##  - geneSets if not missing, it is a list of gene sets that were used to
##    produce dataMatrix.
setMethod("wrapData", signature(container="matrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="dgCMatrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="SVT_SparseMatrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="DelayedMatrix"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              .wrapdata_nonSE(dataMatrix, param, assay, first, last, rem, whdim,
                              dropAssays, geneSets)
          })

setMethod("wrapData", signature(container="ExpressionSet"),
          function(container, dataMatrix, param, assay, first, last, rem, whdim,
                   dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              pdata <- phenoData(container)
              if (!is.na(first) && !is.na(last) && whdim == 2)
                  pdata <- pdata[seq.int(first, last, 1), , drop=FALSE]
              rval <- new("ExpressionSet", exprs=dataMatrix,
                          phenoData=pdata,
                          experimentData=experimentData(container),
                          annotation="")
              omd <- .gsva_output_metadata(dataMatrix, assay)
              attr(rval, "gsvaParam") <- .gsvaParam_as_list(param)
              attr(rval, "assay") <- assay
              attr(rval, "gsvaVersion") <- omd$gsvaVersion
              attr(rval, "ranksNrow") <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  attr(rval, "restrict") <- list(first=first, last=last,
                                                 rem=rem, whdim=whdim)
              if (!missing(geneSets))
                  attr(rval, "geneSets") <- geneSets
              
              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment assayNames
setMethod("wrapData", signature(container="SummarizedExperiment"),
          function(container, dataMatrix, param, assay, first, last, rem,
                   whdim, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## missing geneSets implies adding an assay
                  if (assay %in% assayNames(container))
                      cli_abort(c("x"=paste("Assay {assay} already exists in",
                                            "the input container object.")))
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays) {
                      if (!is.na(first) && !is.na(last) && whdim == 2)
                          adata <- c(assays(container[mask, seq.int(first, last, 1)]), adata)
                      else
                          adata <- c(assays(container[mask, ]), adata)
                  }
                  rdata <- rowData(container)[mask, ]
              }
              cdata <- colData(container)
              if (!is.na(first) && !is.na(last) && whdim == 2)
                  cdata <- cdata[seq.int(first, last, 1), , drop=FALSE]
              rval <- SummarizedExperiment(assays=adata, colData=cdata,
                                           rowData=rdata,
                                           metadata=metadata(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              omd <- .gsva_output_metadata(dataMatrix, assay)
              metadata(rval)$gsvaVersion <- omd$gsvaVersion
              metadata(rval)$ranksNrow <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  metadata(rval)$restrict <- list(first=first, last=last,
                                                  rem=rem, whdim=whdim)
              if (!missing(geneSets)) ## row data has been replaced
                  metadata(rval)$annotation <- NULL

              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims altExps
setMethod("wrapData", signature(container="SingleCellExperiment"),
          function(container, dataMatrix, param, assay, first, last, rem,
                   whdim, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              rdimdata <- reducedDims(container)
              aexpsdata <- altExps(container)
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
                  rdimdata <- aexpsdata <- SimpleList()
              } else { ## missing geneSets implies adding an assay
                  if (assay %in% assayNames(container))
                      cli_abort(c("x"=paste("Assay {assay} already exists in",
                                            "the input container object.")))
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays) {
                      if (!is.na(first) && !is.na(last) && whdim == 2)
                          adata <- c(assays(container[mask, seq.int(first, last, 1)]), adata)
                      else
                          adata <- c(assays(container[mask, ]), adata)
                  }
                  rdata <- rowData(container)[mask, ]
              }
              ## the column data, reduced dimensions and alternative
              ## experiments of a chunk of columns are those of its columns
              cdata <- colData(container)
              if (!is.na(first) && !is.na(last) && whdim == 2) {
                  first_to_last <- seq.int(first, last, 1)
                  cdata <- cdata[first_to_last, , drop=FALSE]
                  rdimdata <- lapply(rdimdata,
                                     function(x) x[first_to_last, , drop=FALSE])
                  aexpsdata <- lapply(aexpsdata,
                                      function(x) x[, first_to_last])
              }
              rval <- SingleCellExperiment(assays=adata, colData=cdata,
                                           rowData=rdata, reducedDims=rdimdata,
                                           altExps=aexpsdata,
                                           metadata=metadata(container))
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              omd <- .gsva_output_metadata(dataMatrix, assay)
              metadata(rval)$gsvaVersion <- omd$gsvaVersion
              metadata(rval)$ranksNrow <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  metadata(rval)$restrict <- list(first=first, last=last,
                                                  rem=rem, whdim=whdim)
              if (!missing(geneSets)) ## row data has been replaced
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })

#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims altExps
#' @importFrom SpatialExperiment SpatialExperiment
setMethod("wrapData", signature(container="SpatialExperiment"),
          function(container, dataMatrix, param, assay, first, last, rem,
                   whdim, dropAssays, geneSets) {
              stopifnot(!missing(param))
              stopifnot(!missing(assay))
              stopifnot(!missing(first))
              stopifnot(!missing(last))
              stopifnot(!missing(rem))
              stopifnot(!missing(whdim))
              stopifnot(!missing(dropAssays))
              rdata <- NULL
              adata <- SimpleList(dataMatrix)
              names(adata) <- assay
              if (!missing(geneSets)) { ## storing enrichment scores only
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## missing geneSets implies adding an assay
                  if (assay %in% assayNames(container))
                      cli_abort(c("x"=paste("Assay {assay} already exists in",
                                            "the input container object.")))
                  stopifnot(all(rownames(dataMatrix) %in% rownames(container)))
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  if (!dropAssays) {
                      if (!is.na(first) && !is.na(last) && whdim == 2)
                          adata <- c(assays(container[mask, seq.int(first, last, 1)]), adata)
                      else
                          adata <- c(assays(container[mask, ]), adata)
                  }
                  rdata <- rowData(container)[mask, ]
              }
              ## the column data, reduced dimensions, alternative experiments
              ## and spatial coordinates of a chunk of columns are those of
              ## its columns
              cdata <- colData(container)
              rdimdata <- reducedDims(container)
              aexpsdata <- altExps(container)
              coords <- spatialCoords(container)
              if (!is.na(first) && !is.na(last) && whdim == 2) {
                  first_to_last <- seq.int(first, last, 1)
                  cdata <- cdata[first_to_last, , drop=FALSE]
                  rdimdata <- lapply(rdimdata,
                                     function(x) x[first_to_last, , drop=FALSE])
                  aexpsdata <- lapply(aexpsdata,
                                      function(x) x[, first_to_last])
                  coords <- coords[first_to_last, , drop=FALSE]
              }
              rval <- SpatialExperiment(
                  assays=adata,
                  colData=cdata,
                  rowData=rdata,
                  reducedDims=rdimdata,
                  altExps=aexpsdata,
                  metadata=metadata(container),
                  imgData=imgData(container),
                  spatialCoords=coords)
              metadata(rval)$gsvaParam <- .gsvaParam_as_list(param)
              omd <- .gsva_output_metadata(dataMatrix, assay)
              metadata(rval)$gsvaVersion <- omd$gsvaVersion
              metadata(rval)$ranksNrow <- omd$ranksNrow
              if (!is.na(first) || !is.na(last))
                  metadata(rval)$restrict <- list(first=first, last=last,
                                                  rem=rem, whdim=whdim)
              if (!missing(geneSets)) ## row data has been replaced
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })

## return direct subclasse of a union class, e.g., 'GsvaExprData'
## or 'GsvaGeneSets'

## obtain direct subclasses of a union class for error reporting purposes
#' @importFrom methods isClass getClass is
#' @importFrom cli cli_abort
.getDirectSubclasses <- function(unionClass) {
    if (!isClass(unionClass))
        cli_abort(c("x"=sprintf("'%s' is not a valid union class name",
                                unionClass)))

    subClasses <- getClass(unionClass)@subclasses
    dist <- vapply(subClasses, function(x) x@distance, numeric(1))
    directSubclasses <- names(dist)[dist == 1]
    directSubclasses
}

.check_input_expr_gene_sets <- function(exprData, geneSets) {
    if (!is(exprData, "GsvaExprData")) {
        subclassnames <- .getDirectSubclasses("GsvaExprData")
        msg <- sprintf(paste("argument 'exprData' must be an object of one of",
                            "the following classes: %s"),
                       paste(subclassnames, collapse=", "))
        cli_abort(c("x"=msg))
    }
    if (!is(geneSets, "GsvaGeneSets")) {
        subclassnames <- .getDirectSubclasses("GsvaGeneSets")
        msg <- sprintf(paste("argument 'geneSets' must be an object of one of",
                            "the following classes: %s"),
                       paste(subclassnames, collapse=", "))
        cli_abort(c("x"=msg))
    }
    if (is(geneSets, "list")) {
        element_types <- vapply(geneSets, typeof, character(1))
        if (length(unique(element_types)) > 1) {
            msg <- paste("All elements of the 'geneSets' list must be of the",
                         "same type (either character or numeric)")
            cli_abort(c("x"=msg))
        }
    }
}

## remove duplicated genes within gene sets given as a list of character or
## numeric vectors and warn about it; 'GeneSetCollection' objects need not be
## checked because 'GSEABase::GeneSet()' does not allow duplicated gene IDs
#' @importFrom cli cli_alert_warning
#' @importFrom methods is
#' @importFrom IRanges IntegerList CharacterList
#' @importFrom BiocGenerics anyDuplicated
.remove_duplicated_genes <- function(geneSets) {
    if (is(geneSets, "GeneSetCollection"))
        return(geneSets)

    geneSets2 <- geneSets
    if (is.numeric(geneSets[[1]]))
        geneSets2 <- IntegerList(geneSets)
    else if (is.character(geneSets[[1]]))
        geneSets2 <- CharacterList(geneSets)
    else {
        msg <- paste("When argument 'geneSets' is a 'list' object it can",
                     "only contain either character or numeric vectors.")
        cli_abort(c("x"=msg))
    }

    dupmask <- anyDuplicated(geneSets2)
    if (any(dupmask)) {
        geneSets2[dupmask] <- unique(geneSets2[dupmask])
        gsnames <- names(geneSets2)[dupmask]
        if (is.null(gsnames))
            gsnames <- as.character(which(dupmask))
        msg <- paste("Duplicated genes removed from {sum(dupmask)} gene",
                     "set{?s}: {.val {gsnames}}")
        cli_alert_warning(msg)
    }

    as.list(geneSets2)
}

## check that 'BPPARAM' is a BiocParallelParam object and, unless a function
## running several steps already did it, see gsva(), whether its workers
## exceed the CPU cores available to this R process, see .check_cores()
.check_bpparam <- function(BPPARAM) {
    if (!is(BPPARAM, "BiocParallelParam")) {
        msg <- paste("Argument 'BPPARAM' must be a 'BiocParallelParam'",
                     "derivative. Please consult the BiocParallel package.")
        cli_abort(c("x"=msg))
    }
    if (gsva_global$check_cores)
        .check_cores(BPPARAM)
}

## warn when the parallel workers of 'BPPARAM' running on this machine, those
## of a MulticoreParam object or of a SnowParam object with a number of
## workers, exceed the 'cores' CPU cores available to this R process, see
## .available_cores(), because workers sharing cores make the calculations
## slower
#' @importFrom BiocParallel bpnworkers bpworkers SnowParam MulticoreParam
#' @importFrom cli cli_warn
.check_cores <- function(BPPARAM, cores=.available_cores()) {
    local <- is(BPPARAM, "MulticoreParam") ||
             (is(BPPARAM, "SnowParam") && is.numeric(bpworkers(BPPARAM)))
    nworkers <- bpnworkers(BPPARAM)
    if (!local || is.na(cores) || nworkers <= cores)
        return(invisible(FALSE))

    cli_warn(c("!"=paste("The {nworkers} parallel workers of {.arg BPPARAM}",
                         "exceed the {cores} CPU core{?s} available to this",
                         "R process, which may make the calculations",
                         "slower."),
               "i"="Consider using at most {cores} worker{?s}."))

    invisible(TRUE)
}

## number of CPU cores available to this R process, or NA when it is not
## known, given by availableCores() of the parallelly package, when it is
## installed, which takes into account the cores allocated by workload
## managers such as SLURM, SGE, PBS or LSF, also when they do not restrict
## the cores on which the process can run, those of containers, and the limits
## set with the option 'mc.cores' or by 'R CMD check', and otherwise by
## .proc_available_cores()
.available_cores <- function() {
    if (requireNamespace("parallelly", quietly=TRUE))
        return(as.integer(parallelly::availableCores()))

    .proc_available_cores()
}

## number of CPU cores available to this R process, or NA when it is not
## known: those of the machine, see detectCores(), restricted on Linux to the
## cores on which the process can run, given by its CPU affinity in the field
## 'Cpus_allowed_list' of 'statusfile', which workload managers such as SLURM
## set to the cores of the job, and to the CPU quota of its control group
## (cgroup), set by container engines, in any of the cgroups along its path
## in 'cgroupfile', using either cgroup version 1 or 2, as .job_memory_limit()
## does with the memory limit. the paths of the files are arguments to allow
## testing this function
#' @importFrom parallel detectCores
.proc_available_cores <- function(statusfile="/proc/self/status",
                                  cgroupfile="/proc/self/cgroup",
                                  cgroupdir="/sys/fs/cgroup") {
    cores <- suppressWarnings(detectCores())
    readfirst <- function(fname) {
        if (!file.exists(fname))
            return(NA_character_)
        tryCatch(readLines(fname, n=1L, warn=FALSE)[1],
                 error=function(e) NA_character_)
    }
    ## the CPU affinity, e.g., 'Cpus_allowed_list:  0-3,8,10-11'
    status <- if (file.exists(statusfile))
                  tryCatch(readLines(statusfile, warn=FALSE),
                           error=function(e) character(0))
              else character(0)
    allowed <- grep("^Cpus_allowed_list:", status, value=TRUE)
    if (length(allowed) == 1L) {
        rngs <- strsplit(trimws(sub("^Cpus_allowed_list:", "", allowed)),
                         ",", fixed=TRUE)[[1]]
        n <- sum(vapply(rngs, function(r) {
            b <- suppressWarnings(as.numeric(strsplit(r, "-",
                                                      fixed=TRUE)[[1]]))
            if (length(b) == 1L) 1 else b[2] - b[1] + 1
        }, numeric(1)))
        if (!is.na(n) && n > 0)
            cores <- min(cores, n, na.rm=TRUE)
    }
    ## the CPU quota of the cgroups, 'quota period' in cpu.max of version 2,
    ## where the quota is 'max' when there is none, or the quota in
    ## cpu.cfs_quota_us, -1 when there is none, and the period in
    ## cpu.cfs_period_us of version 1
    cglines <- if (file.exists(cgroupfile))
                   tryCatch(readLines(cgroupfile, warn=FALSE),
                            error=function(e) character(0))
               else character(0)
    for (cgline in cglines) {
        fields <- strsplit(cgline, ":", fixed=TRUE)[[1]]
        if (length(fields) < 3L)
            next
        v2 <- fields[1] == "0" && fields[2] == ""
        if (!v2 && !"cpu" %in% strsplit(fields[2], ",", fixed=TRUE)[[1]])
            next
        dir <- if (v2) cgroupdir else file.path(cgroupdir, "cpu")
        path <- strsplit(paste(fields[-seq_len(2)], collapse=":"), "/",
                         fixed=TRUE)[[1]]
        path <- path[path != ""]
        for (k in rev(seq_len(length(path) + 1L) - 1L)) {
            d <- do.call(file.path, as.list(c(dir, path[seq_len(k)])))
            if (v2) {
                qp <- strsplit(readfirst(file.path(d, "cpu.max")), " ",
                               fixed=TRUE)[[1]]
                quota <- suppressWarnings(as.numeric(qp[1]))
                period <- suppressWarnings(as.numeric(qp[2]))
            } else {
                quota <- suppressWarnings(as.numeric(
                             readfirst(file.path(d, "cpu.cfs_quota_us"))))
                period <- suppressWarnings(as.numeric(
                              readfirst(file.path(d, "cpu.cfs_period_us"))))
            }
            if (!is.na(quota) && !is.na(period) && quota > 0 && period > 0)
                cores <- min(cores, max(1, ceiling(quota / period)),
                             na.rm=TRUE)
        }
    }

    as.integer(cores)
}

#' @importFrom BiocParallel bpnworkers
#' @importFrom cli cli_alert_info
.check_open_parallelism <- function(expr, BPPARAM, minparrows, minparcols,
                                    verbose) {
    if (bpnworkers(BPPARAM) > 1 && nrow(expr) > minparrows &&
        ncol(expr) > minparcols) {
        if (verbose) {
            msg <- sprintf("Using a %s parallel back-end with %d workers",
                           class(BPPARAM), bpnworkers(BPPARAM))
                      cli_alert_info(msg)
        }
    } else
        BPPARAM <- NULL

    BPPARAM
}

## generate dummy names, e.g. row/col names for object M that knows 'nrow()'
.dummyNames <- function(M, n=nrow(M), prefix="row") {
    fmt <- sprintf("%s%%0%dd", prefix, floor(log10(n)) + 1)
    sprintf(fmt, seq_len(n))
}

## check for presence of valid row/feature names
## and abort or generate dummy names
.check_rowNames <- function(expr, useDummyNames=TRUE, verbose) {
    ## CHECK: is this the right place to check this?
    ## 21/10/24: let's do it at parameter constructor
    if (is.null(rownames(expr))) {
        if (useDummyNames) {
            if (verbose) {
                cli_alert_info("Using dummy rownames for the input assay object.")
                rownames(expr) <- .dummyNames(expr)
            }
        } else {
            cli_abort(c("x"="The input assay object doesn't have rownames"))
        }
    } else if (anyDuplicated(rownames(expr))) {
        cli_abort(c("x"="The input assay object has duplicated rownames"))
    }

    return(expr)
}

## check if the expression data object at hand is a multi-assay container
## required by .check_assayNames() (see below) but may be re-used elsewhere
## *currently* only TRUE for SummarizedExperiment and its subclasses, else FALSE
.isMultiAssayContainer <- function(xd) {
    return(is(xd, "SummarizedExperiment"))
}

## check assay parameter and assay names:
##   abort if selected assay name not found in existing assay name list
##   alert if no assay name selected while there are assay names (and select
##     the first one)
##   alert if assay name selected while there are none
##
##   tricky: no assay names AND no assay name selected is the normal case for
##     non-container data objects, e.g. matrix,  or single-assay containers,
##     e.g. ExpressionSet, and should hence be silently accepted.  Multi-assay
##     containers, e.g. SummarizedExperiment, MAY contain more than one assays
##     BUT no assay names.  In this case we also abort because our parameter
##     objects can only store assay names and we prefer requiring assay names
##     for multi-assay containers (which is not unreasonable!) over making
##     things even more complicated for little practical gain.
##
## 2025-03-12  axel: as an afterthought, if we have a multi-assay container with
##   assay names AND no assay is selected AND one of the assay names happens to
##   be 'logcounts' --> use this one by default rather than the first in list.
#' @importFrom utils head
.check_assayNames <- function(a, xd, verbose) {
    an <- gsvaAssayNames(xd)

    if(length(a) != 1) {
        msg <- "argument 'assay' must be of length 1 (it is {length(a)})"
        cli_abort(c("x"=msg))
    }
    
    if(.isCharNonEmpty(an)) {   # we have assay names
        an <- .omitEmptyChar(an)
        
        if(is.na(a)) {          # but none selected: by default, 
            ## select a common value by default if available -- see afterthought
            ## if unavailable, just select the first available assay name
            def <- grep("logcounts", an, fixed=TRUE, value=TRUE)
            assay <- if(length(def) > 0) head(def, 1) else head(an, 1)
            if (verbose) {
                msg <- "No assay name provided; using default assay '{assay}'"
                cli_alert_info(msg)
            }
        } else {                # check the provided assay name before using it
            if(a %in% an) {
                assay <- a      # found it: OK!
            } else {            # assay name provided but not found: ERROR
                msg <- paste("invalid argument assay='{a}': not part of",
                             "assay names in input argument 'exprData'.")
                cli_abort(c("x"=msg))
            }
        }
    } else {                    # we don't have no assay names at all
        if(.isMultiAssayContainer(xd)) {  # these must have assay names: ERROR
            msg <- "exprData object of class '{class(xd)}' has no assay names."
            cli_abort(c("x"=msg))
        } else {                       # i.e. there is exactly one unnamed assay
            if(verbose && !is.na(a)) { # and the provided name is useless but harmless
                msg <- paste("argument assay='{a}' ignored since input argument",
                             "'exprData' has no assay names.")
                cli_alert_info(msg)
            }

            assay <- NA_character_
        }
    }

    return(assay)
}

## the assays are taken without their dimension names, because assays()
## sets them on every assay, and setting them on a DelayedArray object
## validates it by extracting an empty array from its seeds, which opens
## every HDF5 file bound by gsvaReduce()
#' @importFrom cli cli_abort
.check_unwrapping_assay <- function(container, container_class, assay) {

    anames <- names(assays(container, withDimnames=FALSE))
    if (length(assays(container, withDimnames=FALSE)) == 0L) {
        msg <- "The input {container_class} object has no assay data."
        cli_abort(c("x"=msg))
    }

    if (missing(assay) || is.na(assay)) {
        assay <- anames[1]
    } else {
        if (!is.character(assay))
            cli_abort(c("x"="The 'assay' argument must be a character string."))
        assay <- assay[1]
        if (!assay %in% anames)
            cli_abort(c("x"=paste("assay {assay} not found in the input",
                                  "{container_class} object.")))
    }

    assay
}

## converts a dgCMatrix into a list of its columns, based on
## https://rpubs.com/will_townes/sparse-apply
## it is only slightly more efficient than .sparseToList() below BUT simpler
## and does NOT offer converting to a list of rows which is far less efficient
## on a dgCMatrix object.  if you need lists of rows, simply transpose before
## calling this function, t() is reasonably fast as is calling vapply() on its
## result
#' @importFrom Matrix nnzero
.sparse2columnList <- function(m) {
    return(unname(split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))))
}

## actually, it's not just an apply() but also in-place modification
## ellipsis added for cases such as when FUN=rank where we may need
## to set the parameter 'ties.method' of the 'rank()' function
#' @importFrom BiocParallel SerialParam bplapply
.sparseColumnApplyAndReplace <- function(m, FUN, ...) {
    x <- m@x
    x <- lapply(.sparse2columnList(m), FUN=FUN, ...)
    m@x <- unlist(x, use.names=FALSE)
    if (is.integer(m@x)) ## rank(ties.method="first") returns integers
        mode(m@x) <- "numeric" ## dgCMatrix holds only doubles and logicals
    return(m)
}

#' @importFrom cli cli_abort cli_alert_warning
.check_for_na_values <- function(exprData, assay, checkNA, use) {
    autonaclasseswocheck <- c("matrix", "ExpressionSet",
                              "SummarizedExperiment",
                              "RangedSummarizedExperiment")
    mask <- class(exprData) %in% autonaclasseswocheck
    checkNAyesno <- switch(checkNA, yes="yes", no="no",
                           ifelse(any(mask), "yes", "no"))
    didCheckNA <- any_na <- FALSE
    if (checkNAyesno == "yes") {
        any_na <- anyNA(unwrapData(exprData, assay))
        didCheckNA <- TRUE
        if (any_na) {
            if (use == "all.obs")
                cli_abort(c("x"="Input expression data has NA values."))
            else if (use == "everything")
                cli_alert_warning(paste("Input expression data has NA values,",
                                        "which will be propagated through",
                                        "calculations"))
            else ## na.rm
                cli_alert_warning(paste("Input expression data has NA values,",
                                        "which will be discarded from",
                                        "calculations"))
        }
    }

    list(any_na=any_na, didCheckNA=didCheckNA)
}

.get_NAuse <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam") ||
            inherits(object, "avgParam"))
  return(object@use)
}

.get_checkNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam") ||
            inherits(object, "avgParam"))
  return(object@checkNA)
}

.get_didCheckNA <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam") ||
            inherits(object, "avgParam"))
  return(object@didCheckNA)
}

.get_ondisk <- function(object) {
  stopifnot(inherits(object, "ssgseaParam") ||
            inherits(object, "gsvaParam") ||
            inherits(object, "zscoreParam") ||
            inherits(object, "plageParam") ||
            inherits(object, "avgParam"))
  return(object@ondisk)
}

## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.isPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
  (name %in% loadedNamespaces())
}

.objPkgClass <- function(obj) {
    oc <- class(obj)
    pkg <- attr(oc, "package", exact=TRUE)
    opc <- if(is.null(pkg)) {
               oc[1]
           } else {
               paste(pkg[1], oc[1], sep = "::")
           }
    return(opc)
}

#' @importFrom Biobase selectSome
.showSome <- function(x) {
    paste0(paste(selectSome(x, 4), collapse=", "),
           " (", length(x), " total)")
}

#' @importFrom utils capture.output
.catObj <- function(x, prefix = "  ") {
    if(is.null(x)) {
        cat(paste0(prefix, "none."))
    } else {
        cat(paste0(prefix, capture.output(gsvaShow(x))), sep="\n")
    }
}

.isCharNonEmpty <- function(x) {
    return((!is.null(x)) &&
           (length(x) > 0) &&
           (is.character(x)) &&
           (!all(is.na(x))) &&
           any(nzchar(x)))
}

.omitEmptyChar <- function(x) {
    if(.isCharNonEmpty(x)) {
        return(x[(nzchar(x)) & (!is.na(x))])
    } else {
        return(character(0))
    }
}

.isCharLength1 <- function(x) {
    return((.isCharNonEmpty(x)) && (length(x) == 1))
}

.isNumNonEmpty <- function(x) {
    return((!is.null(x)) &&
           (length(x) > 0) &&
           (is.numeric(x)) &&
           (!all(is.na(x))) &&
           (!all(is.nan(x))))
}

.isNumLength1 <- function(x) {
    return((.isNumNonEmpty(x)) && (length(x) == 1))
}

#' @importFrom cli cli_abort
.check_first_last_values <- function(X, dimfun, dimname, first, last, rmdt) {
    dimfun <- match.fun(dimfun)

    if (all(is.na(first)) && all(is.na(last)))
        return(list(first=NA_real_, last=NA_real_))

    if (!is.null(rmdt))
        cli_abort(c("x"=paste("Input expression already has chunk boundary",
                              "metadata and 'first' and 'last' cannot be set.")))

    if (all(is.na(first)))
        first <- 1
    if (all(is.na(last)))
        last <- dimfun(X)

    if (.isNumLength1(first) && .isNumLength1(last)) {
        if (first < 1 || last < 1 || first != as.integer(first) ||
            last != as.integer(last))
            cli_abort(c("x"="'first' and 'last' must be positive integers."))
        if (first > last)
            cli_abort(c("x"="'first' must be smaller or equal than 'last'."))
    } else
        cli_abort(c("x"="'first' and 'last' must be a single number."))

    if (first > dimfun(X) || last > dimfun(X))
        cli_abort(c("x"=paste("'first' and 'last' must be smaller or equal than",
                              "the number of {dimname} in the input data")))

    return(list(first=first, last=last))
}

## annotation package checks
.isAnnoPkgValid <- function(ap) {
    return(.isCharLength1(ap))
}

#' @importFrom utils installed.packages
.isAnnoPkgInstalled <- function(ap) {
    ap <- c(ap, paste0(ap, ".db"))
    return(any(ap %in% rownames(installed.packages())))
}
