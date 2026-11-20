## This file contains the internal functions of the memory model of GSVA: the
## maximum main memory available to the calculations, the memory of the R
## processes running them, of the data that they hold and of the input data of
## each step, the decision of processing the input data from disk, and the
## estimated memory of each step, with the warnings when it exceeds the maximum
## main memory. The memory factors of each step and the size of the blocks
## that fit in that memory are given in block_processing.R. These functions are
## not exported and are only used internally by GSVA.

#' @importFrom cli cli_abort
.memtext2bytes <- function(x) {
  if (is.numeric(x))
      return(x)

  x <- gsub(",", ".", x)
  pat <- "(\\d*(.\\d+)*)(.*)"
  num  <- as.numeric(sub(pat, "\\1", x))
  unit <- sub(pat, "\\3", x)
  unit[unit==""] <- "1"

  fac <- c("1"=1, "K"=1024, "M"=1024^2, "G"=1024^3, "T"=1024^4)
  if (!toupper(unit) %in% names(fac)) {
      msg <- "Unknown memory unit '{unit}', please use either K, M, G or T."
      cli_abort(c("x"=msg))
  }

  num * unname(fac[toupper(unit)])
}

## memory limit in bytes of the job or container where this R process runs,
## or Inf when there is none. on Linux, workload managers such as SLURM, and
## container engines, enforce this limit through the control group (cgroup)
## of the process, which can be set in any of the cgroups along its path, e.g.,
## at the level of a SLURM job and not of its steps, so that the smallest
## limit along that path is taken, using either cgroup version 1 or 2. the
## limit given by SLURM in its environment variables is also taken into
## account. the paths of the files and the environment variables are
## arguments to allow testing this function
.job_memory_limit <- function(procfile="/proc/self/cgroup",
                              cgroupdir="/sys/fs/cgroup",
                              env=Sys.getenv(c("SLURM_MEM_PER_NODE",
                                               "SLURM_MEM_PER_CPU",
                                               "SLURM_CPUS_ON_NODE"),
                                             unset=NA)) {
    readnum <- function(fname) {
        if (!file.exists(fname))
            return(NA_real_)
        suppressWarnings(tryCatch(as.numeric(readLines(fname, n=1L,
                                                       warn=FALSE)),
                                  error=function(e) NA_real_))
    }
    limit <- Inf
    cglines <- character(0)
    if (file.exists(procfile)) ## only on Linux
        cglines <- tryCatch(readLines(procfile, warn=FALSE),
                            error=function(e) character(0))
    for (cgline in cglines) {
        ## each line has the format 'hierarchy-ID:controllers:path', where
        ## cgroup version 2 has hierarchy-ID 0 and no controllers
        fields <- strsplit(cgline, ":", fixed=TRUE)[[1]]
        if (length(fields) < 3L)
            next
        if (fields[1] == "0" && fields[2] == "") {
            dir <- cgroupdir
            limitfile <- "memory.max" ## 'max' when there is no limit
        } else if ("memory" %in% strsplit(fields[2], ",", fixed=TRUE)[[1]]) {
            dir <- file.path(cgroupdir, "memory")
            limitfile <- "memory.limit_in_bytes"
        } else
            next
        path <- strsplit(paste(fields[-seq_len(2)], collapse=":"), "/",
                         fixed=TRUE)[[1]]
        path <- path[path != ""]
        for (k in rev(seq_len(length(path) + 1L) - 1L)) {
            v <- readnum(do.call(file.path,
                                 as.list(c(dir, path[seq_len(k)],
                                           limitfile))))
            if (!is.na(v) && v > 0)
                limit <- min(limit, v)
        }
    }

    ## SLURM gives memory limits in megabytes, where 0 means no limit
    env <- suppressWarnings(as.numeric(env))
    names(env) <- c("SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU",
                    "SLURM_CPUS_ON_NODE")
    if (!is.na(env["SLURM_MEM_PER_NODE"]) && env["SLURM_MEM_PER_NODE"] > 0)
        limit <- min(limit, env["SLURM_MEM_PER_NODE"] * 1024^2)
    else if (!is.na(env["SLURM_MEM_PER_CPU"]) && env["SLURM_MEM_PER_CPU"] > 0 &&
             !is.na(env["SLURM_CPUS_ON_NODE"]))
        limit <- min(limit, env["SLURM_MEM_PER_CPU"] *
                            env["SLURM_CPUS_ON_NODE"] * 1024^2)

    unname(limit)
}

#' @importFrom cli cli_abort cli_alert_info
#' @importFrom memuse Sys.meminfo mu
.check_maxmem <- function(param, assay=get_assay(param), maxmem, verbose) {
    if (length(maxmem) > 1 || (!is.numeric(maxmem) && !is.character(maxmem))) {
        msg <- paste("'maxmem' should be a vector of length 1 of either a",
                     "number in bytes or a character string formed by a",
                     "number followed by the suffix K, M, G or T.")
        cli_abort(c("x"=msg))
    }

    if (is.character(maxmem) && maxmem == "auto") {
        ## auto takes 90% of RAM or, when smaller, of the memory limit of the
        ## job or container where this R process runs
        totalram <- Sys.meminfo()$totalram
        joblimit <- .job_memory_limit()
        avail <- totalram
        what <- "main memory"
        if (joblimit < as.numeric(totalram)) {
            avail <- mu(joblimit)
            what <- "memory of the job"
        }
        maxmem <- as.numeric(avail * 0.9)
        X <- unwrapData(get_exprData(param), assay)
        if (verbose && is(X, "DelayedArray") &&
            gsva_global$show_start_and_end_messages)
            cli_alert_info(sprintf("Maximum available %s (90%%): %s", what,
                                   as.character(avail * 0.9)))
    }

    maxmem <- .memtext2bytes(maxmem)
    maxmem
}


## fraction of nonzero values of the assay 'assay' of the expression data in
## 'param', from the number of nonzero values in 'param', see
## .estimate_nzcount(), or recomputed with 'recompute_nzcount=TRUE'
#' @importFrom SparseArray nzcount
.input_density <- function(param, assay, recompute_nzcount=FALSE) {
    X <- unwrapData(get_exprData(param), assay)
    tot <- as.numeric(nrow(X)) * as.numeric(ncol(X))
    if (tot == 0)
        return(1)
    if (recompute_nzcount)
        spa <- .estimate_nzcount(get_exprData(param), assay, FALSE) / tot
    else {
        spa <- nzcount(param) / tot
        ## the number of nonzero values in the parameter object
        ## may correspond to a larger data set, e.g., when columns
        ## have been removed, and then it has to be recomputed
        if (spa > 1)
            spa <- .estimate_nzcount(get_exprData(param), assay,
                                     FALSE) / tot
    }

    min(spa, 1)
}

## estimated size in bytes of the input data of a step in main memory, i.e.,
## of the assay 'assay' of the expression data in 'param', restricted to the
## rows (whdim=1) or columns (whdim=2) from 'first' to 'last', when given,
## where sparse data, in an SVT_SparseMatrix object, store each nonzero value
## and its row index, unless 'dense=TRUE', when they are loaded in main memory
## as a dense matrix, as the methods other than GSVA and average do
#' @importFrom cli cli_abort
#' @importFrom BiocGenerics type
#' @importFrom S4Arrays is_sparse
.input_mem_size <- function(param, assay, first, last, whdim,
                            recompute_nzcount=FALSE, dense=FALSE) {
    X <- unwrapData(get_exprData(param), assay)
    nr <- as.numeric(nrow(X))
    nc <- as.numeric(ncol(X))
    spa <- 1
    if (is_sparse(X) && !dense)
        spa <- .input_density(param, assay, recompute_nzcount)
    if (!is.na(first) || !is.na(last)) {
        if (whdim == 1)
            nr <- as.numeric(last - first + 1)
        else if (whdim == 2)
            nc <- as.numeric(last - first + 1)
        else
            cli_abort(c("x"="Invalid internal value for 'whdim' argument."))
    }
    eltbytes <- if (type(X) == "integer") 4 else 8
    if (is_sparse(X) && !dense)
        return(spa * nr * nc * (eltbytes + 4))

    nr * nc * eltbytes
}

## memory in bytes taken by the data of a step processing in blocks a matrix
## with dimensions 'dims', along its rows (whdim=1) or columns (whdim=2), of
## size 'insize' in main memory and 'eltbytes' bytes per value, sparse or not,
## and stored in main memory or not ('inmemory'): the memory that remains
## allocated, see .fixed_mem(), and the working memory of the rows or columns
## of a block of the default size of the DelayedArray package, shared by all
## the workers, see .units_per_block(), which take their size in main memory
## or, when read from disk, the one of their dense form. 'mf' gives the memory
## factors of the step, see .step_mem_factors(). with 'densify=TRUE', the
## calculations give a dense output from a sparse input, e.g., the classical
## algorithm of GSVA on sparse data, and the output and the working memory of
## the blocks take the size of the dense form, while the input keeps its size.
## with 'sendblocks=TRUE', the main R process takes a copy of the block of each
## worker to send it to that worker, see .sends_blocks(). with 'sinkout=TRUE',
## the output is written to disk block by block, see .fixed_mem()
#' @importFrom DelayedArray getAutoBlockSize
.step_data_mem <- function(dims, whdim, insize, eltbytes, sparse, inmemory,
                           mf, heldmem=.held_mem(), densify=FALSE,
                           sendblocks=FALSE, sinkout=FALSE) {
    nunits <- max(1, dims[whdim])
    denseunitbytes <- as.numeric(dims[-whdim]) * eltbytes
    unitbytes <- if (inmemory && !densify) insize / nunits else denseunitbytes
    minunits <- min(nunits, max(1, floor(getAutoBlockSize() / denseunitbytes)))
    .fixed_mem(nunits, insize / nunits, whdim, sparse && !densify, inmemory,
               mf$outfactor, mf$outextra, mf$assembly,
               outunitbytes=if (densify) denseunitbytes else insize / nunits,
               sinkout=sinkout) +
        mf$workfactor * minunits * unitbytes +
        (if (sendblocks && inmemory) minunits * insize / nunits else 0) +
        heldmem
}

## verifies that the 'ondisk' parameter is either 'auto', 'yes' or 'no' and, if
## 'auto', checks whether the calculations fit in the maximum available main
## memory with the input data in it, and sets 'ondisk' to 'yes' or 'no'
## accordingly. the calculations fit when the memory that they require with
## the input data in main memory and the parallel workers of 'BPPARAM', see
## .step_mem_need(), fits in 'maxmem': the one taken by the main R process and
## by each worker by themselves, see .worker_mem(), and the one taken by their
## data, see .step_data_mem(), together with the rest of the input data in main
## memory, see .step_held_input(), from the fraction .mem_fraction_R() of
## 'maxmem' available to the allocations of R, where 'mf' gives the memory
## factors of the step, see .step_mem_factors(), and 'dense=TRUE' indicates
## that the input data is loaded in main memory as a dense matrix, even when
## it is sparse. If
## the input data is a DelayedArray, also reports whether the calculations fit
## in the maximum available main memory.
## If 'ondisk' is set to 'yes', then this function returns TRUE, otherwise it
## returns FALSE, with the estimated size of the input data in main memory and
## 'maxmem' as attributes.

#' @importFrom cli cli_abort cli_alert_info cli_alert_warning
#' @importFrom S4Arrays is_sparse
#' @importFrom BiocGenerics type
.check_ondisk <- function(param, assay=get_assay(param), first, last, whdim,
                          recompute_nzcount=FALSE, maxmem, verbose,
                          mf=list(workfactor=2, outfactor=1, outextra=0),
                          dense=FALSE, densify=FALSE, BPPARAM=NULL) {
    ondisk <- .get_ondisk(param)
    if (ondisk != "auto" && ondisk != "yes" && ondisk != "no")
        cli_abort(c("x"="'ondisk' should be either 'auto', 'yes' or 'no'"))

    insize <- .input_mem_size(param, assay, first, last, whdim,
                              recompute_nzcount, dense)
    if (ondisk == "auto") {
        X <- unwrapData(get_exprData(param), assay)
        dims <- dim(X)
        if (!is.na(first) || !is.na(last))
            dims[whdim] <- last - first + 1
        ## the rest of the input data in main memory remains allocated while
        ## the step runs, see .step_held_input()
        heldmem <- .held_mem() +
                   .step_held_input(X, dims, whdim, FALSE,
                                    .input_held_mem(param))
        eltbytes <- if (type(X) == "integer") 4 else 8
        nworkers <- .n_par_workers(BPPARAM, dims)
        ## each R process takes the memory of this one, besides the matrices
        ## of the assays of the input data that it holds, see .worker_mem()
        permem <- .worker_mem(.held_shared_mem() + .input_shared_mem(param),
                              full=.large_data(dims, eltbytes))
        need <- .step_mem_need(dims, insize, whdim, eltbytes,
                               is_sparse(X) && !dense, TRUE, mf, nworkers,
                               heldmem, densify,
                               .sends_blocks(TRUE, BPPARAM, nworkers),
                               permem=permem, forked=.forked_workers(BPPARAM))
        ondisk <- "no"
        if (need > maxmem) {
            ondisk <- "yes"
            if (is(X, "DelayedArray") && verbose)
                cli_alert_info(paste("Calculations with on-disk input data",
                                     "loaded in main memory do not fit in",
                                     "the maximum available main memory"))
            ## steps giving enrichment scores return them on disk, which
            ## changes the class of their output
            if (mf$outextra > 0)
                cli_alert_warning(paste("The calculations do not fit in the",
                                        "maximum available main memory, and",
                                        "the resulting (dense) matrix of",
                                        "enrichment scores will be returned",
                                        "using an on-disk data structure"))
        } else if (is(X, "DelayedArray") && verbose)
            cli_alert_info(paste("Calculations with on-disk input data",
                                 "loaded in main memory fit in the maximum",
                                 "available main memory"))
    }

    ## the estimated size of the input data and the maximum main memory allow
    ## .check_sparse_load_input_expr() to load the input data into main
    ## memory by blocks that fit in the memory left by it
    structure(ondisk == "yes", insize=insize, maxmem=maxmem)
}

## memory in bytes that each R process running GSVA takes by itself, the main
## one and each of its parallel workers, besides the data of the calculations:
## the memory allocated by R in this process, other than the 'datamem' bytes
## of the data of the calculations that it holds, plus the memory that each R
## process takes outside R, see .proc_extra_mem(). forked workers, such as
## those of a MulticoreParam object, share instead with the main process the
## memory allocated by R, see .forked_worker_mem(). loading a
## SingleCellExperiment object with one million cells adds about 0.2 GB to
## the memory of this process. with 'session=FALSE', it is the memory of a
## new R process loading GSVA, such as the one of a job of gsvaMap() or of a
## socket worker, instead of the one of this process, see 'newprocmem' in
## .memory_model_defaults, together with the memory outside R. the memory
## allocated by R is measured
## with a full collection of garbage with 'full=TRUE', see .r_mem_used(). it can
## be set with the option 'GSVA.workermem'
.worker_mem <- function(datamem=0, session=TRUE, full=FALSE) {
    wm <- getOption("GSVA.workermem")
    if (!is.null(wm))
        return(as.numeric(wm))
    if (!session)
        return(.memory_model("newprocmem") + .proc_extra_mem())

    max(0, .r_mem_used(full) - datamem) + .proc_extra_mem()
}

## memory in bytes allocated by R in this process, as tracked by its garbage
## collector, where each cons cell takes 56 bytes and each vector cell 8 bytes.
## a full collection ('full=TRUE') takes about 0.2 s in a session with GSVA
## loaded, which adds up over the many calls of small calculations, while a
## collection of the youngest generation takes less than 1 ms and gives
## about the same memory, but also the garbage left in the older generations,
## e.g., about 0.1 GB after loading a SingleCellExperiment object with one
## million cells, which each R process would take, see .worker_mem(). within
## a step, which measures it several times, see .open_mem_scope(), a full
## collection is done only when the one of the youngest generation gives more
## than 32 MB of difference with the last full one of the step, giving then
## the memory of the former, which is at least the one in use
.r_mem_used <- function(full=FALSE) {
    g <- gc(full=FALSE)
    used <- sum(g[, 1L] * c(56, 8))
    if (full) {
        lastfull <- gsva_global$rmemfull
        if (is.null(lastfull) || is.na(lastfull) ||
            abs(used - lastfull) > 32 * 1024^2) {
            g <- gc(full=TRUE)
            used <- sum(g[, 1L] * c(56, 8))
            if (!is.null(lastfull))
                gsva_global$rmemfull <- used
        }
    }
    used
}

## release, with a full collection of garbage, the memory that the
## calculations on a block of 'ncols' columns, more than 10000, have allocated
## and no longer need, only in an R process that is not a worker forked from
## the one that loaded GSVA, such as those of a MulticoreParam object: a full
## collection in a forked worker writes to the memory pages of every object
## that it shares with the main R process, which the worker then copies, e.g.,
## about the whole 4.5 GB of an R session holding the ranks of 100,000 cells,
## in each of 4 workers scoring them, while the memory it would release goes
## back to the system anyway when the worker ends
.release_block_mem <- function(ncols) {
    if (ncols > 10000 && identical(Sys.getpid(), gsva_global$pid))
        invisible(gc())
}

## open, in a step, the scope where the memory allocated by R is measured
## with a full collection of garbage only once, see .r_mem_used(), unless
## another step running this one has opened it; it returns whether it opened
## it, to be given to .close_mem_scope() on exit. outside such a scope, such
## as between the steps run by gsva(), whose objects may then be released
## without being collected, each measurement does a full collection
.open_mem_scope <- function() {
    if (!is.null(gsva_global$rmemfull))
        return(FALSE)
    gsva_global$rmemfull <- NA_real_
    TRUE
}

.close_mem_scope <- function(opened) {
    if (opened)
        gsva_global$rmemfull <- NULL
}

## whether a matrix with dimensions 'dims' and 'eltbytes' bytes per value is
## large enough, more than 1 GB in its dense form, for the memory of the R
## processes processing it to be measured with a full collection of garbage,
## see .worker_mem(), whose time is then small compared to the one of its
## calculations
.large_data <- function(dims, eltbytes=8) {
    prod(as.numeric(dims)) * eltbytes > 1024^3
}

## memory in bytes that each R process takes outside the allocations of R,
## such as the memory that R releases but the process keeps, on this
## operating system, see 'procextramem' in .memory_model_defaults
.proc_extra_mem <- function() {
    pe <- .memory_model("procextramem")
    os <- Sys.info()[["sysname"]]

    unname(if (os %in% names(pe)) pe[[os]] else pe[["default"]])
}

## memory in bytes that each worker forked from this R process, such as those
## of a MulticoreParam object, takes by itself besides the blocks of its
## calculations, whose input is in main memory ('inmemory=TRUE') or not:
## with the input in main memory, the memory that it takes outside R, see
## .proc_extra_mem(), because it shares with this process the memory
## allocated by R, except for the memory pages that it copies, see
## .fork_copied_mem(); about 0.04 GB in each worker scoring the ranks of
## 100,000 cells held in main memory, and nothing measurable in each worker
## normalizing their rows, as measured on Linux. with the input read from
## disk, the memory of this process 'permem', see .worker_mem(), because each
## worker then processes several rounds of blocks and keeps memory of its own
## from one round to the next: with less, the steps of one million cells read
## from disk by 10 workers within 'maxmem' exceeded it on macOS. it can be
## set with the option 'GSVA.workermem', as the one of each R process
.forked_worker_mem <- function(inmemory=TRUE, permem=.worker_mem()) {
    wm <- getOption("GSVA.workermem")
    if (!is.null(wm))
        return(as.numeric(wm))
    if (inmemory) .proc_extra_mem() else permem
}

## memory in bytes that each forked worker copies from the 'sessionmem' bytes
## allocated by R in the main R process, including the data of the
## calculations that it holds, in a step whose factor 'forkcopy' gives the
## fraction copied, see .step_mem_factors(): a full collection of garbage in a
## worker writes to the memory page with the header of each object, which the
## worker then copies, so that the data of large vectors remain shared, such
## as the ones of a dense matrix, while objects made of many small vectors,
## such as an SVT_SparseMatrix object, are copied almost completely. workers
## that allocate little memory, as when scoring gene sets, rarely collect
## garbage and copy almost nothing, while workers ranking the columns of
## 100,000 cells held in main memory copied about 0.17 times the memory of the
## main R process, as measured on Linux and macOS, where GSVA does not collect
## garbage itself in forked workers, see .release_block_mem(). workers
## calculating the row ECDFs of 100,000 cells held in main memory as an
## SVT_SparseMatrix object, which take longer, ran full collections at about
## the same time and copied almost the whole memory of the main R process,
## about 3 GB each, on macOS. 'forkcopy' can give two fractions, the first
## for the input in main memory ('inmemory=TRUE') and the second for the input
## read from disk, see .step_mem_factors()
.fork_copied_mem <- function(sessionmem, forkcopy=NULL, inmemory=TRUE) {
    if (is.null(forkcopy))
        return(0)

    forkcopy[if (inmemory) 1L else length(forkcopy)] * sessionmem
}

## memory in bytes that each R process takes outside R while it reads blocks
## of input data from disk, such as the buffers of the HDF5 library to
## decompress the chunks read, see 'diskreadmem' in .memory_model_defaults
.disk_read_mem <- function() {
    .memory_model("diskreadmem")
}

## memory in bytes that other objects keep allocated in the main R process
## while a step processes its input in blocks, which the blocks cannot take,
## see .units_per_block(); gsvaParam's gsva() sets it to the input data in
## main memory while it runs the steps after the row normalization
.held_mem <- function() {
    as.numeric(gsva_global$heldmem)
}

## memory in bytes of the matrices of the assays among the objects held by
## .held_mem(), which forked workers share with the main R process, see
## .worker_mem()
.held_shared_mem <- function() {
    as.numeric(gsva_global$heldshared)
}

## memory in bytes that the matrix 'X' takes in main memory, which is none
## when it is stored on disk
#' @importFrom utils object.size
.inmem_size <- function(X) {
    if (is(X, "DelayedArray"))
        return(0)

    as.numeric(object.size(X))
}

## memory in bytes that the container 'x' of expression data takes in main
## memory. for a SummarizedExperiment object, all of it: all its assays held
## in main memory, including those not used by the calculations, and the rest
## of its data, such as the reduced dimensions and alternative experiments of
## a SingleCellExperiment object, or the spatial coordinates and images of a
## SpatialExperiment object, in its internal slots 'int_colData' and
## 'int_metadata'. for an ExpressionSet object, also the matrices of its assay
## data, which object.size() does not include when they are stored in an
## environment. otherwise, the one of its matrix of the assay 'assay'
#' @importFrom utils object.size
#' @importMethodsFrom Biobase assayData
.container_held_mem <- function(x, assay) {
    if (is(x, "SummarizedExperiment"))
        return(as.numeric(object.size(x)))
    if (is(x, "ExpressionSet")) {
        held <- as.numeric(object.size(x))
        ad <- assayData(x)
        if (is.environment(ad))
            held <- held + sum(vapply(ls(ad), function(nm)
                                          .inmem_size(get(nm, envir=ad)),
                                      numeric(1)))
        return(held)
    }

    .inmem_size(unwrapData(x, assay))
}

## memory in bytes that the input data of the parameter object 'param' takes
## in main memory, including the assays not used by the calculations
.input_held_mem <- function(param) {
    .container_held_mem(get_exprData(param), get_assay(param))
}

## memory in bytes that the matrices of the assays of the container 'x' of
## expression data take in main memory, see .container_held_mem(), whose
## values are in large vectors that forked workers share with the main R
## process, see .fork_copied_mem(), unlike the rest of its data, such as the
## many character strings of the names of the cells of single-cell data,
## which they may copy, see .worker_mem()
#' @importMethodsFrom SummarizedExperiment assays
#' @importMethodsFrom Biobase assayData
.container_shared_mem <- function(x, assay) {
    if (is(x, "SummarizedExperiment"))
        return(sum(vapply(assays(x, withDimnames=FALSE), .inmem_size,
                          numeric(1))))
    if (is(x, "ExpressionSet")) {
        ad <- assayData(x)
        nms <- if (is.environment(ad)) ls(ad) else names(ad)
        return(sum(vapply(nms, function(nm) .inmem_size(ad[[nm]]),
                          numeric(1))))
    }

    .inmem_size(unwrapData(x, assay))
}

.input_shared_mem <- function(param) {
    .container_shared_mem(get_exprData(param), get_assay(param))
}

## memory in bytes held, besides the memory of the step itself, while a step
## processes the input matrix 'X', restricted to dimensions 'dims' along
## 'whdim', from a container whose assays take 'container' bytes in main
## memory, see .container_held_mem(): its assays other than 'X' and, when 'X'
## is in main memory, the whole of 'X' when only part of it is processed, and
## the part processed, which the step copies, when its output is written to
## disk ('ondisk=TRUE'), because otherwise the memory of the step includes it,
## see .step_data_mem(). when 'X' is stored on disk, the step loads in main
## memory only the part processed and only when its output is not written to
## disk, which the memory of the step includes
.step_held_input <- function(X, dims, whdim, ondisk,
                             container=.inmem_size(X)) {
    held <- max(0, container - .inmem_size(X))
    if (is(X, "DelayedArray"))
        return(held)

    frac <- as.numeric(dims[whdim]) / max(1, dim(X)[whdim])
    if (frac < 1)
        held <- held + .inmem_size(X)
    if (ondisk)
        held <- held + frac * .inmem_size(X)

    held
}

## add the memory held by the input data while a step processes its matrix
## 'X', along its rows (whdim=1) or columns (whdim=2), restricted to those
## from 'first' to 'last', when given, see .step_held_input(), to the memory
## held by other objects, see .held_mem(), which the blocks of the
## calculations cannot take, while the calling function runs, restoring it
## when that function exits; 'container' is the memory that the assays of the
## input data take in main memory, see .container_held_mem(), and 'shared'
## the one of its matrices of assays among them, see .container_shared_mem(),
## which is held likewise, see .held_shared_mem(). the part of 'X' processed
## is not included, because the blocks of the calculations include it when
## it is in main memory, see .units_per_block()
.hold_input_mem <- function(X, whdim, first, last, container,
                            shared=.inmem_size(X), envir=parent.frame()) {
    dims <- dim(X)
    if (!is.na(first) || !is.na(last))
        dims[whdim] <- (if (is.na(last)) dims[whdim] else last) -
                       (if (is.na(first)) 1 else first) + 1
    old <- gsva_global$heldmem
    oldshared <- gsva_global$heldshared
    gsva_global$heldmem <- old + .step_held_input(X, dims, whdim, FALSE,
                                                   container)
    gsva_global$heldshared <- oldshared +
                              .step_held_input(X, dims, whdim, FALSE, shared)
    do.call(on.exit, list(substitute({
                              gsva_global$heldmem <- old
                              gsva_global$heldshared <- oldshared
                          }, list(old=old, oldshared=oldshared)),
                          add=TRUE), envir=envir)

    invisible(old)
}

## number of parallel workers of 'BPPARAM' that process a matrix with
## dimensions 'dims', which is 0 when its calculations are not parallelized,
## see .check_open_parallelism()
#' @importFrom BiocParallel bpnworkers
.n_par_workers <- function(BPPARAM, dims, minparrows=100, minparcols=100) {
    if (is.null(BPPARAM) || bpnworkers(BPPARAM) <= 1 ||
        dims[1] <= minparrows || dims[2] <= minparcols)
        return(0L)

    as.integer(bpnworkers(BPPARAM))
}

## whether the parallel workers of 'BPPARAM' are forked from this process, as
## those of a MulticoreParam object, and share its memory until they write to
## it, unlike those of other back-ends, such as SnowParam, which are new R
## processes, see .block_worker_mem()
.forked_workers <- function(BPPARAM) {
    is(BPPARAM, "MulticoreParam")
}

## estimated memory in bytes required by a step processing in blocks a matrix
## with dimensions 'dims', along its rows (whdim=1) or columns (whdim=2), of
## size 'insize' in main memory and 'eltbytes' bytes per value, sparse or not,
## and stored in main memory or not ('inmemory'), with 'nworkers' parallel
## workers, where 'mf' gives the memory factors of the step, see
## .step_mem_factors(). the main R process takes 'permem' bytes, by default
## the ones of .worker_mem() besides the 'datamem' bytes of the data of the
## calculations that this process holds, each worker forked from it the ones
## of .forked_worker_mem(), from 'permem' when the input is read from disk,
## and the data of the
## step take .step_data_mem() from the fraction .mem_fraction_R() of the memory
## available to the allocations of R, including the copies of the blocks sent
## to the workers with 'sendblocks=TRUE', see .sends_blocks(). each R
## process also takes .disk_read_mem() when the input data is read from disk,
## i.e., when it is not loaded in main memory ('inmemory=FALSE'). forked
## workers take a copy of the fraction 'forkcopy' of 'mf' of the memory of
## this process, including the data of the step in main memory, its input and
## the 'heldmem' bytes held by other objects, see .fork_copied_mem(). workers
## not forked from this process ('forked=FALSE'), such as
## those of a SnowParam object, take instead each the memory of a new R
## process, see .worker_mem(), and no copy of those data. blocks of sparse
## data read from disk take their memory by their sparse size, see
## .sparse_disk_factors(), from the ratio 'sparseratio' of the sparse size to
## the dense size of the input, by default the one of 'insize', the estimated
## size of the input in main memory, see .input_mem_size(), and not when it
## is NA. with 'sinkout=TRUE', the output is written to disk block by block,
## see .step_data_mem(). the attribute "workermem" of the result gives the
## memory of each worker, including its copy of the memory of this process
.step_mem_need <- function(dims, insize, whdim, eltbytes, sparse, inmemory,
                           mf, nworkers, heldmem=.held_mem(), densify=FALSE,
                           sendblocks=FALSE, datamem=0,
                           permem=.worker_mem(datamem), forked=TRUE,
                           sparseratio=insize / max(1, prod(as.numeric(dims)) *
                                                       eltbytes),
                           sinkout=FALSE) {
    if (sparse && !densify && !inmemory)
        mf <- .sparse_disk_factors(mf, sparseratio)
    readmem <- if (inmemory) 0 else .disk_read_mem()
    ## forked workers copy a part of the memory of this process, including
    ## the data of the calculations that it holds, see .fork_copied_mem()
    parmem <- .forked_worker_mem(inmemory, permem) + readmem
    copied <- .fork_copied_mem(max(0, permem - .proc_extra_mem()) + heldmem +
                               if (inmemory) insize else 0, mf$forkcopy,
                               inmemory)
    if (!forked) {
        parmem <- .worker_mem(session=FALSE) + readmem
        copied <- 0
    }

    need <- permem + readmem + nworkers * (parmem + copied) +
            .step_data_mem(dims, whdim, insize, eltbytes, sparse, inmemory, mf,
                           heldmem, densify, sendblocks, sinkout) /
            .mem_fraction_R()

    structure(need, workermem=parmem + copied)
}

## warn, before the calculations of a step start, when the estimated memory
## 'need' that they require, with 'nworkers' parallel workers, exceeds the
## maximum main memory 'maxmem', suggesting fewer workers or a larger 'maxmem',
## or what 'hint' says, when given, where the main R process takes 'permem'
## bytes by itself, see .worker_mem(), and each worker forked from it
## 'forkedmem' bytes, see .forked_worker_mem(), together with its copy of the
## memory of this process, see .step_mem_need(), while workers not forked from
## this process ('forked=FALSE') take the memory of a new R process
#' @importFrom cli cli_warn
#' @importFrom memuse mu
.check_mem_need <- function(need, maxmem, nworkers, hint=NULL,
                            permem=.worker_mem(), forked=TRUE,
                            forkedmem=.forked_worker_mem(TRUE, permem)) {
    need <- as.vector(need)
    if (need <= maxmem)
        return(invisible(FALSE))

    needtxt <- as.character(mu(need))
    maxmemtxt <- as.character(mu(maxmem))
    wmtxt <- as.character(mu(permem))
    msg <- c("!"=paste("The memory estimated for these calculations,",
                       "{needtxt}, exceeds the maximum available main memory",
                       "of {maxmemtxt}."))
    if (nworkers > 0 && forked) {
        pwmtxt <- as.character(mu(forkedmem))
        msg <- c(msg,
                 "i"=paste("The main R process takes about {wmtxt} by itself,",
                           "for GSVA, the packages it depends on and the",
                           "other objects of the R session, which its",
                           "{nworkers} parallel workers share, each taking",
                           "about {pwmtxt} besides."))
    } else if (nworkers > 0) {
        pwmtxt <- as.character(mu(.worker_mem(session=FALSE)))
        msg <- c(msg,
                 "i"=paste("The main R process takes about {wmtxt} by itself,",
                           "for GSVA, the packages it depends on and the",
                           "other objects of the R session, and each of its",
                           "{nworkers} parallel workers about {pwmtxt}, for",
                           "GSVA and the packages it depends on."))
    }
    if (!is.null(hint))
        msg <- c(msg, "i"=hint)
    else if (nworkers > 0)
        msg <- c(msg, "i"=paste("Consider using fewer parallel workers in",
                                "{.arg BPPARAM}, or a larger {.arg maxmem}."))
    else
        msg <- c(msg, "i"="Consider using a larger {.arg maxmem}.")
    cli_warn(msg)

    invisible(TRUE)
}

## check the memory required by a step processing in blocks the input data 'X',
## along its rows (whdim=1) or columns (whdim=2), restricted to those from
## 'first' to 'last', when given, with the on-disk decision 'ondisk' given by
## .check_ondisk(), the memory factors 'mf' of the step, the parallel back-end
## 'BPPARAM' and the maximum main memory 'maxmem', see .check_mem_need(). It is
## skipped while gsva() runs the steps, because gsva() checks the memory
## required by all of them before starting, and with the option
## 'GSVA.check_memory=FALSE'. 'dense=TRUE' indicates that the input data is
## processed as a dense matrix, even when it is sparse, and 'densify=TRUE' that
## a sparse input gives a dense output, see .check_ondisk(). 'extramem' bytes
## are taken by the calculations besides those of the blocks, such as the
## products of the genes of each gene set in PLAGE, see .plage_gram_mem().
## the memory of the container of 'X', 'container', see .container_held_mem(),
## and of 'X' itself remains allocated as given by .step_held_input(),
## independently of whether the output is written to disk, where 'shared' is
## the memory of the matrices of its assays, see .container_shared_mem().
## 'ondisksink=TRUE' indicates that the step writes its output to disk block by
## block when it processes its input from disk, see .step_mem_need()
## input data in main memory, which is not loaded from disk, remains
## allocated when the calculations are done on disk
#' @importFrom S4Arrays is_sparse
#' @importFrom BiocGenerics type
.check_step_mem <- function(X, whdim, first, last, ondisk, mf, BPPARAM,
                            maxmem, dense=FALSE, densify=FALSE,
                            extramem=0, container=.inmem_size(X),
                            shared=.inmem_size(X), ondisksink=FALSE) {
    if (!gsva_global$check_memory || !getOption("GSVA.check_memory", TRUE))
        return(invisible(FALSE))

    dims <- dim(X)
    if (!is.na(first) || !is.na(last))
        dims[whdim] <- last - first + 1
    eltbytes <- if (type(X) == "integer") 4 else 8
    insize <- attr(ondisk, "insize")
    ## the size of the input, and of its blocks read from disk, is the one of
    ## its dense form when it is not estimated, e.g., for data in Parquet
    ## files, processed from disk, see .step_mem_need()
    sparseratio <- if (is.null(insize)) NA else
                       insize / max(1, prod(as.numeric(dims)) * eltbytes)
    if (is.null(insize))
        insize <- prod(as.numeric(dims)) * eltbytes
    nworkers <- .n_par_workers(BPPARAM, dims)
    heldmem <- .held_mem() + extramem +
               .step_held_input(X, dims, whdim, ondisk, container)
    ## each R process takes the memory of this one, besides the matrices of
    ## the assays of the input data that it holds, see .worker_mem()
    permem <- .worker_mem(.held_shared_mem() + shared,
                          full=.large_data(dims, eltbytes))
    need <- .step_mem_need(dims, insize, whdim, eltbytes,
                           is_sparse(X) && !dense, !ondisk, mf, nworkers,
                           heldmem, densify,
                           .sends_blocks(!ondisk, BPPARAM, nworkers),
                           permem=permem, forked=.forked_workers(BPPARAM),
                           sparseratio=sparseratio,
                           sinkout=ondisksink && isTRUE(as.vector(ondisk)))

    .check_mem_need(need, maxmem, nworkers, permem=permem,
                    forked=.forked_workers(BPPARAM),
                    forkedmem=attr(need, "workermem"))
}

## check the memory required by each of the three steps that gsva() runs on
## the parameter object 'param', see .check_mem_need(), before running them,
## estimating the size of the input data of each step from the one of the
## previous step, and the on-disk decision of each step as .check_ondisk()
## does
#' @importFrom BiocGenerics type
#' @importFrom S4Arrays is_sparse
.check_gsva_mem <- function(param, BPPARAM, maxmem) {
    if (!getOption("GSVA.check_memory", TRUE))
        return(invisible(FALSE))

    X <- unwrapData(get_exprData(param), get_assay(param))
    dims <- dim(X)
    sparse <- is_sparse(X)
    nworkers <- .n_par_workers(BPPARAM, dims)
    ondisk <- .get_ondisk(param)
    ngs <- length(get_geneSets(param))
    insize <- .input_mem_size(param, get_assay(param), NA, NA, 1L)
    clr <- .get_rowNorm(param) == "clr"
    dgc <- is(X, "dgCMatrix")
    ## the classical algorithm on sparse input normalizes rows into a dense
    ## matrix, which is the input of the next steps
    densify <- sparse && !.get_sparse(param) && .get_rowNorm(param) != "none"
    steps <- list(list(step="rownorm", whdim=1L, int=(type(X) == "integer")),
                  list(step="colranks", whdim=2L, int=FALSE),
                  list(step="scores", whdim=2L, int=TRUE))
    ## the assays of the input data in main memory remain allocated during all
    ## the steps, see gsva(), the other assays during the row normalization,
    ## whose input is the assay used; the input of each step is in main
    ## memory when it is the assay used and this is in main memory, or when
    ## it is the output of the previous step and this was not written to
    ## disk, and it remains allocated when the step writes its own output to
    ## disk, see .step_held_input()
    container <- .input_held_mem(param)
    rownormnone <- .get_rowNorm(param) == "none"
    held <- max(0, container - .inmem_size(X))
    inres <- !is(X, "DelayedArray")
    ## each R process takes the memory of this one, besides the matrices of
    ## the assays of the input data that it holds, see .worker_mem()
    permem <- .worker_mem(.input_shared_mem(param), full=.large_data(dims))
    need <- 0
    forkedmem <- .forked_worker_mem(ondisk == "no", permem)
    for (st in steps) {
        mf <- .step_mem_factors(st$step, ngs=ngs, sparse=sparse && !densify,
                                int=st$int, clr=clr, dgc=dgc,
                                insparse=sparse)
        if (densify && clr) ## dense CLR, see gsvaRowNorm()
            mf$assembly <- 2
        eltbytes <- if (st$int) 4 else 8
        inmemory <- ondisk == "no"
        if (ondisk == "auto")
            inmemory <- .step_mem_need(dims, insize, st$whdim, eltbytes,
                                       sparse, TRUE, mf, nworkers, held,
                                       densify,
                                       .sends_blocks(TRUE, BPPARAM, nworkers),
                                       permem=permem,
                                       forked=.forked_workers(BPPARAM)) <=
                        maxmem
        stepheld <- held + if (inres && !inmemory) insize else 0
        ## the steps processing their input from disk write their output to
        ## disk block by block, see gsvaRowNorm(), gsvaColRanks() and
        ## gsvaColScores()
        stepneed <- .step_mem_need(dims, insize, st$whdim, eltbytes,
                                   sparse, inmemory, mf, nworkers, stepheld,
                                   densify,
                                   .sends_blocks(inmemory, BPPARAM, nworkers),
                                   permem=permem,
                                   forked=.forked_workers(BPPARAM),
                                   sinkout=!inmemory)
        if (stepneed > need) {
            need <- as.vector(stepneed)
            forkedmem <- attr(stepneed, "workermem")
        }
        ## the input of the next step and its size: without row
        ## normalization, it is the assay used, held by the input data
        if (st$step == "rownorm" && rownormnone) {
            held <- max(0, container - .inmem_size(X))
        } else {
            inres <- inmemory
            held <- container
            if (densify)
                insize <- prod(as.numeric(dims)) * 8
            else
                insize <- insize * mf$outfactor
        }
        if (densify)
            sparse <- densify <- FALSE
    }

    .check_mem_need(need, maxmem, nworkers, permem=permem,
                    forked=.forked_workers(BPPARAM), forkedmem=forkedmem)
}
