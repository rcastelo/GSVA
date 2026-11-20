#' @title Memory of the jobs running GSVA
#'
#' @description Estimates the memory that each job running GSVA requires on a
#' workload manager, such as SLURM, to give it, e.g., in the argument `mem`
#' of [`gsvaBatchtoolsSlurmParam`], from the memory model of GSVA, which sizes
#' the blocks of its calculations within the maximum main memory `maxmem`.
#'
#' @param param A parameter object of one of the classes
#' [`gsvaParam`][gsvaParam-class], [`plageParam`][plageParam-class],
#' [`zscoreParam`][zscoreParam-class], [`ssgseaParam`][ssgseaParam-class] or
#' [`avgParam`][avgParam-class], built from the input data of the jobs.
#'
#' @param ncpus_per_task Number of CPU cores of each job, which GSVA uses as
#' parallel workers, as the argument of the same name of
#' [`gsvaBatchtoolsSlurmParam`].
#'
#' @param map When `TRUE`, the default with a [`gsvaParam`][gsvaParam-class]
#' object, the memory of each of the jobs that [`gsvaMap`] sends to run the
#' three steps of the GSVA method on chunks of the input data, see below.
#' When `FALSE`, the default with other parameter objects, the memory of a
#' single job running [`gsva`] on the whole input data, with
#' `ncpus_per_task` parallel workers.
#'
#' @details The memory of the jobs of [`gsvaMap`] determines the size of the
#' chunks of the input data that each one processes: with more memory, each
#' job processes a larger chunk, and fewer jobs are needed. For each of the
#' steps of the GSVA method, row normalization ([`gsvaRowNorm`]), column
#' ranks ([`gsvaColRanks`]) and enrichment scores ([`gsvaColScores`]),
#' `gsvaJobMemory()` gives the smallest memory with which a job fits its own
#' chunk, together with its R process and its parallel workers, and a
#' recommended memory, which, when the input data is stored on disk in chunks
#' of several columns, such as in a dense HDF5 dataset, is the smallest one
#' with which each job takes whole chunks of the data on disk, since chunks
#' read by several jobs are read and decompressed by each of them, which
#' makes the calculations slower. The memory returned is the largest
#' recommended memory of the three steps, so that the same jobs can run all
#' of them. The number of jobs running at the same time, e.g., the argument
#' `nodes` of [`gsvaBatchtoolsSlurmParam`], is not needed, because more jobs
#' only make the chunks smaller.
#'
#' For a single job running [`gsva`], the smallest memory is the one with
#' which the input data, when stored on disk, is processed from disk, and the
#' recommended one follows the argument `ondisk` of the parameter object:
#' with `ondisk="auto"`, the default, or `ondisk="no"`, the memory with which
#' the input data is loaded and processed in main memory, which is faster,
#' and with `ondisk="yes"`, the one with which it is processed from disk,
#' which, with the GSVA method on input data stored in chunks of several
#' columns, is the smallest one with which its blocks take whole chunks, as
#' with the jobs of [`gsvaMap`].
#'
#' The memory model of GSVA was calibrated on single-cell data on Linux and
#' macOS, and its estimates include a margin; they are meant for jobs running
#' on the same operating system as this R session, as the login node of a
#' computer cluster usually does with its compute nodes. For a single job,
#' the memory that this R session takes by itself is the one estimated for
#' the job, so that it should hold the input data and little else. The
#' estimates require no calculations on the input data, besides estimating
#' the fraction of nonzero values of sparse data stored on disk, and, for a
#' single job with methods other than GSVA, filtering and mapping its genes,
#' as [`gsva`] does.
#'
#' @return A character string with the recommended memory, in whole GB or
#' MB, e.g., `"14G"`, which can be given to the argument `mem` of
#' [`gsvaBatchtoolsSlurmParam`]. It has the class `gsvaJobMemory`, and the
#' function `details()` shows, for each step, the smallest and the
#' recommended memory, as well as, for the jobs of [`gsvaMap`], the largest
#' chunk of each job and the number of chunks with the recommended memory.
#'
#' @seealso [`gsvaBatchtoolsSlurmParam`], [`gsvaMap`], [`gsva`]
#'
#' @examples
#' library(GSVA)
#'
#' p <- 1000 ## number of genes
#' n <- 300  ## number of samples
#' set.seed(123)
#' X <- matrix(rnorm(p*n), nrow=p,
#'             dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
#' gs <- as.list(sample(10:100, size=50, replace=TRUE))
#' gs <- lapply(gs, function(n, p) paste0("g", sample(1:p, size=n)), p)
#' names(gs) <- paste0("gs", 1:length(gs))
#'
#' ## memory of each job of gsvaMap() with two CPU cores
#' mem <- gsvaJobMemory(gsvaParam(X, gs, verbose=FALSE), ncpus_per_task=2)
#' mem
#' details(mem)
#'
#' ## memory of a single job running gsva() with the ssGSEA method
#' gsvaJobMemory(ssgseaParam(X, gs, verbose=FALSE))
#'
#' @importFrom cli cli_abort
#' @export
gsvaJobMemory <- function(param, ncpus_per_task=1,
                          map=is(param, "gsvaParam")) {
    if (!is(param, "GsvaMethodParam"))
        cli_abort(c("x"=paste("'param' must be a parameter object built with",
                              "gsvaParam(), plageParam(), zscoreParam(),",
                              "ssgseaParam() or avgParam().")))
    if (length(ncpus_per_task) != 1 || is.na(ncpus_per_task) ||
        ncpus_per_task < 1)
        cli_abort(c("x"="'ncpus_per_task' must be a positive integer."))
    if (!is.logical(map) || length(map) != 1 || is.na(map))
        cli_abort(c("x"="'map' must be either TRUE or FALSE."))
    if (map && !is(param, "gsvaParam"))
        cli_abort(c("x"=paste("gsvaMap() only runs the GSVA method, whose",
                              "parameter object is built with gsvaParam();",
                              "use 'map=FALSE' for a single job running",
                              "gsva().")))
    ncpus <- as.integer(ncpus_per_task)

    det <- if (map) .map_jobs_memory(param, ncpus)
           else .single_job_memory(param, ncpus)

    structure(.maxmem_text(max(det$recommended)),
              class=c("gsvaJobMemory", "character"), details=det,
              ncpus=ncpus, map=map)
}

#' @export
print.gsvaJobMemory <- function(x, ...) {
    print(as.character(x), ...)
    invisible(x)
}

setOldClass(c("gsvaJobMemory", "character"))

#' @aliases details,gsvaJobMemory-method
#' @rdname gsvaJobMemory
#' @param object For `details()`, the output of `gsvaJobMemory()`.
#' @exportMethod details
#' @importFrom GSEABase details
setMethod("details", signature=signature(object="gsvaJobMemory"),
          function(object) {
              det <- attr(object, "details")
              what <- if (attr(object, "map")) "each job of gsvaMap()"
                      else "a single job running gsva()"
              cat(sprintf("memory of %s with %d CPU core%s: %s\n", what,
                          attr(object, "ncpus"),
                          if (attr(object, "ncpus") > 1) "s" else "",
                          as.character(object)))
              tab <- det
              for (col in intersect(c("minimum", "recommended", "inmemory",
                                      "ondisk", "wholechunks"),
                                    colnames(tab)))
                  tab[[col]] <- vapply(det[[col]], function(x)
                                           if (is.na(x)) NA_character_
                                           else .maxmem_text(x),
                                       character(1))
              print(tab, row.names=FALSE)
              invisible(det)
          })

## memory in bytes of each job of gsvaMap() running each of the three steps of
## the GSVA method with the parameter object 'param', using 'ncpus' CPU cores,
## see gsvaJobMemory(), as a data frame with the smallest memory with which a
## job fits its chunk, the recommended one, the largest chunk of each job with
## the latter, in rows (by="rows") or columns (by="columns"), and their number.
## the input of the steps after the row normalization, which gsvaMap() would
## have written to disk, is a matrix of the same dimensions that takes no
## memory, see ConstantArray(), sparse when the row normalization is, with the
## fraction of nonzero values of the input, of doubles with the normalized
## values, and of integers with the ranks, except with the CLR on input
## stored on disk, whose normalized values are a delayed operation on it,
## which the ranks read through its chunks on disk, see .gsva_steps()
.map_jobs_memory <- function(param, ncpus) {
    steps <- .gsva_steps(param)

    do.call(rbind, lapply(steps, function(st)
                              .map_step_memory(st, ncpus,
                                               attr(steps, "clr"))))
}

## the three steps of the GSVA method with the parameter object 'param', as a
## list with, for each step, its name, the function running it, its input data
## and the name of its assay, its matrix, the dimension along which its jobs
## of gsvaMap() take their chunks, the fraction of nonzero values of sparse
## data, and the memory of the input data that each job holds, see
## .map_jobs_memory(), whose attribute "clr" tells whether the row
## normalization is the CLR
#' @importFrom DelayedArray ConstantArray
#' @importFrom S4Arrays is_sparse
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats setNames
.gsva_steps <- function(param) {
    assay <- get_assay(param)
    Xd <- unwrapData(get_exprData(param), assay)
    clr <- .get_rowNorm(param) == "clr"
    density <- .map_density(param, assay, Xd)
    densify <- is_sparse(Xd) && !.get_sparse(param) &&
               .get_rowNorm(param) != "none"
    outsparse <- is_sparse(Xd) && !densify
    dims <- dim(Xd)
    X2 <- if (clr && is(Xd, "DelayedArray")) Xd * 1
          else ConstantArray(dims, value=if (outsparse) 0 else 0.5)
    X3 <- ConstantArray(dims, value=if (outsparse) 0L else 1L)
    stepinput <- function(X, assayname) {
        se <- SummarizedExperiment(assays=setNames(list(X), assayname))
        metadata(se)$gsvaParam <- param
        se
    }
    ## the jobs of the steps after the row normalization receive the rest of
    ## the input data, whose assays in main memory they do not hold
    heldrest <- max(0, .input_held_mem(param) - .inmem_size(Xd))
    outdensity <- if (outsparse) density else NULL
    steps <- list(
        list(step="normalization", funname="gsvaRowNorm", inputData=param,
             assay=assay, X=Xd, whdim=if (clr) 2L else 1L, density=density,
             container=.map_container_mem(param, assay)),
        list(step="ranks", funname="gsvaColRanks",
             inputData=stepinput(X2, "gsvarnorm"), assay="gsvarnorm", X=X2,
             whdim=2L, density=outdensity, container=heldrest),
        list(step="scores", funname="gsvaColScores",
             inputData=stepinput(X3, "gsvaranks"), assay="gsvaranks", X=X3,
             whdim=2L, density=outdensity, container=heldrest))

    structure(steps, clr=clr)
}

## memory in bytes of each job of gsvaMap() running the step described in
## 'st', see .map_jobs_memory(), using 'ncpus' CPU cores: the smallest memory
## with which the memory required by a job, on the largest chunk that
## gsvaMap() gives it, see .map_grid(), fits in it, see .map_job_need(), and
## the recommended one, the smallest memory, not below the former, with which
## the blocks of the calculations of each job on its chunk, which its
## parallel workers share, see .job_block_units(), take whole chunks of the
## input data on disk, when it is stored in chunks of more than one row or
## column along the dimension of the chunks; the jobs calculating the
## statistics of the rows for the CLR take whole chunks anyway, see
## .map_grid(). the chunks are sized for one job at a time, whose chunks are
## not narrowed to give one to each of several jobs
#' @importFrom BiocGenerics dims
.map_step_memory <- function(st, ncpus, clr) {
    stepclr <- clr && st$funname == "gsvaRowNorm"
    grid <- function(mem)
        .map_grid(st$funname, st$inputData, st$assay, st$X, 1L, mem, stepclr,
                  st$density, st$container)
    width <- function(g) {
        w <- max(dims(g)[, st$whdim])
        nrb <- attr(g, "nrowband")
        if (!is.null(nrb)) ## bands of rows of whole chunks, see gsvaMap()
            w <- ceiling(nrb * w / nrow(st$X))
        w
    }
    need <- function(mem)
        .map_job_need(st$funname, st$inputData, st$assay, st$whdim,
                      width(grid(mem)), ncpus, mem, st$density, st$X,
                      st$container)
    baseline <- attr(need(Inf), "baseline")
    minimum <- .smallest_mem(function(mem) need(mem) <= mem, baseline)
    recommended <- minimum
    cd <- .ondisk_chunkdim(st$X)
    if (!stepclr && !is.null(cd) && !is.na(cd[st$whdim]) &&
        cd[st$whdim] > 1) {
        wc <- min(cd[st$whdim], dim(st$X)[st$whdim])
        whole <- .smallest_mem(function(mem) {
            w <- max(dims(grid(mem))[, st$whdim])
            w >= wc && .job_block_units(st, w, ncpus, mem) >= wc
        }, minimum)
        recommended <- max(minimum, whole, na.rm=TRUE)
    }
    g <- grid(recommended)

    data.frame(step=st$step, minimum=minimum, recommended=recommended,
               by=if (st$whdim == 1L) "rows" else "columns",
               chunk=max(dims(g)[, st$whdim]), chunks=length(g))
}

## rows (whdim=1) or columns (whdim=2) of the blocks in which each job of
## gsvaMap() with 'ncpus' CPU cores and 'mem' bytes of memory processes its
## chunk of 'w' rows or columns of the input of the step described in 'st',
## see .map_jobs_memory(), read from disk, as the step sizes them in the job,
## see .units_per_block(), with its parallel workers forked from the R process
## of the job, a new one, see .block_worker_mem(), unless 'session=TRUE', when
## the R process of the job takes the memory of this one, as a single job
## running gsva() on the same input data, and its output written to disk
## block by block
.job_block_units <- function(st, w, ncpus, mem, session=FALSE) {
    X <- if (st$whdim == 1L) st$X[seq_len(w), , drop=FALSE]
         else st$X[, seq_len(w), drop=FALSE]
    mf <- .ondisk_mem_factors(.map_step_factors(st$funname, st$inputData, X),
                              X, st$density)
    heldmem <- max(0, st$container - .inmem_size(X))
    nworkers <- max(1L, .map_job_workers(ncpus, dim(X)))
    if (session) {
        pm <- .block_worker_mem(X, st$whdim, heldmem, forkcopy=mf$forkcopy)
        workermem <- pm$workermem
        parworkermem <- pm$parworkermem
        copiedmem <- pm$copiedmem
    } else {
        readmem <- .disk_read_mem()
        workermem <- parworkermem <- .worker_mem(session=FALSE) + readmem
        copiedmem <- .fork_copied_mem(max(0, workermem - readmem -
                                             .proc_extra_mem()),
                                      mf$forkcopy, inmemory=FALSE)
    }

    .units_per_block(X, st$whdim, nworkers, mem, mf$workfactor, mf$outfactor,
                     mf$outextra, workermem=workermem, dense=mf$densify,
                     heldmem=heldmem, sinkout=TRUE, copiedmem=copiedmem,
                     parworkermem=parworkermem)
}

## smallest memory in bytes, not below 'start', for which 'fits' gives TRUE,
## assuming that it does for any larger memory, found by doubling 'start' and
## then bisecting, to 1% of it, or NA when no memory up to 2^30 times 'start'
## fits
.smallest_mem <- function(fits, start) {
    lo <- hi <- max(start, 1024^2)
    i <- 0L
    while (!fits(hi)) {
        if (i == 30L)
            return(NA_real_)
        lo <- hi
        hi <- 2 * hi
        i <- i + 1L
    }
    if (lo == hi)
        return(hi)
    while (hi - lo > 0.01 * hi) {
        mid <- (lo + hi) / 2
        if (fits(mid))
            hi <- mid
        else
            lo <- mid
    }

    hi
}

## memory in bytes of a single job running gsva() with the parameter object
## 'param', using 'ncpus' CPU cores as parallel workers, see gsvaJobMemory(),
## as a data frame with the smallest memory and the recommended one: the
## memory estimated by gsva() for processing the input data in main memory
## (ondisk="no") and, when it is stored on disk, from disk (ondisk="yes"),
## see .dry_run_mem_need(), the smallest of which is the former, and the
## recommended one the one given by the argument 'ondisk' of 'param', where
## "auto" gives the one in main memory. with the GSVA method processing from
## disk input data stored in chunks of several rows or columns, the
## recommended memory is at least the one with which the blocks of its steps
## take whole chunks, see .whole_chunks_mem()
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
.single_job_memory <- function(param, ncpus) {
    X <- unwrapData(get_exprData(param), get_assay(param))
    bp <- SerialParam()
    if (ncpus > 1)
        bp <- if (.Platform$OS.type == "unix") MulticoreParam(workers=ncpus)
              else SnowParam(workers=ncpus)
    needs <- c(no=.dry_run_mem_need(param, bp, "no"))
    if (is(X, "DelayedArray"))
        needs <- c(needs, yes=.dry_run_mem_need(param, bp, "yes"))
    ondisk <- .get_ondisk(param)
    whole <- NA_real_
    recommended <- needs[["no"]]
    if (ondisk == "yes" && "yes" %in% names(needs)) {
        recommended <- needs[["yes"]]
        if (is(param, "gsvaParam")) {
            whole <- .whole_chunks_mem(param, ncpus, needs[["yes"]])
            recommended <- max(recommended, whole, na.rm=TRUE)
        }
    }

    data.frame(step="gsva()", minimum=min(needs), recommended=recommended,
               inmemory=needs[["no"]],
               ondisk=if ("yes" %in% names(needs)) needs[["yes"]] else NA,
               wholechunks=whole)
}

## smallest memory in bytes, not below 'start', with which the blocks of the
## steps of the GSVA method with the parameter object 'param', processing its
## input data from disk with 'ncpus' parallel workers in a single job, take
## whole chunks of the data on disk, along the dimension of their blocks, with
## the memory of this R process, see .job_block_units(), or NA when the data on disk of no step are stored in
## chunks of more than one row or column along it; the statistics of the rows
## for the CLR take whole chunks anyway, see .rowStats()
.whole_chunks_mem <- function(param, ncpus, start) {
    steps <- .gsva_steps(param)
    clr <- attr(steps, "clr")
    mems <- vapply(steps, function(st) {
        cd <- .ondisk_chunkdim(st$X)
        if ((clr && st$funname == "gsvaRowNorm") || is.null(cd) ||
            is.na(cd[st$whdim]) || cd[st$whdim] <= 1)
            return(NA_real_)
        n <- dim(st$X)[st$whdim]
        wc <- min(cd[st$whdim], n)
        .smallest_mem(function(mem)
                          .job_block_units(st, n, ncpus, mem,
                                           session=TRUE) >= wc, start)
    }, numeric(1))
    if (all(is.na(mems)))
        return(NA_real_)

    max(mems, na.rm=TRUE)
}

## memory in bytes estimated by gsva() with the parameter object 'param' and
## the parallel back-end 'BPPARAM', with its argument 'ondisk' set to
## 'ondisk', which gsva() gives in a memory check before any calculation, see
## .check_mem_need(), where it is stopped, and with a maximum memory large
## enough for the input data to be processed in main memory with
## ondisk="no". the warning of more workers than CPU cores is not given,
## since the job may run on a machine with more cores than this one
.dry_run_mem_need <- function(param, BPPARAM, ondisk) {
    param@ondisk <- ondisk
    olddry <- gsva_global$memdryrun
    oldcores <- gsva_global$check_cores
    oldcheck <- gsva_global$check_memory
    oldopt <- options(GSVA.check_memory=TRUE)
    gsva_global$memdryrun <- TRUE
    gsva_global$check_cores <- FALSE
    gsva_global$check_memory <- TRUE
    on.exit({
        gsva_global$memdryrun <- olddry
        gsva_global$check_cores <- oldcores
        gsva_global$check_memory <- oldcheck
        options(oldopt)
    })
    need <- tryCatch({
        suppressMessages(gsva(param, verbose=FALSE, BPPARAM=BPPARAM,
                              maxmem=2^50))
        NA_real_
    }, gsvaMemDryRun=function(cond) cond$need)
    if (is.na(need))
        cli_abort(c("x"=paste("Internal error, gsva() gave no estimate of",
                              "the memory that it requires.")))

    need
}
