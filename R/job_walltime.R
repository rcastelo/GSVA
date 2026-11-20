#' @title Wall time of the jobs running GSVA
#'
#' @description Estimates the wall time that each job running GSVA requires
#' on a workload manager, such as SLURM, to give it, e.g., in the argument
#' `walltime` of [`gsvaBatchtoolsSlurmParam`], by running a small part of the
#' calculations in this R session, with the same memory and number of parallel
#' workers as the jobs, and extrapolating its time to the whole of each job.
#'
#' @param param A parameter object of one of the classes
#' [`gsvaParam`][gsvaParam-class], [`plageParam`][plageParam-class],
#' [`zscoreParam`][zscoreParam-class], [`ssgseaParam`][ssgseaParam-class] or
#' [`avgParam`][avgParam-class], built from the input data of the jobs.
#'
#' @param mem Memory of each job, as the argument of the same name of
#' [`gsvaBatchtoolsSlurmParam`], e.g., the output of [`gsvaJobMemory`].
#'
#' @param ncpus_per_task Number of CPU cores of each job, which GSVA uses as
#' parallel workers, as the argument of the same name of
#' [`gsvaBatchtoolsSlurmParam`].
#'
#' @param nodes Number of jobs running at the same time, as the argument of
#' the same name of [`gsvaBatchtoolsSlurmParam`], among which [`gsvaMap`]
#' distributes the chunks of the input data.
#'
#' @param arrayjobs As the argument of the same name of
#' [`gsvaBatchtoolsSlurmParam`]: when `TRUE`, the wall time is the one of each
#' chunk, processed by a task of an array job, and otherwise the one of each
#' job, which processes its chunks one after another.
#'
#' @param margin Factor by which the estimated wall time is multiplied, to
#' cover the variability of the time of the calculations, e.g., due to the
#' load of a shared file system, since SLURM ends the jobs that exceed their
#' wall time. Default: 2.
#'
#' @param map When `TRUE`, the default with a [`gsvaParam`][gsvaParam-class]
#' object, the wall time of each of the jobs that [`gsvaMap`] sends to run the
#' three steps of the GSVA method on chunks of the input data. When `FALSE`,
#' the default with other parameter objects, the wall time of a single job
#' running [`gsva`] on the whole input data.
#'
#' @param verbose Gives information about the parts of the calculations that
#' are run. Default: `TRUE`.
#'
#' @details The time of the calculations depends on the CPU and on the file
#' system where they run, so that `gsvaJobWalltime()` should run on a machine
#' like the compute nodes of the jobs, e.g., in an interactive session of
#' SLURM on one of them with as many CPU cores as each job, instead of on the
#' login node of a computer cluster, whose CPUs are often slower and shared.
#'
#' For the jobs of [`gsvaMap`], the chunks of each step are the ones that it
#' gives with the memory `mem` and `nodes` jobs, see [`gsvaJobMemory`], each of
#' which processes its chunk in blocks whose width depends on that memory and
#' on its parallel workers. Each step is run on two parts of the input data,
#' of one and two rounds of such blocks, saving its output to disk as the jobs
#' do, which gives the time of each chunk as a fixed time plus a time
#' proportional to its width; the steps after the row normalization are run
#' on the row-normalized values and the ranks of a part of the input data,
#' calculated beforehand. To the time of the chunks of each job, the time of
#' starting the R process of the job, loading GSVA and its input data, is
#' added, also measured. The wall time returned is the largest one of the
#' three steps, so that the same jobs can run all of them, multiplied by
#' `margin`.
#'
#' For a single job running [`gsva`] with the GSVA method, each of its three
#' steps is run in the same way, on parts of the whole input data, and the
#' wall time is the sum of their times. With the other methods, [`gsva`] is
#' run on two parts of the columns of the input data, whose time is
#' extrapolated to all of them; its blocks may then be wider than the ones of
#' the whole input data, which may give a shorter time.
#'
#' The estimate takes the time of running a part of the calculations, in the
#' order of a minute with millions of cells, and its error is larger than the
#' one of the memory of [`gsvaJobMemory`], hence the `margin`.
#'
#' @return A character string with the estimated wall time, in hours,
#' minutes and seconds rounded up to minutes, e.g., `"1:30:00"`, which can
#' be given to the argument `walltime` of [`gsvaBatchtoolsSlurmParam`]. It has
#' the class `gsvaJobWalltime`, and the function `details()` shows, for each
#' step, the width and number of the chunks, the time of each chunk, the
#' chunks of each job and its time, before multiplying it by `margin`.
#'
#' @seealso [`gsvaJobMemory`], [`gsvaBatchtoolsSlurmParam`], [`gsvaMap`],
#' [`gsva`]
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
#' ## wall time of each job of gsvaMap() with one CPU core and 4 GB
#' wt <- gsvaJobWalltime(gsvaParam(X, gs, verbose=FALSE), mem="4G",
#'                       verbose=FALSE)
#' wt
#' details(wt)
#'
#' @importFrom cli cli_abort cli_warn
#' @export
gsvaJobWalltime <- function(param, mem, ncpus_per_task=1, nodes=2,
                            arrayjobs=FALSE, margin=2,
                            map=is(param, "gsvaParam"), verbose=TRUE) {
    if (!is(param, "GsvaMethodParam"))
        cli_abort(c("x"=paste("'param' must be a parameter object built with",
                              "gsvaParam(), plageParam(), zscoreParam(),",
                              "ssgseaParam() or avgParam().")))
    memb <- .memtext2bytes(as.character(mem))
    if (length(memb) != 1 || is.na(memb) || memb <= 0)
        cli_abort(c("x"="'mem' must be a positive amount of memory."))
    for (a in c("ncpus_per_task", "nodes")) {
        v <- get(a)
        if (length(v) != 1 || is.na(v) || v < 1)
            cli_abort(c("x"="'{a}' must be a positive integer."))
    }
    if (!is.logical(arrayjobs) || length(arrayjobs) != 1 || is.na(arrayjobs))
        cli_abort(c("x"="'arrayjobs' must be either TRUE or FALSE."))
    if (!is.numeric(margin) || length(margin) != 1 || is.na(margin) ||
        margin < 1)
        cli_abort(c("x"="'margin' must be a number not smaller than 1."))
    if (!is.logical(map) || length(map) != 1 || is.na(map))
        cli_abort(c("x"="'map' must be either TRUE or FALSE."))
    if (map && !is(param, "gsvaParam"))
        cli_abort(c("x"=paste("gsvaMap() only runs the GSVA method, whose",
                              "parameter object is built with gsvaParam();",
                              "use 'map=FALSE' for a single job running",
                              "gsva().")))
    ncpus <- as.integer(ncpus_per_task)
    nodes <- as.integer(nodes)

    minmem <- max(attr(gsvaJobMemory(param, ncpus, map), "details")$minimum)
    if (memb < minmem) {
        minmemtxt <- .maxmem_text(minmem)
        cli_warn(c("!"=paste("'mem' is smaller than the memory that each job",
                             "requires, {minmemtxt}, see {.fn gsvaJobMemory},",
                             "and the jobs may fail.")))
    }

    startsecs <- .job_start_secs(param)
    if (map || is(param, "gsvaParam")) {
        ## a single job running gsva() with the GSVA method runs its three
        ## steps on the whole input data, one after another
        det <- .map_jobs_walltime(param, memb, ncpus, nodes, verbose,
                                  single=!map)
    } else
        det <- .single_job_walltime(param, memb, ncpus, verbose)
    ## each task of an array job processes one chunk in its own R process,
    ## and each job of gsvaMap() otherwise its chunks one after another
    det$perjob <- if (arrayjobs || !map) 1 else ceiling(det$chunks / nodes)
    det$jobsecs <- det$perjob * det$chunksecs + startsecs
    secs <- if (map) max(det$jobsecs) else sum(det$chunksecs) + startsecs

    structure(.walltime_text(margin * secs),
              class=c("gsvaJobWalltime", "character"), details=det,
              startsecs=startsecs, margin=margin, ncpus=ncpus, nodes=nodes,
              arrayjobs=arrayjobs, map=map)
}

#' @export
print.gsvaJobWalltime <- function(x, ...) {
    print(as.character(x), ...)
    invisible(x)
}

setOldClass(c("gsvaJobWalltime", "character"))

#' @aliases details,gsvaJobWalltime-method
#' @rdname gsvaJobWalltime
#' @param object For `details()`, the output of `gsvaJobWalltime()`.
#' @exportMethod details
#' @importFrom GSEABase details
setMethod("details", signature=signature(object="gsvaJobWalltime"),
          function(object) {
              det <- attr(object, "details")
              what <- if (!attr(object, "map")) "a single job running gsva()"
                      else if (attr(object, "arrayjobs"))
                          "each task of the array jobs of gsvaMap()"
                      else sprintf("each of the %d jobs of gsvaMap()",
                                   attr(object, "nodes"))
              cat(sprintf(paste("wall time of %s with %d CPU core%s: %s,",
                                "%s times the estimated time\n"), what,
                          attr(object, "ncpus"),
                          if (attr(object, "ncpus") > 1) "s" else "",
                          as.character(object), attr(object, "margin")))
              cat(sprintf(paste("time of starting the R process of a job and",
                                "loading GSVA and its input data: %.1f s\n"),
                          attr(object, "startsecs")))
              tab <- det
              for (col in c("chunksecs", "jobsecs"))
                  tab[[col]] <- round(tab[[col]], 1)
              print(tab, row.names=FALSE)
              invisible(det)
          })

## wall time 'secs' in seconds as hours, minutes and seconds, rounded up to
## whole minutes, as the argument 'walltime' of gsvaBatchtoolsSlurmParam()
## takes it, e.g., "1:30:00"
.walltime_text <- function(secs) {
    mins <- max(1, ceiling(secs / 60))
    sprintf("%d:%02d:00", mins %/% 60, mins %% 60)
}

## time in seconds of starting the R process of a job and loading GSVA, by
## running Rscript once in this R session, and of loading the input data of the
## job, the parameter object 'param', which batchtools sends to the jobs
## serialized
.job_start_secs <- function(param) {
    rsecs <- gsva_global$rscriptsecs
    if (is.null(rsecs)) {
        rscript <- file.path(R.home("bin"), "Rscript")
        rsecs <- system.time(
            system2(rscript,
                    c("-e", shQuote(paste("suppressPackageStartupMessages(",
                                          "library(GSVA))"))),
                    stdout=FALSE, stderr=FALSE))[["elapsed"]]
        gsva_global$rscriptsecs <- rsecs
    }
    ser <- serialize(param, NULL)
    lsecs <- system.time(unserialize(ser))[["elapsed"]]

    rsecs + lsecs
}

## time in seconds of the jobs of gsvaMap() running each of the three steps of
## the GSVA method with the parameter object 'param', using 'ncpus' CPU cores
## and 'mem' bytes of memory each, with 'nodes' jobs at the same time, or,
## with 'single=TRUE', of a single job running each step on the whole input
## data, see gsvaJobWalltime(), as a data frame with, for each step, the
## largest chunk
## of each job, in rows (by="rows") or columns (by="columns"), the number of
## chunks, and the time of the largest chunk, see .chunk_secs(). the steps
## after the row normalization are run on the row-normalized values and the
## ranks of the columns of the input data taken by their parts, calculated
## beforehand and saved to disk, unless the row normalization is the CLR,
## whose values are a delayed operation on the input data, as the jobs of
## gsvaMap() would read them, see .probe_inputs()
#' @importFrom BiocGenerics dims
#' @importFrom IRanges IRanges
.map_jobs_walltime <- function(param, mem, ncpus, nodes, verbose,
                               single=FALSE) {
    steps <- .gsva_steps(param)
    clr <- attr(steps, "clr")
    tmpdir <- tempfile("gsvawalltime")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive=TRUE), add=TRUE)

    ## the chunks of each step and the widths of the parts of the input data
    ## on which it is run
    geom <- lapply(steps, function(st) {
        stepclr <- clr && st$funname == "gsvaRowNorm"
        if (single) {
            w <- dim(st$X)[st$whdim]
            nchunks <- 1L
        } else {
            g <- .map_grid(st$funname, st$inputData, st$assay, st$X, nodes,
                           mem, stepclr, st$density, st$container)
            w <- max(dims(g)[, st$whdim])
            nchunks <- length(g)
        }
        list(chunk=w, chunks=nchunks,
             parts=.probe_widths(st, w, ncpus, mem, stepclr))
    })
    ## the input of the steps after the row normalization, on the columns that
    ## their parts take
    inputs <- .probe_inputs(param, sum(geom[[2]]$parts),
                            sum(geom[[3]]$parts), ncpus, mem, clr, tmpdir,
                            verbose)

    secs <- vapply(seq_along(steps), function(i) {
        st <- steps[[i]]
        if (verbose)
            cli_alert_info(paste("Running the step of {st$step} on",
                                 "{sum(geom[[i]]$parts)} {if (st$whdim == 1L)",
                                 "'rows' else 'columns'}"))
        run <- function(first, last) {
            chunk <- IRanges(start=first, end=last)
            args <- list(output="object", ncpus=ncpus, maxmem=mem)
            if (st$funname == "gsvaRowNorm" && clr)
                args <- c(args, list(WRAPPED_FUN=.map_rowstats_chunk,
                                     param=param))
            else {
                ## the jobs of gsvaMap() save their output to disk, unlike a
                ## single job running gsva()
                if (!single) {
                    chunk <- structure(list(chunk=chunk,
                                            fname=tempfile(tmpdir=tmpdir,
                                                           fileext=".h5")),
                                       class="gsvaMapChunk")
                    args$output <- "HDF5"
                }
                args <- c(args, switch(st$funname,
                    gsvaRowNorm=list(WRAPPED_FUN=gsvaRowNorm, param=param,
                                     dropExistingAssays=TRUE,
                                     errorOnTooFewRows=FALSE),
                    gsvaColRanks=list(WRAPPED_FUN=gsvaColRanks,
                                      rowNormExprData=inputs$rnorm,
                                      dropExistingAssays=TRUE),
                    gsvaColScores=list(WRAPPED_FUN=gsvaColScores,
                                       rankExprData=inputs$ranks,
                                       device="cpu")))
            }
            system.time(suppressMessages(do.call(MAP_FUN_WRAPPER,
                                                 c(list(X=chunk),
                                                   args))))[["elapsed"]]
        }
        .chunk_secs(geom[[i]]$parts, geom[[i]]$chunk, run)
    }, numeric(1))

    data.frame(step=vapply(steps, `[[`, character(1), "step"),
               by=vapply(steps, function(st)
                             if (st$whdim == 1L) "rows" else "columns",
                         character(1)),
               chunk=vapply(geom, `[[`, numeric(1), "chunk"),
               chunks=vapply(geom, `[[`, numeric(1), "chunks"),
               chunksecs=secs)
}

## widths of the two parts of the input of the step described in 'st', see
## .gsva_steps(), on which it is run to estimate the time of its chunks of
## width 'w', see .chunk_secs(), processed by jobs with 'ncpus' CPU cores and
## 'mem' bytes of memory: one and two rounds of the blocks that the parallel
## workers of each job take, see .job_block_units(), aligned to the chunks of
## the input data on disk, when its blocks span them, as in .colgridsize(),
## and at most the width of the chunks, of which two parts take at most the
## whole input data. the parts take at most 2% of the input data, unless one
## round of blocks of whole chunks of the input data on disk is wider, since
## blocks narrower than those chunks would read them several times, which
## makes the calculations slower, while otherwise narrower blocks take about
## the same time per row or column. with 'clr=TRUE', the statistics of the
## rows for the CLR take whole chunks of the input data on disk, when it is
## stored in chunks, in bands of rows, see .colgridsize()
.probe_widths <- function(st, w, ncpus, mem, clr) {
    n <- dim(st$X)[st$whdim]
    nworkers <- max(1L, .map_job_workers(ncpus, dim(st$X)))
    b <- .job_block_units(st, w, ncpus, mem)
    cd <- .ondisk_chunkdim(st$X)
    if (clr && !is.null(cd) && !is.na(cd[st$whdim]))
        b <- max(b, cd[st$whdim])
    b <- .align_to_chunks(min(b, w), st$X, st$whdim)
    cap <- ceiling(n / 50)
    if (!is.null(cd) && !is.na(cd[st$whdim]) && cd[st$whdim] > 1 &&
        b >= cd[st$whdim])
        cap <- max(cap, nworkers * b)
    p1 <- min(w, nworkers * b, cap)
    p2 <- min(w, 2 * p1)
    if (p1 + p2 > n) { ## the whole input data in at most two parts
        p1 <- min(w, n)
        p2 <- min(w, n - p1)
    }

    c(p1, p2)[c(p1, p2) > 0]
}

## time in seconds of a chunk of width 'w', from the times of 'run' on
## consecutive parts of the input data with the widths in 'parts', starting
## at its first row or column: a fixed time plus a time proportional to the
## width, fitted to the times of two parts of different widths, or
## proportional to the width of a single part, either one with no negative
## time
.chunk_secs <- function(parts, w, run) {
    ends <- cumsum(parts)
    starts <- ends - parts + 1
    secs <- mapply(run, starts, ends)
    if (length(parts) < 2 || parts[1] == parts[2])
        return(max(secs) * w / max(parts))

    slope <- (secs[2] - secs[1]) / (parts[2] - parts[1])
    fixed <- secs[1] - slope * parts[1]
    if (slope < 0 || fixed < 0) {
        slope <- max(secs / parts)
        fixed <- 0
    }

    fixed + slope * w
}

## the inputs of the ranks and the scores of the GSVA method with the
## parameter object 'param', calculated with 'ncpus' parallel workers within
## 'mem' bytes of memory, see .map_jobs_walltime(): the row-normalized values
## of its first 'ncolrnorm' columns, saved to disk in the directory 'dir', as
## the jobs of gsvaMap() would read them, unless they are the delayed
## operation on the input data of the CLR, and the ranks of the first
## 'ncolranks' of them, also saved to disk
.probe_inputs <- function(param, ncolrnorm, ncolranks, ncpus, mem, clr, dir,
                          verbose) {
    bp <- .map_job_backend(ncpus)
    if (verbose)
        cli_alert_info(paste("Calculating the row-normalized values of",
                             "{ncolrnorm} columns and the ranks of",
                             "{ncolranks}"))
    ondisk <- function(x) {
        f <- tempfile(tmpdir=dir, fileext=".h5")
        saveHDF5GSVA(x, f)
        loadHDF5GSVA(f, verbose=FALSE)
    }
    rownorm <- function(ncol) {
        p <- .subset_param_cols(param, seq_len(ncol))
        rnorm <- suppressMessages(gsvaRowNorm(p, verbose=FALSE,
                                              dropExistingAssays=TRUE,
                                              BPPARAM=bp, maxmem=mem))
        if (clr) rnorm else ondisk(rnorm)
    }
    rnorm <- rownorm(ncolrnorm)
    ranks <- suppressMessages(gsvaColRanks(rownorm(ncolranks), verbose=FALSE,
                                           dropExistingAssays=TRUE,
                                           BPPARAM=bp, maxmem=mem))

    list(rnorm=rnorm, ranks=ondisk(ranks))
}

## the parameter object 'param' with the columns 'cols' of its input data,
## whose number of nonzero values is the one of 'param' in proportion to
## those columns
.subset_param_cols <- function(param, cols) {
    X <- unwrapData(get_exprData(param), get_assay(param))
    param@exprData <- get_exprData(param)[, cols]
    param@nzcount <- param@nzcount * length(cols) / max(1, ncol(X))
    param
}

## time in seconds of a single job running gsva() with a parameter object
## 'param' of a method other than GSVA, using 'ncpus' CPU cores and 'mem' bytes of memory, see
## gsvaJobWalltime(), as a data frame with the time of gsva() on all the
## columns of the input data, extrapolated from its time on two parts of
## them, see .chunk_secs(), whose widths are the larger of 1,000 columns and
## 2% of them, in whole chunks of the input data on disk, if any
#' @importFrom BiocParallel SerialParam
.single_job_walltime <- function(param, mem, ncpus, verbose) {
    X <- unwrapData(get_exprData(param), get_assay(param))
    n <- ncol(X)
    bp <- if (ncpus > 1) .map_job_backend(ncpus) else SerialParam()
    base <- max(1000, ceiling(n / 50))
    cd <- .ondisk_chunkdim(X)
    if (!is.null(cd) && !is.na(cd[2]) && cd[2] > 1)
        base <- cd[2] * ceiling(base / cd[2])
    p1 <- min(n, base)
    p2 <- min(n - p1, 2 * base)
    parts <- c(p1, p2)[c(p1, p2) > 0]
    if (verbose)
        cli_alert_info("Running gsva() on {sum(parts)} columns")
    run <- function(first, last) {
        p <- .subset_param_cols(param, seq.int(first, last))
        system.time(suppressMessages(gsva(p, verbose=FALSE, BPPARAM=bp,
                                          maxmem=mem)))[["elapsed"]]
    }

    data.frame(step="gsva()", by="columns", chunk=n, chunks=1,
               chunksecs=.chunk_secs(parts, n, run))
}
