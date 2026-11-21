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
#' @param cold When `TRUE`, the default, a small part of each step is also
#' run in a new R process, as each job runs it, to measure what it does the
#' first time, e.g., loading the packages that the step uses, see below. When
#' `FALSE`, this time is not included, and only the time of starting a new R
#' process and loading GSVA is added to the time of each job.
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
#' added, as well as the time of what the R process does the first time that
#' it runs the step, such as loading the packages that the step uses. With
#' `cold=TRUE`, both are measured by running a small part of each step in a
#' new R process, all of them at the same time, which, e.g., loads the
#' annotation package that maps the identifiers of the genes of the gene sets
#' to the ones of the input data only when they are of different kinds, as
#' the jobs would do. The wall time returned is the largest one of the three
#' steps, so that the same jobs can run all of them, multiplied by `margin`.
#' The time that the workload manager and the package batchtools take to start
#' each job, and to read and write its files on a shared file system, is not
#' included, and `margin` should cover it.
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
#' step, the width and number of the chunks, the additional time of the
#' first chunk of a job in a new R process, the time of each chunk, the
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
                            map=is(param, "gsvaParam"), cold=TRUE,
                            verbose=TRUE) {
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
    if (!is.logical(cold) || length(cold) != 1 || is.na(cold))
        cli_abort(c("x"="'cold' must be either TRUE or FALSE."))
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

    if (map || is(param, "gsvaParam")) {
        ## a single job running gsva() with the GSVA method runs its three
        ## steps on the whole input data, one after another
        det <- .map_jobs_walltime(param, memb, ncpus, nodes, cold, verbose,
                                  single=!map)
    } else
        det <- .single_job_walltime(param, memb, ncpus, cold, verbose)
    startsecs <- attr(det, "startsecs")
    ## each task of an array job processes one chunk in its own R process,
    ## and each job of gsvaMap() otherwise its chunks one after another, the
    ## first of which takes in addition the time of what a new R process
    ## does the first time, such as loading the packages that the step uses
    det$perjob <- if (arrayjobs || !map) 1 else ceiling(det$chunks / nodes)
    det$jobsecs <- startsecs + det$firstsecs + det$perjob * det$chunksecs
    secs <- if (map) max(det$jobsecs)
            else startsecs + sum(det$firstsecs + det$chunksecs)

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
              for (col in c("firstsecs", "chunksecs", "jobsecs"))
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
## largest chunk of each job, in rows (by="rows") or columns (by="columns"),
## the number of chunks, the time of the largest chunk in this R session, see
## .chunk_secs(), and the additional time of the first chunk of a job in a new
## R process ("firstsecs"), whose attribute "startsecs" gives the time of
## starting it, see .cold_secs(). the steps after the row normalization are
## run on the row-normalized values and the ranks of the columns of the input
## data taken by their parts, calculated beforehand and saved to disk, unless
## the row normalization is the CLR, whose values are a delayed operation on
## the input data, as the jobs of gsvaMap() would read them, see
## .probe_inputs(). with 'cold=FALSE', no new R process runs a step, and the
## time of starting one is measured by loading GSVA, see .job_start_secs()
#' @importFrom BiocGenerics dims
#' @importFrom IRanges IRanges
.map_jobs_walltime <- function(param, mem, ncpus, nodes, cold, verbose,
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
        parts <- .probe_widths(st, w, ncpus, mem, stepclr)
        list(chunk=w, chunks=nchunks, parts=parts,
             tiny=min(parts[1], if (st$whdim == 1L) 20 else 100))
    })
    ## the input of the steps after the row normalization, on the columns that
    ## their parts take
    inputs <- .probe_inputs(param, sum(geom[[2]]$parts),
                            sum(geom[[3]]$parts), ncpus, mem, clr, tmpdir,
                            verbose)

    ## the arguments of MAP_FUN_WRAPPER() running the i-th step on the rows
    ## or columns from 'first' to 'last', as a job of gsvaMap() does
    stepargs <- function(i, first, last) {
        st <- steps[[i]]
        chunk <- IRanges(start=first, end=last)
        args <- list(output="object", ncpus=ncpus, maxmem=mem)
        if (st$funname == "gsvaRowNorm" && clr)
            return(c(list(X=chunk), args,
                     list(WRAPPED_FUN=.map_rowstats_chunk, param=param)))
        ## the jobs of gsvaMap() save their output to disk, unlike a single
        ## job running gsva()
        if (!single) {
            chunk <- structure(list(chunk=chunk,
                                    fname=tempfile(tmpdir=tmpdir,
                                                   fileext=".h5")),
                               class="gsvaMapChunk")
            args$output <- "HDF5"
        }
        c(list(X=chunk), args, switch(st$funname,
            gsvaRowNorm=list(WRAPPED_FUN=gsvaRowNorm, param=param,
                             dropExistingAssays=TRUE,
                             errorOnTooFewRows=FALSE),
            gsvaColRanks=list(WRAPPED_FUN=gsvaColRanks,
                              rowNormExprData=inputs$rnorm,
                              dropExistingAssays=TRUE),
            gsvaColScores=list(WRAPPED_FUN=gsvaColScores,
                               rankExprData=inputs$ranks, device="cpu")))
    }
    calls <- lapply(seq_along(steps), function(i)
                        list(fun=MAP_FUN_WRAPPER,
                             args=stepargs(i, 1, geom[[i]]$tiny)))
    tiny <- .cold_secs(calls, param, cold, tmpdir, verbose)

    secs <- vapply(seq_along(steps), function(i) {
        st <- steps[[i]]
        if (verbose)
            cli_alert_info(paste("Running the step of {st$step} on",
                                 "{sum(geom[[i]]$parts)} {if (st$whdim == 1L)",
                                 "'rows' else 'columns'}"))
        run <- function(first, last)
            system.time(suppressMessages(do.call(MAP_FUN_WRAPPER,
                                                 stepargs(i, first,
                                                          last))))[["elapsed"]]
        .chunk_secs(geom[[i]]$parts, geom[[i]]$chunk, run)
    }, numeric(1))

    structure(data.frame(step=vapply(steps, `[[`, character(1), "step"),
                         by=vapply(steps, function(st)
                                       if (st$whdim == 1L) "rows"
                                       else "columns", character(1)),
                         chunk=vapply(geom, `[[`, numeric(1), "chunk"),
                         chunks=vapply(geom, `[[`, numeric(1), "chunks"),
                         firstsecs=tiny$firstsecs, chunksecs=secs),
              startsecs=attr(tiny, "startsecs"))
}

## time in seconds of starting a new R process and loading GSVA and the input
## data of a job, and of what the R process of each job does the first time
## it runs each of the 'calls', a list of lists with a function 'fun' and its
## arguments 'args', each running a step on a small part of the input data of
## the parameter object 'param', such as loading the packages that it uses,
## e.g., an annotation package to map the identifiers of the genes of the
## gene sets to the ones of the input data, when they are of different
## kinds: each call is run in a new R process, all of them at the same time,
## see .cold_runs(), and twice in this R session, which then runs them as
## they are run after the first time, and the first-time time of each call
## is the time of its new R process minus the shorter one of this R session.
## the time of starting a new R process is the median of the ones of the
## calls. with 'cold=FALSE', the calls are not run in new R processes, the
## first-time times are zero, and the time of starting a new R process is
## measured by loading GSVA, see .job_start_secs(). 'dir' is a directory for
## the temporary files of the new R processes
#' @importFrom stats median
.cold_secs <- function(calls, param, cold, dir, verbose) {
    startsecs <- NA_real_
    coldsecs <- rep(NA_real_, length(calls))
    if (cold) {
        if (verbose)
            cli_alert_info(paste("Running a small part of each step in a new",
                                 "R process"))
        cr <- .cold_runs(calls, dir)
        startsecs <- median(cr$startsecs, na.rm=TRUE)
        coldsecs <- cr$runsecs
    }
    if (is.na(startsecs))
        startsecs <- .job_start_secs(param)
    ## this R session runs each call after the first time
    warmsecs <- vapply(calls, function(cl)
                           min(replicate(2, system.time(suppressMessages(
                               do.call(cl$fun, cl$args)))[["elapsed"]])),
                       numeric(1))

    structure(data.frame(firstsecs=ifelse(is.na(coldsecs), 0,
                                          pmax(0, coldsecs - warmsecs))),
              startsecs=startsecs)
}

## time in seconds of starting a new R process, loading GSVA and reading the
## arguments of each of the 'calls', a list of lists with a function 'fun'
## and its arguments 'args' ("startsecs"), and of running the call
## ("runsecs"), run each in a new R process with the library paths of this R
## session, all of them at the same time, as a data frame with NA for the
## calls that give an error. 'dir' is a directory for the temporary files
## with the calls and the times of the new R processes, which give an error
## when they do not end within 'timeout' seconds
#' @importFrom cli cli_abort
.cold_runs <- function(calls, dir, timeout=7200) {
    rscript <- file.path(R.home("bin"), "Rscript")
    code <- paste(
        "a <- commandArgs(TRUE); t <- rep(NA_real_, 3); tryCatch({",
        "suppressPackageStartupMessages(library(GSVA));",
        "t[1] <- as.numeric(Sys.time()); x <- readRDS(a[1]);",
        "t[2] <- as.numeric(Sys.time());",
        "invisible(suppressWarnings(suppressMessages(do.call(x$fun,",
        "x$args)))); t[3] <- as.numeric(Sys.time())},",
        "error=function(e) NULL);",
        "writeLines(format(t, digits=15), paste0(a[2], \".tmp\"));",
        "invisible(file.rename(paste0(a[2], \".tmp\"), a[2]))")
    infiles <- file.path(dir, sprintf("coldcall%d.rds", seq_along(calls)))
    outfiles <- file.path(dir, sprintf("coldtimes%d.txt", seq_along(calls)))
    oldlibs <- Sys.getenv("R_LIBS", unset=NA)
    Sys.setenv(R_LIBS=paste(.libPaths(), collapse=.Platform$path.sep))
    on.exit(if (is.na(oldlibs)) Sys.unsetenv("R_LIBS")
            else Sys.setenv(R_LIBS=oldlibs), add=TRUE)
    launched <- numeric(length(calls))
    for (i in seq_along(calls)) {
        saveRDS(calls[[i]], infiles[i])
        launched[i] <- as.numeric(Sys.time())
        system2(rscript, c("-e", shQuote(code), shQuote(infiles[i]),
                           shQuote(outfiles[i])),
                wait=FALSE, stdout=FALSE, stderr=FALSE)
    }
    deadline <- Sys.time() + timeout
    while (!all(file.exists(outfiles))) {
        if (Sys.time() > deadline)
            cli_abort(c("x"=paste("The new R processes running the steps",
                                  "did not end within {timeout} seconds.")))
        Sys.sleep(0.1)
    }
    t <- do.call(rbind, lapply(outfiles, function(f)
                                   as.numeric(readLines(f, warn=FALSE))))

    data.frame(startsecs=t[, 2] - launched, runsecs=t[, 3] - t[, 2])
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
## 'param' of a method other than GSVA, using 'ncpus' CPU cores and 'mem'
## bytes of memory, see gsvaJobWalltime(), as a data frame with the time of
## gsva() on all the columns of the input data, extrapolated from its time on
## two parts of them, see .chunk_secs(), whose widths are the larger of 1,000
## columns and 2% of them, in whole chunks of the input data on disk, if any,
## and the additional time of running it in a new R process, see
## .cold_secs(), on its first 100 columns
#' @importFrom BiocParallel SerialParam
.single_job_walltime <- function(param, mem, ncpus, cold, verbose) {
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
    gsvaargs <- function(first, last)
        list(param=.subset_param_cols(param, seq.int(first, last)),
             verbose=FALSE, BPPARAM=bp, maxmem=mem)
    tmpdir <- tempfile("gsvawalltime")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive=TRUE), add=TRUE)
    tiny <- .cold_secs(list(list(fun=gsva, args=gsvaargs(1, min(n, 100)))),
                       param, cold, tmpdir, verbose)
    if (verbose)
        cli_alert_info("Running gsva() on {sum(parts)} columns")
    run <- function(first, last)
        system.time(suppressMessages(do.call(gsva, gsvaargs(first,
                                                            last))))[["elapsed"]]

    structure(data.frame(step="gsva()", by="columns", chunk=n, chunks=1,
                         firstsecs=tiny$firstsecs,
                         chunksecs=.chunk_secs(parts, n, run)),
              startsecs=attr(tiny, "startsecs"))
}
