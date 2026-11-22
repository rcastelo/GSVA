## the wall time of the jobs running GSVA, see gsvaJobWalltime(): of each job
## of gsvaMap() running the three steps of the GSVA method, by running them on
## parts of the input data, and of a single job running gsva()
test_jobWalltime <- function() {
    message("Running unit tests for the wall time of the jobs running GSVA")

    set.seed(123)
    p <- 200; n <- 300
    x <- matrix(rnorm(p * n), nrow=p,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    gsets <- list(gs1=paste0("g", 1:20), gs2=paste0("g", 21:60))
    par <- gsvaParam(x, gsets, verbose=FALSE)

    ## each job of gsvaMap(), whose time is the largest one of the three
    ## steps, including starting its R process and what it does the first
    ## time it runs the step, measured in new R processes, times the margin
    wt <- gsvaJobWalltime(par, mem="4G", nodes=2, margin=2, verbose=FALSE)
    checkTrue(is(wt, "gsvaJobWalltime"))
    checkTrue(grepl("^[0-9]+:[0-5][0-9]:00$", wt))
    det <- attr(wt, "details")
    checkIdentical(det$step, c("normalization", "ranks", "scores"))
    checkTrue(all(det$chunksecs >= 0) && all(det$firstsecs >= 0))
    checkEqualsNumeric(det$perjob, ceiling(det$chunks / 2))
    checkEqualsNumeric(det$jobsecs, attr(wt, "startsecs") + det$firstsecs +
                                    det$perjob * det$chunksecs)
    checkTrue(attr(wt, "startsecs") > 0)
    secs <- sum(as.numeric(strsplit(as.character(wt), ":")[[1]]) *
                c(3600, 60, 1))
    checkTrue(secs >= 2 * max(det$jobsecs) && secs >= 60)
    out <- capture.output(details(wt))
    checkTrue(any(grepl("each of the 2 jobs of gsvaMap\\(\\)", out)))

    ## with array jobs, the time of each chunk in its own R process, and
    ## without new R processes, no time of doing things the first time
    wta <- gsvaJobWalltime(par, mem="4G", arrayjobs=TRUE, cold=FALSE,
                           verbose=FALSE)
    checkTrue(all(attr(wta, "details")$perjob == 1))
    checkTrue(all(attr(wta, "details")$firstsecs == 0))

    ## a single job running gsva(), with the three steps of the GSVA method
    ## one after another, or gsva() itself with other methods
    wtg <- gsvaJobWalltime(par, mem="4G", map=FALSE, cold=FALSE,
                           verbose=FALSE)
    checkIdentical(nrow(attr(wtg, "details")), 3L)
    wts <- gsvaJobWalltime(ssgseaParam(x, gsets, verbose=FALSE), mem="4G",
                           cold=FALSE, verbose=FALSE)
    checkIdentical(attr(wts, "details")$step, "gsva()")

    ## a memory smaller than the one that the jobs require gives a warning
    w <- tryCatch(gsvaJobWalltime(par, mem="100M", cold=FALSE,
                                  verbose=FALSE),
                  warning=function(w) w)
    checkTrue(is(w, "warning") && grepl("gsvaJobMemory", conditionMessage(w)))

    ## the hours, minutes and seconds of the wall time, rounded up to minutes
    checkIdentical(GSVA:::.walltime_text(1), "0:01:00")
    checkIdentical(GSVA:::.walltime_text(5401), "1:31:00")

    ## a fixed time and a time proportional to the width of the chunks, from
    ## the times of two parts of different widths
    checkEqualsNumeric(GSVA:::.chunk_secs(c(10, 20), 100,
                                          function(first, last)
                                              1 + 0.1 * (last - first + 1)),
                       11)

    checkException(gsvaJobWalltime(ssgseaParam(x, gsets, verbose=FALSE),
                                   mem="4G", map=TRUE), silent=TRUE)
    checkException(gsvaJobWalltime(par, mem="4G", nodes=0), silent=TRUE)
    checkException(gsvaJobWalltime(par, mem="4G", margin=0.5), silent=TRUE)
    ## numbers that as.integer() would truncate or turn into NA, and an
    ## infinite margin, which would give an infinite wall time
    checkException(gsvaJobWalltime(par, mem="4G", ncpus_per_task=1.5),
                   silent=TRUE)
    checkException(gsvaJobWalltime(par, mem="4G", nodes=Inf), silent=TRUE)
    checkException(gsvaJobWalltime(par, mem="4G", margin=Inf), silent=TRUE)
}
