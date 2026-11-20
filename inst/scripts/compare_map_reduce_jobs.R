## compares the runs of measure_map_reduce_jobs.R saved in
## <dir>/map_reduce_jobs_*.rds, where <dir> is the first argument, "results"
## by default: the time of each step and in total, the mean elapsed time of
## the SLURM jobs of each step, from the accounting in
## <dir>/map_reduce_jobs_*_sacct.txt, and whether the scores of every run are equal to
## those of the first one, in the sums of the scores of each column and in
## the scores of every 1000th column. repeated runs of the same configuration,
## whose names end in _r<n>, are also summarized by their mean
args <- commandArgs(trailingOnly=TRUE)
dir <- if (length(args) > 0) args[1] else "results"
files <- list.files(dir, pattern="^map_reduce_jobs_.*\\.rds$",
                    full.names=TRUE)
if (length(files) == 0)
    stop(sprintf("no results of measure_map_reduce_jobs.R in %s/", dir))
runs <- lapply(files, readRDS)
names(runs) <- vapply(runs, function(r) r$tag, character(1))

## mean elapsed time of the jobs of each step, which are told apart by their
## submission times, more than 20 seconds apart between steps; the last three
## groups of jobs are the three steps, after the job checking the library of
## GSVA, if any, see measure_map_reduce_jobs.R
.toS <- function(x) vapply(strsplit(x, ":"), function(v)
    sum(as.numeric(v) * c(3600, 60, 1)[seq_along(v) + 3 - length(v)]),
    numeric(1))
jobtimes <- function(tag) {
    f <- file.path(dir, sprintf("map_reduce_jobs_%s_sacct.txt", tag))
    if (!file.exists(f))
        return(c(norm=NA, rank=NA, score=NA))
    l <- readLines(f)[-(1:2)]
    w <- strsplit(trimws(l), " +")
    w <- w[vapply(w, function(x) length(x) >= 7 && x[3] == "COMPLETED" &&
                                 !grepl("^measure_map", x[2]),
                         logical(1))]
    if (length(w) == 0)
        return(c(norm=NA, rank=NA, score=NA))
    sub <- as.numeric(as.POSIXct(vapply(w, `[`, "", 4),
                                 format="%Y-%m-%dT%H:%M:%S"))
    el <- .toS(vapply(w, `[`, "", 7))
    o <- order(sub)
    grp <- cumsum(c(1, diff(sub[o]) > 20))
    m <- tapply(el[o], grp, mean)
    m <- tail(c(NA, NA, NA, m), 3)
    setNames(as.numeric(m), c("norm", "rank", "score"))
}

times <- t(vapply(runs, function(r) c(r$times, total=r$total),
                  numeric(length(runs[[1]]$times) + 1)))
colnames(times) <- sub("gsvaMap ", "map ", sub("gsvaReduce ", "red ",
                       sub("normalization", "norm", sub("ranking", "rank",
                       sub("scoring", "score", colnames(times))))))
jobs <- t(vapply(names(runs), jobtimes, numeric(3)))
colnames(jobs) <- paste("job", colnames(jobs))
tab <- cbind(times, jobs)
cat("time of each step of the main R session, and mean time of the jobs of",
    "each step, in seconds:\n")
print(round(tab, 1))

config <- sub("_r[0-9]+$", "", rownames(tab))
if (any(duplicated(config))) {
    cat("\nmean of the repeated runs of each configuration:\n")
    print(round(apply(tab, 2, function(x) tapply(x, config, mean)), 1))
}

ref <- runs[[1]]
cat(sprintf("\nscores compared with those of %s:\n", ref$tag))
for (r in runs[-1]) {
    same <- identical(r$dim, ref$dim) && identical(r$rownames, ref$rownames)
    eqsums <- same && isTRUE(all.equal(r$colsums, ref$colsums))
    eqsample <- same && isTRUE(all.equal(r$sample, ref$sample))
    maxdiff <- if (same) max(abs(r$sample - ref$sample), na.rm=TRUE) else NA
    cat(sprintf("%-36s dims %s, column sums %s, sampled scores %s (max diff %.3g)\n",
                r$tag, if (same) "equal" else "DIFFER",
                if (eqsums) "equal" else "DIFFER",
                if (eqsample) "equal" else "DIFFER", maxdiff))
}

cat("\nGSVA library of each run:\n")
for (r in runs)
    cat(sprintf("%-36s %s\n", r$tag,
                if (is.null(r$gsvapath)) "(not recorded)" else r$gsvapath))
cat("\npackage versions:\n")
print(t(vapply(runs, function(r) r$versions, character(length(ref$versions)))))
