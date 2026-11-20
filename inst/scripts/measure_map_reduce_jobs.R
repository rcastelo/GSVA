## runs on SLURM the three-step GSVA pipeline with rowNorm="clr" through
## gsvaMap() and gsvaReduce(), whose CLR row normalization is calculated by
## chunks of columns, on the dense or on the CSC (TENxMatrix) layout of the
## normalized values written by prepare_calibration_data.R, with input read
## from disk by the jobs, whose memory, recorded by the accounting of SLURM,
## see measure_map_reduce_jobs.sh, calibrates the step factors of sparse data
## read from disk ('sparsedisk'), see README.md. usage:
##   Rscript measure_map_reduce_jobs.R <layout> <ncells> <nodes> <ncpus> <mem>
## e.g., 'Rscript measure_map_reduce_jobs.R csc 200000 10 2 6G', run from the
## calibration directory by a job with a time limit long enough for the three
## steps, see measure_map_reduce_jobs.sh. the jobs of gsvaMap() are sent to
## the partition in the environment variable GSVA_PARTITION, "short" by
## default, with the time limit in GSVA_WALLTIME, "01:00:00" by default. the
## timestamps of each step go to stderr, and a summary of the scores goes to
## results/map_reduce_jobs_<layout>_<ncells>_<nodes>x<ncpus>x<mem>.rds, which
## compare_map_reduce_jobs.R compares between runs. the optional environment
## variables GSVA_MAPRUN_LABEL, e.g., "base" or "new", and GSVA_MAPRUN_REP, the
## number of a repeated run, are added to the name of the run, and
## GSVA_MAPRUN_OUT replaces the directory "results", to compare two builds of
## GSVA, each in its library of R_LIBS, also taken by the jobs sent to SLURM
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5)
    stop("usage: Rscript measure_map_reduce_jobs.R <dense|csc> <ncells> <nodes> <ncpus> <mem>")
layout <- match.arg(args[1], c("dense", "csc"))
ncells <- as.integer(args[2])
nodes <- as.integer(args[3])
ncpus <- as.integer(args[4])
mem <- args[5]
tag <- sprintf("%s_%d_%dx%dx%s", layout, ncells, nodes, ncpus, mem)
label <- Sys.getenv("GSVA_MAPRUN_LABEL")
if (nzchar(label))
    tag <- paste(label, tag, sep="_")
if (nzchar(Sys.getenv("GSVA_MAPRUN_REP")))
    tag <- paste0(tag, "_r", Sys.getenv("GSVA_MAPRUN_REP"))
resdir <- Sys.getenv("GSVA_MAPRUN_OUT", "results")
partition <- Sys.getenv("GSVA_PARTITION", "short")
walltime <- Sys.getenv("GSVA_WALLTIME", "01:00:00")

suppressPackageStartupMessages({
    library(HDF5Array); library(SingleCellExperiment)
    library(GSEABase); library(GSVA)
})
cat(sprintf("run: %s\n", tag))
pkgs <- c("GSVA", "HDF5Array", "h5mread", "SparseArray", "DelayedArray",
          "BiocParallel", "batchtools")
versions <- vapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)),
                                              error=function(e) NA_character_),
                   character(1))
print(versions)
gsvapath <- find.package("GSVA")
cat(sprintf("GSVA loaded from: %s\n", gsvapath))
cat(sprintf("R_LIBS: %s\n", Sys.getenv("R_LIBS")))

## print to stderr a timestamp when each of these functions of the main R
## session starts, and when it ends, followed by the cumulative time spent by
## the garbage collector, in seconds, after 'gc'
invisible(gc.time(TRUE))
.tstamp <- function(what) cat(format(Sys.time(), "%H:%M:%OS2"), what,
                              sprintf("gc %.2f", gc.time()[3]), "\n",
                              file=stderr())
## functions missing in the build of GSVA, e.g., .map_clr_stats() in builds
## before the CLR row normalization by chunks of columns, are skipped
.timefuns <- function(funs, ns) for (f in funs)
    if (exists(f, envir=asNamespace(ns), inherits=FALSE))
        suppressMessages(trace(f, where=asNamespace(ns), print=FALSE,
                               tracer=bquote(.tstamp(paste("start", .(f)))),
                               exit=bquote(.tstamp(paste("end  ", .(f))))))
.timefuns(c(".map_fingerprint", ".estimate_nzcount", ".check_map_mem",
            ".write_map_manifest", ".map_chunks", ".map_clr_stats"), "GSVA")
.timefuns(c("submitJobs", "waitForJobs", "reduceResultsList"), "batchtools")

if (layout == "dense") {
    sce <- loadHDF5SummarizedExperiment("data/TENxBrainDataMouseE18_processed_dense.sce")
} else {
    dir <- "data/TENxBrainDataMouseE18_processed_csc.sce"
    sce <- readRDS(file.path(dir, "se.rds"))
    assay(sce, "logcounts", withDimnames=FALSE) <-
        TENxMatrix(file.path(dir, "logcounts.h5"), group="matrix")
}
sce <- sce[, seq_len(min(ncells, ncol(sce)))]
print(seed(assay(sce, "logcounts")))
gsets <- readGMT("data/m8.all.v2025.1.Mm.symbols.gmt")
par <- gsvaParam(sce, gsets, assay="logcounts", rowNorm="clr",
                 minSize=5, maxSize=300, ondisk="yes")

## a separate directory for the intermediate results of each run, deleted at
## the end, so that runs do not find the results of previous ones
outdir <- file.path("GSVAOUTPUT", tag)
if (dir.exists(outdir))
    stop(sprintf("directory %s already exists, delete it first", outdir))
btpar <- gsvaBatchtoolsSlurmParam(dir=outdir, partition=partition,
                                  walltime=walltime, nodes=nodes,
                                  ncpus_per_task=ncpus, mem=mem)

## the jobs sent to SLURM must load GSVA from the same library as this R
## session, which they do when they inherit R_LIBS, as sbatch exports the
## environment by default; this is checked with a single small job, before
## the timed steps, in its own directory
chkpar <- gsvaBatchtoolsSlurmParam(dir=file.path(outdir, "libcheck"),
                                   partition=partition, walltime="00:10:00",
                                   nodes=1, ncpus_per_task=1, mem="2G")
jobpath <- BiocParallel::bplapply(1, function(i) find.package("GSVA"),
                                  BPPARAM=chkpar)[[1]]
cat(sprintf("GSVA loaded by the jobs from: %s\n", jobpath))
if (!identical(normalizePath(jobpath), normalizePath(gsvapath)))
    stop(sprintf("the jobs load GSVA from %s instead of %s", jobpath,
                 gsvapath))
unlink(file.path(outdir, "libcheck"), recursive=TRUE)

times <- numeric(0)
.step <- function(name, expr) {
    .tstamp(paste("start", name))
    t <- system.time(res <- expr)[["elapsed"]]
    .tstamp(paste("end  ", name))
    times[[name]] <<- t
    res
}
t0 <- Sys.time()
m   <- .step("gsvaMap normalization",    gsvaMap(gsvaRowNorm, par, BTPARAM=btpar))
nrm <- .step("gsvaReduce normalization", gsvaReduce(m))
m   <- .step("gsvaMap ranking",          gsvaMap(gsvaColRanks, nrm, BTPARAM=btpar))
cat(sprintf("ranking chunks: %d\n", length(m)))
rnk <- .step("gsvaReduce ranking",       gsvaReduce(m))
m   <- .step("gsvaMap scoring",          gsvaMap(gsvaColScores, rnk, BTPARAM=btpar))
es  <- .step("gsvaReduce scoring",       gsvaReduce(m))
total <- as.numeric(difftime(Sys.time(), t0, units="secs"))
print(round(times, 1))
cat(sprintf("total: %.1f s\n", total))

## the scores, which can be stored in the directory of the intermediate
## results, are summarized before deleting it: the sums of the scores of each
## column, and the scores of every 1000th column
E <- assay(es)
summ <- list(tag=tag, times=times, total=total, dim=dim(E),
             rownames=rownames(E), colsums=as.numeric(colSums(E)),
             sample=as.matrix(E[, seq(1, ncol(E), by=1000), drop=FALSE]),
             versions=versions, gsvapath=gsvapath)
dir.create(resdir, showWarnings=FALSE, recursive=TRUE)
saveRDS(summ, file.path(resdir, sprintf("map_reduce_jobs_%s.rds", tag)))
unlink(outdir, recursive=TRUE)
