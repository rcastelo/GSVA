#!/usr/bin/env Rscript
## runs the three steps of GSVA with rowNorm="clr" on the first <ncells> cells
## of the processed data, in the directory <se_dir>, either in the CSC layout
## written by prepare_calibration_data.R, with the files se.rds and
## logcounts.h5, or in a dense HDF5 dataset written with
## saveHDF5SummarizedExperiment(), with <ncores> forked workers and the memory
## budget <maxmem>, recording the time and the memory allocated by R in each
## step. run it through validate_maxmem_sweep.sh, see README.md:
##   Rscript validate_maxmem_sweep.R <se_dir> <gmt_file> <ncells> <ondisk:auto|yes|no> <maxmem:auto|XG> [ncores]

timestamp1 <- Sys.time()

## warnings are printed as they occur, instead of when the top-level call
## that gives them completes, so that they are also printed by runs killed
## during that call, e.g., for exceeding a hard memory limit
options(warn=1)

args <- commandArgs(trailingOnly=TRUE)

hdf5_se_dirname <- gmt_filename <- ncells <- ondisk <- maxmem <- NA
ncores <- 1

if (length(args) < 5) {
    stop("validate_maxmem_sweep.R <se_dir> <gmt_file> <ncells> <ondisk:auto|yes|no> <maxmem:auto|XG> [ncores]\n")
} else {
    hdf5_se_dirname <- args[1]
    if (!file.exists(hdf5_se_dirname))
        stop(sprintf("<se_dir> (%s) does not exist\n", hdf5_se_dirname))
    if (!file.info(hdf5_se_dirname)$isdir)
        stop(sprintf("<se_dir> (%s) is not a directory\n", hdf5_se_dirname))
    gmt_filename <- args[2]
    if (!file.exists(gmt_filename))
        stop(sprintf("<gmt_file> (%s) does not exist\n", gmt_filename))
    ncells <- as.integer(args[3])
    if (ncells < 1 || ncores < 1)
        stop(sprintf("<ncells> (%) and <ncores> (%d) should be positive integers\n",
		             ncells, ncores))
    ondisk <- args[4]
    if (ondisk != "auto" && ondisk != "yes" && ondisk != "no")
        stop(sprintf("<ondisk> (%s) should be either 'auto', 'yes' or 'no'\n", ondisk))
    maxmem <- args[5]
    if (length(args) > 5)
    ncores <- as.integer(args[6])
}

suppressPackageStartupMessages({
    library(BiocParallel)
    library(ExperimentHub)
    library(HDF5Array)
    library(SingleCellExperiment)
    library(GSEABase)
    library(GSVA)
})

## progress of the phases of the benchmark, written as they begin and end to
## the file given in the environment variable GSVA_PHASEFILE, when set, so that
## the phases completed before the process is killed, e.g., for exceeding a
## memory limit, are kept, see validate_maxmem_sweep.sh. each line has either the form
## 'begin,<phase>,<time>' or 'end,<phase>,<time>,<used>,<maxused>', where times
## are in seconds since the epoch and <used> and <maxused> are the memory
## allocated by R in KiB at the end of the phase and at most during it, or the
## form 'info,<name>,<value>'
phasefile <- Sys.getenv("GSVA_PHASEFILE", unset="")
progress <- function(...) {
    if (nzchar(phasefile))
        cat(paste(..., sep=","), "\n", sep="", file=phasefile, append=TRUE)
}
epoch <- function(t) sprintf("%.3f", as.numeric(t))

## memory in bytes allocated by R in this process, currently ('used') or at
## most since the last reset ('max used'), as tracked by its garbage collector,
## where each cons cell takes 56 bytes and each vector cell 8 bytes; it does
## not include the memory of the parallel workers
rmem <- function(what=c("used", "max used"), reset=FALSE) {
    g <- gc(reset=reset)
    sum(g[, match(match.arg(what), colnames(g))] * c(56, 8))
}

progress("begin", "start", epoch(timestamp1))

## the CSC layout, written by prepare_calibration_data.R, or a dense HDF5
## dataset written by saveHDF5SummarizedExperiment()
if (file.exists(file.path(hdf5_se_dirname, "logcounts.h5"))) {
    sce <- readRDS(file.path(hdf5_se_dirname, "se.rds"))
    assay(sce, "logcounts", withDimnames=FALSE) <-
        TENxMatrix(file.path(hdf5_se_dirname, "logcounts.h5"), group="matrix")
} else
    sce <- loadHDF5SummarizedExperiment(hdf5_se_dirname)
gsets <- readGMT(gmt_filename)

sce <- sce[, 1:ncells]
gsvapar <- gsvaParam(sce, gsets, rowNorm="clr", minSize=5, maxSize=300, verbose=FALSE, ondisk=ondisk)

## memory budget in KiB, as resolved by GSVA, e.g., when maxmem="auto"
maxmemkib <- GSVA:::.check_maxmem(gsvapar, maxmem=maxmem, verbose=FALSE) / 1024
progress("info", "maxmem", sprintf("%.0f", maxmemkib))

## R memory at the start of the steps, from which the peak of each one is
## tracked, resetting it with a garbage collection between steps, which is
## not included in the elapsed time of any of them
rmemstart <- rmem("used", reset=TRUE)
timestamp2 <- Sys.time()
progress("end", "start", epoch(timestamp2), sprintf("%.0f", rmemstart / 1024),
         NA)

## runs a step of GSVA, giving its result, elapsed time and peak of R memory
step <- function(name, expr) {
    t0 <- Sys.time()
    progress("begin", name, epoch(t0))
    res <- expr
    t1 <- Sys.time()
    maxused <- rmem("max used")
    used <- rmem("used", reset=TRUE)
    progress("end", name, epoch(t1), sprintf("%.0f", used / 1024),
             sprintf("%.0f", maxused / 1024))
    list(res=res, secs=as.numeric(difftime(t1, t0, units="secs")),
         maxused=maxused)
}

norm <- step("norm", gsvaRowNorm(gsvapar, verbose=FALSE, maxmem=maxmem,
                                 BPPARAM=MulticoreParam(workers=ncores)))
ranks <- step("ranks", gsvaColRanks(norm$res, verbose=FALSE, maxmem=maxmem,
                                    BPPARAM=MulticoreParam(workers=ncores)))
scores <- step("scores", gsvaColScores(ranks$res, verbose=FALSE, maxmem=maxmem,
                                       BPPARAM=MulticoreParam(workers=ncores)))

ngsets <- length(geneSets(scores$res))
progress("info", "ngsets", ngsets)
rmemhigh <- max(norm$maxused, ranks$maxused, scores$maxused)

cat(sprintf("%d,%d,%d,%.0f,%.3f,%.3f,%.3f,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,",
            ncells, ngsets, ncores, maxmemkib,
            difftime(timestamp2, timestamp1, units="secs"), norm$secs,
            ranks$secs, scores$secs, rmemstart / 1024, norm$maxused / 1024,
            ranks$maxused / 1024, scores$maxused / 1024, rmemhigh / 1024))
