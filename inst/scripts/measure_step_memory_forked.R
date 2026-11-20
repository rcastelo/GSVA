#!/usr/bin/env Rscript
## runs one step of GSVA on the first <ncells> cells of the CSC data held in
## main memory (ondisk="no"), with <workers> forked workers and the memory
## budget <maxmem> (default "auto"), e.g., 16G, printing GSVA's memory
## estimate for the step and the times at which the step starts and ends,
## from which measure_step_memory_forked.sh takes the peak memory of the step;
## with step "prep", it saves the row normalization and ranks that the steps
## "ranks" and "scores" take as input. run it through
## measure_step_memory_forked.sh, from the calibration directory, see
## README.md:
##   Rscript measure_step_memory_forked.R <prep|norm_ecdf|norm_clr|ranks|scores> <ncells> <workers> [maxmem]
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) %in% 3:4)
    stop("usage: Rscript measure_step_memory_forked.R <prep|norm_ecdf|norm_clr|ranks|scores> <ncells> <workers> [maxmem]")
step <- match.arg(args[1], c("prep", "norm_ecdf", "norm_clr", "ranks", "scores"))
ncells <- as.integer(args[2])
w <- as.integer(args[3])
maxmem <- if (length(args) == 4) args[4] else "auto"

suppressPackageStartupMessages({
    library(GSVA)
    library(BiocParallel)
    library(HDF5Array)
    library(SingleCellExperiment)
})
cat("GSVA", as.character(packageVersion("GSVA")), "from", find.package("GSVA"),
    "\n")

prepfile <- file.path("data", sprintf("forked_ranks_%d.rds", ncells))
gib <- function(x) sprintf("%.3f", x / 1024^3)
param <- function(rowNorm) {
    dir <- "data/TENxBrainDataMouseE18_processed_csc.sce"
    sce <- readRDS(file.path(dir, "se.rds"))
    assay(sce, "logcounts", withDimnames=FALSE) <-
        TENxMatrix(file.path(dir, "logcounts.h5"), group="matrix")
    sce <- sce[, seq_len(ncells)]
    gsets <- readGMT("data/m8.all.v2025.1.Mm.symbols.gmt")
    gsvaParam(sce, gsets, assay="logcounts", rowNorm=rowNorm, minSize=5,
              maxSize=300, ondisk="no", verbose=FALSE)
}
bp <- if (w == 1) SerialParam() else MulticoreParam(workers=w, progressbar=FALSE)

if (step == "prep") {
    dir.create("data", showWarnings=FALSE)
    rk <- gsvaColRanks(gsvaRowNorm(param("clr"), verbose=FALSE, BPPARAM=bp),
                       verbose=FALSE, BPPARAM=bp)
    ## the object keeps the original assay on disk, which base::saveRDS()
    ## stores as a reference to its file
    base::saveRDS(rk, prepfile)
    quit(save="no")
}

input <- switch(step,
                norm_ecdf=param("ecdf"),
                norm_clr=param("clr"),
                readRDS(prepfile))
## the ranks step adds the assay of the ranks
if (step == "ranks")
    assay(input, "gsvaranks") <- NULL
invisible(gc())
cat("RMEM", gib(sum(gc()[, 1] * c(56, 8))), "\n")

## GSVA's estimate of the memory of the step, in GiB, and the memory it
## charges each R process, in each of its checks
trace(GSVA:::.check_mem_need, where=asNamespace("GSVA"), print=FALSE,
      tracer=quote(cat("ESTIMATE", gib(need), gib(permem), nworkers, "\n")))

cat("PHASE start", format(as.numeric(Sys.time()), nsmall=3), "\n")
t <- system.time(res <- switch(step,
    norm_ecdf=,
    norm_clr=gsvaRowNorm(input, verbose=FALSE, BPPARAM=bp, maxmem=maxmem),
    ranks=gsvaColRanks(input, verbose=FALSE, BPPARAM=bp, maxmem=maxmem),
    scores=gsvaColScores(input, verbose=FALSE, BPPARAM=bp, device="cpu",
                         maxmem=maxmem)))
cat("PHASE end", format(as.numeric(Sys.time()), nsmall=3), "\n")
cat("ELAPSED", t[["elapsed"]], "\n")
