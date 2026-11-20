#!/usr/bin/env Rscript
## time of loading into main memory the first <ncells> cells of the processed
## data, as GSVA does with ondisk="no": as sparse, by GSVA's loader of sparse
## on-disk data, see .load_sparse_by_blocks(), in this process and with
## <workers> forked workers, and as dense, with as.matrix() in this process
## and, when the installed GSVA has it, by GSVA's loader of dense data with
## <workers> forked workers, see .load_dense(), which uses them only for data
## stored in chunks of several columns, from the CSC layout (TENxMatrix) and
## from the dense HDF5 dataset. each measurement is repeated <reps> times,
## alternating the serial and parallel loads, so that the first one shows the
## effect of reading the files cold. it calibrates 'maxparloaddensity', the
## largest fraction of nonzero values with which reading the dense dataset in
## parallel as sparse and converting it into a dense matrix in this process is
## faster than reading it in this process, see README.md. it needs both
## layouts written by prepare_calibration_data.R, with the option --dense.
## run it from the calibration directory on a machine, or in a job, with at
## least <workers> cores, e.g.:
##   sbatch -p short -c 6 --mem=48G -o results/input_load_time.out \
##          --wrap "Rscript measure_input_load_time.R 100000 6 3"
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 3)
    stop("usage: Rscript measure_input_load_time.R [ncells] [workers] [reps]")
ncells <- if (length(args) > 0) as.integer(args[1]) else 100000L
w <- if (length(args) > 1) as.integer(args[2]) else 6L
reps <- if (length(args) > 2) as.integer(args[3]) else 3L

suppressPackageStartupMessages({
    library(GSVA)
    library(BiocParallel)
    library(HDF5Array)
})
cat("GSVA", as.character(packageVersion("GSVA")), "from", find.package("GSVA"),
    "\nhost", Sys.info()[["nodename"]], "ncells", ncells, "workers", w,
    "\n")

csc <- TENxMatrix("data/TENxBrainDataMouseE18_processed_csc.sce/logcounts.h5",
                  group="matrix")[, seq_len(ncells)]
dense <- HDF5Array("data/TENxBrainDataMouseE18_processed_dense.sce/assays.h5",
                   "assay001", as.sparse=TRUE)[, seq_len(ncells)]
elapsed <- function(expr) {
    t <- system.time(res <- expr)[["elapsed"]]
    rm(res)
    invisible(gc())
    t
}

hasdense <- exists(".load_dense", envir=asNamespace("GSVA"))
res <- NULL
for (r in seq_len(reps)) {
    for (src in c("csc", "dense")) {
        X <- if (src == "csc") csc else dense
        res <- rbind(res, data.frame(rep=r, source=src, load="sparse",
                                     workers=1L,
                                     elapsed=elapsed(
                                         GSVA:::.load_sparse_by_blocks(X))))
        res <- rbind(res, data.frame(rep=r, source=src, load="sparse",
                                     workers=w,
                                     elapsed=elapsed(
                                         GSVA:::.load_sparse_by_blocks(X,
                                             BPPARAM=MulticoreParam(w)))))
        res <- rbind(res, data.frame(rep=r, source=src, load="dense",
                                     workers=1L,
                                     elapsed=elapsed(as.matrix(X))))
        if (hasdense)
            res <- rbind(res, data.frame(rep=r, source=src, load="dense",
                                         workers=w,
                                         elapsed=elapsed(
                                             GSVA:::.load_dense(X,
                                                 BPPARAM=MulticoreParam(w)))))
        print(tail(res, if (hasdense) 4 else 3), row.names=FALSE)
    }
}
cat("\nmedian elapsed time in seconds over the repetitions\n")
print(aggregate(elapsed ~ source + load + workers, data=res, FUN=median),
      row.names=FALSE)
