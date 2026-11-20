#!/usr/bin/env Rscript
## plots the time and the peak memory of each step and in total of the runs of
## validate_maxmem_sweep.sh against their number of cells, see README.md:
##   Rscript plot_memory_time.R <input_csv_file> <output_png_file> [ncores]

args <- commandArgs(trailingOnly=TRUE)

inputfname <- outputfname <- ncores <- NA

if (length(args) < 2) {
    stop("plot_memory_time.R <input_csv_file> <output_png_file> [ncores]")
} else {
    inputfname <- args[1]
    if (!file.exists(inputfname))
        stop(paste("Input file", inputfname, "does not exist"))
    outputfname <- args[2]
    if (length(args) > 2)
        ncores <- as.integer(args[3])
}

dat <- read.csv(inputfname)
## only the runs that completed, leaving out those that did not, e.g., because
## they were killed for exceeding a hard memory limit, whose times are NA
done <- !is.na(dat$ENORM) & !is.na(dat$ERANKS) & !is.na(dat$ESCORES)
if ("STATUS" %in% colnames(dat))
    done <- done & dat$STATUS == 0
if (!all(done))
    message(sprintf("Leaving out %d run(s) that did not complete: NCELLS %s",
                    sum(!done), paste(dat$NCELLS[!done], collapse=", ")))
if (!any(done))
    stop(paste("No run in", inputfname, "completed"))
dat <- dat[done, ]

## the number of cores from the runs, unless given
if (is.na(ncores))
    ncores <- paste(unique(dat$NCORES), collapse=",")

png(file=outputfname, height=5, width=10, units="in", res=300)
par(mfrow=c(1, 2), mar=c(4, 5, 2, 1))
tot <- dat$ENORM + dat$ERANKS + dat$ESCORES
plot(dat$NCELLS, tot/60, type="b", lwd=2, xlab=paste("# Cells (x 10^5)"),
     ylab=sprintf("Time (min.) using %s cores", ncores), xaxt="n", las=1,
     cex.axis=1.2, cex.lab=1.5, ylim=c(0, ceiling(max(tot/60))),
     panel.first=grid(), main="")
axis(1, at=seq(1e5, 1e6, by=1e5), labels=1:10, cex.axis=1.2)
lines(dat$NCELLS, dat$ENORM/60, type="b", lwd=2, lty=2)
lines(dat$NCELLS, dat$ERANKS/60, type="b", lwd=2, lty=3)
lines(dat$NCELLS, dat$ESCORES/60, type="b", lwd=2, lty=4)
legend("topleft", c("Total", "Norm.", "Ranks", "Scores"), lty=c(1, 2, 3, 4),
       lwd=2, inset=0.05, bg="white")

## peak memory of the R process and its parallel workers, overall and in each
## phase, from the samples of the memory profiler, excluding the page cache of
## the files that they read, i.e., the anonymous memory of the cgroup on Linux
## and the physical footprint of the processes on macOS, so that it has the same
## meaning on both, and is the memory to compare with the budget 'maxmem'
mem <- dat$GRP_ANON_HIGH / 1024^2
plot(dat$NCELLS, mem, type="b", lwd=2, xlab="# Cells (x 10^5)",
     ylab=sprintf("Peak memory (GiB) using %s cores", ncores), xaxt="n",
     las=1, cex.axis=1.2, cex.lab=1.5, ylim=c(0, ceiling(max(mem))),
     panel.first=grid(), main="")
axis(1, at=seq(1e5, 1e6, by=1e5), labels=1:10, cex.axis=1.2)
lines(dat$NCELLS, dat$GRP_ANON_NORM / 1024^2, type="b", lwd=2, lty=2)
lines(dat$NCELLS, dat$GRP_ANON_RANKS / 1024^2, type="b", lwd=2, lty=3)
lines(dat$NCELLS, dat$GRP_ANON_SCORES / 1024^2, type="b", lwd=2, lty=4)
legend("topleft", c("Overall", "Norm.", "Ranks", "Scores"), lty=c(1, 2, 3, 4),
       lwd=2, inset=0.05, bg="white")
dev.off()
