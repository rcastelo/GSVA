#!/usr/bin/env Rscript
## plots the peak memory and the time of the runs of validate_maxmem_sweep.sh
## against their memory budget, for each number of cells, marking the runs
## killed for exceeding it, see README.md:
##   Rscript plot_maxmem_sweep.R <input_csv_file> <output_png_file>

args <- commandArgs(trailingOnly=TRUE)

inputfname <- outputfname <- NA

if (length(args) < 2) {
    stop("plot_maxmem_sweep.R <input_csv_file> <output_png_file>")
} else {
    inputfname <- args[1]
    if (!file.exists(inputfname))
        stop(paste("Input file", inputfname, "does not exist"))
    outputfname <- args[2]
}

dat <- read.csv(inputfname)
dat <- dat[!is.na(dat$MAXMEM), ]
dat <- dat[order(dat$NCELLS, dat$MAXMEM), ]
ncores <- unique(dat$NCORES)
cells <- sort(unique(dat$NCELLS))
killed <- dat$STATUS != 0

## memory in GiB; the peak memory of the R process and its parallel workers,
## from the samples of the memory profiler, excluding the page cache of the
## files that they read, i.e., the anonymous memory of the cgroup on Linux and
## the physical footprint of the processes on macOS, so that it has the same
## meaning on both, and is the memory to compare with the budget 'maxmem'
dat$MAXMEM <- dat$MAXMEM / 1024^2
dat$GRP_ANON_HIGH <- dat$GRP_ANON_HIGH / 1024^2
## memory allocated by R during the steps above the one at their start, only
## with one core, because it does not include the memory of the parallel workers
rmem <- all(dat$NCORES == 1)
if (rmem)
    dat$R_MEM_STEPS <- (dat$R_MEM_HIGH - dat$R_MEM_START) / 1024^2
tot <- (dat$ENORM + dat$ERANKS + dat$ESCORES) / 60

png(file=outputfname, height=5, width=10, units="in", res=300)
par(mfrow=c(1, 2), mar=c(4, 5, 2, 1))
lim <- c(0, max(c(dat$MAXMEM, dat$GRP_ANON_HIGH), na.rm=TRUE))
## headroom above the data for the legend
plot(dat$MAXMEM, dat$GRP_ANON_HIGH, type="n", xlim=lim, ylim=lim * c(1, 1.6), las=1,
     xlab="maxmem (GiB)", ylab=sprintf("Peak memory (GiB) using %s cores",
                                     paste(ncores, collapse=",")),
     cex.axis=1.2, cex.lab=1.5, panel.first=grid(), main="")
abline(0, 1, col="grey40")
abline(0, 0.7, col="grey40", lty=3)
for (i in seq_along(cells)) {
    mask <- dat$NCELLS == cells[i]
    lines(dat$MAXMEM[mask], dat$GRP_ANON_HIGH[mask], type="b", lwd=2, lty=i,
          pch=ifelse(killed[mask], 4, 1), cex=ifelse(killed[mask], 1.5, 1))
    if (rmem)
        lines(dat$MAXMEM[mask], dat$R_MEM_STEPS[mask], type="b", lwd=2, lty=i,
              pch=2, col="grey50")
}
leg <- c(sprintf("%d cells", cells), "maxmem", "70% maxmem", "killed")
legpch <- c(rep(1, length(cells)), NA, NA, 4)
leglty <- c(seq_along(cells), 1, 3, NA)
leglwd <- c(rep(2, length(cells)), 1, 1, 2)
legcol <- c(rep("black", length(cells)), "grey40", "grey40", "black")
if (rmem) {
    leg <- c("all processes", "R in steps", leg)
    legpch <- c(1, 2, legpch)
    leglty <- c(NA, NA, leglty)
    leglwd <- c(2, 2, leglwd)
    legcol <- c("black", "grey50", legcol)
}
legend("topleft", leg, lty=leglty, pch=legpch, lwd=leglwd, col=legcol,
       inset=0.02, bg="white", cex=0.8, ncol=2)

## the time of the runs that completed, if any
ok <- !killed
plot(dat$MAXMEM[ok], tot[ok], type="n", xlim=lim,
     ylim=c(0, if (any(ok)) max(tot[ok]) else 1),
     xlab="maxmem (GiB)", ylab="Time (min.)", las=1, cex.axis=1.2, cex.lab=1.5,
     panel.first=grid(), main="")
for (i in seq_along(cells)) {
    mask <- ok & dat$NCELLS == cells[i]
    lines(dat$MAXMEM[mask], tot[mask], type="b", lwd=2, lty=i)
}
legend("topright", sprintf("%d cells", cells), lty=seq_along(cells), lwd=2,
       inset=0.05, bg="white")
dev.off()
