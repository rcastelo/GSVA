## the memory of the jobs running GSVA, see gsvaJobMemory(): of each job of
## gsvaMap() running the three steps of the GSVA method, and of a single job
## running gsva(), with which the memory checks of gsvaMap() and gsva() pass
test_jobMemory <- function() {
    message("Running unit tests for the memory of the jobs running GSVA")

    suppressPackageStartupMessages({
        library(HDF5Array)
        library(SummarizedExperiment)
    })

    set.seed(123)
    p <- 200; n <- 400
    x <- matrix(rnorm(p * n), nrow=p,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    gsets <- list(gs1=paste0("g", 1:20), gs2=paste0("g", 21:60))
    gb <- 1024^3

    ## each job of gsvaMap(), whose memory fits the largest chunk of each
    ## step, together with its R process and its parallel workers
    par <- gsvaParam(x, gsets, verbose=FALSE)
    mem <- gsvaJobMemory(par, ncpus_per_task=2)
    checkTrue(is(mem, "gsvaJobMemory"))
    checkTrue(grepl("^[0-9]+[MG]$", mem))
    det <- attr(mem, "details")
    checkIdentical(det$step, c("normalization", "ranks", "scores"))
    checkTrue(all(det$minimum <= det$recommended))
    checkEqualsNumeric(GSVA:::.memtext2bytes(as.character(mem)),
                       ceiling(max(det$recommended) / gb) * gb)
    memb <- GSVA:::.memtext2bytes(as.character(mem))
    need <- GSVA:::.map_job_need("gsvaRowNorm", par, NA, 1L, det$chunk[1], 2L,
                                 memb)
    checkTrue(need <= memb)
    checkTrue(attr(need, "baseline") <= memb)
    out <- capture.output(details(mem))
    checkTrue(any(grepl("each job of gsvaMap\\(\\) with 2 CPU cores", out)))

    ## with the input data stored on disk in chunks of several columns, the
    ## recommended memory of the ranks, at least the smallest one, gives
    ## blocks of whole chunks to the jobs, which their workers share
    h <- writeHDF5Array(x, chunkdim=c(50L, 100L))
    se <- SummarizedExperiment(assays=list(counts=h))
    parh <- gsvaParam(se, gsets, assay="counts", verbose=FALSE)
    deth <- attr(gsvaJobMemory(parh, ncpus_per_task=2), "details")
    checkTrue(all(deth$minimum <= deth$recommended))
    checkTrue(deth$chunk[deth$step == "ranks"] >= 100)

    ## a single job running gsva(), with any method, whose estimate is the one
    ## of gsva() itself, while its calculations are not run
    memss <- gsvaJobMemory(ssgseaParam(x, gsets, verbose=FALSE))
    detss <- attr(memss, "details")
    checkIdentical(nrow(detss), 1L)
    checkTrue(detss$minimum > 0 && detss$minimum <= detss$recommended)
    checkTrue(!GSVA:::gsva_global$memdryrun)
    checkTrue(GSVA:::gsva_global$check_cores)
    checkTrue(GSVA:::gsva_global$check_memory)
    memg <- gsvaJobMemory(par, map=FALSE)
    checkIdentical(nrow(attr(memg, "details")), 1L)
    ## with the input data on disk, it can also be processed from disk
    detgh <- attr(gsvaJobMemory(parh, map=FALSE), "details")
    checkTrue(!is.na(detgh$ondisk))
    checkEqualsNumeric(detgh$minimum, min(detgh$inmemory, detgh$ondisk))

    ## gsvaMap() runs only the GSVA method
    checkException(gsvaJobMemory(ssgseaParam(x, gsets, verbose=FALSE),
                                 map=TRUE), silent=TRUE)
    checkException(gsvaJobMemory(par, ncpus_per_task=0), silent=TRUE)
}
