test_ondisk <- function() {

    message("Running unit tests for ondisk input.")

    suppressPackageStartupMessages({
        library(Matrix)
        library(HDF5Array)
        library(SparseArray)
    })

    p <- 50 ## number of genes
    n <- 100 ## number of samples
    m <- 5 ## number of gene sets

    ## create random gene sets of size between 3 and 10
    set.seed(123)
    gssizes <- sample(3:10, size=m, replace=TRUE)
    gsets <- lapply(gssizes, function(x) sample(paste0("g", 1:p), size=x, replace=FALSE))
    names(gsets) <- paste0("gs", 1:m)

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    s <- ceiling(0.15 * p * n)
    sam <- sample(1:(p * n), size=s, replace=FALSE)
    x <- numeric(p * n)
    x[sam] <- runif(s)
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    M <- Matrix(y, sparse=TRUE)
    H5 <- as(M, "HDF5Matrix")

    ## estimate GSVA enrichment scores with and without HDF5 input and check that they are identical
    es_noh5 <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE)
    es_h5 <- gsva(gsvaParam(H5, gsets, verbose=FALSE), verbose=FALSE)
    checkIdentical(es_noh5, es_h5)

    ## estimate GSVA enrichment scores with HDF5 input and output and check that they are identical
    es_h5ondisk <- gsva(gsvaParam(H5, gsets, ondisk="yes", verbose=FALSE), verbose=FALSE)
    es_h5ondiskmat <- as.matrix(es_h5ondisk)
    checkEqualsNumeric(es_noh5, es_h5ondiskmat)

    ## estimate ssGSEA enrichment scores with and without HDF5 input and check that they are identical
    es_noh5 <- gsva(ssgseaParam(M, gsets, verbose=FALSE), verbose=FALSE)
    es_h5 <- gsva(ssgseaParam(H5, gsets, verbose=FALSE), verbose=FALSE)
    checkIdentical(es_noh5, es_h5)

    ## estimate ssGSEA enrichment scores with HDF5 input and output and check that they are identical
    es_h5ondisk <- gsva(ssgseaParam(H5, gsets, ondisk="yes", verbose=FALSE),
                        verbose=FALSE)
    es_h5ondiskmat <- as.matrix(es_h5ondisk)
    checkEqualsNumeric(es_noh5, es_h5ondiskmat)

    ## estimate Z-scores enrichment scores with and without HDF5 input and check that they are identical
    es_noh5 <- gsva(zscoreParam(M, gsets, verbose=FALSE), verbose=TRUE)
    es_h5 <- gsva(zscoreParam(H5, gsets, verbose=FALSE), verbose=FALSE)
    ## not identical also due to the rowSds() vs sd() differences
    checkEqualsNumeric(es_noh5, es_h5)

    ## estimate Z-scores enrichment scores with HDF5 input and output and check that they are identical
    es_h5ondisk <- gsva(zscoreParam(H5, gsets, ondisk="yes", verbose=FALSE),
                        verbose=FALSE)
    es_h5ondiskmat <- as.matrix(es_h5ondisk)
    ## not identical also due to the rowSds() vs sd() differences
    checkEqualsNumeric(es_noh5, es_h5ondiskmat)

    ## estimate average enrichment scores with and without HDF5 input and check that they are identical
    es_noh5 <- gsva(avgParam(M, gsets, verbose=FALSE), verbose=TRUE)
    es_h5 <- gsva(avgParam(H5, gsets, verbose=FALSE), verbose=FALSE)
    checkIdentical(es_noh5, es_h5)

    ## estimate average enrichment scores with HDF5 input and output and check that they are identical
    es_h5ondisk <- gsva(avgParam(H5, gsets, ondisk="yes", verbose=FALSE),
                        verbose=FALSE)
    es_h5ondiskmat <- as.matrix(es_h5ondisk)
    attr(es_noh5, "gsvaParam") <- attr(es_noh5, "assay") <- attr(es_noh5, "geneSets") <- NULL
    attr(es_noh5, "gsvaVersion") <- NULL
    checkIdentical(es_noh5, es_h5ondiskmat)

    ## test the block processing of a small toy HDF5 input and output by
    ## setting a small block size and maximum available memory
    oldautoblocksize <- getAutoBlockSize()
    setAutoBlockSize(1024)
    es_noh5 <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=FALSE)
    es_chunks <- gsva(gsvaParam(M, gsets, verbose=FALSE), verbose=TRUE, maxmem="25K")
    checkIdentical(es_noh5, es_chunks)
    setAutoBlockSize(oldautoblocksize)
}

test_job_memory_limit <- function() {

    message("Running unit tests for the memory limit of the job")

    noslurm <- c(NA, NA, NA)
    root <- tempfile("cgroups")
    on.exit(unlink(root, recursive=TRUE), add=TRUE)
    writefile <- function(path, lines) {
        dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
        writeLines(lines, path)
    }
    procfile <- file.path(root, "proc_self_cgroup")
    gb <- 1024^3

    ## no cgroup information and no SLURM job: no limit
    checkIdentical(GSVA:::.job_memory_limit(procfile=file.path(root, "none"),
                                            cgroupdir=root, env=noslurm), Inf)

    ## cgroup version 2, where SLURM sets the limit of a job in an ancestor of
    ## the cgroup of the process, whose own limit is 'max', i.e., no limit
    cgdir <- file.path(root, "v2")
    writefile(procfile, "0::/system.slice/slurmstepd.scope/job_1/step_0/task_0")
    writefile(file.path(cgdir, "system.slice", "slurmstepd.scope", "job_1",
                        "memory.max"), format(10 * gb, scientific=FALSE))
    writefile(file.path(cgdir, "system.slice", "slurmstepd.scope", "job_1",
                        "step_0", "task_0", "memory.max"), "max")
    checkEqualsNumeric(GSVA:::.job_memory_limit(procfile=procfile,
                                                cgroupdir=cgdir,
                                                env=noslurm), 10 * gb)

    ## cgroup version 1, with the memory controller in its own hierarchy
    cgdir <- file.path(root, "v1")
    writefile(procfile, c("12:cpu,cpuacct:/slurm/uid_1/job_2",
                          "11:memory:/slurm/uid_1/job_2/step_batch"))
    writefile(file.path(cgdir, "memory", "slurm", "uid_1", "job_2",
                        "step_batch", "memory.limit_in_bytes"),
              format(4 * gb, scientific=FALSE))
    checkEqualsNumeric(GSVA:::.job_memory_limit(procfile=procfile,
                                                cgroupdir=cgdir,
                                                env=noslurm), 4 * gb)

    ## SLURM environment variables, in megabytes, where the smallest of all
    ## limits is taken, and 0 means no limit
    nocg <- file.path(root, "none")
    checkEqualsNumeric(GSVA:::.job_memory_limit(procfile=nocg, cgroupdir=root,
                                                env=c("2048", NA, NA)),
                       2 * gb)
    checkEqualsNumeric(GSVA:::.job_memory_limit(procfile=nocg, cgroupdir=root,
                                                env=c(NA, "1024", "4")),
                       4 * gb)
    checkIdentical(GSVA:::.job_memory_limit(procfile=nocg, cgroupdir=root,
                                            env=c("0", NA, NA)), Inf)
    checkEqualsNumeric(GSVA:::.job_memory_limit(procfile=procfile,
                                                cgroupdir=cgdir,
                                                env=c("2048", NA, NA)),
                       2 * gb)
}

test_load_sparse_by_blocks <- function() {

    message("Running unit tests for loading sparse on-disk data by blocks")

    suppressPackageStartupMessages({
        library(Matrix)
        library(HDF5Array)
    })

    ## sparse on-disk data, in a budget that fits it but only blocks of a few
    ## columns, is loaded by blocks into the same object as with as()
    set.seed(123)
    m <- rsparsematrix(200, 90, density=0.1)
    dimnames(m) <- list(paste0("g", 1:200), paste0("s", 1:90))
    for (type in c("double", "integer")) {
        if (type == "integer")
            m@x <- round(abs(m@x) * 10) + 1
        x <- as(writeHDF5Array(m, as.sparse=TRUE), "DelayedMatrix")
        if (type == "integer")
            x <- DelayedArray::DelayedArray(x)
        type(x) <- type
        insize <- as.numeric(object.size(as(x, "SVT_SparseMatrix")))
        ## budget for blocks of 5 columns
        bytespercol <- 200 * if (type == "integer") 4 else 8
        maxmem <- insize + 2 * 2 * 5 * bytespercol
        res <- GSVA:::.load_sparse_by_blocks(x, maxmem=maxmem, insize=insize)
        checkTrue(is(res, "SVT_SparseMatrix"))
        checkIdentical(type(res), type)
        checkIdentical(res, as(x, "SVT_SparseMatrix"))
        ## without a budget, the default block size of DelayedArray
        checkIdentical(GSVA:::.load_sparse_by_blocks(x),
                       as(x, "SVT_SparseMatrix"))
    }
}
