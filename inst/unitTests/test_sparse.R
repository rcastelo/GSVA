source(system.file("unitTests", "test_helpers.R", package = "GSVA"))

test_sparseMethods <- function(){
    message("Running unit tests for sparse methods")
    
    suppressPackageStartupMessages({
        library(Matrix)
        library(cli) ## for cli_fmt()
    })

    p <- 50 ## number of genes
    n <- 100 ## number of samples/cells
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
    m <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    M <- Matrix(m, sparse=TRUE)
    
    checkException(mg <- gsva(gsvaParam(x, gsets), verbose=TRUE))
    out <- cli_fmt(mg <- gsva(gsvaParam(m, gsets), verbose=TRUE, device="cpu"))
    out <- cli_fmt(Mg <- gsva(gsvaParam(M, gsets, sparse=FALSE), verbose=TRUE, device="cpu"))
    checkEqualsNumeric(mg, Mg)
    .run_with_gpu({
        mg_gpu <- gsva(gsvaParam(m, gsets), verbose=FALSE, device="gpu")
        .check_gpu_equiv(mg, mg_gpu)

        Mg_gpu <- gsva(gsvaParam(M, gsets, sparse=FALSE), verbose=FALSE, device="gpu")
        .check_gpu_equiv(Mg, Mg_gpu)
    })

    out <- cli_fmt(mp <- gsva(plageParam(m, gsets), verbose=TRUE))
    out <- cli_fmt(Mp <- gsva(plageParam(M, gsets), verbose=TRUE))
    checkEqualsNumeric(mp, Mp)
    
    out <- cli_fmt(mz <- gsva(zscoreParam(m, gsets), verbose=TRUE))
    out <- cli_fmt(Mz <- gsva(zscoreParam(M, gsets), verbose=TRUE))
    checkEqualsNumeric(mz, Mz)
    
    out <- cli_fmt(ms <- gsva(ssgseaParam(m, gsets), verbose=TRUE))
    out <- cli_fmt(Ms <- gsva(ssgseaParam(M, gsets), verbose=TRUE))
    checkEqualsNumeric(ms, Ms)

    M[1, 1] <- NA
    checkException(Mg <- gsva(gsvaParam(M, gsets), verbose=TRUE))
    
}

test_ecdfvals <- function() {
    message("Running unit tests for ECDF values calculations")

    ecdfvals_dense <- function(X) t(apply(X, 1, function(rx) ecdf(rx)(rx)))
    ecdfvals_sparse_to_sparse <- function(X) {
        stopifnot(is(X, "dgCMatrix")) ## QC
        for (i in 1:nrow(X)) {
            rx <- X[i, , drop=FALSE]
            vals <- sort(unique(rx@x))
            mt <- match(rx@x, vals)
            tab <- tabulate(mt, nbins=length(vals))
            ecdfvals <- cumsum(tab) / nnzero(rx)
            X[i, which(diff(rx@p) > 0)] <- ecdfvals[mt]
        }
        X
    }

    suppressPackageStartupMessages({
        library(Matrix)
	library(HDF5Array)
	library(SparseArray)
    })

    n <- 100
    p <- 100
    z <- numeric(p * n)
    nnz <- ceiling(0.05 * p * n) ## 5% nonzero values
    z[sample(1:(p*n), size=nnz, replace=FALSE)] <- rnorm(nnz)
    zz <- matrix(z, nrow=p, ncol=n)
    zzs <- Matrix(zz, sparse=TRUE)
    res_R_dense <- ecdfvals_dense(zz)
    res_C_dense <- GSVA:::compute.gene.cdf(zz, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_dense)

    res_C_sparse_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_sparse_to_dense)

    res_R_sparse_to_sparse <- ecdfvals_sparse_to_sparse(zzs)
    res_C_sparse_to_sparse <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparse_to_sparse)

    zzs <- SparseArray(zzs)
    res_C_svt_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=FALSE, kernel=FALSE, sparse=FALSE)
    checkEqualsNumeric(res_R_dense, res_C_svt_to_dense)
    res_C_svt_to_svt <- GSVA:::.ecdfvals_svt_to_svt(zzs, FALSE)
    checkEqualsNumeric(SparseArray(res_C_sparse_to_sparse), res_C_svt_to_svt)

    ## on-disk input is processed by blocks of rows, read into main memory, and
    ## its output written on disk, see GSVA:::.ondisk_blocks()
    zz <- as(zz, "HDF5Array")
    res_C_denseh5_to_denseh5 <- GSVA:::.processMatrixRows(zz, GSVA:::compute.gene.cdf, Gaussk=FALSE, kernel=FALSE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_dense, res_C_denseh5_to_denseh5)
    zzs <- as(zzs, "HDF5Array")
    res_C_sparseh5_to_denseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=FALSE, kernel=FALSE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_dense, res_C_sparseh5_to_denseh5)
    res_C_sparseh5_to_sparseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=FALSE, kernel=FALSE, sparse=TRUE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparseh5_to_sparseh5)
}

test_kcdfvals <- function() {
    message("Running unit tests for KCDF values calculations")

    kcdfegaussianvals_dense_to_dense <- function(x) {
        x <- as.matrix(x)
        t(apply(x, 1, function(rx, n) {
           bw <- sd(rx) / 4
	   if (is.na(bw) || bw == 0)
               bw = 0.001
           kecdf <- rowSums(pnorm(outer(rx, rx, FUN="-")/bw)/n)
	   -log((1-kecdf)/kecdf)
        }, ncol(x)))
    }

    kcdfegaussianvals_sparse_to_dense <- function(x) {
        x <- as.matrix(x)
        t(apply(x, 1, function(rx, n) {
           bw <- sd(rx) / 4
	   if (is.na(bw) || bw == 0)
               bw = 0.001
           rowSums(pnorm(outer(rx, rx, FUN="-")/bw)/n)
        }, ncol(x)))
    }

    kcdfgaussianvals_sparse_to_sparse <- function(X) {
        stopifnot(is(X, "dgCMatrix")) ## QC
        for (i in 1:nrow(X)) {
            rx <- X[i, , drop=FALSE]
            bw <- sd(rx@x) / 4
	    if (is.na(bw) || bw == 0)
                bw = 0.001
            kcdfvals <- rowSums(pnorm(outer(rx@x, rx@x, FUN="-")/bw)/nnzero(rx))
            X[i, which(diff(rx@p) > 0)] <- kcdfvals
        }
        X
    }

    suppressPackageStartupMessages({
        library(Matrix)
	library(HDF5Array)
	library(SparseArray)
    })

    n <- 100
    p <- 100
    z <- numeric(p * n)
    nnz <- ceiling(0.05 * p * n) ## 5% nonzero values
    z[sample(1:(p*n), size=nnz, replace=FALSE)] <- rnorm(nnz)
    zz <- matrix(z, nrow=p, ncol=n)
    zzs <- Matrix(zz, sparse=TRUE)

    res_R_sparse_to_dense <- kcdfegaussianvals_sparse_to_dense(zzs)
    res_C_sparse_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=FALSE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_sparse_to_dense, tolerance=0.00001)

    res_R_sparse_to_sparse <- kcdfgaussianvals_sparse_to_sparse(zzs)
    res_C_sparse_to_sparse <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparse_to_sparse, tolerance=0.0001)

    zzs <- SparseArray(zzs)
    res_C_svt_to_dense <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=FALSE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_svt_to_dense, tolerance=0.00001)
    res_C_svt_to_svt <- GSVA:::compute.gene.cdf(zzs, Gaussk=TRUE, kernel=TRUE, sparse=TRUE)
    checkEqualsNumeric(SparseArray(res_C_sparse_to_sparse), res_C_svt_to_svt)

    res_R_dense_to_dense <- kcdfegaussianvals_dense_to_dense(zz)
    ## on-disk input is processed by blocks of rows, read into main memory, and
    ## its output written on disk, see GSVA:::.ondisk_blocks()
    zz <- as(zz, "HDF5Array")
    res_C_denseh5_to_denseh5 <- GSVA:::.processMatrixRows(zz, GSVA:::compute.gene.cdf, Gaussk=TRUE, kernel=TRUE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_dense_to_dense, res_C_denseh5_to_denseh5, tolerance=0.0001)
    zzs <- as(zzs, "HDF5Array")
    res_C_sparseh5_to_denseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=TRUE, kernel=TRUE, sparse=FALSE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_sparse_to_dense, res_C_sparseh5_to_denseh5, tolerance=0.00001)
    res_C_sparseh5_to_sparseh5 <- GSVA:::.processMatrixRows(zzs, GSVA:::compute.gene.cdf, Gaussk=TRUE, kernel=TRUE, sparse=TRUE, verbose=FALSE, sinkout=TRUE)
    checkEqualsNumeric(res_R_sparse_to_sparse, res_C_sparseh5_to_sparseh5, tolerance=0.0001)
}

test_sparseWalkInterpolation <- function() {
    message("Running unit tests for the sparse walk interpolation")

    suppressPackageStartupMessages(library(Matrix))

    set.seed(42)
    p <- 2000L                     ## genes
    n <- 2000L                     ## cells
    m <- 100L                      ## gene set size
    n_gsets <- 200L
    nnz <- as.integer(round(0.03 * p)) ## expressed genes per cell

    ## ranks in sparse format; absent entries have rank zero
    ii <- integer(0); jj <- integer(0); xx <- numeric(0)
    for (j in seq_len(n)) {
        idx <- sample.int(p, nnz)
        ii <- c(ii, idx)
        jj <- c(jj, rep(j, nnz))
        xx <- c(xx, sample.int(nnz))
    }
    ranks <- as(Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(p, n),
                                     dimnames = list(paste0("g", seq_len(p)),
                                                     paste0("c", seq_len(n)))),
                "dgCMatrix")

    gene_names <- paste0("g", seq_len(p))
    gsetsidx <- lapply(seq_len(n_gsets),
                       function(s) match(paste0("g", sample.int(p, m)),
                                         gene_names))

    ## the set must be small enough to be dispatched to the GPU kernel, and
    ## must contain zero-expressed genes, otherwise this test would pass
    ## without ever exercising the interpolated zero tail
    checkTrue(m <= GSVA:::.gsva_cuda_max_gset_size())
    checkTrue(all(vapply(gsetsidx[seq_len(5)],
                         function(g) sum(ranks[g, 1L] == 0) > 1L,
                         logical(1))))

    run_walk <- function(device) {
        wna_env <- new.env()
        assign("w", FALSE, envir = wna_env)
        assign("u", 0, envir = wna_env)
        do.call(GSVA:::.gsva_score_genesets,
                list(R = ranks, geneSetsIdx = gsetsidx, intrnks = FALSE,
                     sparse = TRUE, maxDiff = TRUE, absRanking = FALSE,
                     tau = 1, any_na = FALSE, na_use = "everything",
                     minSize = 1L, wna_env = wna_env, verbose = FALSE,
                     device = device))
    }

    cpu <- run_walk("cpu")
    checkTrue(!any(is.na(cpu)))

    .run_with_gpu({
        gpu <- run_walk("gpu")
        maxabs <- max(abs(as.vector(cpu) - as.vector(gpu)))
        checkTrue(maxabs < 1e-6,
                  msg = sprintf(paste("largest CPU/GPU enrichment-score",
                                      "discrepancy %.3e exceeds the 1e-6",
                                      "per-element bound; a discrepancy means",
                                      "that the zero-position interpolation of",
                                      "gsea_sparse_walk_kernel disagrees with",
                                      "src/rnd_walk.c on some event position"),
                                maxabs))
    })
}

test_svt_lacunar <- function() {
    message("Running unit tests for ECDF and KCDF values on lacunar SVT leaves")

    suppressPackageStartupMessages({
        library(Matrix)
        library(SparseArray)
    })

    ## sparse count data, where the first column has only ones as nonzero
    ## values, which SparseArray stores in a 'lacunar' leaf without values,
    ## and the fourth column has only zeros
    m <- matrix(0L, nrow=6, ncol=5)
    m[c(1, 3, 5), 1] <- 1L
    m[c(1, 2, 4), 2] <- c(2L, 3L, 1L)
    m[c(1, 3, 5, 6), 3] <- c(1L, 4L, 2L, 7L)
    m[c(2, 3, 4, 6), 5] <- c(5L, 1L, 2L, 3L)

    for (type in c("integer", "double")) {
        storage.mode(m) <- type
        svt <- as(m, "SVT_SparseMatrix")
        ## check that the first column is stored in a lacunar leaf
        checkTrue(is.null(svt@SVT[[1]][[1]]))
        dgc <- as(m, "dgCMatrix")

        res_svt <- GSVA:::.ecdfvals_svt_to_svt(svt, FALSE)
        checkEqualsNumeric(as.matrix(GSVA:::.ecdfvals_sparse_to_sparse(dgc,
                                                                       FALSE)),
                           as.matrix(res_svt))

        for (Gaussk in c(TRUE, FALSE)) {
            res_svt <- GSVA:::.kcdfvals_svt_to_svt(svt, Gaussk, FALSE)
            res_dgc <- GSVA:::.kcdfvals_sparse_to_sparse(dgc, Gaussk, FALSE)
            checkEqualsNumeric(as.matrix(res_dgc), as.matrix(res_svt))
        }
    }

    ## a matrix with only zeros, whose SVT is NULL
    svt <- as(matrix(0L, nrow=3, ncol=2), "SVT_SparseMatrix")
    checkEqualsNumeric(matrix(0, nrow=3, ncol=2),
                       as.matrix(GSVA:::.ecdfvals_svt_to_svt(svt, FALSE)))
    checkEqualsNumeric(matrix(0, nrow=3, ncol=2),
                       as.matrix(GSVA:::.kcdfvals_svt_to_svt(svt, TRUE, FALSE)))
}

test_undefinedWalks <- function() {
    message("Running unit tests for undefined random walks on the GPU path")

    runWalk <- function(R, gsetsidx, sparse, device) {
        wna_env <- new.env()
        assign("w", FALSE, envir = wna_env)
        assign("u", 0, envir = wna_env)
        sco <- do.call(GSVA:::.gsva_score_genesets,
                       list(R = R, geneSetsIdx = gsetsidx, intrnks = FALSE,
                            sparse = sparse, maxDiff = TRUE,
                            absRanking = FALSE, tau = 1, any_na = FALSE,
                            na_use = "everything", minSize = 1L,
                            wna_env = wna_env, verbose = FALSE,
                            device = device))
        list(sco = sco, nundef = get("u", envir = wna_env))
    }

    set.seed(123)
    p <- 10L
    ranks <- matrix(as.numeric(sample.int(p)), ncol = 1L,
                    dimnames = list(paste0("g", seq_len(p)), "c1"))

    ## a gene set with all the genes leaves no out-of-gene-set steps, a set
    ## with only the gene at the middle rank has zero in-gene-set steps; both
    ## make the walk undefined
    midRank <- which(ranks[, 1L] == p / 2L)
    checkTrue(identical(length(midRank), 1L))
    gsetsidx <- list(allGenes = seq_len(p),
                     midRank = as.integer(midRank),
                     control = as.integer(c(1L, 3L, 7L)))

    cpu <- runWalk(ranks, gsetsidx, FALSE, "cpu")
    checkTrue(all(is.na(cpu$sco[1:2, , drop = FALSE])))
    checkTrue(!any(is.na(cpu$sco[3, , drop = FALSE])))
    checkEquals(cpu$nundef, 2)

    .run_with_gpu({
        gpu <- runWalk(ranks, gsetsidx, FALSE, "gpu")
        checkEquals(is.na(gpu$sco), is.na(cpu$sco))
        checkEqualsNumeric(cpu$sco[3, ], gpu$sco[3, ], tolerance = 1e-6)
        checkEquals(gpu$nundef, cpu$nundef)
    })

    ## hybrid dispatch: sets larger than GSVA_MAX_GPU_GSET_SIZE run on the CPU
    ## while the rest run on the GPU, and both counts have to add up
    P <- 5000L
    checkTrue(P > GSVA:::.gsva_cuda_max_gset_size())
    set.seed(123)
    ranks_big <- matrix(as.numeric(sample.int(P)), ncol = 1L,
                        dimnames = list(paste0("g", seq_len(P)), "c1"))
    mid_big <- which(ranks_big[, 1L] == P / 2L)
    checkTrue(identical(length(mid_big), 1L))
    gsetsidx_big <- list(gpuUndefined = as.integer(mid_big),
                         cpuUndefined = seq_len(P))

    cpu_big <- runWalk(ranks_big, gsetsidx_big, FALSE, "cpu")
    checkEquals(cpu_big$nundef, 2)

    .run_with_gpu({
        gpu_big <- runWalk(ranks_big, gsetsidx_big, FALSE, "gpu")
        checkEquals(is.na(gpu_big$sco), is.na(cpu_big$sco))
        checkEquals(gpu_big$nundef, cpu_big$nundef)
    })
}
