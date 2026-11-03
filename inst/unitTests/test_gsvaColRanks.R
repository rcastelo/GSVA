source(system.file("unitTests", "test_helpers.R", package = "GSVA"))

test_gsvaColRanks <- function() {
    message("Running unit tests for GSVA column ranks")

    p <- 10 ## number of genes
    n <- 30 ## number of samples
    nGrp1 <- 15 ## number of samples in group 1
    nGrp2 <- n - nGrp1 ## number of samples in group 2

    ## consider three disjoint gene sets
    gsets <- list(set1=paste("g", 1:3, sep=""),
                  set2=paste("g", 4:6, sep=""),
                  set3=paste("g", 7:10, sep=""))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(y, gsets)

    ## calculate GSVA scores in one step
    gsva_es1 <- gsva(gsvapar, verbose=FALSE, device="cpu")

    ## calculate GSVA scores in three steps
    ## first calculate row-normalized expression values
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## second calculate GSVA column ranks
    gsvaranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)

    ## third calculate GSVA scores from column ranks
    gsva_es2 <- gsvaColScores(gsvaranks, device = "cpu", verbose = FALSE)

    ## both approaches to calculate GSVA scores must give
    ## the same result with the same input gene sets
    checkEqualsNumeric(gsva_es1, gsva_es2)

    ## CPU vs GPU
    .run_with_gpu({
      gsva_es3 <- gsva(gsvapar, device = "gpu", verbose = FALSE)
      .check_gpu_equiv(gsva_es1, gsva_es3)
      gsva_es2_gpu <- gsvaColScores(gsvaranks, device = "gpu", verbose = FALSE)
      .check_gpu_equiv(gsva_es2, gsva_es2_gpu)
    })

    ## check that gsvaEnrichment() works
    gsvaenrich <- gsvaEnrichment(gsvaranks, plot="no")
    checkEqualsNumeric(gsva_es1[1, 1], gsvaenrich$score, tolerance = 1e-6)
    gsvaenrich2 <- gsvaEnrichment(gsvaranks,
                                  geneSet=c("g1", "g4", "g7"),
                                  plot="no")
    checkTrue(!is.na(gsvaenrich2$score))

    ## test the ggplotting from gsvaEnrichment()
    ggp <- gsvaEnrichment(gsvaranks, plot="ggplot")
    checkTrue(is(ggp, "ggplot"))
    checkTrue(identical(gsvaenrich$stats, ggp@data))

    ## calculate again row-normalized expression values,
    ## but this time in chunks
    gsvarownorm1 <- gsvaRowNorm(gsvapar, first=1, last=5, verbose=FALSE)
    gsvarownorm2 <- gsvaRowNorm(gsvapar, first=6, last=10, verbose=FALSE)
    checkEqualsNumeric(gsvarownorm, rbind(gsvarownorm1, gsvarownorm2))
    checkException(gsvaRowNorm(gsvapar, first=10, last=6, verbose=FALSE))
    checkException(gsvaRowNorm(gsvapar, first=11, last=20, verbose=FALSE))

    ## calculate again GSVA column ranks, but this time in chunks
    gsvaranks1 <- gsvaColRanks(gsvarownorm, first=1, last=10, verbose=FALSE)
    gsvaranks2 <- gsvaColRanks(gsvarownorm, first=11, last=20, verbose=FALSE)
    gsvaranks3 <- gsvaColRanks(gsvarownorm, first=21, last=30, verbose=FALSE)
    checkEqualsNumeric(gsvaranks, cbind(gsvaranks1, gsvaranks2, gsvaranks3))

    ## calculate again GSVA scores from column ranks, but this time in chunks
    gsva_es_c1 <- gsvaColScores(gsvaranks, first=1, last=10, device="cpu", verbose=FALSE)
    gsva_es_c2 <- gsvaColScores(gsvaranks, first=11, last=20, device="cpu", verbose=FALSE)
    gsva_es_c3 <- gsvaColScores(gsvaranks, first=21, last=30, device="cpu", verbose=FALSE)
    checkEqualsNumeric(gsva_es1, cbind(gsva_es_c1, gsva_es_c2, gsva_es_c3))

    ## CPU vs GPU
    .run_with_gpu({
      gsva_es_c1_gpu <- gsvaColScores(gsvaranks, first=1, last=10, device="gpu", verbose=FALSE)
      gsva_es_c2_gpu <- gsvaColScores(gsvaranks, first=11, last=20, device="gpu", verbose=FALSE)
      gsva_es_c3_gpu <- gsvaColScores(gsvaranks, first=21, last=30, device="gpu", verbose=FALSE)
      .check_gpu_equiv(gsva_es_c1, gsva_es_c1_gpu)
      .check_gpu_equiv(gsva_es_c2, gsva_es_c2_gpu)
      .check_gpu_equiv(gsva_es_c3, gsva_es_c3_gpu)
    })
}

test_gsvaColRanksPermutation <- function() {
    message("Running unit tests for GSVA column-rank permutation property")

    ## the random walk in both the CPU and the CUDA code treats each
    ## gene in a gene set as one up-step at a distinct position; this is
    ## only sound if the rank matrix produced by gsvaColRanks is a
    ## permutation of the stored values in every column

    set.seed(123)
    p <- 40; n <- 20
    cnt <- matrix(rpois(p * n, 0.8), p, n,
                  dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    gsets <- list(set1=paste0("g", 1:5))

    ## dense ranks: each column a permutation of 1..p
    gpar <- gsvaParam(cnt, gsets)
    R <- as.matrix(gsvaColRanks(gsvaRowNorm(gpar, verbose=FALSE),
                                verbose=FALSE))
    for (j in seq_len(ncol(R)))
        checkIdentical(seq_len(p), as.integer(sort.int(R[, j])))

    ## sparse ranks: the stored values of each column are 1..nnz without
    ## repeats, the implicitly-zero genes carry rank zero.
    gsp <- gsvaParam(as(cnt, "dgCMatrix"), gsets, sparse=TRUE)
    gsvarnorm <- gsvaRowNorm(gsp, verbose=FALSE)
    Rs <- gsvaColRanks(gsvarnorm, verbose=FALSE)
    checkTrue(is(Rs, "sparseMatrix"))

    ## one stored entry per expressed gene, none of them stored as zero
    checkIdentical(nrow(as.data.frame(summary(Rs))), Matrix::nnzero(gsvarnorm))
    Rs <- as.matrix(Rs)
    Zd <- as.matrix(gsvarnorm)
    for (j in seq_len(ncol(Rs))) {
        whnz <- Zd[, j] != 0 ## expressed genes of the input column
        checkIdentical(seq_len(sum(whnz)), as.integer(sort.int(Rs[whnz, j])))
        checkTrue(all(Rs[!whnz, j] == 0))
    }
}

test_gsvaColRanksSparseGpu <- function() {
    message("Running unit tests for GSVA sparse-rank scoring on GPU.")

    p <- 10
    n <- 8
    genes <- paste0("g", 1:p)
    set.seed(123)
    y <- matrix(rnorm(p*n), nrow=p, ncol=n,
                dimnames=list(genes, paste0("s", 1:n)))
    y[y < 0] <- 0                   ## roughly half of the entries are zero
    y[, 1] <- 0                     ## make the first sample a single nonzero
    y["g1", 1] <- 1
    ys <- Matrix(y, sparse=TRUE)

    ## in the first sample only g1 is nonzero; the members of 'set2' and
    ## 'zeroset' are all zero there, their rank statistic is
    ## (nnz+1)/2-1 = 0 and their walk weights vanish, which makes those
    ## walks degenerate and their scores NA on both devices
    gsets <- list(set1=genes[c(1, 3, 5)], set2=genes[4:6], zeroset="g7")
    spar <- gsvaParam(ys, gsets, kcdf="none", tau=1)

    es_cpu <- gsva(spar, device="cpu", verbose=FALSE)
    checkTrue(is.na(es_cpu["set2", 1]))
    checkTrue(is.na(es_cpu["zeroset", 1]))
    checkTrue(all(is.finite(es_cpu["set1", ])))
    checkTrue(all(is.finite(es_cpu[-1, -1])))

    if (!GSVA:::.gsva_cuda_available()) {
        message("CUDA not available; skipping GPU assertions.")
        return(invisible(NULL))
    }

    es_gpu <- gsva(spar, device="gpu", verbose=FALSE)
    checkTrue(is.na(es_gpu["set2", 1]))
    checkTrue(is.na(es_gpu["zeroset", 1]))
    .check_gpu_equiv(es_cpu, es_gpu)
}
