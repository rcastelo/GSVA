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

test_gsvaWalkDivert <- function() {
    message("Running unit tests for the CPU diversion of overflow-prone walks.")

    safe <- GSVA:::.gsva_cuda_walk_safe
    checkTrue(safe(20000L, 1500L, 6))
    ## the weights, not the walk total, set the bound: at tau=9 the largest
    ## weight is ~1e36 and fits gsva_float_t, the total over the 1500 genes is
    ## ~1.5e39 and does not; the device side of this case is the tau=9 run below
    checkTrue(safe(20000L, 1500L, 9))
    checkTrue(safe(20000L, 1500L, 0))
    checkTrue(safe(2L, 2L, 1e6))
    checkTrue(!safe(10L, 1L, Inf))
    checkTrue(!safe(10L, 1L, NaN))
    checkTrue(!safe(10L, 1L, -1))

    ## that safety bound is a property of the float type the walk was compiled
    ## with. tau=100 at G=10 overflows float but not double, so it separates
    ## the builds
    narrowWalk <- !safe(10L, 1L, 100L)
    if (narrowWalk) {
        ## tau=10 at G=20000 is over the limit: tau * log(G/2) is 92.1, above
        ## log(FLT_MAX) = 88.72
        checkIdentical(c(FALSE, FALSE), safe(20000L, c(1500L, 2L), 10))
        ## the boundary of that limit: tau * log(G/2) < log(FLT_MAX) holds up to
        ## tau=55 at G=10, since 5^55 < FLT_MAX < 5^56
        checkTrue(safe(10L, 1L, 55))
        checkTrue(!safe(10L, 1L, 56))
    } else {
        message("Walk weights are doubles; FLT_MAX-bound assertions skipped.")
    }

    if (!GSVA:::.gsva_cuda_available()) {
        message("CUDA not available; skipping GPU diversion assertions.")
        return(invisible(NULL))
    }

    p <- 20000
    n <- 10
    genes <- paste0("g", seq_len(p))
    set.seed(123)
    y <- matrix(rnorm(p*n), nrow=p, ncol=n,
                dimnames=list(genes, paste0("s", seq_len(n))))
    set.seed(456)
    gsets <- list(bigset=sample(genes, 1500))

    ## tau=9 at G=20000 is the case the double accumulator rescues: the
    ## largest weight is (G/2)^9 ~ 1e36, which fits gsva_float_t, while the
    ## 1500-term total is ~1.5e39, which does not.  The bound tracks the type
    ## the weights saturate in and charges the walk no extra log(k) headroom,
    ## so this set stays on the GPU and must match the CPU walk.  Weights are
    ## computed with GSVA_POW into gsva_float_t and the running peak and valley
    ## are held in gsva_float_t (struct PeakValley), so terms and extremes carry
    ## a float rounding even though the sum is double: the relative error
    ## against the CPU walk is ~1e-5 at p=20000, looser than
    ## checkEqualsNumeric's default tolerance.
    risky_par <- gsvaParam(y, gsets, tau=9)
    checkTrue(safe(p, lengths(gsets), 9))
    msgs9 <- character()
    es_gpu9 <- withCallingHandlers(gsva(risky_par, device="gpu", verbose=TRUE),
                                   message=function(m) {
                                     msgs9 <<- c(msgs9, conditionMessage(m))
                                     invokeRestart("muffleMessage") })
    checkTrue(all(is.finite(es_gpu9)))
    checkTrue(!any(grepl("at risk of float overflow", msgs9)))
    es_cpu9 <- gsva(risky_par, device="cpu", verbose=FALSE)
    .check_gpu_equiv(es_cpu9, es_gpu9, tol=1e-4)

    ## a set whose weights would saturate the GPU walk is diverted to the
    ## CPU and must score identically to the CPU-only run, without NA or
    ## sentinel values leaking from the kernel.  Such a tau exists only while
    ## the walk is narrower than the CPU's own accumulation; in a double build
    ## GPU and CPU saturate at the same point, so there is nothing to divert.
    if (!narrowWalk) {
        message("Walk weights are doubles; saturation-diversion assertions skipped.")
    } else {
        bigpar <- gsvaParam(y, gsets, tau=10)
        es_cpu <- gsva(bigpar, device="cpu", verbose=FALSE)
        es_gpu <- gsva(bigpar, device="gpu", verbose=FALSE)
        checkTrue(all(is.finite(es_gpu)))
        .check_gpu_equiv(es_cpu, es_gpu)

        ## just below the overflow risk
        okpar <- gsvaParam(y, gsets, tau=3)
        checkTrue(GSVA:::.gsva_cuda_walk_safe(p, lengths(gsets), 3))
        es_cpu_ok <- gsva(okpar, device="cpu", verbose=FALSE)
        es_gpu_ok <- gsva(okpar, device="gpu", verbose=FALSE)
        checkTrue(all(is.finite(es_gpu_ok)))
        .check_gpu_equiv(es_cpu_ok, es_gpu_ok, tol=1e-4)
    }

    ## a gene set covering all genes is diverted to the CPU
    p2 <- 100
    n2 <- 6
    genes2 <- paste0("h", seq_len(p2))
    set.seed(789)
    y2 <- matrix(rnorm(p2*n2), nrow=p2, ncol=n2,
                 dimnames=list(genes2, paste0("c", seq_len(n2))))
    set.seed(987)
    fsets <- list(fullset=genes2, ctrl=sample(genes2, 30))

    ## neither set is diverted by size or by the saturation predicate:
    ## the full-gene-set term must be what sends `fullset` to the CPU
    checkTrue(all(lengths(fsets) <= GSVA:::.gsva_cuda_max_gset_size()))
    checkTrue(all(GSVA:::.gsva_cuda_walk_safe(p2, lengths(fsets), 1)))

    fullpar <- gsvaParam(y2, fsets, tau=1)
    es_cpu_full <- suppressWarnings(gsva(fullpar, device="cpu", verbose=FALSE))
    checkTrue(all(is.na(es_cpu_full["fullset", ])))
    checkTrue(all(is.finite(es_cpu_full["ctrl", ])))

    es_gpu_full <- suppressWarnings(gsva(fullpar, device="gpu", verbose=FALSE))
    checkTrue(all(is.na(es_gpu_full["fullset", ])))
    checkTrue(all(is.finite(es_gpu_full["ctrl", ])))
    .check_gpu_equiv(es_cpu_full, es_gpu_full)

    ## the diversion of full gene sets is announced in verbose mode
    msgs <- character()
    outs <- capture.output(
        suppressWarnings(
            withCallingHandlers(
                gsva(fullpar, device="gpu", verbose=TRUE),
                message=function(m) {
                    msgs <<- c(msgs, conditionMessage(m))
                    invokeRestart("muffleMessage")
                })))
    checkTrue(grepl("containing all genes",
                    paste(c(outs, msgs), collapse="\n"), fixed=TRUE))

    ## the same diversion applies to the sparse layout
    p3 <- 200
    n3 <- 6
    genes3 <- paste0("q", seq_len(p3))
    set.seed(654)
    cnt3 <- matrix(rpois(p3*n3, 0.8), nrow=p3, ncol=n3,
                   dimnames=list(genes3, paste0("d", seq_len(n3))))
    ssets <- list(fullset=genes3, ctrl=sample(genes3, 40))
    spar <- gsvaParam(as(cnt3, "dgCMatrix"), ssets, kcdf="none", tau=1)
    es_cpu_sp <- suppressWarnings(gsva(spar, device="cpu", verbose=FALSE))
    es_gpu_sp <- suppressWarnings(gsva(spar, device="gpu", verbose=FALSE))
    checkTrue(all(is.na(es_cpu_sp["fullset", ])))
    checkTrue(all(is.na(es_gpu_sp["fullset", ])))
    checkTrue(all(is.finite(es_gpu_sp["ctrl", ])))
    .check_gpu_equiv(es_cpu_sp, es_gpu_sp)
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
