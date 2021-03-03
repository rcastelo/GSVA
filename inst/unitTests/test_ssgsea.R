test_ssgsea <- function() {

  p <- 10 ## number of genes
  n <- 30 ## number of samples
  nGrp1 <- 15 ## number of samples in group 1
  nGrp2 <- n - nGrp1 ## number of samples in group 2

  ## consider three disjoint gene sets
  geneSets <- list(set1=paste("g", 1:3, sep=""),
                   set2=paste("g", 4:6, sep=""),
                   set3=paste("g", 7:10, sep=""))

  ## sample data from a normal distribution with mean 0 and st.dev. 1
  ## seeding the random number generator for the purpose of this test
  set.seed(123)
  y <- matrix(rnorm(n*p), nrow=p, ncol=n,
              dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

  ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
  y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

  ## estimate GSVA enrichment scores for the three sets
  es <- gsva(y, geneSets, method="ssgsea", verbose=FALSE)

  checkTrue(max(abs(rowMeans(es) - c(0.22893323, -0.04400744, -0.08289233))) < 1e-08)
  checkTrue(max(abs(apply(es, 1,sd) - c(0.2562903, 0.2260589, 0.2268853))) < 1e-07)

  gset.idx.list <- lapply(geneSets,
                          function(x, y) na.omit(match(x, y)),
                          rownames(y))
  fast.gset.idx.list <- lapply(geneSets,
                               function(x, y) na.omit(match(x, y)),
                               rownames(y))
  checkIdentical(gset.idx.list, fast.gset.idx.list)

  R <- apply(y, 2, function(x ,p) as.integer(rank(x)), p)
  alpha <- 0.25
  Ra <- abs(R)^alpha

  for (i in 1:n) {
    geneRanking <- order(R[, i], decreasing=TRUE)
    frw <- GSVA:::.fastRndWalk(gset.idx.list[[1]], geneRanking, i, Ra)
    rw <- GSVA:::.rndWalk(gset.idx.list[[1]], geneRanking, i, R, alpha)
    checkEqualsNumeric(rw, frw)
  }
}

