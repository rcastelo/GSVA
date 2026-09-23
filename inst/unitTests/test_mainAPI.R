
###
### unit tests for the main API
###

test_mainAPI <- function() {
    message("Running unit tests for the main API")
    
    p <- 100; n <- 30; ngs <- 5 # number of genes, samples, gene sets
    xf <- matrix(rnorm(n * p), nrow=p, ncol=n,
                 dimnames=list(paste0("g", seq.int(p)),
                               paste0("s", seq.int(n))))
    xi <- round(abs(xf * 1e6))
    gs <- replicate(ngs, sample(rownames(xf), 25, replace=FALSE), simplify=FALSE)
    names(gs) <- paste0("gs", seq_len(ngs))

    checkException(g <- gsvaRowNorm(gsvaParam(exprData=xf, geneSets=gs),
                                    verbose=FALSE, maxmem=c(1, 2)))

    g <- gsvaRowNorm(gsvaParam(xf, gs), verbose=FALSE, maxmem="1M")

    ## check discarding rows with constant values
    library(cli)
    xf2 <- rbind(rep(1, 30), xf)
    gsvapar <- gsvaParam(xf2, gs)
    out <- cli_fmt(g1 <- gsva(gsvapar, verbose=TRUE))
    checkTrue(grepl("1 rows with constant values throughout the columns", out[3]))

    g1 <- gsva(gsvaParam(exprData=xf, geneSets=gs), verbose=FALSE)
    checkIdentical(class(xf), class(g1))
    checkEquals(names(gs), rownames(g1))
    checkEquals(colnames(xf), colnames(g1))
    checkTrue(min(g1) >= -1 && (max(g1) <= 1))
    checkTrue(!any(is.na(g1)))

    g2 <- gsva(gsvaParam(exprData=xi, geneSets=gs, kcdf="Poisson"), verbose=FALSE)
    checkIdentical(class(xi), class(g2))
    checkEquals(names(gs), rownames(g2))
    checkEquals(colnames(xi), colnames(g2))
    checkTrue(min(g2) >= -1 && (max(g2) <= 1))
    checkTrue(!any(is.na(g2)))

    g3 <- gsva(gsvaParam(exprData=xf, geneSets=gs, kcdf="none"), verbose=FALSE)
    checkIdentical(class(xf), class(g3))
    checkEquals(names(gs), rownames(g3))
    checkEquals(colnames(xf), colnames(g3))
    checkTrue(min(g3) >= -1 && (max(g3) <= 1))
    checkTrue(!any(is.na(g3)))

    ## CLR cannot be applied to nonpositive values
    checkException(g4 <- gsva(gsvaParam(exprData=xf, geneSets=gs,
                                        rowNorm="clr"), verbose=FALSE))
    g4 <- gsva(gsvaParam(exprData=abs(xf), geneSets=gs,
                         rowNorm="clr"), verbose=FALSE)
    checkIdentical(class(xf), class(g4))
    checkEquals(names(gs), rownames(g4))
    checkEquals(colnames(xf), colnames(g4))
    checkTrue(min(g4) >= -1 && (max(g4) <= 1))
    checkTrue(!any(is.na(g4)))

    g5 <- gsva(gsvaParam(exprData=abs(xf), geneSets=gs,
                         rowNorm="none"), verbose=FALSE)
    checkIdentical(class(xf), class(g5))
    checkEquals(names(gs), rownames(g5))
    checkEquals(colnames(xf), colnames(g5))
    checkTrue(min(g5) >= -1 && (max(g5) <= 1))
    checkTrue(!any(is.na(g5)))

    p1 <- gsva(plageParam(exprData=xf, geneSets=gs), verbose=FALSE)
    checkIdentical(class(xf), class(p1))
    checkEquals(names(gs), rownames(p1))
    checkEquals(colnames(xf), colnames(p1))
    checkTrue(min(p1) >= -1 && (max(p1) <= 1))
    checkTrue(!any(is.na(p1)))

    s1 <- gsva(ssgseaParam(exprData=xf, geneSets=gs), verbose=FALSE)
    checkIdentical(class(xf), class(s1))
    checkEquals(names(gs), rownames(s1))
    checkEquals(colnames(xf), colnames(s1))
    checkTrue(min(s1) >= -1 && (max(s1) <= 1))
    checkTrue(!any(is.na(s1)))

    z1 <- gsva(zscoreParam(exprData=xf, geneSets=gs), verbose=FALSE)
    checkIdentical(class(xf), class(z1))
    checkEquals(names(gs), rownames(z1))
    checkEquals(colnames(xf), colnames(z1))
    checkTrue((min(z1) >= -4 && (max(z1) <= 4)))
    checkTrue(!any(is.na(z1)))

    ngs <- 40 ## to check the Z-score pipeline that iterates over gene sets
    gs <- replicate(ngs, sample(rownames(xf), 25, replace=FALSE), simplify=FALSE)
    names(gs) <- paste0("gs", seq_len(ngs))

    z2 <- gsva(zscoreParam(exprData=xf, geneSets=gs), verbose=FALSE)
    checkIdentical(class(xf), class(z2))
    checkEquals(names(gs), rownames(z2))
    checkEquals(colnames(xf), colnames(z2))
    checkTrue((min(z2) >= -4 && (max(z2) <= 4)))
    checkTrue(!any(is.na(z2)))

    a1 <- gsva(avgParam(exprData=xf, geneSets=gs), verbose=FALSE)
    checkIdentical(class(xf), class(a1))
    checkEquals(names(gs), rownames(a1))
    checkEquals(colnames(xf), colnames(a1))
    checkTrue((min(a1) >= -1 && (max(a1) <= 1)))
    checkTrue(!any(is.na(a1)))

    ngs <- 40 ## to check the average pipeline that iterates over gene sets
    gs <- replicate(ngs, sample(rownames(xf), 25, replace=FALSE), simplify=FALSE)
    names(gs) <- paste0("gs", seq_len(ngs))

    a2 <- gsva(avgParam(exprData=xf, geneSets=gs), verbose=FALSE)
    checkIdentical(class(xf), class(a2))
    checkEquals(names(gs), rownames(a2))
    checkEquals(colnames(xf), colnames(a2))
    checkTrue((min(a2) >= -1 && (max(a2) <= 1)))
    checkTrue(!any(is.na(a2)))
}
