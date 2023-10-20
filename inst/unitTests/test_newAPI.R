
###
### unit tests for GSVA's new API
###

test_newAPI <- function() {
    p <- 10; n <- 30; ngs <- 5 # number of genes, samples, gene sets
    xf <- matrix(rnorm(n * p), nrow=p, ncol=n,
                 dimnames=list(paste0("g", seq.int(p)),
                               paste0("s", seq.int(n))))
    xi <- round(abs(xf * 1e6))
    gs <- replicate(ngs, sample(rownames(xf), 25, replace=TRUE), simplify=FALSE)

    ## library("Biobase")
    ## es1 <- ExpressionSet(xf)
    ## es2 <- ExpressionSet(xi)

    ## library("SummarizedExperiment")
    ## se <- SummarizedExperiment(assays=list(expr=xf, counts=xi))

    g1 <- gsva(gsvaParam(exprData=xf, geneSets=gs))
    checkIdentical(class(xf), class(g1))
    checkEquals(names(gs), rownames(g1))
    checkEquals(colnames(xf), colnames(g1))
    ## checkTrue((min(g1) >= -1 && (max(g1) <= 1)))
    checkTrue(!any(is.na(g1)))

    g2 <- gsva(gsvaParam(exprData=xi, geneSets=gs, kcdf="Poisson"))
    checkIdentical(class(xi), class(g2))
    checkEquals(names(gs), rownames(g2))
    checkEquals(colnames(xi), colnames(g2))
    ## checkTrue((min(g2) >= -1 && (max(g2) <= 1)))
    checkTrue(!any(is.na(g2)))

    g3 <- gsva(gsvaParam(exprData=xf, geneSets=gs, kcdf="none"))
    checkIdentical(class(xf), class(g3))
    checkEquals(names(gs), rownames(g3))
    checkEquals(colnames(xf), colnames(g3))
    ## checkTrue((min(g3) >= -1 && (max(g3) <= 1)))
    checkTrue(!any(is.na(g3)))

    g4 <- gsva(gsvaParam(xf, gs))
    checkIdentical(class(xf), class(g4))
    checkEquals(names(gs), rownames(g4))
    checkEquals(colnames(xf), colnames(g4))
    ## checkTrue((min(g4) >= -1 && (max(g4) <= 1)))
    checkTrue(!any(is.na(g4)))

    g5 <- gsva(gsvaParam(xi, gs, kcdf = "Poisson"))
    checkIdentical(class(xf), class(g5))
    checkEquals(names(gs), rownames(g5))
    checkEquals(colnames(xf), colnames(g5))
    ## checkTrue((min(g5) >= -1 && (max(g5) <= 1)))
    checkTrue(!any(is.na(g5)))

    g6 <- gsva(gsvaParam(xf, gs, kcdf = "none"))
    checkIdentical(class(xf), class(g6))
    checkEquals(names(gs), rownames(g6))
    checkEquals(colnames(xf), colnames(g6))
    ## checkTrue((min(g6) >= -1 && (max(g6) <= 1)))
    checkTrue(!any(is.na(g6)))

    checkIdentical(g1, g4)
    checkIdentical(g2, g5)
    checkIdentical(g3, g6)

    p1 <- gsva(plageParam(exprData=xf, geneSets=gs))
    checkIdentical(class(xf), class(p1))
    checkEquals(names(gs), rownames(p1))
    checkEquals(colnames(xf), colnames(p1))
    ## checkTrue((min(p1) >= -1 && (max(p1) <= 1)))
    checkTrue(!any(is.na(p1)))

    p2 <- gsva(plageParam(xf, gs))
    checkIdentical(class(xf), class(p2))
    checkEquals(names(gs), rownames(p2))
    checkEquals(colnames(xf), colnames(p2))
    ## checkTrue((min(p2) >= -1 && (max(p2) <= 1)))
    checkTrue(!any(is.na(p2)))

    checkIdentical(p1, p2)

    z1 <- gsva(zscoreParam(exprData=xf, geneSets=gs))
    checkIdentical(class(xf), class(z1))
    checkEquals(names(gs), rownames(z1))
    checkEquals(colnames(xf), colnames(z1))
    ## checkTrue((min(z1) >= -1 && (max(z1) <= 1)))
    checkTrue(!any(is.na(z1)))

    z2 <- gsva(zscoreParam(xf, gs))
    checkIdentical(class(xf), class(z2))
    checkEquals(names(gs), rownames(z2))
    checkEquals(colnames(xf), colnames(z2))
    ## checkTrue((min(z2) >= -1 && (max(z2) <= 1)))
    checkTrue(!any(is.na(z2)))

    checkIdentical(z1, z2)

    s1 <- gsva(ssgseaParam(exprData=xf, geneSets=gs))
    checkIdentical(class(xf), class(s1))
    checkEquals(names(gs), rownames(s1))
    checkEquals(colnames(xf), colnames(s1))
    ## checkTrue((min(s1) >= -1 && (max(s1) <= 1)))
    checkTrue(!any(is.na(s1)))

    s2 <- gsva(ssgseaParam(xf, gs))
    checkIdentical(class(xf), class(s2))
    checkEquals(names(gs), rownames(s2))
    checkEquals(colnames(xf), colnames(s2))
    ## checkTrue((min(s2) >= -1 && (max(s2) <= 1)))
    checkTrue(!any(is.na(s2)))

    checkIdentical(s1, s2)
}
