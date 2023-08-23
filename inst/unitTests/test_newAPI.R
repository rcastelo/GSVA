
###
### 2023-08-14  axel: unit tests for GSVA's new API
###

test_newAPI <- function() {
    tid <- GSVA:::generateTestInputData()

    gs <- tid[["geneSets"]]
    xf <- tid[["mxFloat"]]
    xi <- tid[["mxCount"]]

    g1 <- gsva(expr = xf, gset.idx.list = gs)
    checkIdentical(class(xf), class(g1))
    checkEquals(names(gs), rownames(g1))
    checkEquals(colnames(xf), colnames(g1))
    
    g2 <- gsva(expr = xi, gset.idx.list = gs, kcdf = "Poisson")
    checkIdentical(class(xi), class(g2))
    checkEquals(names(gs), rownames(g2))
    checkEquals(colnames(xi), colnames(g2))

    g3 <- gsva(expr = xf, gset.idx.list = gs, kcdf = "none")
    checkIdentical(class(xf), class(g3))
    checkEquals(names(gs), rownames(g3))
    checkEquals(colnames(xf), colnames(g3))

    g4 <- gsva(expr = xf, gset.idx.list = gs, method = "gsva")
    checkIdentical(class(xf), class(g4))
    checkEquals(names(gs), rownames(g4))
    checkEquals(colnames(xf), colnames(g4))

    g5 <- gsva(expr = xi, gset.idx.list = gs, method = "gsva", kcdf = "Poisson")
    checkIdentical(class(xi), class(g5))
    checkEquals(names(gs), rownames(g5))
    checkEquals(colnames(xi), colnames(g5))

    g6 <- gsva(expr = xf, gset.idx.list = gs, method = "gsva", kcdf = "none")
    checkIdentical(class(xf), class(g6))
    checkEquals(names(gs), rownames(g6))
    checkEquals(colnames(xf), colnames(g6))

    checkIdentical(g1, g4)
    checkIdentical(g2, g5)
    checkIdentical(g3, g6)

    p1 <- gsva(expr = xf, gset.idx.list = gs, method = "plage")
    checkIdentical(class(xf), class(p1))
    checkEquals(names(gs), rownames(p1))
    checkEquals(colnames(xf), colnames(p1))

    p2 <- gsva(expr = xf, gset.idx.list = gs, param = plageParam())
    checkIdentical(class(xf), class(p2))
    checkEquals(names(gs), rownames(p2))
    checkEquals(colnames(xf), colnames(p2))

    checkIdentical(p1, p2)
    
    z1 <- gsva(expr = xf, gset.idx.list = gs, method = "zscore")

    s1 <- gsva(expr = xf, gset.idx.list = gs, method = "ssgsea")
}
