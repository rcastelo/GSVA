
###
### 2023-08-14  axel: unit tests for GSVA's new API
###

test_newAPI <- function() {
    tid <- GSVA:::generateTestInputData()

    gs <- tid[["geneSets"]]
    xf <- tid[["continuousData"]]
    xi <- tid[["countData"]]

    ## check: container class, dimnames, values; with default method 'gsva'; Gaussian kernel
    g1 <- gsva(expr = xf, gset.idx.list = gs)
    checkIdentical(class(xf), class(g1))
    checkEquals(names(gs), rownames(g1))
    checkEquals(colnames(xf), colnames(g1))
    checkTrue((min(g1) >= -1 && (max(g1) <= 1)))
    checkTrue(!any(is.na(g1)))
    
    ## check: container class, dimnames, values; with default method 'gsva'; Poisson kernel
    g2 <- gsva(expr = xi, gset.idx.list = gs, kcdf = "Poisson")
    checkIdentical(class(xi), class(g2))
    checkEquals(names(gs), rownames(g2))
    checkEquals(colnames(xi), colnames(g2))
    checkTrue((min(g2) >= -1 && (max(g2) <= 1)))
    checkTrue(!any(is.na(g2)))

    ## check: container class, dimnames, values; with default method 'gsva'; no kernel
    g3 <- gsva(expr = xf, gset.idx.list = gs, kcdf = "none")
    checkIdentical(class(xf), class(g3))
    checkEquals(names(gs), rownames(g3))
    checkEquals(colnames(xf), colnames(g3))
    checkTrue((min(g3) >= -1 && (max(g3) <= 1)))
    checkTrue(!any(is.na(g3)))

    ## check: container class, dimnames, values; with explicit method 'gsva'; Gaussian kernel
    g4 <- gsva(expr = xf, gset.idx.list = gs, method = "gsva")
    checkIdentical(class(xf), class(g4))
    checkEquals(names(gs), rownames(g4))
    checkEquals(colnames(xf), colnames(g4))
    checkTrue((min(g4) >= -1 && (max(g4) <= 1)))
    checkTrue(!any(is.na(g4)))

    ## check: container class, dimnames, values; with explicit method 'gsva'; Poisson kernel
    g5 <- gsva(expr = xi, gset.idx.list = gs, method = "gsva", kcdf = "Poisson")
    checkIdentical(class(xi), class(g5))
    checkEquals(names(gs), rownames(g5))
    checkEquals(colnames(xi), colnames(g5))
    checkTrue((min(g5) >= -1 && (max(g5) <= 1)))
    checkTrue(!any(is.na(g5)))

    ## check: container class, dimnames, values; with explicit method 'gsva'; no kernel
    g6 <- gsva(expr = xf, gset.idx.list = gs, method = "gsva", kcdf = "none")
    checkIdentical(class(xf), class(g6))
    checkEquals(names(gs), rownames(g6))
    checkEquals(colnames(xf), colnames(g6))
    checkTrue((min(g6) >= -1 && (max(g6) <= 1)))
    checkTrue(!any(is.na(g6)))

    ## check: container class, dimnames, values; method 'gsva' by param object; Gaussian kernel
    g7 <- gsva(expr = xf, gset.idx.list = gs, param = gsvaParam())
    checkIdentical(class(xf), class(g7))
    checkEquals(names(gs), rownames(g7))
    checkEquals(colnames(xf), colnames(g7))
    checkTrue((min(g7) >= -1 && (max(g7) <= 1)))
    checkTrue(!any(is.na(g7)))

    ## check: container class, dimnames, values; method 'gsva' by param object; Poisson kernel
    g8 <- gsva(expr = xi, gset.idx.list = gs, param = gsvaParam(kcdf = "Poisson"))
    checkIdentical(class(xf), class(g8))
    checkEquals(names(gs), rownames(g8))
    checkEquals(colnames(xf), colnames(g8))
    checkTrue((min(g8) >= -1 && (max(g8) <= 1)))
    checkTrue(!any(is.na(g8)))

    ## check: container class, dimnames, values; method 'gsva' by param object; no kernel
    g9 <- gsva(expr = xf, gset.idx.list = gs, param = gsvaParam(kcdf = "none"))
    checkIdentical(class(xf), class(g9))
    checkEquals(names(gs), rownames(g9))
    checkEquals(colnames(xf), colnames(g9))
    checkTrue((min(g9) >= -1 && (max(g9) <= 1)))
    checkTrue(!any(is.na(g9)))

    ## check: container class, dimnames, values; method 'gsva' by param object; Gaussian kernel; unnamed arguments
    g10 <- gsva(xf, gs, gsvaParam())
    checkIdentical(class(xf), class(g10))
    checkEquals(names(gs), rownames(g10))
    checkEquals(colnames(xf), colnames(g10))
    checkTrue((min(g10) >= -1 && (max(g10) <= 1)))
    checkTrue(!any(is.na(g10)))

    ## check: container class, dimnames, values; method 'gsva' by param object; Poisson kernel; unnamed arguments
    g11 <- gsva(xi, gs, gsvaParam(kcdf = "Poisson"))
    checkIdentical(class(xf), class(g11))
    checkEquals(names(gs), rownames(g11))
    checkEquals(colnames(xf), colnames(g11))
    checkTrue((min(g11) >= -1 && (max(g11) <= 1)))
    checkTrue(!any(is.na(g11)))

    ## check: container class, dimnames, values; method 'gsva' by param object; no kernel; unnamed arguments
    g12 <- gsva(xf, gs, gsvaParam(kcdf = "none"))
    checkIdentical(class(xf), class(g12))
    checkEquals(names(gs), rownames(g12))
    checkEquals(colnames(xf), colnames(g12))
    checkTrue((min(g12) >= -1 && (max(g12) <= 1)))
    checkTrue(!any(is.na(g12)))

    ## compare results using default method vs explicit method
    checkIdentical(g1, g4)
    checkIdentical(g2, g5)
    checkIdentical(g3, g6)
    checkIdentical(g1, g7)
    checkIdentical(g2, g8)
    checkIdentical(g3, g9)
    checkIdentical(g1, g10)
    checkIdentical(g2, g11)
    checkIdentical(g3, g12)


    ## check: container class, dimnames, values; method 'plage' by name
    p1 <- gsva(expr = xf, gset.idx.list = gs, method = "plage")
    checkIdentical(class(xf), class(p1))
    checkEquals(names(gs), rownames(p1))
    checkEquals(colnames(xf), colnames(p1))
    checkTrue((min(p1) >= -1 && (max(p1) <= 1)))
    checkTrue(!any(is.na(p1)))

    ## check: container class, dimnames, values; method 'plage' by param object
    p2 <- gsva(expr = xf, gset.idx.list = gs, param = plageParam())
    checkIdentical(class(xf), class(p2))
    checkEquals(names(gs), rownames(p2))
    checkEquals(colnames(xf), colnames(p2))
    checkTrue((min(p2) >= -1 && (max(p2) <= 1)))
    checkTrue(!any(is.na(p2)))

    ## check: container class, dimnames, values; method 'plage' by param object; unnamed arguments
    p3 <- gsva(xf, gs, plageParam())
    checkIdentical(class(xf), class(p3))
    checkEquals(names(gs), rownames(p3))
    checkEquals(colnames(xf), colnames(p3))
    checkTrue((min(p3) >= -1 && (max(p3) <= 1)))
    checkTrue(!any(is.na(p3)))

    ## compare results using method name vs param object
    checkIdentical(p1, p2)
    checkIdentical(p1, p3)
    

    ## check: container class, dimnames, values; method 'zscore' by name
    z1 <- gsva(expr = xf, gset.idx.list = gs, method = "zscore")
    checkIdentical(class(xf), class(z1))
    checkEquals(names(gs), rownames(z1))
    checkEquals(colnames(xf), colnames(z1))
    ## checkTrue((min(z1) >= -1 && (max(z1) <= 1)))
    checkTrue(!any(is.na(z1)))

    ## check: container class, dimnames, values; method 'zscore' by param object
    z2 <- gsva(expr = xf, gset.idx.list = gs, param = zscoreParam())
    checkIdentical(class(xf), class(z2))
    checkEquals(names(gs), rownames(z2))
    checkEquals(colnames(xf), colnames(z2))
    ## checkTrue((min(z2) >= -1 && (max(z2) <= 1)))
    checkTrue(!any(is.na(z2)))

    ## check: container class, dimnames, values; method 'zscore' by param object; unnamed arguments
    z3 <- gsva(xf, gs, zscoreParam())
    checkIdentical(class(xf), class(z3))
    checkEquals(names(gs), rownames(z3))
    checkEquals(colnames(xf), colnames(z3))
    ## checkTrue((min(z3) >= -1 && (max(z3) <= 1)))
    checkTrue(!any(is.na(z3)))

    ## compare results of using method name vs param object
    checkIdentical(z1, z2)
    checkIdentical(z1, z3)
    

    ## check: container class, dimnames, values; method ssgsea' by name; default parameters
    s1 <- gsva(expr = xf, gset.idx.list = gs, method = "ssgsea")
    checkIdentical(class(xf), class(s1))
    checkEquals(names(gs), rownames(s1))
    checkEquals(colnames(xf), colnames(s1))
    checkTrue((min(s1) >= -1 && (max(s1) <= 1)))
    checkTrue(!any(is.na(s1)))

    ## check: container class, dimnames, values; method ssgsea' by param object; default parameters
    s2 <- gsva(expr = xf, gset.idx.list = gs, param = ssgseaParam())
    checkIdentical(class(xf), class(s2))
    checkEquals(names(gs), rownames(s2))
    checkEquals(colnames(xf), colnames(s2))
    checkTrue((min(s2) >= -1 && (max(s2) <= 1)))
    checkTrue(!any(is.na(s2)))

    ## check: container class, dimnames, values; method ssgsea' by name; without normalization
    s3 <- gsva(expr = xf, gset.idx.list = gs, method = "ssgsea", ssgsea.norm = FALSE)
    checkIdentical(class(xf), class(s3))
    checkEquals(names(gs), rownames(s3))
    checkEquals(colnames(xf), colnames(s3))
    ## checkTrue((min(s3) >= -1 && (max(s3) <= 1)))
    checkTrue(!any(is.na(s3)))

    ## check: container class, dimnames, values; method ssgsea' by param object; without normalization
    s4 <- gsva(expr = xf, gset.idx.list = gs, param = ssgseaParam(normalize = FALSE))
    checkIdentical(class(xf), class(s4))
    checkEquals(names(gs), rownames(s4))
    checkEquals(colnames(xf), colnames(s4))
    ## checkTrue((min(s4) >= -1 && (max(s4) <= 1)))
    checkTrue(!any(is.na(s4)))

    ## check: container class, dimnames, values; method ssgsea' by param object; default parameters; unnamed arguments
    s5 <- gsva(xf, gs, ssgseaParam())
    checkIdentical(class(xf), class(s5))
    checkEquals(names(gs), rownames(s5))
    checkEquals(colnames(xf), colnames(s5))
    checkTrue((min(s5) >= -1 && (max(s5) <= 1)))
    checkTrue(!any(is.na(s5)))

    ## check: container class, dimnames, values; method ssgsea' by param object; without normalization; unnames arguments
    s6 <- gsva(xf, gs, ssgseaParam(normalize = FALSE))
    checkIdentical(class(xf), class(s6))
    checkEquals(names(gs), rownames(s6))
    checkEquals(colnames(xf), colnames(s6))
    ## checkTrue((min(s6) >= -1 && (max(s6) <= 1)))
    checkTrue(!any(is.na(s6)))

    ## compare results of using method name vs param object
    checkIdentical(s1, s2)
    checkIdentical(s3, s4)
    checkIdentical(s1, s5)
    checkIdentical(s3, s6)
}
