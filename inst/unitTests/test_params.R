
# 2023-08-11  axel: test functions for various parameter classes and their constructors
test_plageParam <- function() {
    plgPrm <- plageParam()

    checkTrue(inherits(plgPrm, "PlageParam"))
    checkTrue(inherits(plgPrm, "EmptyParam"))
    checkTrue(!inherits(plgPrm, "GsvaParam"))

    checkException(plageParam(42), msg = "unexpected argument to plageParam()")
}

test_zscoreParam <- function() {
    zscPrm <- zscoreParam()

    checkTrue(inherits(zscPrm, "ZScoreParam"))
    checkTrue(inherits(zscPrm, "EmptyParam"))
    checkTrue(!inherits(zscPrm, "GsvaParam"))

    checkException(zscoreParam(42), msg = "unexpected argument to zscoreParam()")
}

test_ssgseaParam <- function() {
    ssgPrm <- ssgseaParam()

    checkTrue(inherits(ssgPrm, "SsGseaParam"))
    checkTrue(inherits(ssgPrm, "EmptyParam"))
    checkTrue(!inherits(ssgPrm, "GsvaParam"))

    checkEqualsNumeric(ssgPrm@alpha, 0.25)
    checkTrue(ssgPrm@normalize)

    ssgPrm <- ssgseaParam(0.5, FALSE)

    checkTrue(inherits(ssgPrm, "SsGseaParam"))
    checkTrue(inherits(ssgPrm, "EmptyParam"))
    checkTrue(!inherits(ssgPrm, "GsvaParam"))

    checkEqualsNumeric(ssgPrm@alpha, 0.5)
    checkTrue(!ssgPrm@normalize)

    checkException(ssgseaParam(0.5, FALSE, "gaga"), msg = "unexpected argument to ssgseaParam()")
    checkException(ssgseaParam(bla = "gaga"), msg = "unexpected argument to ssgseaParam()")
}

test_gsvaParam = function() {
    gsvaPrm <- gsvaParam()

    checkTrue(inherits(gsvaPrm, "GsvaParam"))
    checkTrue(inherits(gsvaPrm, "EmptyParam"))
    checkTrue(!inherits(gsvaPrm, "SsGseaParam"))

    checkEquals(gsvaPrm@kcdf, "Gaussian")
    checkEqualsNumeric(gsvaPrm@tau, 1)
    checkTrue(gsvaPrm@mx.diff)
    checkTrue(!gsvaPrm@abs.ranking)

    gsvaPrm <- gsvaParam(kcdf = "Poisson", tau = 0.5,
                         mx.diff = FALSE, abs.ranking = TRUE)

    checkTrue(inherits(gsvaPrm, "GsvaParam"))
    checkTrue(inherits(gsvaPrm, "EmptyParam"))
    checkTrue(!inherits(gsvaPrm, "SsGseaParam"))

    checkEquals(gsvaPrm@kcdf, "Poisson")
    checkEqualsNumeric(gsvaPrm@tau, 0.5)
    checkTrue(!gsvaPrm@mx.diff)
    checkTrue(gsvaPrm@abs.ranking)

    gsvaPrm <- gsvaParam(kcdf = "none")

    checkTrue(inherits(gsvaPrm, "GsvaParam"))
    checkTrue(inherits(gsvaPrm, "EmptyParam"))
    checkTrue(!inherits(gsvaPrm, "SsGseaParam"))

    checkEquals(gsvaPrm@kcdf, "none")

    checkException(gsvaParam(kcdf = "Poison"))
    checkException(gsvaParam(kcdf = "Poisson", bla= "gaga"))
    checkException(gsvaParam(kcdf = "Poison", tau = 0.42,
                             mx.diff = TRUE, abs.ranking = FALSE, 42))
}
