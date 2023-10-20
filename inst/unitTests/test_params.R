
###
### test functions for GSVA method parameter objects and their constructors
###

test_plageParam <- function() {
    p <- 10; n <- 30; ngs <- 5 # number of genes, samples, gene sets
    y <- matrix(rnorm(n * p), nrow=p, ncol=n,
                dimnames=list(paste0("g", seq.int(p)),
                              paste0("s", seq.int(n))))
    gs <- replicate(ngs, sample(rownames(y), 25, replace=TRUE), simplify=FALSE)

    pp <- plageParam(y, gs)

    checkTrue(inherits(pp, "plageParam"))
    checkTrue(inherits(pp, "GsvaMethodParam"))

    checkException(plageParam(42))
}

test_zscoreParam <- function() {
    p <- 10; n <- 30; ngs <- 5 # number of genes, samples, gene sets
    y <- matrix(rnorm(n * p), nrow=p, ncol=n,
                dimnames=list(paste0("g", seq.int(p)),
                              paste0("s", seq.int(n))))
    gs <- replicate(ngs, sample(rownames(y), 25, replace=TRUE), simplify=FALSE)
    
    zp <- zscoreParam(y, gs)

    checkTrue(inherits(zp, "zscoreParam"))
    checkTrue(inherits(zp, "GsvaMethodParam"))

    checkException(zscoreParam(42))
}

test_ssgseaParam <- function() {
    p <- 10; n <- 30; ngs <- 5 # number of genes, samples, gene sets
    y <- matrix(rnorm(n * p), nrow=p, ncol=n,
                dimnames=list(paste0("g", seq.int(p)),
                              paste0("s", seq.int(n))))
    gs <- replicate(ngs, sample(rownames(y), 25, replace=TRUE), simplify=FALSE)
    
    sp <- ssgseaParam(y, gs)

    checkTrue(inherits(sp, "ssgseaParam"))
    checkTrue(inherits(sp, "GsvaMethodParam"))

    checkEqualsNumeric(sp@alpha, 0.25)
    checkTrue(sp@normalize)

    sp <- ssgseaParam(y, gs, alpha=0.5, normalize=FALSE)

    checkTrue(inherits(sp, "ssgseaParam"))
    checkTrue(inherits(sp, "GsvaMethodParam"))

    checkEqualsNumeric(sp@alpha, 0.5)
    checkTrue(!sp@normalize)

    checkException(ssgseaParam(y, gs, alpha=0.5, normalize=FALSE, bla="gaga"))
}

test_gsvaParam <- function() {
    p <- 10; n <- 30; ngs <- 5 # number of genes, samples, gene sets
    y <- matrix(rnorm(n * p), nrow=p, ncol=n,
                dimnames=list(paste0("g", seq.int(p)),
                              paste0("s", seq.int(n))))
    gs <- replicate(ngs, sample(rownames(y), 25, replace=TRUE), simplify=FALSE)
    
    gp <- gsvaParam(y, gs)

    checkTrue(inherits(gp, "gsvaParam"))
    checkTrue(inherits(gp, "GsvaMethodParam"))

    checkEquals(gp@kcdf, "Gaussian")
    checkEqualsNumeric(gp@tau, 1)
    checkTrue(gp@maxDiff)
    checkTrue(!gp@absRanking)

    gp <- gsvaParam(y, gs,
                    kcdf="Poisson", tau=0.5, maxDiff=FALSE, absRanking=TRUE)

    checkTrue(inherits(gp, "gsvaParam"))
    checkTrue(inherits(gp, "GsvaMethodParam"))

    checkEquals(gp@kcdf, "Poisson")
    checkEqualsNumeric(gp@tau, 0.5)
    checkTrue(!gp@maxDiff)
    checkTrue(gp@absRanking)

    gp <- gsvaParam(y, gs, kcdf="none")

    checkTrue(inherits(gp, "gsvaParam"))
    checkTrue(inherits(gp, "GsvaMethodParam"))

    checkEquals(gp@kcdf, "none")

    checkException(gsvaParam(kcdf="Poison"))
    checkException(gsvaParam(kcdf="Poisson", bla= "gaga"))
    checkException(gsvaParam(kcdf="Poison", tau=0.42,
                             maxDiff=TRUE, absRanking=FALSE, 42))
}
