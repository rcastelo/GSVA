test_mapReduce <- function() {

    message("Running unit tests for map reduce")

    suppressPackageStartupMessages({
        library(Matrix)
        library(GSEABase)
        library(SingleCellExperiment)
        library(scrapper)
    })

    p <- 10 ## number of genes
    n <- 30 ## number of samples

    ## consider three disjoint gene sets
    gsets <- list(gset1=paste0("g", 1:3),
                  gset2=paste0("g", 4:6),
                  gset3=paste0("g", 7:10))

    ## build a random sparse count matrix with 85% sparsity
    set.seed(123)
    cnt <- integer(n*p)
    idx <- sample(length(cnt), size=round(length(cnt)*0.15)) ## 85% sparsity
    cnt[idx] <- rpois(length(idx), lambda=2)+1
    cnt <- matrix(cnt, nrow=p, ncol=n,
                  dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
    cnt <- Matrix(cnt, sparse=TRUE)

    sce <- SingleCellExperiment(assays=list(counts=cnt))

    ## process as if it were single-cell RNA-seq data
    sce <- quickRnaQc.se(sce, subsets=list(mito=rep(FALSE, nrow(sce))))
    sce <- sce[, sce$keep]
    sce <- normalizeRnaCounts.se(sce, size.factors=sce$sum)

    ## build GSVA parameter object
    gsvapar <- gsvaParam(sce, gsets, verbose=FALSE)

    ## calculate row-normalized expression values without map-reduce
    gsvarownorm <- gsvaRowNorm(gsvapar, errorOnTooFewRows=FALSE, verbose=TRUE)

    ## calculate row-normalized expression values with map-reduce
    gsvarownorm2 <- gsvaReduce(gsvaMap(gsvapar, verbose=FALSE), verbose=FALSE)

    ## check that both approaches yield the same row-normalized expression values
    checkEqualsNumeric(assay(gsvarownorm, "gsvarnorm"),
                       assay(gsvarownorm2, "gsvarnorm"))
}
