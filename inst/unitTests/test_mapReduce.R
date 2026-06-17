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

    ## process it as if it were single-cell RNA-seq data
    sce <- quickRnaQc.se(sce, subsets=list(mito=rep(FALSE, nrow(sce))))
    sce <- sce[, sce$keep]
    sce <- normalizeRnaCounts.se(sce, size.factors=sce$sum)

    ## build GSVA parameter object
    gsvapar <- gsvaParam(sce, gsets, verbose=FALSE)

    ## calculate row-normalized expression values without map-reduce
    gsvarnorm <- gsvaRowNorm(gsvapar, verbose=FALSE)

    ## calculate row-normalized expression values with map-reduce
    gsvarnorm2 <- gsvaReduce(gsvaMap(gsvaRowNorm, gsvapar, verbose=FALSE), verbose=FALSE)

    ## check that both approaches yield the same row-normalized expression values
    checkEqualsNumeric(assay(gsvarnorm, "gsvarnorm"),
                       assay(gsvarnorm2, "gsvarnorm"))

    ## calculate column rank values without map-reduce
    gsvaranks <- gsvaColRanks(gsvarnorm, verbose=FALSE)

    ## calculate column rank values with map-reduce
    gsvaranks2 <- gsvaReduce(gsvaMap(gsvaColRanks, gsvarnorm, verbose=FALSE), verbose=FALSE)

    ## check that both approaches yield the same column rank values
    checkEqualsNumeric(assay(gsvaranks, "gsvaranks"),
                       assay(gsvaranks2, "gsvaranks"))

    ## calculate column GSVA scores without map-reduce
    gsvaes <- gsvaColScores(gsvaranks, verbose=FALSE)

    ## calculate column GSVA scores with map-reduce
    gsvaes2 <- gsvaReduce(gsvaMap(gsvaColScores, gsvaranks, verbose=FALSE), verbose=FALSE)

    ## check that both approaches yield the same column GSVA scores
    checkEqualsNumeric(assay(gsvaes, "es"),
                       assay(gsvaes2, "es"))

}
