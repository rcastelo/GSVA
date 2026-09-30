## GSVA ranks out of range, e.g., after removing rows from the output of
## gsvaColRanks(), give an error instead of reading and writing out of the
## bounds of the arrays used to calculate GSVA scores
test_rankbounds <- function() {

    message("Running unit tests for GSVA ranks out of range")

    suppressPackageStartupMessages({
        library(Matrix)
        library(SparseArray)
        library(SummarizedExperiment)
    })

    scores <- function(R, any_na=FALSE)
        GSVA:::.gsva_score_genesets(R, list(1:5, 6:12),
                                    intrnks=is.integer(R[1, 1]),
                                    sparse=is(R, "sparseMatrix") ||
                                           is(R, "SVT_SparseMatrix"),
                                    maxDiff=TRUE, absRanking=FALSE, tau=1,
                                    any_na=any_na, na_use="na.rm", minSize=1L,
                                    wna_env=new.env(), verbose=FALSE)

    p <- 20 ## number of genes
    n <- 8 ## number of samples
    set.seed(123)

    ## dense ranks, as integer and double values
    R <- vapply(seq_len(n), function(j) sample.int(p), integer(p))
    checkTrue(is.matrix(scores(R)))
    checkTrue(is.matrix(scores(1 * R)))
    Rbad <- R
    Rbad[3, 2] <- p + 5L
    checkException(scores(Rbad), silent=TRUE)
    checkException(scores(1 * Rbad), silent=TRUE)
    Rbad <- R
    Rbad[3, 2] <- 0L ## dense ranks cannot be zero
    checkException(scores(Rbad), silent=TRUE)

    ## sparse ranks, whose nonzero values range between 1 and the number of
    ## nonzero values in each column
    S <- matrix(0L, nrow=p, ncol=n)
    for (j in seq_len(n)) {
        nzidx <- sort(sample.int(p, 10))
        S[nzidx, j] <- sample.int(10)
    }
    Sbad <- S
    Sbad[which(S[, 2] != 0)[1], 2] <- 15L ## only 10 nonzero values
    for (f in list(function(x) as(1 * x, "dgCMatrix"),
                   function(x) as(x, "SVT_SparseMatrix"),
                   function(x) as(1 * x, "SVT_SparseMatrix"))) {
        checkTrue(is.matrix(scores(f(S))))
        checkException(scores(f(Sbad)), silent=TRUE)
    }

    ## ranks with missing values, whose nonmissing values range between 1
    ## and the number of nonmissing values in each column
    Rna <- R
    Rna[c(4, 9), 3] <- NA
    Rna[!is.na(Rna[, 3]), 3] <- sample.int(p - 2)
    checkTrue(is.matrix(scores(Rna, any_na=TRUE)))
    Rnabad <- Rna
    Rnabad[1, 3] <- p - 1L ## only p - 2 nonmissing values
    checkException(scores(Rnabad, any_na=TRUE), silent=TRUE)

    ## ranks after removing rows, without the number of rows on which they
    ## were calculated in the metadata, as in objects produced by previous
    ## versions of GSVA, give an error in both gsvaColScores(), which calls
    ## C code, and gsvaEnrichment(), which calls R code
    gsets <- list(gset1=paste0("g", 1:5), gset2=paste0("g", 6:12))
    y <- matrix(rnorm(n*p*5), nrow=p, ncol=n*5,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:(n*5))))
    for (input in list(y, as(Matrix(pmax(y, 0), sparse=TRUE), "dgCMatrix"))) {
        se <- SummarizedExperiment(assays=list(exprs=input))
        gsvacolranks <- gsvaColRanks(gsvaRowNorm(gsvaParam(se, gsets,
                                                           verbose=FALSE),
                                                 verbose=FALSE),
                                     verbose=FALSE)
        subranks <- gsvacolranks[1:15, ]
        metadata(subranks)$ranksNrow <- NULL
        msg <- tryCatch(gsvaColScores(subranks, verbose=FALSE),
                        error=conditionMessage)
        checkTrue(is.character(msg) && grepl("out of range", msg))
        ## a column with ranks out of range
        outofrange <- function(r) {
            nz <- r != 0
            any(r[nz] > sum(nz))
        }
        badcol <- which(apply(as.matrix(assay(subranks, "gsvaranks")), 2,
                              outofrange))[1]
        msg <- tryCatch(gsvaEnrichment(subranks, column=badcol, geneSet=1,
                                       plot="no"),
                        error=conditionMessage)
        checkTrue(is.character(msg) && grepl("out of range", msg))
    }
}
