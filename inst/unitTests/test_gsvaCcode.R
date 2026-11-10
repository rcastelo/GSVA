test_gsvaCcode <- function() {
    message("Running unit tests for GSVA C code without missing data.")

    ## geneids: gene identifiers
    ## k: number of genesets
    ## m: size of gene sets
    sampleGeneSets <- function(geneids, k, m) {
        stopifnot(m < length(geneids)) ## QC

        gsets <- replicate(k, sample(geneids, m, replace=FALSE), simplify=FALSE)
        names(gsets) <- paste0("gs", seq_len(k))

        gsets
    }

    p <- 10 ## number of genes
    n <- 30 ## number of samples
    k <- 5  ## number of gene sets
    m <- 3  ## size of every gene set

    ## sample k-1 gene sets at random using a seed for reproducibility
    set.seed(123)
    gsets <- sampleGeneSets(paste0("g", 1:p), k=k-1, m=m)
    ## a fifth gene set is planted to check for nonNA propagation below
    gsets[[5]] <- paste0("g", 1:3)

    ## sample data from a standard normal distribution
    ## seeding the random number generator for reproducibility
    set.seed(123)
    x <- rnorm(p*n)
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(y, gsets)

    ## calculate GSVA ranks
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    param <- GSVA:::.pull_param(gsvacolranks)
    exprData <- GSVA:::get_exprData(param)
    R <- GSVA:::unwrapData(exprData, get_assay(param))
    gsetsidx <- GSVA:::.filterAndMapGeneSets(param=param,
                                             filteredDataMatrix=R,
                                             verbose=FALSE)

    sco_R <- lapply(as.list(1:ncol(R)), function(j, R) {
        rnkstats <- GSVA:::.ranks2stats(R[, j], sparse=FALSE)
        ## calculate GSVA scores using the R implementation
        GSVA:::.gsva_score_genesets_Rimp(gsetsidx,
                                         decOrdStat=rnkstats$dos,
                                         symRnkStat=rnkstats$srs,
                                         whz=rnkstats$whz,
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param))
    }, R=R)
    sco_R <- do.call("cbind", sco_R)

    ## calculate GSVA scores using the C implementation
    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)
    sco_C <- GSVA:::.gsva_score_genesets(R, gsetsidx, is.integer(R[1, 1]),
					                     sparse=GSVA:::.get_sparse(param), 
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param),
                                         wna_env=wna_env, verbose=FALSE)

    ## both approaches to calculate GSVA scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    ## with missing data
    message("Running unit tests for GSVA C code with missing data.")

    set.seed(123)
    x[sample(1:(n*p), size=floor(0.05*n*p), replace=FALSE)] <- NA ## 5% missing data
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    ## remove some NAs to check for correct nonNA value propagation
    y[1, is.na(y[1, ])] <- rnorm(sum(is.na(y[1, ])))
    y[2, is.na(y[2, ])] <- rnorm(sum(is.na(y[2, ])))
    y[3, is.na(y[3, ])] <- rnorm(sum(is.na(y[3, ])))

    ## build GSVA parameter object discarding NA values
    gsvapar <- gsvaParam(y, gsets, use="na.rm")

    ## calculate GSVA ranks
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    param <- GSVA:::.pull_param(gsvacolranks)
    exprData <- GSVA:::get_exprData(param)
    R <- GSVA:::unwrapData(exprData, get_assay(param))
    gsetsidx <- GSVA:::.filterAndMapGeneSets(param=param,
                                             filteredDataMatrix=R,
                                             verbose=FALSE)

    sco_R <- lapply(as.list(1:ncol(R)), function(j, R) {
        if (anyNA(param))
            rnkstats <- GSVA:::.ranks2stats_nas(R[, j], sparse=FALSE)
        else
            rnkstats <- GSVA:::.ranks2stats(R[, j], sparse=FALSE)

        ## calculate GSVA scores using the R implementation
        GSVA:::.gsva_score_genesets_Rimp(gsetsidx,
                                         decOrdStat=rnkstats$dos,
                                         symRnkStat=rnkstats$srs,
                                         whz=rnkstats$whz,
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param))
    }, R=R)
    sco_R <- do.call("cbind", sco_R)

    ## calculate GSVA scores using the C implementation
    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)
    sco_C <- GSVA:::.gsva_score_genesets(R, gsetsidx, is.integer(R[1, 1]),
                                         sparse=GSVA:::.get_sparse(param), 
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param),
                                         wna_env=wna_env, verbose=FALSE)

    ## both approaches to calculate GSVA scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    ## build GSVA parameter object propagating NA values
    gsvapar <- gsvaParam(y, gsets, use="everything")

    ## calculate GSVA ranks
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    param <- GSVA:::.pull_param(gsvacolranks)
    exprData <- GSVA:::get_exprData(param)
    R <- GSVA:::unwrapData(exprData, get_assay(param))
    gsetsidx <- GSVA:::.filterAndMapGeneSets(param=param,
                                             filteredDataMatrix=R,
                                             verbose=FALSE)

    sco_R <- lapply(as.list(1:ncol(R)), function(j, R) {
        if (anyNA(param))
            rnkstats <- GSVA:::.ranks2stats_nas(R[, j], sparse=FALSE)
        else
            rnkstats <- GSVA:::.ranks2stats(R[, j], sparse=FALSE)

        ## calculate GSVA scores using the R implementation
        GSVA:::.gsva_score_genesets_Rimp(gsetsidx,
                                         decOrdStat=rnkstats$dos,
                                         symRnkStat=rnkstats$srs,
                                         whz=rnkstats$whz,
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param))
    }, R=R)
    sco_R <- do.call("cbind", sco_R)

    ## calculate GSVA scores using the C implementation
    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)
    assign("u", 0, envir=wna_env)
    sco_C <- GSVA:::.gsva_score_genesets(R, gsetsidx, is.integer(R[1, 1]),
                                         sparse=GSVA:::.get_sparse(param), 
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param),
                                         wna_env=wna_env, verbose=FALSE)

    ## both approaches to calculate GSVA scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    message("Running unit tests for GSVA C code with sparse data.")

    ## check now the C code on sparse input data
    library(Matrix)
    library(SparseArray)

    p <- 100 ## number of genes
    n <- 50  ## number of samples
    k <- 15  ## number of gene sets
    m <- 5   ## size of gene sets

    ## sample k gene sets at random using a seed for reproducibility
    set.seed(123)
    gsets <- sampleGeneSets(paste0("g", 1:p), k=k, m=m)

    ## build a random sparse matrix with 85% sparsity
    set.seed(123)
    x <- numeric(p*n)
    idx <- sample(length(x), size=round(length(x)*0.15)) ## 85% sparsity
    x[idx] <- rnorm(length(idx))
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p) , paste0("s", 1:n)))
    y <- Matrix(y, sparse=TRUE)

    ## build parameter object
    gsvapar <- gsvaParam(y, gsets, verbose=FALSE)

    ## calculate GSVA ranks
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    param <- GSVA:::.pull_param(gsvacolranks)
    exprData <- GSVA:::get_exprData(param)
    R <- GSVA:::unwrapData(exprData, get_assay(param))
    gsetsidx <- GSVA:::.filterAndMapGeneSets(param=param,
                                             filteredDataMatrix=R,
                                             verbose=FALSE)

    sco_R <- lapply(as.list(1:ncol(R)), function(j, R) {
        rnkstats <- GSVA:::.ranks2stats(R[, j], sparse=GSVA:::.get_sparse(param))

        ## calculate GSVA scores using the R implementation
        GSVA:::.gsva_score_genesets_Rimp(gsetsidx,
                                         decOrdStat=rnkstats$dos,
                                         symRnkStat=rnkstats$srs,
                                         whz=rnkstats$whz,
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param))
    }, R=R)
    sco_R <- do.call("cbind", sco_R)

    ## calculate GSVA scores using the C implementation
    sco_C <- GSVA:::.gsva_score_genesets(R, gsetsidx, is.integer(R[1, 1]),
                                         sparse=GSVA:::.get_sparse(param), 
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param),
                                         wna_env=wna_env, verbose=FALSE)

    ## both approaches to calculate GSVA scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    ## check now with an SVT_SparseMatrix object
    y <- SparseArray(y)

    ## build parameter object
    gsvapar <- gsvaParam(y, gsets, verbose=FALSE)

    ## calculate GSVA ranks
    gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
    gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
    param <- GSVA:::.pull_param(gsvacolranks)
    exprData <- GSVA:::get_exprData(param)
    R <- GSVA:::unwrapData(exprData, get_assay(param))
    gsetsidx <- GSVA:::.filterAndMapGeneSets(param=param,
                                             filteredDataMatrix=R,
                                             verbose=FALSE)

    sco_R <- lapply(as.list(1:ncol(R)), function(j, R) {
        rnkstats <- GSVA:::.ranks2stats(R[, j], sparse=GSVA:::.get_sparse(param))

        ## calculate GSVA scores using the R implementation
        GSVA:::.gsva_score_genesets_Rimp(gsetsidx,
                                         decOrdStat=rnkstats$dos,
                                         symRnkStat=rnkstats$srs,
                                         whz=rnkstats$whz,
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param))
    }, R=R)
    sco_R <- do.call("cbind", sco_R)

    ## calculate GSVA scores using the C implementation
    sco_C <- GSVA:::.gsva_score_genesets(R, gsetsidx, is.integer(R[1, 1]),
                                         sparse=GSVA:::.get_sparse(param), 
                                         maxDiff=GSVA:::.get_maxDiff(param),
                                         absRanking=GSVA:::.get_absRanking(param),
                                         tau=GSVA:::.get_tau(param),
                                         any_na=anyNA(param),
                                         na_use=GSVA:::.get_NAuse(param),
                                         minSize=GSVA:::get_minSize(param),
                                         wna_env=wna_env, verbose=FALSE)

    ## both approaches to calculate GSVA scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)
}

test_gsvaCcode_radix <- function() {
    message("Running unit tests for GSVA C code with large gene sets.")

    ## the C code sorts the positions of the genes in a gene set with radix
    ## sort when there are 20 or more of them (20 or more nonzero ones in the
    ## sparse random walk), and with two radix passes when p > 256 genes
    library(Matrix)

    ## calculates GSVA scores with the R and the C implementations
    scoreRandC <- function(y, gsets, use="everything") {
        gsvapar <- gsvaParam(y, gsets, minSize=1, use=use)
        gsvarownorm <- gsvaRowNorm(gsvapar, verbose=FALSE)
        gsvacolranks <- gsvaColRanks(gsvarownorm, verbose=FALSE)
        param <- GSVA:::.pull_param(gsvacolranks)
        exprData <- GSVA:::get_exprData(param)
        R <- GSVA:::unwrapData(exprData, get_assay(param))
        gsetsidx <- GSVA:::.filterAndMapGeneSets(param=param,
                                                 filteredDataMatrix=R,
                                                 verbose=FALSE)
        sparse <- GSVA:::.get_sparse(param)

        sco_R <- lapply(as.list(1:ncol(R)), function(j, R) {
            if (anyNA(param))
                rnkstats <- GSVA:::.ranks2stats_nas(R[, j], sparse=FALSE)
            else
                rnkstats <- GSVA:::.ranks2stats(R[, j], sparse=sparse)

            GSVA:::.gsva_score_genesets_Rimp(gsetsidx,
                                             decOrdStat=rnkstats$dos,
                                             symRnkStat=rnkstats$srs,
                                             whz=rnkstats$whz,
                                             maxDiff=GSVA:::.get_maxDiff(param),
                                             absRanking=GSVA:::.get_absRanking(param),
                                             tau=GSVA:::.get_tau(param),
                                             any_na=anyNA(param),
                                             na_use=GSVA:::.get_NAuse(param),
                                             minSize=GSVA:::get_minSize(param))
        }, R=R)
        sco_R <- do.call("cbind", sco_R)

        wna_env <- new.env()
        assign("w", FALSE, envir=wna_env)
        assign("u", 0, envir=wna_env)
        sco_C <- GSVA:::.gsva_score_genesets(R, gsetsidx, is.integer(R[1, 1]),
                                             sparse=sparse,
                                             maxDiff=GSVA:::.get_maxDiff(param),
                                             absRanking=GSVA:::.get_absRanking(param),
                                             tau=GSVA:::.get_tau(param),
                                             any_na=anyNA(param),
                                             na_use=GSVA:::.get_NAuse(param),
                                             minSize=GSVA:::get_minSize(param),
                                             wna_env=wna_env, verbose=FALSE)

        list(R=sco_R, C=sco_C)
    }

    p <- 1000 ## number of genes
    n <- 8    ## number of samples
    genes <- paste0("g", 1:p)

    ## gene sets with sizes below, at and above the radix sort threshold
    set.seed(123)
    gsets <- lapply(c(19, 20, 21, 50, 300, p-1), function(s) sample(genes, s))
    names(gsets) <- paste0("gs", seq_along(gsets))

    ## dense data
    set.seed(123)
    x <- rnorm(p*n)
    y <- matrix(x, nrow=p, ncol=n, dimnames=list(genes, paste0("s", 1:n)))
    sco <- scoreRandC(y, gsets)
    checkEqualsNumeric(sco$R, sco$C)

    ## dense data with missing values
    set.seed(123)
    y[sample(p*n, size=floor(0.05*p*n))] <- NA ## 5% missing data
    sco <- scoreRandC(y, gsets, use="na.rm")
    checkEqualsNumeric(sco$R, sco$C)

    ## sparse data with 50% zeros, where gene sets with 50 or more genes have
    ## at least 20 nonzero genes in every column
    set.seed(123)
    x <- numeric(p*n)
    idx <- sample(length(x), size=round(length(x)*0.5))
    x[idx] <- rexp(length(idx))
    y <- matrix(x, nrow=p, ncol=n, dimnames=list(genes, paste0("s", 1:n)))
    y <- Matrix(y, sparse=TRUE)
    sco <- scoreRandC(y, gsets)
    checkEqualsNumeric(sco$R, sco$C)

    ## sparse data with few genes, where redistributing the zero genes of
    ## small gene sets among all zeros leads to positions ending in .5, which
    ## the R and the C code must round in the same way
    p <- 60
    n <- 20
    genes <- paste0("g", 1:p)
    set.seed(123)
    x <- numeric(p*n)
    idx <- sample(length(x), size=round(length(x)*0.5))
    x[idx] <- rexp(length(idx))
    y <- matrix(x, nrow=p, ncol=n, dimnames=list(genes, paste0("s", 1:n)))
    y <- Matrix(y, sparse=TRUE)
    set.seed(123)
    gsets <- lapply(rep(3:8, length.out=40), function(s) sample(genes, s))
    names(gsets) <- paste0("gs", seq_along(gsets))
    sco <- scoreRandC(y, gsets)
    checkEqualsNumeric(sco$R, sco$C)
}

test_gsvaCcodeUndefWalkOverflow <- function() {
    message("Running unit tests for random walks with overflowing steps.")

    ## the walk accumulates the symmetric rank statistics raised to tau.
    ## when such exponent is large, the sum might overflow to infinity, and
    ## the normalized walk would silently collapse to zeros and NaNs.
    ## the C entry point is called directly on hand-crafted rank matrices,
    ## exposing this issue.
    tau <- 200
    p <- 200

    ## dense columns with controlled ranks
    Rd <- cbind(1:p, p:1)
    storage.mode(Rd) <- "integer"
    gsets <- list(ovf=1:100, ctl=96:105, bnd=68:77)
    sco <- .Call("gsva_score_genesets_R", Rd, gsets, TRUE, FALSE, TRUE,
                 FALSE, as.double(tau), FALSE, 1L, 1L, FALSE,
                 PACKAGE="GSVA")
    checkTrue(all(is.na(sco[1, ])))                ## overflow -> undefined walks
    checkTrue(all(is.finite(sco[2, ])))            ## small weights are unaffected
    checkTrue(all(is.finite(sco[3, ])))            ## large finite sums stay defined
    checkIdentical(attr(sco, "attrUndefWalks"), 2) ## one per overflowing column

    ## without overflowing gene sets no walk is reported as undefined
    sco <- .Call("gsva_score_genesets_R", Rd, list(ctl=96:105), TRUE, FALSE,
                 TRUE, FALSE, as.double(tau), FALSE, 1L, 1L, FALSE,
                 PACKAGE="GSVA")
    checkTrue(is.null(attr(sco, "attrUndefWalks")))

    ## sparse column with 150 nonzero genes and 50 zero genes
    library(Matrix)
    Rs <- new("dgCMatrix", i=as.integer(0:149), p=c(0L, 150L),
              x=as.double(1:150), Dim=c(200L, 1L))
    gsets <- list(ovf=151:200, ctl=61:90)
    sco <- .Call("gsva_score_genesets_R", Rs, gsets, FALSE, TRUE, TRUE,
                 FALSE, as.double(tau), FALSE, 1L, 1L, FALSE,
                 PACKAGE="GSVA")
    checkTrue(is.na(sco[1, 1]))
    checkTrue(is.finite(sco[2, 1]))
    checkIdentical(attr(sco, "attrUndefWalks"), 1)

    ## column with 20 missing genes
    Rn <- matrix(as.double(1:200), ncol=1)
    Rn[181:200, 1] <- NA
    gsets <- list(ovf=1:60, ctl=81:100)
    sco <- .Call("gsva_score_genesets_R", Rn, gsets, FALSE, FALSE, TRUE,
                 FALSE, as.double(tau), TRUE, 3L, 1L, FALSE,
                 PACKAGE="GSVA")
    checkTrue(is.na(sco[1, 1]))
    checkTrue(is.finite(sco[2, 1]))
    checkIdentical(attr(sco, "attrUndefWalks"), 1)
}
