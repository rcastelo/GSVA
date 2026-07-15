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
}
