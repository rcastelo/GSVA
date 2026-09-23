test_avgCcode <- function() {
    ## here we use the acronym AVG to refer to the average method for
    ## calculating enrichment scores
    message("Running unit tests for AVG C code without missing data")

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
    gsets$gs5 <- paste0("g", 1:3)

    ## sample data from a standard normal distribution
    ## seeding the random number generator for reproducibility
    set.seed(123)
    x <- rnorm(p*n)
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))

    ## build AVG parameter object
    avgpar <- avgParam(y, gsets)

    famGaGS <- GSVA:::.filterAndMapGenesAndGeneSets(avgpar)

    filtDataMatrix <- famGaGS[["filteredDataMatrix"]]
    filtMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    ## calculate AVG scores using the R implementation
    sco_R <- GSVA:::.compute_average_scores_R(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)

    ## calculate AVG scores using the C implementation
    sco_C <- GSVA:::.compute_average_scores_block(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)
    dimnames(sco_C) <- dimnames(sco_R)

    ## both approaches to calculate AVG scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    ##
    ## with missing data
    ##
    message("Running unit tests for AVG C code with missing data")

    set.seed(123)
    x[sample(1:(n*p), size=floor(0.05*n*p), replace=FALSE)] <- NA ## 5% missing data
    y <- matrix(x, nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
    ## remove some NAs to check for correct nonNA value propagation
    y[1, is.na(y[1, ])] <- rnorm(sum(is.na(y[1, ])))
    y[2, is.na(y[2, ])] <- rnorm(sum(is.na(y[2, ])))
    y[3, is.na(y[3, ])] <- rnorm(sum(is.na(y[3, ])))

    ## build GSVA parameter object discarding NA values
    avgpar <- avgParam(y, gsets, use="na.rm")

    famGaGS <- GSVA:::.filterAndMapGenesAndGeneSets(avgpar)

    filtDataMatrix <- famGaGS[["filteredDataMatrix"]]
    filtMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    ## calculate AVG scores using the R implementation
    sco_R <- GSVA:::.compute_average_scores_R(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", TRUE, "na.rm",
                                              1, wna_env, FALSE)

    ## calculate AVG scores using the C implementation
    sco_C <- GSVA:::.compute_average_scores_block(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", TRUE, "na.rm",
                                              1, wna_env, FALSE)
    dimnames(sco_C) <- dimnames(sco_R)

    ## both approaches to calculate AVG scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    ## build AVG parameter object propagating NA values
    avgpar <- avgParam(y, gsets, use="everything")

    ## calculate AVG scores using the R implementation
    sco_R <- GSVA:::.compute_average_scores_R(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)

    ## calculate AVG scores using the C implementation
    sco_C <- GSVA:::.compute_average_scores_block(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)
    dimnames(sco_C) <- dimnames(sco_R)

    ## both approaches to calculate AVG scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    message("Running unit tests for GSVA C code with sparse data")

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

    ## build AVG parameter object
    avgpar <- avgParam(y, gsets)

    famGaGS <- GSVA:::.filterAndMapGenesAndGeneSets(avgpar)

    filtDataMatrix <- famGaGS[["filteredDataMatrix"]]
    filtMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    ## calculate AVG scores using the R implementation
    sco_R <- GSVA:::.compute_average_scores_R(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)

    ## calculate AVG scores using the C implementation
    sco_C <- GSVA:::.compute_average_scores_block(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)
    dimnames(sco_C) <- dimnames(sco_R)

    ## both approaches to calculate AVG scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)

    ## check it now with an SVT_SparseMatrix object
    y <- SparseArray(y)

    ## build AVG parameter object
    avgpar <- avgParam(y, gsets)

    famGaGS <- GSVA:::.filterAndMapGenesAndGeneSets(avgpar)

    filtDataMatrix <- famGaGS[["filteredDataMatrix"]]
    filtMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

    wna_env <- new.env()
    assign("w", FALSE, envir=wna_env)

    ## calculate AVG scores using the R implementation
    sco_R <- GSVA:::.compute_average_scores_R(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)

    ## calculate AVG scores using the C implementation
    sco_C <- GSVA:::.compute_average_scores_block(filtDataMatrix,
                                              filtMappedGeneSets,
                                              "mean", FALSE, "everything",
                                              1, wna_env, FALSE)
    dimnames(sco_C) <- dimnames(sco_R)

    ## both approaches to calculate AVG scores must give the same result
    checkEqualsNumeric(sco_R, sco_C)
}
