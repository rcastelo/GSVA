test_gsvaCcode <- function() {
    message("Running unit tests for GSVA C code without missing data.")

    p <- 10 ## number of genes
    n <- 30 ## number of samples
    nGrp1 <- 15 ## number of samples in group 1
    nGrp2 <- n - nGrp1 ## number of samples in group 2

    ## consider three disjoint gene sets
    geneSets <- list(set1=paste("g", 1:3, sep=""),
                     set2=paste("g", 4:6, sep=""),
                     set3=paste("g", 7:10, sep=""))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    x <- rnorm(n*p)
    y <- matrix(x, nrow=p, ncol=n, dimnames=list(paste("g", 1:p, sep=""),
                                                 paste("s", 1:n, sep="")))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(y, geneSets)

    ## calculate GSVA ranks
    gsvarankspar <- gsvaRanks(gsvapar, verbose=FALSE)
    exprData <- GSVA:::get_exprData(gsvarankspar)
    R <- GSVA:::unwrapData(exprData, get_assay(gsvarankspar))
    geneSetsIdx <- GSVA:::.filterAndMapGeneSets(param=gsvarankspar,
                                                filteredDataMatrix=R,
                                                verbose=FALSE)

    rnkstats <- GSVA:::.ranks2stats(R[, 1], sparse=FALSE)

    lapply(as.list(1:ncol(R)), function(j, R) {
        ## calculate GSVA scores using the R implementation
        sco_R <- GSVA:::.gsva_score_genesets_Rimp(geneSetsIdx,
                                                  decOrdStat=rnkstats$dos,
                                                  symRnkStat=rnkstats$srs,
                                                  maxDiff=GSVA:::get_maxDiff(gsvarankspar),
                                                  absRanking=GSVA:::get_absRanking(gsvarankspar),
                                                  tau=GSVA:::get_tau(gsvarankspar),
                                                  any_na=anyNA(gsvarankspar),
                                                  na_use=GSVA:::get_NAuse(gsvarankspar),
                                                  minSize=GSVA:::get_minSize(gsvarankspar))

        ## calculate GSVA scores using the C implementation
        wna_env <- new.env()
        assign("w", FALSE, envir=wna_env)
        sco_C <- GSVA:::.gsva_score_genesets(geneSetsIdx, decOrdStat=rnkstats$dos,
                                             symRnkStat=rnkstats$srs,
                                             maxDiff=GSVA:::get_maxDiff(gsvarankspar),
                                             absRanking=GSVA:::get_absRanking(gsvarankspar),
                                             tau=GSVA:::get_tau(gsvarankspar),
                                             any_na=anyNA(gsvarankspar),
                                             na_use=GSVA:::get_NAuse(gsvarankspar),
                                             minSize=GSVA:::get_minSize(gsvarankspar),
                                             wna_env=wna_env)
        ## both approaches to calculate GSVA scores must give
        ## the same result

        checkEqualsNumeric(sco_R, sco_C)
    }, R=R)

    ## with missing data
    message("Running unit tests for GSVA C code with missing data.")

    x[sample(1:(n*p), size=floor(0.1*n*p), replace=FALSE)] <- NA ## 10% missing data
    y <- matrix(x, nrow=p, ncol=n, dimnames=list(paste("g", 1:p, sep=""),
                                                 paste("s", 1:n, sep="")))

    ## build GSVA parameter object
    gsvapar <- gsvaParam(y, geneSets, use="na.rm")

    ## calculate GSVA ranks
    gsvarankspar <- gsvaRanks(gsvapar, verbose=FALSE)
    exprData <- GSVA:::get_exprData(gsvarankspar)
    R <- GSVA:::unwrapData(exprData, get_assay(gsvarankspar))
    geneSetsIdx <- GSVA:::.filterAndMapGeneSets(param=gsvarankspar,
                                                filteredDataMatrix=R,
                                                verbose=FALSE)

    lapply(as.list(1:ncol(R)), function(j, R) {
        if (anyNA(gsvarankspar))
            rnkstats <- GSVA:::.ranks2stats_nas(R[, j], sparse=FALSE)
        else
            rnkstats <- GSVA:::.ranks2stats(R[, j], sparse=FALSE)

        ## calculate GSVA scores using the R implementation
        sco_R <- GSVA:::.gsva_score_genesets_Rimp(geneSetsIdx,
                                                  decOrdStat=rnkstats$dos,
                                                  symRnkStat=rnkstats$srs,
                                                  maxDiff=GSVA:::get_maxDiff(gsvarankspar),
                                                  absRanking=GSVA:::get_absRanking(gsvarankspar),
                                                  tau=GSVA:::get_tau(gsvarankspar),
                                                  any_na=anyNA(gsvarankspar),
                                                  na_use=GSVA:::get_NAuse(gsvarankspar),
                                                  minSize=GSVA:::get_minSize(gsvarankspar))

        ## calculate GSVA scores using the C implementation
        wna_env <- new.env()
        assign("w", FALSE, envir=wna_env)
        sco_C <- GSVA:::.gsva_score_genesets(geneSetsIdx, decOrdStat=rnkstats$dos,
                                             symRnkStat=rnkstats$srs,
                                             maxDiff=GSVA:::get_maxDiff(gsvarankspar),
                                             absRanking=GSVA:::get_absRanking(gsvarankspar),
                                             tau=GSVA:::get_tau(gsvarankspar),
                                             any_na=anyNA(gsvarankspar),
                                             na_use=GSVA:::get_NAuse(gsvarankspar),
                                             minSize=GSVA:::get_minSize(gsvarankspar),
                                             wna_env=wna_env)
        ## both approaches to calculate GSVA scores must give
        ## the same result

        checkEqualsNumeric(sco_R, sco_C)
    }, R=R)
}
