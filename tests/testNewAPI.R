
#
# 2023-08-14  axel: some utilities for developing and testing unit tests
#

# ----- generate test input data -----
generateTestInputData <- function(p = 10, nGS = 3,
                                  n = 30, nGrp1 = round(n / 2), nGrp2 = n - nGrp1) {
    # ff <- formals(generateTestInputData)
    # for(i in seq.int(length(ff))) assign(names(ff[i]), eval(ff[[i]]))

    # we need to be able to reproduce this
    set.seed(2023-08-18)
    tid <- list()

    geneSets <- split(paste0("g", seq.int(p)),
                      cut(seq.int(p), breaks = nGS, labels = paste0("GeneSet", seq.int(nGS))))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p) , paste0("s", 1:n)))
    ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
    y[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] <- y[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] + 2
    tid[["matrix_continuous"]] <- y

    # same structure with count data
    d <- matrix(rpois(n*p, 42), nrow=p, ncol=n,
                dimnames=list(paste0("g", 1:p) , paste0("s", 1:n)))
    d[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] <- d[geneSets[["GeneSet1"]], seq.int(nGrp1+1, n)] + 23
    tid[["matrix_counts"]] <- d

    return(tid)
}


# ---- generate and save test input data if necessary ----
# gti <- file.path(BiocGenerics:::packageRoot(getwd()), "tests", "GSVA_test_input.RData")
gti <- file.path(BiocGenerics:::packageRoot("/home/axel/dev/R/pkg/GSVA"), "tests", "GSVA_test_input.RData")


if(!file.exists(gti)) {
    testInputData <- generateTestInputData()
    save(testInputData, file = gti)
}

