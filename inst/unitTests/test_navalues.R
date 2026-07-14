test_navalues <- function() {
    message("Running unit tests for handling of NA values")

    suppressPackageStartupMessages({
        library(cli) ## for cli_fmt()
    })

    p <- 10 ## number of genes
    n <- 30 ## number of samples

    ## consider three disjoint gene sets
    gsets <- list(set1=paste("g", 1:3, sep=""),
                  set2=paste("g", 4:6, sep=""),
                  set3=paste("g", 7:10, sep=""))

    ## sample data from a normal distribution with mean 0 and st.dev. 1
    ## seeding the random number generator for the purpose of this test
    set.seed(123)
    y <- matrix(rnorm(n*p), nrow=p, ncol=n,
                dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

    ## set some values to NA
    y[1, 1] <- y[3, 1] <- NA
    y[2, 3] <- y[4, 3] <- NA

    ## calculate GSVA scores with default use="everything"
    out <- cli_fmt(gsvapar <- gsvaParam(y, gsets))
    checkTrue(grepl("Input expression data has NA values, which will be propagated through calculations", out))
    es <- gsva(gsvapar, verbose=FALSE)

    ## calculate GSVA scores with use="na.rm"
    out <- cli_fmt(gsvapar <- gsvaParam(y, gsets, use="na.rm"))
    checkTrue(grepl("Input expression data has NA values, which will be discarded from calculations", out))
    es <- gsva(gsvapar, verbose=FALSE)
    out <- cli_fmt(gsvapar <- gsvaParam(y, gsets, kcdf="none", use="na.rm"))
    checkTrue(grepl("Input expression data has NA values, which will be discarded from calculations", out))
    es <- gsva(gsvapar, verbose=FALSE)

    ## calculate GSVA scores with use="all.obs", which
    ## should prompt an error when building the parameter object
    checkException(gsvapar <- gsvaParam(y, gsets, use="all.obs"))

    ## calculate ssGSEA scores with default use="everything"
    out <- cli_fmt(ssgseapar <- ssgseaParam(y, gsets))
    checkTrue(grepl("Input expression data has NA values, which will be propagated through calculations", out))
    es <- gsva(ssgseapar, verbose=FALSE)

    ## calculate ssGSEA scores with use="na.rm"
    out <- cli_fmt(ssgseapar <- ssgseaParam(y, gsets, use="na.rm"))
    checkTrue(grepl("Input expression data has NA values, which will be discarded from calculations", out))
    es <- gsva(ssgseapar, verbose=FALSE)

    ## calculate ssGSEA scores with use="all.obs", which
    ## should prompt an error when building the parameter object
    checkException(ssgseapar <- ssgseaParam(y, gsets, use="all.obs"))
}
