
## 2025-11-12 regression test script for method PLAGE

suppressPackageStartupMessages(library(GSVAdata))
suppressPackageStartupMessages(library(GSVA))

data(c2BroadSets)
data(geneprotExpCostaEtAl2021)

set.seed(2025-11-12)
options(width=1024)
nRowsToPrint <- 25

p <- 10000 ## number of genes
n <- 30    ## number of samples
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
X[1:5, 1:5]
gs <- as.list(sample(10:100, size=100, replace=TRUE))
gs <- lapply(gs, function(n, p)
                   paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))
plagePar <- plageParam(X, gs)
plage.es <- gsva(plagePar, verbose=FALSE)
dim(plage.es)
plage.es[seq.int(min(nRowsToPrint, nrow(plage.es))),]

c2BroadSets <- c2BroadSets[c(grep("_UP$", names(c2BroadSets)),
                             grep("_DN$", names(c2BroadSets)))]
firPar <- plageParam(geneExpCostaEtAl2021, c2BroadSets,
                     minSize=10, maxSize=500, verbose=FALSE)
fir_es <- gsva(firPar, verbose=FALSE)
assay(fir_es)[seq.int(min(nRowsToPrint, nrow(fir_es))),]
