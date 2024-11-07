
## 2024-10-25 regression test script for method GSVA

suppressPackageStartupMessages(library(GSVAdata))
suppressPackageStartupMessages(library(GSVA))

data(c2BroadSets)
data(commonPickrellHuang)
data(genderGenesEntrez)
data(gbm_VerhaakEtAl)
data(brainTxDbSets)
data(leukemia)

set.seed(2024-10-25)
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
gsvaPar <- gsvaParam(X, gs)
gsva.es <- gsva(gsvaPar, verbose=FALSE)
dim(gsva.es)
gsva.es[seq.int(min(nRowsToPrint, nrow(gsva.es))),]

stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset),
                    featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))
## until the updated GSVAdata goes through the build system
## remove duplicated rows
fnames <- featureNames(huangArrayRMAnoBatchCommon_eset)
mask <- duplicated(fnames)
huangArrayRMAnoBatchCommon_eset <- huangArrayRMAnoBatchCommon_eset[!mask, ]
pickrellCountsArgonneCQNcommon_eset <- pickrellCountsArgonneCQNcommon_eset[!mask, ]

canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"),
               setName="MSY")
XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"),
               setName="XiE")
canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))
huangPar <- gsvaParam(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets,
                      minSize=5, maxSize=500)
esmicro <- gsva(huangPar, verbose=FALSE)
exprs(esmicro)[seq.int(min(nRowsToPrint, nrow(esmicro))),]
pickrellPar <- gsvaParam(pickrellCountsArgonneCQNcommon_eset,
                         canonicalC2BroadSets, minSize=5, maxSize=500,
                         kcdf="Poisson")
esrnaseq <- gsva(pickrellPar, verbose=FALSE)
exprs(esrnaseq)[seq.int(min(nRowsToPrint, nrow(esrnaseq))),]

gbmPar <- gsvaParam(gbm_eset, brainTxDbSets, maxDiff=FALSE)
gbm_es <- gsva(gbmPar, verbose=FALSE)
exprs(gbm_es)[seq.int(min(nRowsToPrint, nrow(gbm_es))),]

cgpC2BroadSets <- c2BroadSets[c(grep("_UP$", names(c2BroadSets)),
                                grep("_DN$", names(c2BroadSets)))]
leukPar <- gsvaParam(leukemia_eset, cgpC2BroadSets,
                     minSize=10, maxSize=500)
leukemia_es <- gsva(leukPar, verbose=FALSE)
exprs(leukemia_es)[seq.int(min(nRowsToPrint, nrow(leukemia_es))),]
