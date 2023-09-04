
###
### 2023-09-03  axele
### some useful information about package tests found in WRE 1.1.5, as well as
### WRE 1.3.1, item 19, WRE 2.1.1, section \examples{...}
###

## ----setup, include=FALSE---------------------------------------------------------------
set.seed(2023-09-03)
options(warn = 1)
## options(width=80)
## knitr::opts_chunk$set(collapse=TRUE,
##                       message=FALSE,
##                       comment="")


## ----library_install, message=FALSE, cache=FALSE, eval=FALSE----------------------------
## install.packages("BiocManager")
## BiocManager::install("GSVA")


## ----load_library, message=FALSE, warning=FALSE, cache=FALSE----------------------------
suppressPackageStartupMessages(library(GSVA))


## ---------------------------------------------------------------------------------------
p <- 10000 ## number of genes
n <- 30    ## number of samples
## simulate expression values from a standard Gaussian distribution
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
X[1:5, 1:5]


## ---------------------------------------------------------------------------------------
## sample gene set sizes
gs <- as.list(sample(10:100, size=100, replace=TRUE))
## sample gene sets
gs <- lapply(gs, function(n, p)
                   paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))


## ---------------------------------------------------------------------------------------
gsva.es <- gsva(X, gs, verbose=FALSE)
dim(gsva.es)
gsva.es[1:5, 1:5]


## ---------------------------------------------------------------------------------------
class(gs)
length(gs)
head(lapply(gs, head))


## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(org.Hs.eg.db))

goannot <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
head(goannot)
genesbygo <- split(goannot$ENTREZID, goannot$GO)
length(genesbygo)
head(genesbygo)


## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(GSVAdata))

data(c2BroadSets)
class(c2BroadSets)
c2BroadSets


## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(Biobase))

data(commonPickrellHuang)

stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset),
                    featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))


## ---------------------------------------------------------------------------------------
canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
canonicalC2BroadSets


## ---------------------------------------------------------------------------------------
data(genderGenesEntrez)

MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),
                              collectionType=BroadCollection(category="c2"), setName="MSY")
MSY
XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),
                              collectionType=BroadCollection(category="c2"), setName="XiE")
XiE

canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))
canonicalC2BroadSets


## ---- results="hide"--------------------------------------------------------------------
esmicro <- gsva(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500)
esrnaseq <- gsva(pickrellCountsArgonneCQNcommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500,
                 kcdf="Poisson")


## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(edgeR))

lcpms <- cpm(exprs(pickrellCountsArgonneCQNcommon_eset), log=TRUE)


## ---------------------------------------------------------------------------------------
genecorrs <- sapply(1:nrow(lcpms),
                    function(i, expmicro, exprnaseq) cor(expmicro[i, ], exprnaseq[i, ],
                                                         method="spearman"),
                    exprs(huangArrayRMAnoBatchCommon_eset), lcpms)
names(genecorrs) <- rownames(lcpms)


## ---------------------------------------------------------------------------------------
pwycorrs <- sapply(1:nrow(esmicro),
                   function(i, esmicro, esrnaseq) cor(esmicro[i, ], esrnaseq[i, ],
                                                      method="spearman"),
                   exprs(esmicro), exprs(esrnaseq))
names(pwycorrs) <- rownames(esmicro)


## ----compcorrs, height=500, width=1000, fig.cap="Comparison of correlation values of gene and pathway expression profiles derived from microarray and RNA-seq data."----
## par(mfrow=c(1, 2), mar=c(4, 5, 3, 2))
## hist(genecorrs, xlab="Spearman correlation", main="Gene level\n(RNA-seq log-CPMs vs microarray RMA)",
##      xlim=c(-1, 1), col="grey", las=1)
## hist(pwycorrs, xlab="Spearman correlation", main="Pathway level\n(GSVA enrichment scores)",
##      xlim=c(-1, 1), col="grey", las=1)


## ----compsexgenesets, height=500, width=1000, fig.cap="Comparison of GSVA enrichment scores obtained from microarray and RNA-seq data for two gene sets formed by genes with sex-specific expression."----
rmsy <- cor(exprs(esrnaseq)["MSY", ], exprs(esmicro)["MSY", ])
## par(mfrow=c(1, 2))
## plot(exprs(esrnaseq)["MSY", ], exprs(esmicro)["MSY", ], xlab="GSVA scores RNA-seq",
##      ylab="GSVA scores microarray", main=sprintf("MSY R=%.2f", rmsy), las=1, type="n")
## abline(lm(exprs(esmicro)["MSY", ] ~ exprs(esrnaseq)["MSY", ]), lwd=2, lty=2, col="grey")
## points(exprs(esrnaseq["MSY", pickrellCountsArgonneCQNcommon_eset$Gender == "Female"]),
##        exprs(esmicro)["MSY", huangArrayRMAnoBatchCommon_eset$Gender == "Female"],
##        col="red", pch=21, bg="red", cex=1)
## points(exprs(esrnaseq)["MSY", pickrellCountsArgonneCQNcommon_eset$Gender == "Male"],
##        exprs(esmicro)["MSY", huangArrayRMAnoBatchCommon_eset$Gender == "Male"],
##        col="blue", pch=21, bg="blue", cex=1)
## legend("topleft", c("female", "male"), pch=21, col=c("red", "blue"), pt.bg=c("red", "blue"), inset=0.01)
## rxie <- cor(exprs(esrnaseq)["XiE", ], exprs(esmicro)["XiE", ])
## plot(exprs(esrnaseq)["XiE", ], exprs(esmicro)["XiE", ], xlab="GSVA scores RNA-seq",
##      ylab="GSVA scores microarray", main=sprintf("XiE R=%.2f", rxie), las=1, type="n")
## abline(lm(exprs(esmicro)["XiE", ] ~ exprs(esrnaseq)["XiE", ]), lwd=2, lty=2, col="grey")
## points(exprs(esrnaseq["XiE", pickrellCountsArgonneCQNcommon_eset$Gender == "Female"]),
##        exprs(esmicro)["XiE", huangArrayRMAnoBatchCommon_eset$Gender == "Female"],
##        col="red", pch=21, bg="red", cex=1)
## points(exprs(esrnaseq)["XiE", pickrellCountsArgonneCQNcommon_eset$Gender == "Male"],
##        exprs(esmicro)["XiE", huangArrayRMAnoBatchCommon_eset$Gender == "Male"],
##        col="blue", pch=21, bg="blue", cex=1)
## legend("topleft", c("female", "male"), pch=21, col=c("red", "blue"), pt.bg=c("red", "blue"), inset=0.01)


## ---------------------------------------------------------------------------------------
data(gbm_VerhaakEtAl)
gbm_eset
head(featureNames(gbm_eset))
table(gbm_eset$subtype)
data(brainTxDbSets)
lengths(brainTxDbSets)
lapply(brainTxDbSets, head)


## ---- results="hide"--------------------------------------------------------------------
gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff=FALSE)


## ----gbmSignature, height=500, width=700, fig.cap="Heatmap of GSVA scores for cell-type brain signatures from murine models (y-axis) across GBM samples grouped by GBM subtype."----
suppressPackageStartupMessages(library(RColorBrewer))
subtypeOrder <- c("Proneural", "Neural", "Classical", "Mesenchymal")
sampleOrderBySubtype <- sort(match(gbm_es$subtype, subtypeOrder),
                             index.return=TRUE)$ix
subtypeXtable <- table(gbm_es$subtype)
subtypeColorLegend <- c(Proneural="red", Neural="green",
                        Classical="blue", Mesenchymal="orange")
geneSetOrder <- c("astroglia_up", "astrocytic_up", "neuronal_up",
                  "oligodendrocytic_up")
geneSetLabels <- gsub("_", " ", geneSetOrder)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]

## heatmap(exprs(gbm_es)[geneSetOrder, sampleOrderBySubtype], Rowv=NA,
##         Colv=NA, scale="row", margins=c(3,5), col=hmcol,
##         ColSideColors=rep(subtypeColorLegend[subtypeOrder],
##                           times=subtypeXtable[subtypeOrder]),
##         labCol="", gbm_es$subtype[sampleOrderBySubtype],
##         labRow=paste(toupper(substring(geneSetLabels, 1,1)),
##                      substring(geneSetLabels, 2), sep=""),
##         cexRow=2, main=" \n ")
## par(xpd=TRUE)
## text(0.23,1.21, "Proneural", col="red", cex=1.2)
## text(0.36,1.21, "Neural", col="green", cex=1.2)
## text(0.47,1.21, "Classical", col="blue", cex=1.2)
## text(0.62,1.21, "Mesenchymal", col="orange", cex=1.2)
## mtext("Gene sets", side=4, line=0, cex=1.5)
## mtext("Samples          ", side=1, line=4, cex=1.5)


## ---------------------------------------------------------------------------------------
data(leukemia)
leukemia_eset


## ---- results="hide"--------------------------------------------------------------------
leukemia_es <- gsva(leukemia_eset, c2BroadSets, min.sz=10, max.sz=500)


## ---------------------------------------------------------------------------------------
class(leukemia_es)
leukemia_es
head(featureNames(leukemia_es))


## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(limma))

mod <- model.matrix(~ factor(leukemia_es$subtype))
colnames(mod) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)
summary(res)


## ----leukemiavolcano, height=700, width=500, fig.cap="Volcano plot for the differential expression analysis at pathway level between two leukemia subtypes."----
tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]
## plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
##      main="", xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~Raw~P-value))
## abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), col=grey(0.5), lwd=1, lty=2)
## points(tt$logFC[match(DEpwys, rownames(tt))],
##        -log10(tt$P.Value[match(DEpwys, rownames(tt))]), pch=".", cex=5, col="darkred")
## text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), "1% FDR", pos=3)


## ----leukemiaheatmap, height=500, width=1200, fig.cap="Heatmap of GSVA enrichment scores for the differentially expressed pathways between two leukemia subtypes."----
DEpwys_es <- exprs(leukemia_es[DEpwys, ])
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("ALL", "MLL")
sample.color.map <- colorLegend[pData(leukemia_es)[, "subtype"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")
## heatmap(DEpwys_es, ColSideColors=sample.color.map, xlab="samples",
##         ylab="Pathways", margins=c(2, 20),
##         labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "",
##                                           rownames(DEpwys_es))), 1, 35),
##         labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
##         Rowv=as.dendrogram(geneSetClustering))
## legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")


## ---- eval=FALSE------------------------------------------------------------------------
## res <- igsva()


## ----session_info, cache=FALSE----------------------------------------------------------
## sessionInfo()

