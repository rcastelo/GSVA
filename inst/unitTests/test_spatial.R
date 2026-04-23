test_spatial <- function() {

    message("Running unit tests for spatial input")

    suppressPackageStartupMessages({
        library(Matrix)
        library(GSEABase)
        library(SpatialExperiment)
        library(GSVAdata)
        library(cli)
    })

    spe <- HumanCerebellumNormSubset()
    gsvaAnnotation(spe) <- ENSEMBLIdentifier("org.Hs.eg.db")
 
    set.seed(123) ## for reproducibility of the random gene sets
    ## build two gene sets with 4 randomly chosen genes and one
    ## third gene set with a few microglia marker genes
    gsets <- list(gset1=sample(rownames(spe), size=4, replace=FALSE),
                  gset2=sample(rownames(spe), size=4, replace=FALSE),
                  microglia=c("ENSG00000078808", "ENSG00000116251",
                              "ENSG00000142583", "ENSG00000173372"))

    ## calculate GSVA enrichment scores and check output
    gsvapar <- gsvaParam(spe, gsets, verbose=FALSE)
    es <- gsva(gsvapar, verbose=FALSE)
    checkTrue(is(es, "SpatialExperiment"))
    checkTrue(all(dim(es) == c(length(gsets), ncol(spe))))
    checkTrue(all(colnames(es) == colnames(spe)))
    checkTrue(is(geneSets(es), "list"))
    out <- cli_fmt(es <- gsva(gsvapar, verbose=FALSE, maxmem="100K"))
    checkTrue(is(assay(es), "DelayedMatrix"))
    checkTrue(grepl("on-disk", out))

    ## calculate spatial autocorrelation on the GSVA enrichment scores
    r <- spatCor(es, verbose=FALSE)
    checkTrue(all(r$observed[r$gene_id == "microglia"] > r$observed[r$gene_id != "microglia"]))
}
