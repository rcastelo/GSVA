test_spatial <- function() {

    message("Running unit tests for spatial input")

    suppressPackageStartupMessages({
        library(Matrix)
        library(GSEABase)
        library(SpatialExperiment)
        library(GSVAdata)
        library(cli)
        library(org.Hs.eg.db)
    })

    spe <- HumanCerebellumNormSubset()
 
    set.seed(123) ## for reproducibility of the random gene sets
    ## build two gene sets with 4 randomly chosen genes and one third gene set
    ## with a few markers for granule cells, the most abundant cell type in the
    ## cerebellum: SLC17A7, RBFOX3, PAX6, KCND2
    granulecellmarkers <- c("SLC17A7", "RBFOX3", "PAX6", "KCND2")
    granulecellmarkers <- mapIds(org.Hs.eg.db, granulecellmarkers,
				 "ENSEMBL", "SYMBOL")
    ## remove when GSVAdata 1.49.1 is available on the Bioconductor build system
    microgliamarkers <- c("ENSG00000078808", "ENSG00000116251",
                          "ENSG00000142583", "ENSG00000173372")

    gsets <- list(gset1=sample(rownames(spe), size=4, replace=FALSE),
                  gset2=sample(rownames(spe), size=4, replace=FALSE),
                  microglia=microgliamarkers)         ## remove when GSVAdata 1.49.1 is available on the Bioconductor build system
                  ## granulecells=granulecellmarkers) ## and uncomment this line

    ## calculate GSVA enrichment scores and check output
    gsvapar <- gsvaParam(spe, gsets, verbose=TRUE)
    es <- gsva(gsvapar, verbose=TRUE)
    checkTrue(is(es, "SpatialExperiment"))
    checkTrue(all(dim(es) == c(length(gsets), ncol(spe))))
    checkTrue(all(colnames(es) == colnames(spe)))
    checkTrue(is(geneSets(es), "list"))
    out <- cli_fmt(es <- gsva(gsvapar, verbose=FALSE, maxmem="100K"))
    checkTrue(is(assay(es), "DelayedMatrix"))
    checkTrue(grepl("on-disk", out))

    ## comment until GSVAdata 1.49.1 is available on the Bioconductor build system
    ## calculate spatial autocorrelation on the GSVA enrichment scores
    ## r <- spatCor(es, verbose=FALSE)
    ## checkTrue(all(r$observed[r$gene_id == "granulecells"] > r$observed[r$gene_id != "granulecells"]))
}
