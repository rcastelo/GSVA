test_spatial <- function() {
    message("Running unit tests for spatial input")

    suppressPackageStartupMessages({
        library(Matrix)
        library(SpatialExperiment)
    })

    ## build a SpatialExperiment object
    syspath <- system.file("extdata", package="GSVA")
    fname <- "human_cerebellum_norm_logcounts_250x4816.mtx.gz"
    logcounts <- as(readMM(gzfile(file.path(syspath, fname))), "CsparseMatrix")
    fname <- "human_cerebellum_rowdata_250x4816.csv.gz"
    rowdata <- read.csv(gzfile(file.path(syspath, fname)), row.names=1)
    fname <- "human_cerebellum_coldata_250x4816.csv.gz"
    coldata <- read.csv(gzfile(file.path(syspath, fname)), row.names=1)
    fname <- "human_cerebellum_spatialcoords_250x4816.csv.gz"
    spatialcoords <- as.matrix(read.csv(gzfile(file.path(syspath, fname)),
                                        row.names=1))
    spe <- SpatialExperiment(assays=list(logcounts=logcounts),
                             rowData=rowdata,
                             colData=coldata,
                             spatialCoords=spatialcoords,
                             sample_id="HumanCerebellum_WholeTranscriptome")
    spe <- addImg(spe, sample_id="HumanCerebellum_WholeTranscriptome",
                  image_id="lowres",
                  imageSource=file.path(syspath, "human_cerebellum_lowres.png"),
                  scaleFactor=0.0450045, load=TRUE)
 
    set.seed(123) ## for reproducibility of the random gene sets
    ## build two gene sets with 4 randomly chosen genes and one
    ## third gene set with a few microglia marker genes
    gsets <- list(gset1=sample(rownames(spe), size=4, replace=FALSE),
                  gset2=sample(rownames(spe), size=4, replace=FALSE),
                  microglia=c("ENSG00000078808", "ENSG00000116251",
                              "ENSG00000142583", "ENSG00000173372"))

    ## calculate GSVA enrichment scores
    gsvapar <- gsvaParam(spe, gsets, verbose=FALSE)
    es <- gsva(gsvapar, verbose=FALSE)
    checkTrue(is(es, "SpatialExperiment"))
    checkTrue(all(dim(es) == c(length(gsets), ncol(spe))))
    checkTrue(all(colnames(es) == colnames(spe)))

    ## calculate spatial autocorrelation on the GSVA enrichment scores
    r <- spatCor(es, verbose=FALSE)
    checkTrue(all(r$observed[r$gene_id == "microglia"] > r$observed[r$gene_id != "microglia"]))
}
