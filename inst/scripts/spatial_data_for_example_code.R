library(Matrix)
library(TENxVisiumData)
library(scuttle)
library(R.utils)
library(EBImage)

spe <- HumanCerebellum()
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$symbol)
spe <- addPerCellQC(spe, subsets=list(Mito=is_mito))
discardmask <- spe$sum < 250 | spe$detected < 200 | spe$subsets_mito_percent > 40
spe <- spe[, !discardmask]
spe <- spe[rowSums(assay(spe)) > 100, ]
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

## subset to 250 random genes, including a few microglia markers
## ("ENSG00000078808", "ENSG00000116251", "ENSG00000142583", "ENSG00000173372")
microgliamarkers <- c("ENSG00000078808", "ENSG00000116251", "ENSG00000142583",
                      "ENSG00000173372")
samplegenes <- setdiff(rownames(spe), microgliamarkers)
set.seed(12345)
samplegenes <- sample(samplegenes, 250 - length(microgliamarkers))
samplegenes <- c(samplegenes, microgliamarkers)
spe <- spe[samplegenes, ]

## remove columns with less than 5 genes detected (~2% of 250)
spe <- spe[, colSums(assay(spe) > 0) >= 5]

fname <- sprintf("human_cerebellum_norm_logcounts_250x%d.mtx", ncol(spe))
writeMM(assays(spe)$logcounts, fname)
compressFile(fname, ext="gz", FUN=gzfile)

fname <- sprintf("human_cerebellum_rowdata_250x%d.csv", ncol(spe))
write.csv(rowData(spe), fname, row.names=TRUE)
compressFile(fname, ext="gz", FUN=gzfile)

fname <- sprintf("human_cerebellum_coldata_250x%d.csv", ncol(spe))
write.csv(colData(spe)[, "sizeFactor", drop=FALSE], fname, row.names=TRUE)
compressFile(fname, ext="gz", FUN=gzfile)

fname <- sprintf("human_cerebellum_spatialcoords_250x%d.csv", ncol(spe))
write.csv(spatialCoords(spe), fname, row.names=TRUE)
compressFile(fname, ext="gz", FUN=gzfile)

raw_img <- Image(imgRaster(get_img(spe))) ## convert to EBImage format
fname <- sprintf("human_cerebellum_raw_image_250x%d.png", ncol(spe))
writeImage(raw_img, fname, type="png") ## use EBImage::writeImage() to save the image


## how to build back the SpatialExperiment object from the files,
## but without the counts assay

fname <- "human_cerebellum_norm_logcounts_250x4816.mtx.gz"
logcounts <- as(readMM(gzfile(fname)), "CsparseMatrix")
fname <- "human_cerebellum_rowdata_250x4816.csv.gz"
rowdata <- read.csv(gzfile(fname), row.names=1)
fname <- "human_cerebellum_coldata_250x4816.csv.gz"
coldata <- read.csv(gzfile(fname), row.names=1)
fname <- "human_cerebellum_spatialcoords_250x4816.csv.gz"
spatialcoords <- as.matrix(read.csv(gzfile(fname), row.names=1))

spe <- SpatialExperiment(assays=list(logcounts=logcounts),
                         rowData=rowdata,
                         colData=coldata,
                         spatialCoords=spatialcoords,
                         sample_id="HumanCerebellum_WholeTranscriptome")
spe <- addImg(spe, sample_id="HumanCerebellum_WholeTranscriptome",
              image_id="lowres",
              imageSource="human_cerebellum_lowres.png",
              scaleFactor=0.0450045, load=TRUE)
