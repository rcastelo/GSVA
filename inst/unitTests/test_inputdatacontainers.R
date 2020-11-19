test_inputdatacontainers <- function() {

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
  y <- matrix(rnorm(n*p), nrow=p, ncol=n,
              dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

  ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
  y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

  ## estimate GSVA enrichment scores with input as a matrix
  es.mat <- gsva(y, geneSets, verbose=FALSE)

  ## estimate GSVA enrichment scores with input as an ExpressionSet object
  y2 <- y
  rownames(y2) <- NULL
  eset <- Biobase::ExpressionSet(assayData=y2,
                                 phenoData=as(data.frame(dummy=1:ncol(y),
                                                         row.names=colnames(y)),
                                              "AnnotatedDataFrame"),
                                 featureData=as(data.frame(dummy=1:nrow(y),
                                                           row.names=rownames(y)),
                                                "AnnotatedDataFrame"))
  es.eset <- gsva(eset, geneSets, verbose=FALSE)

  checkTrue(identical(es.mat, Biobase::exprs(es.eset)))

  ## estimate GSVA enrichment scores with input as a SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(assay=list(counts=y2),
                             rowData=S4Vectors::DataFrame(data.frame(dummy=1:nrow(y),
                                                          row.names=rownames(y))),
                             colData=S4Vectors::DataFrame(data.frame(dummy=1:ncol(y),
                                                          row.names=colnames(y))))
  es.se <- gsva(se, geneSets, verbose=FALSE)

  checkTrue(identical(es.mat, SummarizedExperiment::assays(es.se)[[1]]))

  ## estimate GSVA enrichment scores with input as a dgCMatrix object
  yMat <- Matrix::Matrix(y, sparse=TRUE)
  es.dgCMat <- gsva(yMat, geneSets, verbose=FALSE)

  checkTrue(identical(es.mat, es.dgCMat))
}
