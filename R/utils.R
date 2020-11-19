.filterFeatures <- function(expr, method) {

  ## filter out genes with constant expression values
  sdGenes <- apply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  } 

  expr
}

## transforms a dgCMatrix into a list of its
## non-zero values by MARGIN (1 for row, 2 for column)
.sparseToList <-function(dgCMat, MARGIN){
  MARGIN <- as.integer(MARGIN)
  J <- rep(1:ncol(dgCMat), diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    result <- split(x, I)
    names(result) <- rownames(dgCMat)[as.numeric(names(result))]
  } else if (MARGIN == 2L) {
    result <- split(x, J)
    names(result) <- colnames(dgCMat)[as.numeric(names(result))]
  }
  else {
    warning("invalid MARGIN; return NULL")
    result <- NULL
  }
  result
}

.dgCapply<-function(m,f, MARGIN){
  x <- lapply(.sparseToList(m, MARGIN), f)
  m@x <- unlist(x, use.names=FALSE)
  m
}

## filter out genes which non-zero values have
## constant expression values
.filterFeaturesSparse <- function(expr, method) {
  
  sdGenes <- sapply(.sparseToList(expr, 1), sd)
  
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  }
  
  expr
}