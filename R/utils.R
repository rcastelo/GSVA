.filterFeatures <- function(expr, method) {

  ## filter out genes with constant expression values
  ## DelayedMatrixStats::rowSds() works for both base and 
  ## DelayedArray matrices
  sdGenes <- DelayedMatrixStats::rowSds(expr)
  ## the following fixes this bug, see issues
  ## https://github.com/rcastelo/GSVA/issues/54
  ## https://github.com/HenrikBengtsson/matrixStats/issues/204
  sdGenes[sdGenes < 1e-10] <- 0
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  } 

  if (nrow(expr) < 2)
    stop("Less than two genes in the input assay object\n")
  
  if(is.null(rownames(expr)))
    stop("The input assay object doesn't have rownames\n")
  
  expr
}

## maps gene sets content in 'gsets' to 'features', where 'gsets'
## is a 'list' object with character string vectors as elements,
## and 'features' is a character string vector object. it assumes
## features in both input objects follow the same nomenclature,
.mapGeneSetsToFeatures <- function(gsets, features) {

  ## Aaron Lun's suggestion at
  ## https://github.com/rcastelo/GSVA/issues/39#issuecomment-765549620
  gsets2 <- CharacterList(gsets)
  mt <- match(gsets2, features)
  mapdgenesets <- as.list(mt[!is.na(mt)])

  if (length(unlist(mapdgenesets, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")

  mapdgenesets
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
