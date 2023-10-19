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
            " genes with constant expression values throughout the samples.")
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


.guessIfCountData <- function(x, tolerance = sqrt(.Machine$double.eps)) {
    return(typeof(x) == "integer" ||
           (all(x >= 0) && all(x - round(x) < tolerance)))
}


.objPkgClass <- function(obj) {
    oc <- class(obj)
    pkg <- attr(oc, "package", exact=TRUE)
    opc <- if(is.null(pkg)) {
               oc[1]
           } else {
               paste(pkg[1], oc[1], sep = "::")
           }
    return(opc)
}

.showSome <- function(x) {
    paste0(paste(Biobase::selectSome(x, 4), collapse=", "),
           " (", length(x), " total)")
}

.catObj <- function(x, prefix = "  ") {
    cat(paste0(prefix, capture.output(gsvaShow(x))), sep="\n")
}

.isCharNonEmpty <- function(x) {
    return((!is.null(x)) &&
           (length(x) > 0) &&
           (is.character(x)) &&
           (!all(is.na(x))) &&
           (any(nchar(x) > 0)))
}

.isCharLength1 <- function(x) {
    return((.isCharNonEmpty(x)) && (length(x) == 1))
}

## annotation package checks
.isAnnoPkgValid <- function(ap) {
    return(.isCharLength1(ap))
}

.isAnnoPkgInstalled <- function(ap) {
    ap <- c(ap, paste0(ap, ".db"))
    return(any(ap %in% rownames(installed.packages())))
}

