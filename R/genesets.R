##
## functions and methods to deal with gene sets
##


#' @title Filter gene sets
#' 
#' @description Filters gene sets through a given minimum and maximum set size.
#' 
#' This function filters the input gene sets according to a given minimum and
#' maximum set size.
#' 
#' @aliases filterGeneSets
#'
#' @name filterGeneSets
#' 
#' @rdname filterGeneSets
#'
#' @param gSets Gene sets given either as a `list` or a
#' `GeneSetCollection` object.
#' 
#' @param minSize Minimum size.
#' 
#' @param maxSize Maximum size.
#' 
#' @return A collection of gene sets that meet the given minimum and maximum
#' set size.
#' 
#' @author J. Guinney
#' 
#' @seealso [`computeGeneSetsOverlap`]
#' 
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' 
#' @keywords Gene set
#' 
#' @examples
#' geneSets <- list(set1=as.character(1:4), set2=as.character(4:10))
#' filterGeneSets(geneSets, minSize=5)
NULL

#' @aliases filterGeneSets,list-method
#' @rdname filterGeneSets
#' @exportMethod filterGeneSets
setMethod("filterGeneSets", signature(gSets="list"),
          function(gSets, minSize=1, maxSize=Inf) {
	gSetsLen <- lengths(gSets)
	return (gSets[gSetsLen >= minSize & gSetsLen <= maxSize])	
})

#' @aliases filterGeneSets,GeneSetCollection-method
#' @rdname filterGeneSets
#' @exportMethod filterGeneSets
setMethod("filterGeneSets", signature(gSets="GeneSetCollection"),
          function(gSets, minSize=1, maxSize=Inf) {
  filterGeneSets(geneIds(gSets), minSize, maxSize)
})


#' @title Compute gene-sets overlap
#' 
#' @description Calculates the overlap among every pair of gene-sets given as
#' input.
#' 
#' This function calculates the overlap between every pair of gene sets of the
#' input argument `gSets`. Before this calculation takes place, the gene
#' sets in `gSets` are firstly filtered to discard genes that do not match
#' to the identifiers in `uniqGenes`. Secondly, they are further filtered
#' to meet the minimum and/or maximum size specified with the arguments
#' `minSize` and `maxSize`. The overlap between two gene sets is
#' calculated as the number of common genes between the two gene sets divided
#' by the smallest size of the two gene sets.
#' 
#' @aliases computeGeneSetsOverlap
#'
#' @name computeGeneSetsOverlap
#'
#' @rdname computeGeneSetsOverlap
#' 
#' @param gSets Gene sets given either as a `list` or a
#' `GeneSetCollection` object.
#' 
#' @param uniqGenes Vector of unique genes to be considered when calculating
#' the overlaps.
#' 
#' @param minSize Minimum size.
#' 
#' @param maxSize Maximum size.
#' 
#' @return A gene-set by gene-set matrix of the overlap among every pair of
#' gene sets.
#' 
#' @author J. Guinney
#' 
#' @seealso [`filterGeneSets`]
#' 
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' 
#' @keywords Gene set
#' 
#' @examples
#' geneSets <- list(set1=as.character(1:4), set2=as.character(4:10))
#' computeGeneSetsOverlap(geneSets, unique(unlist(geneSets)))
NULL

#' @aliases computeGeneSetsOverlap,list,character-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="character"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
  totalGenes <- length(uniqGenes)

  ## map to the actual features for which expression data is available
  gSets <- .mapGeneSetsToFeatures(gSets, uniqGenes)

  lenGsets <- lengths(gSets)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,list,ExpressionSet-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
  uniqGenes <- featureNames(uniqGenes)
  totalGenes <- length(uniqGenes)

  ## map to the actual features for which expression data is available
  gSets <- .mapGeneSetsToFeatures(gSets, uniqGenes)

  lenGsets <- lengths(gSets)
  totalGsets <- length(gSets)

  gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members <- cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] <- 1

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,GeneSetCollection,character-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="character"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,GeneSetCollection,ExpressionSet-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap", signature(gSets="GeneSetCollection", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
  ## map gene identifiers of the gene sets to the features in the chip
  ## Biobase::annotation() is necessary to disambiguate from the
  ## 'annotation' argument
  gSets <- mapIdentifiers(gSets, AnnoOrEntrezIdentifier(Biobase::annotation(uniqGenes)))
  
  uniqGenes <- featureNames(uniqGenes)

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

.computeGeneSetsOverlap <- function(gSetsMembershipMatrix, minSize=1, maxSize=Inf) {
  ## gSetsMembershipMatrix should be a (genes x gene-sets) incidence matrix

  lenGsets <- colSums(gSetsMembershipMatrix)

  szFilterMask <- lenGsets >= max(1, minSize) & lenGsets <= maxSize
  if (!any(szFilterMask))
    stop("No gene set meets the minimum and maximum size filter\n")

  gSetsMembershipMatrix <- gSetsMembershipMatrix[, szFilterMask]
  lenGsets <- lenGsets[szFilterMask]

  totalGsets <- ncol(gSetsMembershipMatrix)

  M <- t(gSetsMembershipMatrix) %*% gSetsMembershipMatrix

  M1 <- matrix(lenGsets, nrow=totalGsets, ncol=totalGsets,
               dimnames=list(colnames(gSetsMembershipMatrix), colnames(gSetsMembershipMatrix)))
  M2 <- t(M1)
  M.min <- matrix(0, nrow=totalGsets, ncol=totalGsets)
  M.min[M1 < M2] <- M1[M1 < M2]
  M.min[M2 <= M1] <- M2[M2 <= M1]
  overlapMatrix <- M / M.min

  return (overlapMatrix)
}
