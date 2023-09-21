
#' @export 
setGeneric("gsva",
           function(expr, gset.idx.list, ...) standardGeneric("gsva"))

#' @export 
setGeneric("filterGeneSets",
           function(gSets, ...) standardGeneric("filterGeneSets"))

#' @export 
setGeneric("computeGeneSetsOverlap",
           function(gSets, uniqGenes=unique(unlist(gSets, use.names=FALSE)), ...) standardGeneric("computeGeneSetsOverlap"))


#' @export 
setGeneric("unwrapData",
           function(container, ...) standardGeneric("unwrapData"))

#' @export 
setGeneric("wrapData",
           function(dataMatrix, container) standardGeneric("wrapData"))

#' @export 
setGeneric("mapGeneSetsToAnno",
           function(geneSets, ...) standardGeneric("mapGeneSetsToAnno"))

#' @export 
setGeneric("get_exprData", function(object) standardGeneric("get_exprData"))

#' @export 
setGeneric("get_geneSets", function(object) standardGeneric("get_geneSets"))


