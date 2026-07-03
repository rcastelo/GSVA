
#' @export 
setGeneric("gsva",
           function(param, ...) standardGeneric("gsva"))

#' @export 
setGeneric("gsvaRanks",
           function(param, ...) standardGeneric("gsvaRanks"))

#' @export 
setGeneric("gsvaScores",
           function(param, ...) standardGeneric("gsvaScores"))

#' @export 
setGeneric("filterGeneSets",
           function(gSets, ...) standardGeneric("filterGeneSets"))

#' @export 
setGeneric("computeGeneSetsOverlap",
           function(gSets, uniqGenes=unique(unlist(gSets, use.names=FALSE)), ...) standardGeneric("computeGeneSetsOverlap"))

#' @export
setGeneric("geneSets",
           function(obj, ...) standardGeneric("geneSets"))

#' @export
setGeneric("geneSets<-",
           function(object, value) standardGeneric("geneSets<-"))

#' @export
setGeneric("geneSetSizes",
           function(obj, ...) standardGeneric("geneSetSizes"))

#' @export
setGeneric("gsvaAnnotation",
           function(object) standardGeneric("gsvaAnnotation"))

#' @export
setGeneric("gsvaAnnotation<-",
           function(object, value) standardGeneric("gsvaAnnotation<-"))


## for now, these should be private methods

setGeneric("unwrapData",
           function(container, ...) standardGeneric("unwrapData"))

setGeneric("wrapData",
           function(container, ...) standardGeneric("wrapData"))

setGeneric("mapGeneSetsToAnno",
           function(geneSets, anno, ...) standardGeneric("mapGeneSetsToAnno"))

setGeneric("get_geneSets", function(object) standardGeneric("get_geneSets"))

setGeneric("get_assay", function(object) standardGeneric("get_assay"))

setGeneric("get_annotation", function(object) standardGeneric("get_annotation"))

setGeneric("get_minSize", function(object) standardGeneric("get_minSize"))

setGeneric("get_maxSize", function(object) standardGeneric("get_maxSize"))

setGeneric("gsvaShow", function(object) standardGeneric("gsvaShow"))

## spatial methods

#' @export
setGeneric("spatCor",
           function(spe, ...) standardGeneric("spatCor"))




