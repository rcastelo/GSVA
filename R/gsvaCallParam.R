
### get methods for the common slots of all GSVA call parameter classes

setGeneric("get_dataSet", function(object) standardGeneric("get_dataSet"))
setGeneric("get_geneSets", function(object) standardGeneric("get_geneSets"))

setMethod("get_dataSet", signature("gsvaCallParam"),
          function(object) {
              return(object@dataSet)
          })


setMethod("get_geneSets", signature("gsvaCallParam"),
          function(object) {
              return(object@geneSets)
          })
