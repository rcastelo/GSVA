
## methods for the slots common to all GSVA method parameter classes

## ----- getters -----

setGeneric("get_exprData", function(object) standardGeneric("get_exprData"))
setGeneric("get_geneSets", function(object) standardGeneric("get_geneSets"))

setMethod("get_exprData", signature("GsvaMethodParam"),
          function(object) {
              return(object@exprData)
          })


setMethod("get_geneSets", signature("GsvaMethodParam"),
          function(object) {
              return(object@geneSets)
          })
