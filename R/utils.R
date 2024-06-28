## 2024-02-06  axel: function .filterGenes() is intended to detect genes (rows)
##  with constant expression (and, hence, no information), warn about them and
##  optionally remove them (in particular, ssGSEA's choice is to keep them).
##  the original approach tried to identify genes with a standard deviation of
##  exactly 0 but failed in certain cases of all identical values due to the use
##  of floating point arithmetic, see issues:
## https://github.com/rcastelo/GSVA/issues/54
## https://github.com/HenrikBengtsson/matrixStats/issues/204
##  an improvement in matrixStats::rowSds() fixes the original issue but cannot
##  guarantee that there won't be other problematic cases.
##
## We propose to detect cases of constant gene expression by comparing genewise
##  min and max values rather than computing the SD, which *should* avoid using
##  floating point arithmetic in favour of comparisons and scale linearly with
##  the number of samples (columns) -- we'll of course have to check that. ;-)
##
## A related but different issue has recently surfaced when methods PLAGE and
##  z-scores are applied to sparse matrices and attempt to scale the non-zero
##  values of genes: genes that are constant in their non-zero values will have
##  an SD of 0 and therefore scaling them will result in division by 0.

.filterGenes <- function(expr, removeConstant=TRUE, removeNzConstant=TRUE) {
    geneRanges <- rowRanges(expr, na.rm=TRUE, useNames=FALSE)
    constantGenes <- (geneRanges[, 1] == geneRanges[, 2])

    if(any(constantGenes) || anyNA(constantGenes)) {
        invalidGenes <- (constantGenes | is.na(constantGenes))
        warning(sum(invalidGenes),
                " genes with constant values throughout the samples.")
        if(removeConstant) {
            warning("Genes with constant values are discarded.")
            expr <- expr[!invalidGenes, ]
        }
    }

    if(is(expr, "dgCMatrix")) {
        nzGeneList <- .sparse2columnList(t(expr))
        nzGeneRanges <- vapply(nzGeneList, FUN=range, FUN.VALUE=double(2))
        constantNzGenes <- (nzGeneRanges[1,] == nzGeneRanges[2,])

        if(any(constantNzGenes) || anyNA(constantNzGenes)) {
            invalidNzGenes <- (constantNzGenes | is.na(constantNzGenes))
            warning(sum(invalidNzGenes),
                    " genes with constant non-zero values throughout the sample.")
            if(removeNzConstant) {
                warning("Genes with constant non-zero values are discarded.")
                expr <- expr[!invalidNzGenes, ]
            }
        }
    }

    if(nrow(expr) < 2)
        stop("Less than two genes in the input assay object\n")
    
    ## CHECK: is this the right place to check this?
    if(is.null(rownames(expr)))
        stop("The input assay object doesn't have rownames\n")
    
    return(expr)
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


.filterAndMapGenesAndGeneSets <- function(param,
                                          removeConstant=TRUE,
                                          removeNzConstant=TRUE) {
    exprData <- get_exprData(param)
    dataMatrix <- unwrapData(exprData, get_assay(param))
    
    ## filter genes according to various criteria,
    ## e.g., constant expression
    filteredDataMatrix <- .filterGenes(dataMatrix,
                                       removeConstant=removeConstant,
                                       removeNzConstant=removeNzConstant)

    ## note that the method for 'GeneSetCollection' calls geneIds(), i.e., 
    ## whatever the input, from here on we have a list of character vectors
    geneSets <- mapGeneSetsToAnno(geneSets=get_geneSets(param),
                                  anno=gsvaAnnotation(exprData))
    
    ## map to the actual features for which expression data is available
    ## note that the result is a list of integer vectors (indices to rownames)
    ## and not a list of character vector any longer
    mappedGeneSets <- .mapGeneSetsToFeatures(geneSets, rownames(filteredDataMatrix))
    
    ## remove gene sets from the analysis for which no features are available
    ## and meet the minimum and maximum gene-set size specified by the user
    filteredMappedGeneSets <- filterGeneSets(mappedGeneSets,
                                             minSize=get_minSize(param),
                                             maxSize=get_maxSize(param))

    if(length(filteredMappedGeneSets) == 0)
        stop("The gene set list is empty! Filter may be too stringent.")

    ## this should NEVER happen -- just to make sure it doesn't...
    if(anyDuplicated(names(filteredMappedGeneSets)) > 0)
        stop("The gene set list contains duplicated gene set names.")

    if(any(lengths(filteredMappedGeneSets) == 1))
        warning("Some gene sets have size one. Consider setting 'minSize' to a value > 1.")

    return(list(filteredDataMatrix=filteredDataMatrix,
                filteredMappedGeneSets=filteredMappedGeneSets))
}


## (re-)extract a list of gene names from a list of indices
## (indices resulting from the matching above)
.geneSetsIndices2Names <- function(indices, names) {
    return(lapply(indices, function(i, n) n[i], n=names))
}


## access to gene set attribute without explicit use of attributes
.geneSets <- function(obj) {
    gs <- attr(obj, "geneSets", exact=TRUE)

    if(is.null(gs))
        stop("The object does not contain information about gene sets.")

    return(gs)
}


## converts a dgCMatrix into a list of its columns, based on
## https://rpubs.com/will_townes/sparse-apply
## it is only slightly more efficient than .sparseToList() below BUT simpler
## and does NOT offer converting to a list of rows which is far less efficient
## on a dgCMatrix object.  if you need lists of rows, simply transpose before
## calling this function, t() is reasonably fast as is calling vapply() on its
## result
.sparse2columnList <- function(m) {
    return(split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE)))
}

## actually, it's not just an apply() but also in-place modification
.sparseColumnApplyAndReplace <- function(m, FUN) {
    x <- lapply(.sparse2columnList(m), FUN=FUN)
    m@x <- unlist(x, use.names=FALSE)
    return(m)
}

.sparseScaleMessage <- function() {
    message("Please bear in mind that this method first scales the values of ",
            "the gene expression data. In order to take advantage of the ",
            "sparse Matrix type, the scaling will only be applied to the ",
            "non-zero values of the data. This is a provisional solution in ",
            "order to give support to the dgCMatrix format.")
}

## transforms a dgCMatrix into a list of its
## non-zero values by MARGIN (1 for row, 2 for column)
##
## currently unused because replaced by .sparse2columnList() and kept for the
## time being, just in case...
## .sparseToList <-function(dgCMat, MARGIN){
##   MARGIN <- as.integer(MARGIN)
##   J <- rep(1:ncol(dgCMat), diff(dgCMat@p))
##   I <- dgCMat@i + 1
##   x <- dgCMat@x
##   if (MARGIN == 1L) {
##     result <- split(x, I)
##     names(result) <- rownames(dgCMat)[as.numeric(names(result))]
##   } else if (MARGIN == 2L) {
##     result <- split(x, J)
##     names(result) <- colnames(dgCMat)[as.numeric(names(result))]
##   }
##   else {
##     warning("invalid MARGIN; return NULL")
##     result <- NULL
##   }
##   result
## }

## .dgCapply<-function(m, MARGIN, FUN){
##   x <- lapply(.sparseToList(m, MARGIN), FUN)
##   m@x <- unlist(x, use.names=FALSE)
##   m
## }


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

## utility function to make sure matrix to dgCMatrix coercion is uniform
## (since direct coercion to dgCMatrix is deprecated (!) by Matrix pkg)
.matrix2dgCMatrix <- function(m) {
    return(as(as(as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix"))
}
