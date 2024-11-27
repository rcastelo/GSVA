## check for presence of valid row/feature names
#' @importFrom Biobase featureNames
.check_rownames <- function(expr) {
    ## CHECK: is this the right place to check this?
    ## 21/10/24: let's do it at parameter constructor
    if (is.null(rownames(expr)))
        cli_abort(c("x"="The input assay object doesn't have rownames"))
    else if (any(duplicated(rownames(expr)))) {
        cli_abort(c("x"="The input assay object has duplicated rownames"))
    } 
}


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

#' @importFrom sparseMatrixStats rowRanges
#' @importFrom cli cli_alert_warning cli_abort
.filterGenes <- function(expr, removeConstant=TRUE, removeNzConstant=TRUE) {
    geneRanges <- rowRanges(expr, na.rm=TRUE, useNames=FALSE)
    constantGenes <- (geneRanges[, 1] == geneRanges[, 2])

    if (any(constantGenes) || anyNA(constantGenes)) {
        invalidGenes <- (constantGenes | is.na(constantGenes))
        msg <- sprintf("%d genes with constant values throughout the samples",
                       sum(invalidGenes))
        cli_alert_warning(msg)
        if (removeConstant) {
            cli_alert_warning("Genes with constant values are discarded")
            expr <- expr[!invalidGenes, ]
        }
    }

    if (is(expr, "dgCMatrix")) {
        nzGeneList <- .sparse2columnList(t(expr))
        nzGeneRanges <- vapply(nzGeneList, FUN=range, FUN.VALUE=double(2))
        constantNzGenes <- (nzGeneRanges[1,] == nzGeneRanges[2,])

        if (any(constantNzGenes) || anyNA(constantNzGenes)) {
            invalidNzGenes <- (constantNzGenes | is.na(constantNzGenes))
            msg <- sprintf("%d genes with constant non-zero values throughout the samples",
                           sum(invalidNzGenes))
            cli_alert_warning(msg)
            if (removeNzConstant) {
                cli_alert_warning("Genes with constant non-zero values are discarded")
                expr <- expr[!invalidNzGenes, ]
            }
        }
    }

    if (nrow(expr) < 2)
        cli_abort(c("x"="Less than two genes in the input assay object"))
    
    return(expr)
}


## maps gene sets content in 'gsets' to 'features', where 'gsets'
## is a 'list' object with character string vectors as elements,
## and 'features' is a character string vector object. it assumes
## features in both input objects follow the same nomenclature,

#' @importFrom cli cli_abort
.mapGeneSetsToFeatures <- function(gsets, features) {

    ## Aaron Lun's suggestion at
    ## https://github.com/rcastelo/GSVA/issues/39#issuecomment-765549620
    gsets2 <- CharacterList(gsets)
    mt <- match(gsets2, features)
    mapdgenesets <- as.list(mt[!is.na(mt)])

    if (length(unlist(mapdgenesets, use.names=FALSE)) == 0) {
      msg <- paste("No identifiers in the gene sets could be matched to the",
                   "identifiers in the expression data.")
      cli_abort(c("x"=msg))
    }

    mapdgenesets
}

## it assumes that all arguments have been already checked for correctness
#' @importFrom cli cli_abort
.filterAndMapGeneSets <- function(param, wgset=NA, filteredDataMatrix, verbose) {

    geneSets <- get_geneSets(param)
    if (!is.na(wgset))
        geneSets <- geneSets[wgset]

    minSize <- get_minSize(param)
    maxSize <- get_maxSize(param)

    ## note that the method for 'GeneSetCollection' calls geneIds(), i.e., 
    ## whatever the input, from here on we have a list of character vectors
    geneSets <- mapGeneSetsToAnno(geneSets=geneSets,
                                  anno=get_annotation(param),
                                  verbose=verbose)
    
    ## map to the actual features for which expression data is available
    ## note that the result is a list of integer vectors (indices to rownames)
    ## and not a list of character vector any longer
    mappedGeneSets <- .mapGeneSetsToFeatures(geneSets, rownames(filteredDataMatrix))
    
    ## remove gene sets from the analysis for which no features are available
    ## and meet the minimum and maximum gene-set size specified by the user
    filteredMappedGeneSets <- filterGeneSets(mappedGeneSets,
                                             minSize=minSize,
                                             maxSize=maxSize)
    
    if (length(filteredMappedGeneSets) == 0) {
        msg <- "No gene set left after mapping and filtering."
        cli_abort(c("x"=msg))
    }

    ## this should NEVER happen -- just to make sure it doesn't...
    if (anyDuplicated(names(filteredMappedGeneSets)) > 0) {
        msg <- "The gene set list contains duplicated gene set names."
        cli_abort(c("x"=msg))
    }

    if (any(lengths(filteredMappedGeneSets) == 1)) {
        msg <- "Some gene sets have size one. Consider setting minSize > 1"
        cli_alert_warning(msg)
    }

    return(filteredMappedGeneSets)
}

#' @importFrom cli cli_alert_warning
.filterAndMapGenesAndGeneSets <- function(param,
                                          removeConstant=TRUE,
                                          removeNzConstant=TRUE,
                                          verbose=FALSE) {
    exprData <- get_exprData(param)
    dataMatrix <- unwrapData(exprData, get_assay(param))
    
    ## filter genes according to various criteria,
    ## e.g., constant expression
    filteredDataMatrix <- .filterGenes(dataMatrix,
                                       removeConstant=removeConstant,
                                       removeNzConstant=removeNzConstant)

    filteredMappedGeneSets <- .filterAndMapGeneSets(param=param,
                                                    filteredDataMatrix=filteredDataMatrix,
                                                    verbose=verbose)

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

    if (is.null(gs))
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
## ellipsis added for cases such as when FUN=rank where we may need
## to set the parameter 'ties.method' of the 'rank()' function
.sparseColumnApplyAndReplace <- function(m, FUN, ...) {
    x <- lapply(.sparse2columnList(m), FUN=FUN, ...)
    m@x <- unlist(x, use.names=FALSE)
    if (is.integer(m@x)) ## rank(ties.method="first") returns integers
        mode(m@x) <- "numeric" ## dgCMatrix holds only doubles and logicals
    return(m)
}

#' @importFrom cli cli_abort cli_alert_danger
.check_for_na_values <- function(exprData, assay, checkNA, use) {
    autonaclasseswocheck <- c("matrix", "ExpressionSet",
                              "SummarizedExperiment",
                              "RangedSummarizedExperiment")
    mask <- class(exprData) %in% autonaclasseswocheck
    checkNAyesno <- switch(checkNA, yes="yes", no="no",
                           ifelse(any(mask), "yes", "no"))
    didCheckNA <- any_na <- FALSE
    if (checkNAyesno == "yes") {
        any_na <- anyNA(unwrapData(exprData, assay))
        didCheckNA <- TRUE
        if (any_na) {
            if (use == "all.obs")
                cli_abort(c("x"="Input expression data has NA values."))
            else if (use == "everything")
                cli_alert_warning(paste("Input expression data has NA values,",
                                       "which will be propagated through",
                                       "calculations."))
            else ## na.rm
                cli_alert_warning(paste("Input expression data has NA values,",
                                       "which will be discarded from",
                                       "calculations."))
        }
    }

    list(any_na=any_na, didCheckNA=didCheckNA)
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


## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.isPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
  (name %in% loadedNamespaces())
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
    if(is.null(x)) {
        cat(paste0(prefix, "none."))
    } else {
        cat(paste0(prefix, capture.output(gsvaShow(x))), sep="\n")
    }
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


### 2024-08-02  axel: the following three functions have been copied from
### GSEABase/R/utilities.R (v. 1.66.0) as our implementation of readGMT()
### is mostly based on (a copy of) GSEABase::getGmt() which is making use
### of these utility functions.  Since we decided that GSVA::readGMT() may
### return a list of gene sets as well as a GeneSetCollection, it should work
### if a user doesn't have GSEABase installed at all.

## Placeholder 'till something appropriate decided
.uniqueIdentifier <- local({
    node <- NULL
    pid <- NULL
    uid <- 0L
    function() {
        if (is.null(node)) {
            node <<- Sys.info()['nodename']
            pid <<- Sys.getpid()
        }
        uid <<- uid + 1L
        base::paste(node, pid, date(), uid, sep=":")
    }
})

.stopf <- function(...) {
    call <- match.call(call=sys.call(sys.parent(1)))
    msg <- paste(sprintf(...), collapse="\n")
    stop(simpleError(msg, call=call))
}

.warningf <- function(...) {
    call <- match.call(call=sys.call(sys.parent(1)))
    msg <- paste(sprintf(...), collapse="\n")
    warning(simpleWarning(msg, call=call))
}
### end of copy from GSEABase/R/utilities.R
