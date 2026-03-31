##
## functions and methods to deal with gene sets
##

## ----- methods for retrieving gene sets -----

#' @title Retrieve or Determine Gene Sets
#' 
#' @description Retrieves or determines the gene sets that have been used
#' or would be used in a `gsva()` gene set analysis.  These are not necessarily
#' the same as the input gene sets.  See Details.
#' 
#' @param obj An object of one of the following classes:
#' * An expression data object of one of the classes described in
#' [`GsvaExprData-class`] that is the return value of a call to `gsva()`.
#' * A parameter object of one of the classes described in
#' [`GsvaMethodParam-class`] that could be used in a call to `gsva()`.
#'
#' @return The `geneSets()` methods return a named list of character vectors
#' where each character vector contains the gene IDs of a gene set.
#' The `geneSetSizes()` methods return a named integer vector of gene set sizes.
#' 
#' @details The gene sets used in a `gsva()` gene set analysis, or just their
#' sizes, may be a valuable input to subsequent analyses.  However, they are not
#' necessarily the same as the original input gene sets, or their sizes: based
#' on user choices, the gene annotation used, or presence/absence of genes in
#' gene sets and expression data set, `gsva()` may have to modify them during
#' the preparation of an analysis run.
#' In order to make use of these gene sets or their sizes, you can either
#' * retrieve them from the object returned by `gsva()` by passing this object
#' to `geneSets()` or `geneSetSizes()`, or
#' * predict them by calling `geneSets()` or `geneSetSizes()` on the parameter
#' object that would also be passed to `gsva()`.  This is much slower and should
#' only be done if you do not intend to run an actual gene set analysis.
#'
#' `geneSetSizes()` is a convenience wrapper running `lengths()` on the list of
#' gene sets returned by `geneSets()`.
#'
#' @aliases geneSets geneSetSizes
#'
#' @examples
#'
#' library(GSVA)
#'
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#'
#' gsets <- list(set1=paste0("g", 1:3),
#'               set2=paste0("g", 4:6),
#'               set3=paste0("g", 7:10),
#'               set4=paste0("g", 10:13)) ## genes not in the expression data
#' gsets
#'
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' gsvapar <- gsvaParam(y, gsets)
#' geneSets(gsvapar)
#'
#' es <- gsva(gsvapar)
#'
#' geneSets(es)
#'
#' @name geneSets
#' @rdname geneSets
#' 
NULL

#' @aliases geneSets,GsvaMethodParam-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("GsvaMethodParam"),
          function(obj) {
              famGaGS <- .filterAndMapGenesAndGeneSets(obj)

              return(.geneSetsIndices2Names(
                  indices=famGaGS[["filteredMappedGeneSets"]],
                  names=rownames(famGaGS[["filteredDataMatrix"]])
              ))
          })

#' @aliases geneSets,SummarizedExperiment-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("SummarizedExperiment"),
          function(obj) {
              return(as(rowData(obj)$gs, "list"))
          })

#' @aliases geneSets,SingleCellExperiment-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("SingleCellExperiment"),
          function(obj) {
              return(as(rowData(obj)$gs, "list"))
          })

#' @aliases geneSets,SpatialExperiment-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("SpatialExperiment"),
          function(obj) {
              return(as(rowData(obj)$gs, "list"))
          })

#' @aliases geneSets,GsvaExprData-method
#' @rdname geneSets
#' @exportMethod geneSets
setMethod("geneSets", signature("GsvaExprData"),
          function(obj) {
              return(.geneSets(obj))
          })


#' @aliases geneSetSizes,GsvaMethodParam-method
#' @rdname geneSets
#' @exportMethod geneSetSizes
setMethod("geneSetSizes", signature("GsvaMethodParam"),
          function(obj) {
              return(lengths(geneSets(obj)))
          })

#' @aliases geneSetSizes,GsvaExprData-method
#' @rdname geneSets
#' @exportMethod geneSetSizes
setMethod("geneSetSizes", signature("GsvaExprData"),
          function(obj) {
              return(lengths(geneSets(obj)))
          })


## ----- helper functions for gene set I/O and preprocessing -----

#' @title Handling of Duplicated Gene Set Names
#' 
#' @description Offers a choice of ways for handling duplicated gene set names
#' that may not be suitable as input to other gene set analysis functions.
#' 
#' @param geneSets A named list of gene sets represented as character vectors
#' of gene IDs as e.g. returned by [`readGMT`].
#'
#' @param deduplUse A character vector of length 1 specifying one of several
#' methods to handle duplicated gene set names.
#' Duplicated gene set names are explicitly forbidden by the
#' [GMT file format specification](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
#' but can nevertheless be encountered in the wild.
#' The available choices are:
#' * `first` (the default): drops all gene sets whose names are duplicated
#' according to the base R function and retains only the first occurence of a
#' gene set name.
#' * `drop`:  removes *all* gene sets that have a duplicated name, including its
#' first occurrence.
#' * `union`: replaces gene sets with duplicated names by a single gene set
#' containing the union of all their gene IDs.
#' * `smallest`: drops gene sets with duplicated names and retains only the
#' smallest of them, i.e. the one with the fewest gene IDs.  If there are
#' several smallest gene sets, the first will be selected.
#' * `largest`: drops gene sets with duplicated names and retains only the
#' largest of them, i.e. the one with the most gene IDs.  If there are
#' several largest gene sets, the first will be selected.
#'
#' @return A named list of gene sets represented as character vectors of
#' gene IDs.
#' 
#' @examples
#'
#' library(GSVA)
#'
#' gsets <- list(gs1=LETTERS[1:3], gs2=LETTERS[4:6], gs2=LETTERS[5:8])
#' gsets
#'
#' deduplicateGeneSets(gsets)
#' deduplicateGeneSets(gsets, deduplUse="drop")
#' deduplicateGeneSets(gsets, deduplUse="union")
#' deduplicateGeneSets(gsets, deduplUse="smallest")
#' deduplicateGeneSets(gsets, deduplUse="largest")
#'
#' fname <- system.file("extdata", "c2.subsetdups.v7.5.symbols.gmt.gz",
#'                      package="GSVAdata")
#'
#' ## readGMT() calls internally deduplicateGeneSets() and it takes the
#' ## parameter 'deduplUse' which is passed to the internal call
#' c2.dupgenesets <- readGMT(fname, deduplUse="union")
#' c2.dupgenesets
#' any(duplicated(names(c2.dupgenesets)))
#'
#' @aliases deduplicateGeneSets
#' @name deduplicateGeneSets
#' @rdname deduplicateGeneSets
#' @export
#' 
deduplicateGeneSets <- function(geneSets,
                                deduplUse = c("first", "drop", "union",
                                              "smallest", "largest")) {
    ddUse <- match.arg(deduplUse)
    isNameDuplicated <- duplicated(names(geneSets))
    duplicatedNames <- unique(names(geneSets[isNameDuplicated]))

    ## a nested list containing sublists of duplicated gene sets
    duplicatedGeneSets <- lapply(duplicatedNames,
                                 function(dn, gs) unname(gs[dn == names(gs)]),
                                 gs = geneSets)

    ## transformation function operating on sublists of such nested lists,
    ## returning a single deduplicated gene set, i.e. character vector
    ddFunc <- switch(ddUse,
                     union=function(dgs) Reduce(union, dgs),
                     smallest=function(dgs) dgs[[which.min(lengths(dgs))]],
                     largest=function(dgs) dgs[[which.max(lengths(dgs))]])

    ## apply transformation function to deduplicate gene sets (if requested)
    if(!is.null(ddFunc))
        dedupl <- lapply(duplicatedGeneSets, FUN=ddFunc)

    ## drop all duplicate gene sets (sufficient for default of "first")
    geneSets[isNameDuplicated] <- NULL

    ## remove or replace non-duplicated with deduplicated gene sets
    if(ddUse == "drop") {
        geneSets[duplicatedNames] <- NULL
    } else if(!is.null(ddFunc)) {
        geneSets[duplicatedNames] <- dedupl
    }

    return(geneSets)
}


deduplicateGmtLines <- function(geneSets,
                                deduplUse = c("first", "drop", "union",
                                              "smallest", "largest")) {
    ddUse <- match.arg(deduplUse)
    gsName <- vapply(X=geneSets, FUN="[", FUN.VALUE=character(1), 1)
    isNameDuplicated <- which(duplicated(gsName))

    if(length(isNameDuplicated) > 0) {
        warning("GMT contains duplicated gene set names; deduplicated",
                " using method: ", ddUse)
        duplicatedNames <- unique(gsName[isNameDuplicated])
        lIdxDuplGS <- lapply(duplicatedNames,
                             function(DN, GSN) which(DN == GSN),
                             GSN = gsName)
        idxReplace <- vapply(X=lIdxDuplGS, FUN="[", FUN.VALUE=integer(1), 1)
        idxRemove <- unique(unlist(lapply(lIdxDuplGS, tail, -1)))

        ddFunc <- switch(ddUse,
                         union = function(idxGS, lGS) {
                             ld <- lGS[idxGS]
                             gsn <- lapply(ld, head, 1)[[1]]
                             gsd <- do.call("paste", c(lapply(ld, "[", 2), sep = " | "))
                             gsg <- Reduce(union, lapply(ld, tail, -2))
                             c(gsn, gsd, gsg)
                         },
                         smallest = function(idxGS, lGS) {
                             lGS[idxGS][[which.min(lengths(lGS[idxGS]))]]
                         },
                         largest = function(idxGS, lGS) {
                             lGS[idxGS][[which.max(lengths(lGS[idxGS]))]]
                         })

        if(!is.null(ddFunc)) {
            dedupl <- lapply(lIdxDuplGS, FUN = ddFunc, lGS = geneSets)
            geneSets[idxReplace] <- dedupl
        } else if(ddUse == "drop") {
            idxRemove <- union(idxRemove, idxReplace)
        }
        
        geneSets[idxRemove] <- NULL
    }

    return(geneSets)
}


#' @title Guess the gene identifier type from a list of character vectors
#' 
#' @description This function tries to derive the type of gene IDs used in a
#' named list of `character` vectors provided as input.
#' 
#' @param geneIdsList A named list of character vectors like the ones returned
#' by `geneIds()`.
#'
#' @return An object of a subclass of `GeneIdentifierType` derived from the
#' input.
#'
#' @details In order to make this function useful and keep it as simple as
#' possible, we limit ourselves to the most common types of gene identifiers:
#' "Gene IDs" consisting of digits only are considered ENTREZ IDs, anything
#' starting with 'ENS' an ENSEMBL identifier and anything else a HuGO gene
#' symbol.
#' 
#' @seealso [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class]
#'
#' @aliases guessGeneIdType
#'
#' @examples
#'
#' library(GSVA)
#'
#' gsets <- list(INNATE_RESPONSE=c("AIM2", "ALPK1", "AP3B1"),
#'               ADAPTIVE_RESPONSE=c("CD27", "CD70", "EBAG9"))
#'
#' idtype <- guessGeneIdType(gsets)
#' idtype
#' class(idtype)
#'
#' @name guessGeneIdType
#' @rdname guessGeneIdType
#' @importFrom GSEABase EntrezIdentifier ENSEMBLIdentifier SymbolIdentifier
#' @export
#' 
guessGeneIdType <- function(geneIdsList) {
    allIds <- unlist(geneIdsList)
    
    if(all(grepl("^[[:digit:]]+$", allIds))) {
        retVal <- EntrezIdentifier()
    } else if(all(grepl("^ENS", allIds))) {
        retVal <- ENSEMBLIdentifier()
    } else {
        retVal <- SymbolIdentifier()
    }

    return(retVal)
}


#' @title Construct a GeneSetCollection object from a list of character vectors
#' 
#' @description This function is essentially the reverse of
#' `GSEABase::geneIds()`, i.e., it takes as input a named list of `character`
#' vectors representing gene sets and returns the corresponding
#' GeneSetCollection object.
#' 
#' @param geneIdsList A named list of character vectors like the ones returned
#' by `geneIds()`.  Names must be unique; otherwise see `deduplicateGeneSets()`
#' for a number of strategies to resolve this issue.  
#'
#' @param geneIdType By default a character vector of length 1 with the special
#' value `"auto"` or an object of a subclass of `GeneIdentifierType`.  If set
#' to `"auto"`, the function will try to derive the gene ID type from argument
#' `geneIdsList` using [`guessGeneIdType`].
#' Other values, including `NULL`, will be ignored with a warning and
#' `geneIdType=NullIdentifier()` will be used instead.
#' The gene ID type of all `GeneSet` objects in the resulting
#' `GeneSetCollection` will be set to this value.
#' 
#' @param collectionType An object of class `CollectionType`.  The collection
#' type of all `GeneSet` objects in the resulting `GeneSetCollection` will be
#' set to this value but can afterwards be modified for individual `GeneSet`s
#' if necessary.
#'
#' @return An object of class `GeneSetCollection` with all its `GeneSet`
#' objects using the gene ID and collection types specified by the corresponding
#' arguments.  Applying function `geneIds()` to this object should return a list
#' identical to the `geneIdsList` argument.
#' 
#' @seealso [`GeneSetCollection`][GSEABase::GeneSetCollection-class],
#' [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class],
#' [`geneIds`][GSEABase::geneIds],
#' [`deduplicateGeneSets`],
#' [`guessGeneIdType`],
#' [`GeneSet`][GSEABase::GeneSet-class]
#'
#' @aliases geneIdsToGeneSetCollection
#'
#' @examples
#'
#' library(GSVA)
#'
#' gsets <- list(INNATE_RESPONSE=c("AIM2", "ALPK1", "AP3B1"),
#'               ADAPTIVE_RESPONSE=c("CD27", "CD70", "EBAG9"))
#' gsets
#' geneIdsToGeneSetCollection(gsets)
#'
#' @name geneIdsToGeneSetCollection
#' @rdname geneIdsToGeneSetCollection
#' @importFrom GSEABase NullIdentifier NullCollection
#' @export
#' 
geneIdsToGeneSetCollection <- function(geneIdsList,
                                       geneIdType="auto",
                                       collectionType=NullCollection()) {
    if (inherits(geneIdType, "character") && (geneIdType == "auto")) {
        git <- gsvaAnnotation(geneIdsList)
        if (is.null(git)) {
            git <- guessGeneIdType(geneIdsList)
        }
    } else if (inherits(geneIdType, "GeneIdentifierType")) {
        git <- geneIdType
    } else {
        git <- NullIdentifier()
        cli_alert_warning(paste0("Invalid value of argument `geneIdType` ",
                                 "ignored, using `NullIdentifier()` instead."))
    }
    
    return(GeneSetCollection(mapply(function(gn, gs) {
        if(anyDuplicated(gs) > 0) {
            gs <- unique(gs)
            msg <- sprintf("Duplicated gene IDs removed from gene set %s", gn)
            cli_alert_warning(msg)
        }
        
        GeneSet(gs,
                geneIdType=git,
                collectionType=collectionType,
                setName=gn)
    }, gn=names(geneIdsList), gs=geneIdsList)))
}


#' @title Import Gene Sets from a GMT File
#' 
#' @description Imports a list of gene sets from a GMT (Gene Matrix Transposed)
#' format file, offering a choice of ways to handle duplicated gene set names.
#' 
#' @param con A connection object or a non-empty character string of length 1
#' containing e.g. the filename or URL of a (possibly compressed) GMT file. 
#'
#' @param sep The character string separating members of each gene set in the
#' GMT file.
#'
#' @param geneIdType By default a character vector of length 1 with the special
#' value `"auto"` or an object of a subclass of `GeneIdentifierType`.  If set
#' to `"auto"`, the function will try to derive the gene ID type from argument
#' `geneIdsList` using [`guessGeneIdType`].
#' Other values, including `NULL`, will be ignored with a warning and
#' `geneIdType=NullIdentifier()` will be used instead.
#' Depending on the value of argument `valueType`, the gene ID type of the
#' resulting list or of all `GeneSet` objects in the resulting
#' `GeneSetCollection` will be set to this value.
#' 
#' @param collectionType Only used when `valueType == "GeneSetCollection"`. See
#' `getGmt` for more information.
#'
#' @param valueType A character vector of length 1 specifying the desired type
#' of return value.  It must be one of:
#' * `GeneSetCollection` (the default): a `GeneSetCollection` object as defined
#' and described by package `GSEABase`.
#' * `list`: a named list of gene sets represented as character vectors of gene IDs.
#' This format is much simpler and cannot store the metadata required for automatic
#' mapping of gene IDs.
#'
#' @param deduplUse A character vector of length 1 specifying one of several
#' methods to handle duplicated gene set names.
#' Duplicated gene set names are explicitly forbidden by the
#' [GMT file format specification](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
#' but can nevertheless be encountered in the wild.
#' The available choices are:
#' * `first` (the default): drops all gene sets whose names are duplicated
#' according to the base R function and retains only the first occurence of a
#' gene set name.
#' * `drop`:  removes *all* gene sets that have a duplicated name, including its
#' first occurrence.
#' * `union`: replaces gene sets with duplicated names by a single gene set
#' containing the union of all their gene IDs.
#' * `smallest`: drops gene sets with duplicated names and retains only the
#' smallest of them, i.e. the one with the fewest gene IDs.  If there are
#' several smallest gene sets, the first will be selected.
#' * `largest`: drops gene sets with duplicated names and retains only the
#' largest of them, i.e. the one with the most gene IDs.  If there are
#' several largest gene sets, the first will be selected.
#'
#' @param ... Further arguments passed on to `readLines()`
#' 
#' @return The gene sets imported from the GMT file, with duplicate gene sets
#' resolved according to argument `deduplUse` and in the format determined by
#' argument `valueType`.
#' 
#' @seealso [`deduplicateGeneSets`], [`readLines`],
#' [`GeneSetCollection`][GSEABase::GeneSetCollection-class],
#' [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class],
#' \code{\link[GSEABase]{getGmt}},
### we are using the plain Rd above because
###  #' [`getGmt`][GSEABase::getGmt]
### results in the following R CMD check NOTE:
### Non-topic package-anchored link(s) in Rd file 'readGMT.Rd':
###   ‘[GSEABase:getObjects]{getGmt}’
#'
#' @examples
#' library(GSVA)
#' suppressPackageStartupMessages(library(GSVAdata))
#'
#' fname <- file.path(system.file("extdata", package="GSVAdata"),
#'    "c2.subsetdups.v7.5.symbols.gmt.gz")
#' 
#' ## by default, guess geneIdType from content and return a GeneSetCollection
#' genesets <- readGMT(fname)
#' genesets
#'
#' ## how to manually override the geneIdType
#' genesets <- readGMT(fname, geneIdType=NullIdentifier())
#' genesets
#' 
#' ## how to drop *all* gene sets with duplicated names (instead of ignoring
#' ## only the duplicated one)
#' genesets <- readGMT(fname, deduplUse="drop")
#' genesets
#' 
#' ## return a simple list instead of a GeneSetCollection
#' genesets <- readGMT(fname, valueType="list")
#' head(genesets, 2)
#'
#' ## the list has a geneIdType, too
#' gsvaAnnotation(genesets)
#'
#' @aliases readGMT
#' @name readGMT
#' @rdname readGMT
#' @importFrom Biobase selectSome
#' @export
#' 
readGMT <- function (con,
                     sep = "\t",
                     geneIdType = "auto",
                     collectionType = NullCollection(), 
                     valueType = c("GeneSetCollection", "list"),
                     deduplUse = c("first", "drop", "union", "smallest", "largest"),
                     ...) {
    valueType <- match.arg(valueType)

    if ((!.isCharLength1(con)) && (!inherits(con, "connection"))) {
        msg <- paste("Argument 'con' is not a valid filename, URL",
                     "or connection.")
        cli_abort(c("x"=msg))
    }
    
    ## from GSEABase::getGmt()
    lines <- strsplit(readLines(con, ...), sep)
    if (any(lengths(lines) < 2)) {
        txt <- paste("all records in the GMT file must have >= 2 fields", 
                     "\n  first invalid line:  %s\n", collapse = "")
        .stopf(txt, lines[lengths(lines) < 2][[1]])
    }
    dups <- new.env(parent = emptyenv())
    lines <- lapply(lines, function(elt, dups) {
        dmask <- duplicated(elt[-(1:2)])
        if (any(dmask)) {
            dups[[elt[[1]]]] <- unique(elt[-(1:2)][dmask])
            elt <- c(elt[1:2], unique(elt[-(1:2)]))
        }
        elt
    }, dups)
    if (length(dups)) 
        .warningf("%d record(s) contain duplicate ids: %s", length(dups), 
                  paste(selectSome(sort(ls(dups))), collapse = ", "))

    ## our small addition to tolerate duplicate gene set names
    lines <- deduplicateGmtLines(lines, deduplUse)

    if(inherits(geneIdType, "character") && (geneIdType == "auto")) {
        geneIdType <- guessGeneIdType(lapply(lines, tail, -2))
    } else if(!inherits(geneIdType, "GeneIdentifierType")) {
        geneIdType <- NullIdentifier()
        cli_alert_warning(paste0("Invalid value of argument `geneIdType` ",
                                 "ignored, using `NullIdentifier()` instead."))
    } ## else: fine, no?
    
    ## on second thoughts, another small addition: let the user choose the return type
    if(valueType == "GeneSetCollection") {
        ## from GSEABase::getGmt()
        template <- GeneSet(geneIdType = geneIdType, collectionType = collectionType)
        return(GeneSetCollection(lapply(lines, function(line) {
            initialize(template, geneIds = unlist(line[-(1:2)]), 
                       setName = line[[1]], shortDescription = line[[2]], 
                       setIdentifier = .uniqueIdentifier())
        })))
    } else if(valueType == "list") {
        gs <- lapply(lines, tail, -2)
        names(gs) <- vapply(X=lines, FUN=head, FUN.VALUE=character(1), 1)

        ## even more thoughts, now that we make use of gene ID metadata in lists
        if(!is.null(geneIdType)) {
            gsvaAnnotation(gs) <- geneIdType
        }
        
        return(gs)
    }
}


## ----- methods for retrieving and setting annotation metadata -----

#' @title Store and Retrieve Annotation Metadata
#' 
#' @description Methods for storing and retrieving annotation metadata in
#' expression data objects that support it.  If gene sets and expression data
#' are using different but known gene identifier types and an appropriate
#' annotation database is available, gene set identifiers can be mapped to
#' expression data identifiers without manual user intervention, e.g. from
#' an MSigDb gene set using ENTREZ IDs or gene symbols to an expression data
#' set using ENSEMBL IDs.
#' 
#' @param object An expression data object of one of the classes described in
#' [`GsvaExprData-class`].  Simple `matrix` and `dgCMatrix` objects are not
#' capable of storing annotation metadata and will return `NULL`.
#'
#' @param value For the replacement methods, the annotation metadata to be
#' stored in the object.  For `ExpressionSet` objects, this must be a
#' character of length 1 specifying the name of the annotation database to be
#' used.  For `SummarizedExperiment` and its subclasses, this must be
#' a `GeneIdentifierType` created by one of the constructors from package
#' `GSEABase` where the `annotation` argument is typically the name of an
#' organism or annotation database, e.g. `org.Hs.eg.db`.  Simple `matrix` and
#' `dgCMatrix` objects are not capable of storing annotation metadata and the
#' attempt to do so will result in an error.
#'
#' @return For the retrieval methods, the annotation metadata stored in the
#' object or `NULL`.  For the replacement methods, the updated object.
#'
#' @examples
#'
#' library(GSEABase)
#' library(GSVA)
#' library(GSVAdata)
#'
#' data(geneprotExpCostaEtAl2021)
#' se <- geneExpCostaEtAl2021
#' se
#'
#' gsvaAnnotation(se)
#' gsvaAnnotation(se) <- EntrezIdentifier("org.Hs.eg.db")
#' gsvaAnnotation(se)
#'
#' @seealso
#' \code{\link[Biobase]{ExpressionSet}},
### we are using the plain Rd above because
###  #' [`ExpressionSet`][Biobase::ExpressionSet-class],
### results in the following R CMD check NOTE:
### Non-topic package-anchored link(s) in Rd file 'gsvaAnnotation.Rd':
###  ‘[Biobase:class.ExpressionSet]{ExpressionSet}’
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment-class],
#' [`GeneIdentifierType`][GSEABase::GeneIdentifierType-class],
#' [`dgCMatrix`][Matrix::dgCMatrix-class]
#' 
#' @aliases gsvaAnnotation gsvaAnnotation<-
#' @name gsvaAnnotation
#' @rdname gsvaAnnotation
#' 
NULL


#' @aliases gsvaAnnotation,GsvaExprData-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="GsvaExprData"),
          function(object) {
              return(attr(object, which="geneIdType", exact=TRUE))
          })

#' @aliases gsvaAnnotation<-,GsvaExprData,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                     object="GsvaExprData",
                     value="GeneIdentifierType"),
          function(object, value) {
              attr(object, which="geneIdType") <- value
              object
          })

#' @aliases gsvaAnnotation,ExpressionSet-method
#' @rdname gsvaAnnotation
#' @importFrom BiocGenerics annotation
#' @importFrom GSEABase AnnoOrEntrezIdentifier
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="ExpressionSet"),
          function(object) {
              ## always a character giving the db pkg, potentially empty ("")
              ## unfortunately, as it turns out, even character(0) sometimes
              ao <- annotation(object)
              if(.isCharLength1(ao)) {
                  return(AnnoOrEntrezIdentifier(ao))
              } else {
                  return(NULL)
              }
          })

#' @aliases gsvaAnnotation<-,ExpressionSet,character-method
#' @rdname gsvaAnnotation
#' @importFrom BiocGenerics "annotation<-"
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="ExpressionSet",
                   value="character"),
                 function(object, value) {
                     annotation(object) <- value
                     object
                 })

#' @aliases gsvaAnnotation<-,ExpressionSet,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @importFrom BiocGenerics annotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="ExpressionSet",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     gsvaAnnotation(object) <- annotation(value)
                     object
                 })

#' @aliases gsvaAnnotation,SummarizedExperiment-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation", signature("SummarizedExperiment"),
          function(object) {
              ## NULL if unset; otherwise anything but we *expect* and handle
              ## a GSEABase::GeneIdentifierType with or without annotation(),
              ## i.e., db pkg, available.  Same for subclasses below.
              return(metadata(object)$annotation)
          })

#' @aliases gsvaAnnotation<-,SummarizedExperiment,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="SummarizedExperiment",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     metadata(object)$annotation <- value
                     object
                 })

#' @aliases gsvaAnnotation,SingleCellExperiment-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation", signature("SingleCellExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })

#' @aliases gsvaAnnotation<-,SingleCellExperiment,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="SingleCellExperiment",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     metadata(object)$annotation <- value
                     object
                 })

#' @aliases gsvaAnnotation,SpatialExperiment-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation", signature("SpatialExperiment"),
          function(object) {
              return(metadata(object)$annotation)
          })

#' @aliases gsvaAnnotation<-,SpatialExperiment,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                   object="SpatialExperiment",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     metadata(object)$annotation <- value
                     object
                 })


#' @aliases gsvaAnnotation,list-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="list"),
          function(object) {
              return(attr(object, which="geneIdType", exact=TRUE))
          })

#' @aliases gsvaAnnotation<-,list,GeneIdentifierType-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setReplaceMethod("gsvaAnnotation",
                 signature=signature(
                     object="list",
                     value="GeneIdentifierType"),
          function(object, value) {
              attr(object, which="geneIdType") <- value
              object
          })

#' @aliases gsvaAnnotation,GeneSetCollection-method
#' @rdname gsvaAnnotation
#' @exportMethod gsvaAnnotation
setMethod("gsvaAnnotation",
          signature=signature(object="GeneSetCollection"),
          function(object) {
              lgit <- unique(lapply(object, geneIdType))
              return(if(length(lgit) == 1) lgit[[1]] else NULL)
          })


## these are internal helper methods for mapping gene/feature identifiers
## in gene sets from one type to another using the annotation metadata
## provided in the input argument 'anno'.

setMethod("mapGeneSetsToAnno", signature(geneSets="list", anno="NULL"),
          function(geneSets, anno, verbose=FALSE) {
              return(geneSets)
          })

setMethod("mapGeneSetsToAnno", signature(geneSets="list", anno="character"),
          function(geneSets, anno, verbose=FALSE) {
              gsc <- geneIdsToGeneSetCollection(geneIdsList=geneSets)
              return(mapGeneSetsToAnno(gsc, anno))
          })

setMethod("mapGeneSetsToAnno",
          signature(geneSets="list", anno="GeneIdentifierType"),
          function(geneSets, anno, verbose=FALSE) {
              gsc <- geneIdsToGeneSetCollection(geneIdsList=geneSets)
              return(mapGeneSetsToAnno(gsc, anno))
          })

setMethod("mapGeneSetsToAnno",
          signature(geneSets="GeneSetCollection", anno="NULL"),
          function(geneSets, anno, verbose=FALSE) {
              return(geneIds(geneSets))
          })

#' @importFrom cli cli_alert_info cli_alert_warning
#' @importFrom GSEABase AnnoOrEntrezIdentifier mapIdentifiers
setMethod("mapGeneSetsToAnno",
          signature(geneSets="GeneSetCollection", anno="character"),
          function(geneSets, anno, verbose=FALSE) {
              if(.isAnnoPkgValid(anno)) {
                  if(!.isAnnoPkgInstalled(anno)) {
                      msg <- "Please install the annotation package %s"
                      stop(sprintf(msg, anno))
                  }

                  if (verbose)
                      cli_alert_info("Mapping identifiers")

                  mappedGeneSets <- mapIdentifiers(geneSets,
                                                   AnnoOrEntrezIdentifier(anno))
                  rval <- geneIds(mappedGeneSets)

              } else {
                  if (verbose) {
                      msg <- paste("No annotation metadata available in the",
                                   "input expression data object")
                      cli_alert_warning(msg)
                      msg <- paste("Attempting to directly match identifiers",
                                   "in expression data to gene sets")
                      cli_alert_warning(msg)
                  }

                  rval <- geneIds(geneSets)
              }

              return(rval)
          })

#' @importFrom cli cli_alert_info cli_alert_warning
#' @importFrom GSEABase mapIdentifiers geneIds
#' @importFrom BiocGenerics annotation
setMethod("mapGeneSetsToAnno",
          signature(geneSets="GeneSetCollection",
                    anno="GeneIdentifierType"),
          function(geneSets, anno, verbose=FALSE) {
              annoDb <- annotation(anno)

              if(.isAnnoPkgValid(annoDb)) {
                  if(!.isAnnoPkgInstalled(annoDb)) {
                      msg <- "Please install the annotation package %s"
                      stop(sprintf(msg, annoDb))
                  }

                  if (verbose)
                      cli_alert_info("Mapping identifiers")

                  mappedGeneSets <- mapIdentifiers(geneSets, anno)
                  rval <- geneIds(mappedGeneSets)

              } else {
                  if (verbose) {
                      msg <- paste("No annotation metadata available in the",
                                   "input expression data object")
                      cli_alert_warning(msg)
                      msg <- paste("Attempting to directly match identifiers",
                                   "in expression data to gene sets")
                      cli_alert_warning(msg)
                  }

                  rval <- geneIds(geneSets)
              }

              return(rval)
          })


## ----- methods for filtering gene sets -----

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

        return(gSets[gSetsLen >= minSize & gSetsLen <= maxSize])	
})

#' @aliases filterGeneSets,GeneSetCollection-method
#' @rdname filterGeneSets
#' @exportMethod filterGeneSets
setMethod("filterGeneSets", signature(gSets="GeneSetCollection"),
    function(gSets, minSize=1, maxSize=Inf) {
        filterGeneSets(geneIds(gSets), minSize, maxSize)
})



## ----- methods for computing gene sets overlap -----

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
setMethod("computeGeneSetsOverlap",
          signature(gSets="list", uniqGenes="character"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
    totalGenes <- length(uniqGenes)

    ## map to the actual features for which expression data is available
    gSets <- .mapGeneSetsToFeatures(gSets, uniqGenes)

    lenGsets <- lengths(gSets)
    totalGsets <- length(gSets)

    gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                    dimnames=list(uniqGenes, names(gSets)))
    members <- cbind(unlist(gSets, use.names=FALSE),
                     rep(seq_len(totalGsets), times=lenGsets))
    gSetsMembershipMatrix[members] <- 1

    .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,list,ExpressionSet-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap",
          signature(gSets="list", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
    uniqGenes <- featureNames(uniqGenes)
    totalGenes <- length(uniqGenes)

    ## map to the actual features for which expression data is available
    gSets <- .mapGeneSetsToFeatures(gSets, uniqGenes)

    lenGsets <- lengths(gSets)
    totalGsets <- length(gSets)

    gSetsMembershipMatrix <- matrix(0, nrow=totalGenes, ncol=totalGsets,
                                    dimnames=list(uniqGenes, names(gSets)))
    members <- cbind(unlist(gSets, use.names=FALSE),
                     rep(seq_len(totalGsets), times=lenGsets))
    gSetsMembershipMatrix[members] <- 1

    .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,GeneSetCollection,character-method
#' @rdname computeGeneSetsOverlap
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap",
          signature(gSets="GeneSetCollection", uniqGenes="character"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {

    gSetsMembershipMatrix <- incidence(gSets)
    mask <- colnames(gSetsMembershipMatrix) %in% uniqGenes
    gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, mask])

    .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

#' @aliases computeGeneSetsOverlap,GeneSetCollection,ExpressionSet-method
#' @rdname computeGeneSetsOverlap
#' @importFrom GSEABase mapIdentifiers AnnoOrEntrezIdentifier
#' @importFrom Biobase featureNames
#' @importFrom BiocGenerics annotation
#' @exportMethod computeGeneSetsOverlap
setMethod("computeGeneSetsOverlap",
          signature(gSets="GeneSetCollection", uniqGenes="ExpressionSet"),
          function(gSets, uniqGenes, minSize=1, maxSize=Inf) {
    ## map gene identifiers of the gene sets to the features in the chip
    gSets <- mapIdentifiers(gSets,
                            AnnoOrEntrezIdentifier(annotation(uniqGenes)))
  
    uniqGenes <- featureNames(uniqGenes)

    gSetsMembershipMatrix <- incidence(gSets)
    mask <- colnames(gSetsMembershipMatrix) %in% uniqGenes
    gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, mask])

    .computeGeneSetsOverlap(gSetsMembershipMatrix, minSize, maxSize)
})

.computeGeneSetsOverlap <- function(gSetsMembershipMatrix,
                                    minSize=1, maxSize=Inf) {
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
                 dimnames=list(colnames(gSetsMembershipMatrix),
                               colnames(gSetsMembershipMatrix)))
    M2 <- t(M1)
    M.min <- matrix(0, nrow=totalGsets, ncol=totalGsets)
    M.min[M1 < M2] <- M1[M1 < M2]
    M.min[M2 <= M1] <- M2[M2 <= M1]
    overlapMatrix <- M / M.min

    return(overlapMatrix)
}
