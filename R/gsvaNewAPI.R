#' @title Gene Set Variation Analysis
#' 
#' @description Estimates GSVA enrichment scores.
#' 
#' @param param A parameter object of one of the following classes:
#' * A [`gsvaParam`] object built using the constructor function
#' [`gsvaParam`].
#'   This object will trigger `gsva()` to use the GSVA algorithm by
#'   Hänzelmann et al. (2013).
#' * A [`plageParam`] object built using the constructor function
#' [`plageParam`].
#'   This object will trigger `gsva()` to use the PLAGE algorithm by
#'   Tomfohr et al. (2005).
#' * A [`zscoreParam`] object built using the constructor function
#' [`zscoreParam`].
#'   This object will trigger `gsva()` to use the combined z-score algorithm by
#'   Lee et al. (2008).
#' * A [`ssgseaParam`] object built using the constructor function
#' [`ssgseaParam`].
#'   This object will trigger `gsva()` to use the ssGSEA algorithm by
#'   Barbie et al. (2009).
#'
#' @param verbose Gives information about each calculation step. Default: `TRUE`.
#' 
#' @param BPPARAM An object of class `BiocParallelParam` specifying parameters
#'   related to the parallel execution of some of the tasks and calculations
#'   within this function.
#' 
#' @return A gene-set by sample matrix of GSVA enrichment scores stored in a
#' container object of the same type as the input expression data container. If
#' the input was a base matrix or a `dgCMatrix` object, then the output will
#' be a base matrix object with the gene sets employed in the calculations
#' stored in an attribute called `geneSets`. If the input was an
#' `ExpressionSet` object, then the output will be also an `ExpressionSet`
#' object with the gene sets employed in the calculations stored in an
#' attribute called `geneSets`. If the input was an object of one of the
#' classes described in [`GsvaExprData`], such as a `SingleCellExperiment`,
#' then the output will be of the same class, where enrichment scores will be
#' stored in an assay called `es` and the gene sets employed in the
#' calculations will be stored in the `rowData` slot of the object under the
#' column name `gs`.
#' 
#' @seealso [`plageParam`], [`zscoreParam`], [`ssgseaParam`], [`gsvaParam`],
#' [`BiocParallelParam`][BiocParallel::BiocParallelParam-class],
#' [`dgCMatrix`][Matrix::dgCMatrix-class],
#' \code{\link[Biobase]{ExpressionSet}},
### we are using the plain Rd above because
###  #' [`ExpressionSet`][Biobase::ExpressionSet-class],
### results in the following R CMD check NOTE:
### Non-topic package-anchored link(s) in Rd file 'gsva.Rd':
###  ‘[Biobase:class.ExpressionSet]{ExpressionSet}’
#' [`SingleCellExperiment`][SingleCellExperiment::SingleCellExperiment-class]
#'
#' @aliases gsva
#' @name gsva
#' @rdname gsva
#' 
#' @references Barbie, D.A. et al. Systematic RNA interference reveals that
#' oncogenic KRAS-driven cancers require TBK1.
#' *Nature*, 462(5):108-112, 2009.
#' \doi{10.1038/nature08460}
#'
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' \doi{10.1186/1471-2105-14-7}
#'
#' @references Lee, E. et al. Inferring pathway activity toward precise
#' disease classification.
#' *PLoS Comp Biol*, 4(11):e1000217, 2008.
#' \doi{10.1371/journal.pcbi.1000217}
#'
#' @references Tomfohr, J. et al. Pathway level analysis of gene expression
#' using singular value decomposition.
#' *BMC Bioinformatics*, 6:225, 2005.
#' \doi{10.1186/1471-2105-6-225}
#'
#' @examples
#' library(GSVA)
#' library(limma)
#' 
#' p <- 10 ## number of genes
#' n <- 30 ## number of samples
#' nGrp1 <- 15 ## number of samples in group 1
#' nGrp2 <- n - nGrp1 ## number of samples in group 2
#' 
#' ## consider three disjoint gene sets
#' geneSets <- list(set1=paste("g", 1:3, sep=""),
#'                  set2=paste("g", 4:6, sep=""),
#'                  set3=paste("g", 7:10, sep=""))
#'
#' ## sample data from a normal distribution with mean 0 and st.dev. 1
#' y <- matrix(rnorm(n*p), nrow=p, ncol=n,
#'             dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#'
#' ## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
#' y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2
#' 
#' ## build design matrix
#' design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nGrp1), rep(1, nGrp2)))
#' 
#' ## fit linear model
#' fit <- lmFit(y, design)
#' 
#' ## estimate moderated t-statistics
#' fit <- eBayes(fit)
#' 
#' ## genes in set1 are differentially expressed
#' topTable(fit, coef="sampleGroup2vs1")
#' 
#' ## build GSVA parameter object
#' gsvapar <- gsvaParam(y, geneSets)
#' 
#' ## estimate GSVA enrichment scores for the three sets
#' gsva_es <- gsva(gsvapar)
#' 
#' ## fit the same linear model now to the GSVA enrichment scores
#' fit <- lmFit(gsva_es, design)
#' 
#' ## estimate moderated t-statistics
#' fit <- eBayes(fit)
#' 
#' ## set1 is differentially expressed
#' topTable(fit, coef="sampleGroup2vs1")
NULL

#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @aliases gsva,plageParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="plageParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       verbose=verbose,
                                                       BPPARAM=BPPARAM)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if(verbose)
                  cli_alert_info(sprintf("Calculating PLAGE scores for %d gene sets",
                                         length(filteredMappedGeneSets)))

              plageScores <- plage(X=filteredDataMatrix,
                                   geneSets=filteredMappedGeneSets,
                                   verbose=verbose,
                                   BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), plageScores, gs)

              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @aliases gsva,zscoreParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="zscoreParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=TRUE,
                                                       removeNzConstant=TRUE,
                                                       verbose=verbose,
                                                       BPPARAM=BPPARAM)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              ## if (rnaseq)
              ##     stop("rnaseq=TRUE does not work with method='zscore'.")

              if(verbose)
                  cli_alert_info(sprintf("Calculating Z-scores for %d gene sets",
                                         length(filteredMappedGeneSets)))

              zScores <- zscore(X=filteredDataMatrix,
                                geneSets=filteredMappedGeneSets,
                                verbose=verbose,
                                BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), zScores, gs)
              
              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })

#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom utils packageDescription
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @aliases gsva,ssgseaParam-method
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="ssgseaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose)
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))

              famGaGS <- .filterAndMapGenesAndGeneSets(param,
                                                       removeConstant=FALSE,
                                                       removeNzConstant=FALSE,
                                                       verbose=verbose,
                                                       BPPARAM=BPPARAM)
              filteredDataMatrix <- famGaGS[["filteredDataMatrix"]]
              filteredMappedGeneSets <- famGaGS[["filteredMappedGeneSets"]]

              if (!inherits(BPPARAM, "SerialParam") && verbose) {
                  msg <- sprintf("Using a %s parallel back-end with %d workers",
                                 class(BPPARAM), bpnworkers(BPPARAM))
                  cli_alert_info(msg)
              }

              if(verbose)
                  cli_alert_info(sprintf("Calculating ssGSEA scores for %d gene sets",
                                         length(filteredMappedGeneSets)))

              ssgsea_sco <- ssgsea(X=filteredDataMatrix,
                                   geneSets=filteredMappedGeneSets,
                                   alpha=get_alpha(param), 
                                   normalization=get_normalize(param),
                                   any_na=anyNA(param),
                                   na_use=get_NAuse(param),
                                   minSize=get_minSize(param),
                                   verbose=verbose,
                                   BPPARAM=BPPARAM)

              gs <- .geneSetsIndices2Names(
                  indices=filteredMappedGeneSets,
                  names=rownames(filteredDataMatrix))
              rval <- wrapData(get_exprData(param), ssgsea_sco, gs)
              
              if (verbose)
                  cli_alert_success("Calculations finished")
              
              return(rval)
          })


#' @aliases gsva,gsvaParam-method
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils packageDescription
#' @rdname gsva
#' @exportMethod gsva
setMethod("gsva", signature(param="gsvaParam"),
          function(param,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose))
          {
              if (verbose) {
                  cli_alert_info(sprintf("GSVA version %s",
                                         packageDescription("GSVA")[["Version"]]))
                  gsva_global$show_start_and_end_messages <- FALSE
              }

              rankspar <- gsvaRanks(param=param, verbose=verbose,
                                    BPPARAM=BPPARAM)

              es <- gsvaScores(param=rankspar, verbose=verbose,
                               BPPARAM=BPPARAM)

              if (verbose) {
                  cli_alert_success("Calculations finished")
                  gsva_global$show_start_and_end_messages <- TRUE
              }
              
              return(es)
          })



### ----- methods for data pre-/post-processing -----

## unwrapData: extract a data matrix from a container object
setMethod("unwrapData", signature("matrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("dgCMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("SVT_SparseArray"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("DelayedMatrix"),
          function(container, assay) {
              return(container)
          })

setMethod("unwrapData", signature("ExpressionSet"),
          function(container, assay) {
              return(exprs(container))
          })

setMethod("unwrapData", signature("SummarizedExperiment"),
          function(container, assay) {
              if (length(assays(container)) == 0L)
                  stop("The input SummarizedExperiment object has no assay data.")

              if (missing(assay) || is.na(assay)) {
                  assay <- names(assays(container))[1]
              } else {
                  if (!is.character(assay))
                      stop("The 'assay' argument must contain a character string.")

                  assay <- assay[1]

                  if (!assay %in% names(assays(container)))
                      stop(sprintf("Assay %s not found in the input SummarizedExperiment object.", assay))
              }

              return(assays(container)[[assay]])
          })

setMethod("unwrapData", signature("SingleCellExperiment"),
          function(container, assay) {
              if (length(assays(container)) == 0L)
                  stop("The input SingleCellExperiment object has no assay data.")

              if (missing(assay) || is.na(assay)) {
                  assay <- names(assays(container))[1]
              } else {
                  if (!is.character(assay))
                      stop("The 'assay' argument must contain a character string.")

                  assay <- assay[1]

                  if (!assay %in% names(assays(container)))
                      stop(sprintf("Assay %s not found in the input SingleCellExperiment object.", assay))
              }

              return(assays(container)[[assay]])
          })

setMethod("unwrapData", signature("SpatialExperiment"),
          function(container, assay) {
            if (length(assays(container)) == 0L)
              stop("The input SpatialExperiment object has no assay data.")
            
            if (missing(assay) || is.na(assay)) {
              assay <- names(assays(container))[1]
            } else {
              if (!is.character(assay))
                stop("The 'assay' argument must contain a character string.")
              
              assay <- assay[1]
              
              if (!assay %in% names(assays(container)))
                stop(sprintf("Assay %s not found in the input SpatialExperiment object.", assay))
            }
            
            return(assays(container)[[assay]])
          })


## wrapData: put the resulting data and gene sets into the original data container type
setMethod("wrapData", signature(container="matrix"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="dgCMatrix"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="SVT_SparseArray"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="DelayedMatrix"),
          function(container, dataMatrix, geneSets) {
              if (!missing(geneSets))
                  attr(dataMatrix, "geneSets") <- geneSets
              return(dataMatrix)
          })

setMethod("wrapData", signature(container="ExpressionSet"),
          function(container, dataMatrix, geneSets) {
              rval <- new("ExpressionSet", exprs=dataMatrix,
                          phenoData=phenoData(container),
                          experimentData=experimentData(container),
                          annotation="")
              if (!missing(geneSets))
                  attr(rval, "geneSets") <- geneSets
              
              return(rval)
          })

setMethod("wrapData", signature(container="SummarizedExperiment"),
          function(container, dataMatrix, geneSets) {
              rdata <- adata <- NULL
              if (!missing(geneSets)) {
                  adata <- SimpleList(es=dataMatrix)
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## assume missing geneSets imples dataMatrix are ranks
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  adata <- c(assays(container[mask, ]),
                             SimpleList(gsvaranks=dataMatrix))
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SummarizedExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container))
              if (!missing(geneSets))
                  metadata(rval)$annotation <- NULL

              return(rval)
          })

setMethod("wrapData", signature(container="SingleCellExperiment"),
          function(container, dataMatrix, geneSets) {
              rdata <- adata <- NULL
              if (!missing(geneSets)) {
                  adata <- SimpleList(es=dataMatrix)
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## assume missing geneSets imples dataMatrix are ranks
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  adata <- c(assays(container[mask, ]),
                             SimpleList(gsvaranks=dataMatrix))
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SingleCellExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container))
              if (!missing(geneSets))
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })

setMethod("wrapData", signature(container="SpatialExperiment"),
          function(container, dataMatrix, geneSets) {
              rdata <- adata <- NULL
              if (!missing(geneSets)) {
                  adata <- SimpleList(es=dataMatrix)
                  rdata <- DataFrame(gs=CharacterList(geneSets))
              } else { ## assume missing geneSets imples dataMatrix are ranks
                  mask <- rownames(container) %in% rownames(dataMatrix)
                  adata <- c(assays(container[mask, ]),
                             SimpleList(gsvaranks=dataMatrix))
                  rdata <- rowData(container)[mask, ]
              }
              rval <- SpatialExperiment(
                  assays=adata,
                  colData=colData(container),
                  rowData=rdata,
                  metadata=metadata(container),
                  imgData=imgData(container),
                  spatialCoords=spatialCoords(container))
              if (!missing(geneSets))
                  metadata(rval)$annotation <- NULL
              
              return(rval)
          })


