#' @title Compute Spatial Autocorrelation for SpatialExperiment objects
#' 
#' @description Computes spatial autocorrelation using Moran's I statistic
#' for a \code{SpatialExperiment} object, using an inverse squared distance
#' weight matrix as default, or an inverse distance weight matrix as an
#' alternative. It also tests for spatial autocorrelation assuming normality.
#' 
#' @param spe An object of \code{SpatialExperiment} class.
#' @param assay Character vector of length 1, specifying the name of the assay
#' to use. By default, an assay called 'logcounts' will be used if present,
#' otherwise the first assay is used.
#' @param alternative A character string specifying the alternative hypothesis
#' tested against the null hypothesis of no spatial autocorrelation; must be
#' one of "two.sided", "less", or "greater", or any unambiguous abbreviation of
#' these.
#' @param na.rm A logical indicating whether missing values should be removed.
#' @param squared A logical indicating whether the inverse distance weight
#' matrix should be squared or not.
#' @param verbose Gives information about each calculation step. Default:
#' `TRUE`.
#' @param BPPARAM An object of class `BiocParallelParam` specifying parameters
#' related to the parallel execution of some of the tasks and calculations
#' within this function.
#' 
#' @return A \code{data.frame} with the same row names as the original
#' `SpatialExperiment` object. Columns include the observed Moran's I statistic,
#' the expected Moran's I statistic under no spatial autocorrelation, the
#' expected standard deviation under no spatial autocorrelation, and the
#' p-value of the test.
#'
#' @seealso [`BiocParallelParam`][BiocParallel::BiocParallelParam-class]
#' 
#' @aliases spatCor spatCor,SpatialExperiment-method
#'
#' @name spatCor
#' @rdname spatCor
#'
#' @examples
#'
#' library(Matrix)
#' library(SpatialExperiment)
#'
#' ## build a SpatialExperiment object
#' syspath <- system.file("extdata", package="GSVA")
#' fname <- "human_cerebellum_norm_logcounts_250x4816.mtx.gz"
#' logcounts <- as(readMM(gzfile(file.path(syspath, fname))), "CsparseMatrix")
#' fname <- "human_cerebellum_rowdata_250x4816.csv.gz"
#' rowdata <- read.csv(gzfile(file.path(syspath, fname)), row.names=1)
#' fname <- "human_cerebellum_coldata_250x4816.csv.gz"
#' coldata <- read.csv(gzfile(file.path(syspath, fname)), row.names=1)
#' fname <- "human_cerebellum_spatialcoords_250x4816.csv.gz"
#' spatialcoords <- as.matrix(read.csv(gzfile(file.path(syspath, fname)),
#'                                     row.names=1))
#'
#' spe <- SpatialExperiment(assays=list(logcounts=logcounts),
#'                          rowData=rowdata,
#'                          colData=coldata,
#'                          spatialCoords=spatialcoords,
#'                          sample_id="HumanCerebellum_WholeTranscriptome")
#' spe <- addImg(spe, sample_id="HumanCerebellum_WholeTranscriptome",
#'               image_id="lowres",
#'               imageSource=file.path(syspath, "human_cerebellum_lowres.png"),
#'               scaleFactor=0.0450045, load=TRUE)
#'
#' set.seed(123) ## for reproducibility of the random gene sets
#' ## build two gene sets with 4 randomly chosen genes and one
#' ## third gene set with a few microglia marker genes
#' gsets <- list(gset1=sample(rownames(spe), size=4, replace=FALSE),
#'               gset2=sample(rownames(spe), size=4, replace=FALSE),
#'               microglia=c("ENSG00000078808", "ENSG00000116251",
#'                           "ENSG00000142583", "ENSG00000173372"))
#'
#' ## calculate GSVA enrichment scores
#' gsvapar <- gsvaParam(spe, gsets, verbose=FALSE)
#' es <- gsva(gsvapar, verbose=FALSE)
#'
#' ## calculate spatial autocorrelation on the GSVA enrichment scores
#' spatCor(es, verbose=FALSE)
#'
#' @importFrom stats dist pnorm
#' @importFrom SummarizedExperiment assay
#' @importFrom cli cli_abort cli_progress_bar cli_progress_update
#' @importFrom BiocParallel bpnworkers bplapply bpprogressbar SerialParam
#' @exportMethod spatCor
#' @export

setMethod("spatCor", signature("SpatialExperiment"),
    function(spe, assay= NA_character_, na.rm=FALSE,
             alternative="two.sided", squared=TRUE, verbose=TRUE,
             BPPARAM=SerialParam(progressbar = verbose)) {

        .check_bpparam(BPPARAM)
        nworkers <- bpnworkers(BPPARAM)

        assay <- .check_assayNames(assay, spe, verbose)
        weight_list <- .spe_dist_weight_matrix(spe,squared)
        rowns <- rownames(spe)
        if (is.null(rowns))
            cli_abort(c("x"="Input rows have no names"))
        df_res <- data.frame(observed=numeric(), expected=numeric(),
                             sd=numeric(), p.value=numeric(),
                             sample_id=character())
        for (sample in unique(colData(spe)$sample_id)) {
            spe_Moran <- list()
            logc <- assay(spe[,colData(spe)$sample_id == sample], assay)
            logc <- .filterGenes(logc, verbose=verbose, BPPARAM=BPPARAM)
            if (is.null(BPPARAM) || nworkers == 1L) {
                env <- NULL
                if (verbose) {
                    env <- new.env(parent=globalenv())
                    progressmsg <- sprintf("Sample %s", sample)
                    assign("idpb", cli_progress_bar(progressmsg,
                                                    total=length(rowns)),
                           envir=env)
                }
                spe_Moran <- lapply(rowns, function(x, idpbe) {
                    res <- list(observed=NA, expected=NA, sd=NA, p.value=NA)
                    if (x %in% rownames(logc)) {
                        res <- .internal_moran(logc[x, ], weight_list,
                                               na.rm=na.rm,
                                               alternative=alternative)
                        if (verbose && is(idpbe, "environment"))
                            cli_progress_update(id=get("idpb", envir=idpbe), 1)
                    }
                    res
                }, idpb=env)
                if (verbose)
                    cli_progress_done(get("idpb", envir=env))
            } else {                               ## parallel execution
                if (verbose)
                    bpprogressbar(BPPARAM) <- TRUE ## reporting progress wo/ cli
                spe_Moran <- bplapply(rowns, function(x) {
                    res <- list(observed=NA, expected=NA, sd=NA, p.value=NA)
                    if (x %in% rownames(logc))
                        res <- .internal_moran(logc[x, ], weight_list, na.rm=na.rm,
                                               alternative=alternative)	       
                    res
                }, BPPARAM = BPPARAM)
            }
            df_sample <- do.call(rbind, lapply(spe_Moran, as.data.frame))
            df_sample <- data.frame(gene_id = rownames(spe), df_sample) 
            df_sample$sample_id <- sample
            df_res <- rbind(df_res,df_sample)
        }
        return(df_res)
})

.internal_moran <- function(x, weight_list, na.rm=FALSE,
                            alternative=c("two.sided", "less", "greater")) {

    alternative <- match.arg(alternative)
    weight <- weight_list$weight[names(x),names(x),drop=FALSE]
    n <- length(x)
    ei <- -1/(n - 1)
    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        } else {
            warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
            return(list(observed=NA, expected=ei, sd=NA, 
                        p.value=NA))
        }
    }
    m <- mean(x, na.rm = na.rm)
    y <- x - m
    cv <- sum(weight * y %o% y, na.rm = na.rm)
    v <- sum(y^2, na.rm = na.rm)
    s <- weight_list$s
    obs <- (n/s) * (cv/v)
    k <- (sum(y^4)/n)/(v/n)^2
    s.sq <- weight_list$s.sq
    S1 <- weight_list$S1
    S2 <- weight_list$S2
    sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) - 
                 k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq)) /
                 ((n - 1) * (n - 2) * (n - 3) * s.sq) - 1 / ((n - 1)^2))
    pv <- pnorm(obs, mean = ei, sd = sdi)
    if (alternative == "two.sided") 
        pv <- if (obs <= ei) 
                  2 * pv
              else
                  2 * (1 - pv)
    if (alternative == "greater") 
        pv <- 1 - pv

    list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}



.spe_dist_weight_matrix <- function(spe, squared=TRUE) {
    xy <- spatialCoords(spe)[unique(rownames(spatialCoords(spe))), ]
    weight <- as.matrix(dist(xy))
    if(squared == TRUE)
        weight <- 1 / weight^2
    else
        weight <- 1 / weight
    diag(weight) <- 0
    s <- sum(weight)
    ROWSUM <- rowSums(weight)
    ROWSUM[ROWSUM == 0] <- 1
    return(list(weight=weight,
                ROWSUM=ROWSUM,
                s=s,
                S1=0.5 * sum((2*weight)^2),
                S2=sum((ROWSUM + colSums(weight))^2),
                s.sq=s^2))
}
