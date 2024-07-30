#' @importFrom stats dist pnorm
#' @importFrom SummarizedExperiment assay

### ----- Method for computing Spatial Autocorrelation for SpatialExperiment object -----

#' @title Compute Spatial Autocorrelation for SpatialExperiment objects
#' 
#' @description Computes spatial autocorrelation using Moran's I statistic
#' for a \code{SpatialExperiment} object, using an inverse squared distance weight matrix as default,
#' or an inverse distance weight matrix as an alternative. It also tests for spatial autocorrelation
#' assuming normality.
#' 
#' @param spe An object of \code{SpatialExperiment} class.
#' @param alternative A character string specifying the alternative hypothesis tested against the null hypothesis of no spatial autocorrelation;
#'  must be one of "two.sided", "less", or "greater", or any unambiguous abbreviation of these.
#' @param na.rm A logical indicating whether missing values should be removed.
#' @param squared A logical indicating whether the inverse distance weight matrix should be squared or not.
#' 
#' @return A \code{data.frame} with the same row names as the original \code{SpatialExperiment} object.
#' Columns include the observed Moran's I statistic, the expected Moran's I statistic under no spatial autocorrelation, the expected
#' standard deviation under no spatial autocorrelation, and the p-value of the test.
#' 
#' @aliases spatCor spatCor,SpatialExperiment-method
#' @name spatCor
#' @rdname spatCor
#' @exportMethod spatCor
#' @export

setMethod("spatCor", signature("SpatialExperiment"),
          function(spe, na.rm = FALSE, alternative = "two.sided", squared = TRUE) {
            weight_list <- .spe_dist_weight_matrix(spe, squared)
            logc <- assay(spe)
            
            spe_Moran <- list()
            rowns <- rownames(spe)
            spe_Moran <- lapply(rowns, function(x){
              .internal_moran(logc[rowns == x, ], weight_list, na.rm = na.rm, alternative = alternative)
            })
            names(spe_Moran) <- rownames(spe)
            df_res <- do.call(rbind, lapply(spe_Moran, as.data.frame))
            return(df_res)
          })

.internal_moran <- function (x, weight_list, na.rm = FALSE, alternative = "two.sided") 
{
  weight <- weight_list$weight
  n <- length(x)
  ei <- -1/(n - 1)
  nas <- is.na(x)
  if (any(nas)) {
    if (na.rm) {
      x <- x[!nas]
      n <- length(x)
      weight <- weight[!nas, !nas]
    }
    else {
      warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
      return(list(observed = NA, expected = ei, sd = NA, 
                  p.value = NA))
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
                 k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n - 
                                                                     1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
  alternative <- match.arg(alternative, c("two.sided", "less", 
                                          "greater"))
  pv <- pnorm(obs, mean = ei, sd = sdi)
  if (alternative == "two.sided") 
    pv <- if (obs <= ei) 
      2 * pv
  else 2 * (1 - pv)
  if (alternative == "greater") 
    pv <- 1 - pv
  list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}



.spe_dist_weight_matrix<-function(spe, squared = TRUE){
  xy <- spatialCoords(spe)
  weight<-as.matrix(dist(xy))
  if(squared == TRUE)
    weight=1/weight^2
  else weight = 1/weight
  diag(weight)<-0
  s <- sum(weight)
  ROWSUM <- rowSums(weight)
  ROWSUM[ROWSUM == 0] <- 1
  return(list(weight=weight,
              ROWSUM = ROWSUM,
              s = s,
              S1 = 0.5 * sum((2*weight)^2),
              S2 = sum((ROWSUM + colSums(weight))^2),
              s.sq = s^2))
}