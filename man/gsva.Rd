\name{gsva}
\alias{gsva}
\alias{gsva,ExpressionSet,list-method}
\alias{gsva,ExpressionSet,GeneSetCollection-method}
\alias{gsva,matrix,list-method}

\title{
Gene Set Variation Analysis
}
\description{
Estimates GSVA enrichment scores.
}
\usage{
\S4method{gsva}{ExpressionSet,list}(expr, gset.idx.list,
    abs.ranking=FALSE,
    min.sz=1,
    max.sz=Inf,
    no.bootstraps=0,
    bootstrap.percent = .632,
    parallel.sz=0,
    parallel.type="SOCK",
    verbose=TRUE,
    mx.diff=TRUE)
\S4method{gsva}{ExpressionSet,GeneSetCollection}(expr, gset.idx.list,
    abs.ranking=FALSE,
    min.sz=1,
    max.sz=Inf,
    no.bootstraps=0,
    bootstrap.percent = .632,
    parallel.sz=0,
    parallel.type="SOCK",
    verbose=TRUE,
    mx.diff=TRUE)
\S4method{gsva}{matrix,list}(expr, gset.idx.list,
    abs.ranking=FALSE,
    min.sz=1,
    max.sz=Inf,
    no.bootstraps=0,
    bootstrap.percent = .632,
    parallel.sz=0,
    parallel.type="SOCK",
    verbose=TRUE,
    mx.diff=TRUE)

}
\arguments{
  \item{expr}{Gene expression data which can be given either as an \code{ExpressionSet}
              object or as a matrix of expression values where rows correspond
              to genes and columns correspond to samples.}
  \item{gset.idx.list}{Gene sets provided either as a \code{list} object or as a
                       \code{GeneSetCollection} object.}
  \item{abs.ranking}{Flag to determine whether genes should be ranked according to 
  					their sign (flag=FALSE) or by absolute value (flag=TRUE). 
  					In the latter, pathways with genes enriched on either extreme
  					(high or low) will be regarded as 'highly' activated. }
  \item{min.sz}{Minimum size of the resulting gene sets.}
  \item{max.sz}{Maximum size of the resulting gene sets.}
  \item{no.bootstraps}{Number of bootstrap iterations to perform.}
  \item{bootstrap.percent}{.632 is the ideal percent samples bootstrapped.}
  \item{parallel.sz}{Number of processors to use when doing the calculations in parallel.
                     This requires to previously load either the \code{multicore} or the
                     \code{snow} library. If \code{multicore} is loaded and this argument
                     is left with its default value (\code{parallel.sz=0}) then it will use
                     all available core processors unless we set this argument with a
                     smaller number. If \code{snow} is loaded then we must set this argument
                     to a positive integer number that specifies the number of processors to
                     employ in the parallel calculation.}
  \item{parallel.type}{Type of cluster architecture when using \code{snow}.}
  \item{verbose}{Gives information about each calculation step. Default: \code{FALSE}.}
  \item{mx.diff}{Offers two approaches to calculate the enrichment statistic (ES)
                 from the KS random walk statistic. \code{mx.diff=FALSE}: ES is calculated as
                 the maximum distance of the random walk from 0. \code{mx.diff=TRUE} (default): ES
                 is calculated as the magnitude difference between the largest positive
                 and negative random walk deviations.
}
}

\details{
GSVA assesses the relative enrichment of gene sets across samples using
a non-parametric approach.  Conceptually, GSVA transforms a p-gene by n-sample
gene expression matrix into a g-geneset by n-sample pathway enrichment matrix.
This facilitates many forms of statistical analysis in the 'space' of pathways
rather than genes, providing a higher level of interpretability.

The \code{gsva()} function first maps the identifiers in the gene sets to the
identifiers in the input expression data leading to a filtered collection of
gene sets. This collection can be further filtered to require a minimun and/or
maximum size of the gene sets for which we want to calculate GSVA enrichment
scores, by using the arguments \code{min.sz} and \code{max.sz}.
}
\value{
A gene-set by sample matrix of GSVA enrichment scores.
}
\references{
H\"anzelmann, S., Castelo, R. and Guinney, J.
GSVA: Gene Set Variation Analysis, \emph{submitted}
}
\author{J. Guinney}
\seealso{
  \code{\link{filterGeneSets}}
  \code{\link{computeGeneSetsOverlap}}
}
\examples{

library(limma)

p <- 10 ## number of genes
n <- 30 ## number of samples
nGrp1 <- 15 ## number of samples in group 1
nGrp2 <- n - nGrp1 ## number of samples in group 2

## consider three disjoint gene sets
geneSets <- list(set1=paste("g", 1:3, sep=""),
                 set2=paste("g", 4:6, sep=""),
                 set3=paste("g", 7:10, sep=""))

## sample data from a normal distribution with mean 0 and st.dev. 1
y <- matrix(rnorm(n*p), nrow=p, ncol=n,
            dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

## genes in set1 are expressed at higher levels in the last 10 samples
y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

## build design matrix
design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nGrp1), rep(1, nGrp2)))

## fit linear model
fit <- lmFit(y, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

## genes in set1 are differentially expressed
topTable(fit, coef="sampleGroup2vs1")

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(y, geneSets, mx.diff=1)$es.obs

## fit the same linear model now to the GSVA enrichment scores
fit <- lmFit(gsva_es, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

## set1 is differentially expressed
topTable(fit, coef="sampleGroup2vs1")

}
\keyword{Pathway variation}
