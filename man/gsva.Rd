\name{gsva}
\alias{gsva}
\alias{gsva,ExpressionSet,list,missing-method}
\alias{gsva,ExpressionSet,GeneSetCollection,missing-method}
\alias{gsva,matrix,GeneSetCollection,character-method}
\alias{gsva,matrix,list,missing-method}

\title{
Gene Set Variation Analysis
}
\description{
Estimates GSVA enrichment scores.
}
\usage{
\S4method{gsva}{ExpressionSet,list,missing}(expr, gset.idx.list, annotation=NULL,
    method=c("gsva", "ssgsea", "zscore", "plage"),
    rnaseq=FALSE,
    abs.ranking=FALSE,
    min.sz=1,
    max.sz=Inf,
    no.bootstraps=0,
    bootstrap.percent = .632,
    parallel.sz=0,
    parallel.type="SOCK",
    mx.diff=TRUE,
    tau=switch(method, gsva=1, ssgsea=0.25, NA),
    kernel=TRUE,
    verbose=TRUE)
\S4method{gsva}{ExpressionSet,GeneSetCollection,missing}(expr, gset.idx.list, annotation=NULL,
    method=c("gsva", "ssgsea", "zscore", "plage"),
    rnaseq=FALSE,
    abs.ranking=FALSE,
    min.sz=1,
    max.sz=Inf,
    no.bootstraps=0,
    bootstrap.percent = .632,
    parallel.sz=0,
    parallel.type="SOCK",
    mx.diff=TRUE,
    tau=switch(method, gsva=1, ssgsea=0.25, NA),
    kernel=TRUE,
    verbose=TRUE)
\S4method{gsva}{matrix,GeneSetCollection,character}(expr, gset.idx.list, annotation=NA,
    method=c("gsva", "ssgsea", "zscore", "plage"),
    rnaseq=FALSE,
    abs.ranking=FALSE,
    min.sz=1,
    max.sz=Inf,
    no.bootstraps=0,
    bootstrap.percent = .632,
    parallel.sz=0,
    parallel.type="SOCK",
    mx.diff=TRUE,
    tau=switch(method, gsva=1, ssgsea=0.25, NA),
    kernel=TRUE,
    verbose=TRUE)
\S4method{gsva}{matrix,list,missing}(expr, gset.idx.list, annotation=NULL,
    method=c("gsva", "ssgsea", "zscore", "plage"),
    rnaseq=FALSE,
    abs.ranking=FALSE,
    min.sz=1,
    max.sz=Inf,
    no.bootstraps=0,
    bootstrap.percent = .632,
    parallel.sz=0,
    parallel.type="SOCK",
    mx.diff=TRUE,
    tau=switch(method, gsva=1, ssgsea=0.25, NA),
    kernel=TRUE,
    verbose=TRUE)
}
\arguments{
  \item{expr}{Gene expression data which can be given either as an \code{ExpressionSet}
              object or as a matrix of expression values where rows correspond
              to genes and columns correspond to samples.}
  \item{gset.idx.list}{Gene sets provided either as a \code{list} object or as a
                       \code{GeneSetCollection} object.}
  \item{annotation}{In the case of calling \code{gsva()} with expression data in a \code{matrix}
                    and gene sets as a \code{GeneSetCollection} object, the \code{annotation} argument
                    can be used to supply the name of the Bioconductor package that contains
                    annotations for the class of gene identifiers occurring in the row names of
                    the expression data matrix. By default \code{annotation=NULL} which implies that
                    \code{gsva()} will try to match the identifiers in \code{expr} to the identifiers
                    in \code{gset.idx.list} just as they are.}
  \item{method}{Method to employ in the estimation of gene-set enrichment scores per sample. By default
                this is set to \code{gsva} (H\"anzelmann et al, submitted) and other options are
                \code{ssgsea} (Barbie et al, 2009), \code{zscore} (Lee et al, 2008) or \code{plage}
                (Tomfohr et al, 2005). The latter two standardize first expression profiles into z-scores
                over the samples and, in the case of \code{zscore}, it combines them together as their sum
                divided by the square-root of the size of the gene set,
                while in the case of \code{plage} they are used to calculate the singular value decomposition
                (SVD) over the genes in the gene set and use the coefficients of the first right-singular vector
                as pathway activity profile.}
  \item{rnaseq}{Flag to inform whether the input gene expression data comes from microarray
                (\code{rnaseq=FALSE}, default) or RNA-Seq (\code{rnaseq=TRUE}) experiments.}
  \item{abs.ranking}{Flag to determine whether genes should be ranked according to 
  					their sign (\code{abs.ranking=FALSE}) or by absolute value (\code{abs.ranking=TRUE}). 
  					In the latter, pathways with genes enriched on either extreme
  					(high or low) will be regarded as 'highly' activated. }
  \item{min.sz}{Minimum size of the resulting gene sets.}
  \item{max.sz}{Maximum size of the resulting gene sets.}
  \item{no.bootstraps}{Number of bootstrap iterations to perform.}
  \item{bootstrap.percent}{.632 is the ideal percent samples bootstrapped.}
  \item{parallel.sz}{Number of processors to use when doing the calculations in parallel.
                     This requires to previously load either the \code{parallel} or the
                     \code{snow} library. If \code{parallel} is loaded and this argument
                     is left with its default value (\code{parallel.sz=0}) then it will use
                     all available core processors unless we set this argument with a
                     smaller number. If \code{snow} is loaded then we must set this argument
                     to a positive integer number that specifies the number of processors to
                     employ in the parallel calculation.}
  \item{parallel.type}{Type of cluster architecture when using \code{snow}.}
  \item{mx.diff}{Offers two approaches to calculate the enrichment statistic (ES)
                 from the KS random walk statistic. \code{mx.diff=FALSE}: ES is calculated as
                 the maximum distance of the random walk from 0. \code{mx.diff=TRUE} (default): ES
                 is calculated as the magnitude difference between the largest positive
                 and negative random walk deviations.}
  \item{tau}{Exponent defining the weight of the tail in the random walk performed by both the \code{gsva}
             (H\"anzelmann et al., submitted) and the \code{ssgsea} (Barbie et al., 2009) methods. By default,
             this \code{tau=1} when \code{method="gsva"} and \code{tau=0.25} when \code{method="ssgsea"} just
             as specified by Barbie et al. (2009) where this parameter is called \code{alpha}.}
  \item{kernel}{Logical, set to \code{TRUE} when the GSVA method employes a kernel non-parametric
                estimation of the empirical cumulative distribution function (default) and \code{FALSE}
                when this function is directly estimated from the observed data. This last option is
                justified in the limit of the size of the sample by the so-called Glivenko-Cantelli theorem.}
  \item{verbose}{Gives information about each calculation step. Default: \code{FALSE}.}
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
Barbie, D.A. et al. Systematic RNA interference reveals that oncogenic KRAS-driven
cancers require TBK1. \emph{Nature}, 462(5):108-112, 2009.

H\"anzelmann, S., Castelo, R. and Guinney, J.
GSVA: Gene Set Variation Analysis, \emph{submitted}

Lee, E. et al. Inferring pathway activity toward precise disease classification.
\emph{PLoS Comp Biol}, 4(11):e1000217, 2008.

Tomfohr, J. et al. Pathway level analysis of gene expression using singular value decomposition.
\emph{BMC Bioinformatics}, 6:225, 2005.
}
\author{J. Guinney and R. Castelo}
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
