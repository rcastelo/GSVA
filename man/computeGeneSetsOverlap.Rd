\name{computeGeneSetsOverlap}
\alias{computeGeneSetsOverlap}
\alias{computeGeneSetsOverlap,list,character-method}
\alias{computeGeneSetsOverlap,list,ExpressionSet-method}
\alias{computeGeneSetsOverlap,GeneSetCollection,character-method}
\alias{computeGeneSetsOverlap,GeneSetCollection,ExpressionSet-method}

\encoding{latin1}

\title{
Compute gene-sets overlap
}
\description{
Calculates the overlap among every pair of gene-sets given as input.
}
\usage{
\S4method{computeGeneSetsOverlap}{list,character}(gSets, uniqGenes, min.sz=1, max.sz=Inf)
\S4method{computeGeneSetsOverlap}{list,ExpressionSet}(gSets, uniqGenes, min.sz=1, max.sz=Inf)
\S4method{computeGeneSetsOverlap}{GeneSetCollection,character}(gSets, uniqGenes, min.sz=1, max.sz=Inf)
\S4method{computeGeneSetsOverlap}{GeneSetCollection,ExpressionSet}(gSets, uniqGenes, min.sz=1, max.sz=Inf)
}
\arguments{
  \item{gSets}{Gene sets given either as a \code{list} or a \code{GeneSetCollection} object.}
  \item{uniqGenes}{Vector of unique genes to be considered when calculating the overlaps.}
  \item{min.sz}{Minimum size.}
  \item{max.sz}{Maximum size.}
}

\details{
This function calculates the overlap between every pair of gene sets
of the input argument \code{gSets}. Before this calculation takes place,
the gene sets in \code{gSets} are firstly filtered to discard genes that do not
match to the identifiers in \code{uniqGenes}. Secondly, they are further
filtered to meet the minimum and/or maximum size specified with the arguments
\code{min.sz} and \code{max.sz}. The overlap between two gene sets is calculated
as the number of common genes between the two gene sets divided by the smallest
size of the two gene sets.
}
\value{
A gene-set by gene-set matrix of the overlap among every pair of gene sets.
}
\references{
\enc{Hänzelmann}{Hanzelmann}, S., Castelo, R. and Guinney, J.
GSVA: Gene set variation analysis for microarray and RNA-Seq data.
\emph{BMC Bioinformatics}, 14:7, 2013.
}
\author{J. Guinney}
\seealso{
  \code{\link{filterGeneSets}}
}
\examples{

geneSets <- list(set1=as.character(1:4), set2=as.character(4:10))

computeGeneSetsOverlap(geneSets, unique(unlist(geneSets)))

}
\keyword{Gene set}
