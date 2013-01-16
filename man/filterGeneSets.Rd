\name{filterGeneSets}
\alias{filterGeneSets}
\alias{filterGeneSets,list-method}
\alias{filterGeneSets,GeneSetCollection-method}

\encoding{latin1}

\title{
Filter gene sets
}
\description{
Filters gene sets through a given minimum and maximum set size.
}
\usage{
\S4method{filterGeneSets}{list}(gSets, min.sz=1, max.sz=Inf)
\S4method{filterGeneSets}{GeneSetCollection}(gSets, min.sz=1, max.sz=Inf)
}
\arguments{
  \item{gSets}{Gene sets given either as a \code{list} or a \code{GeneSetCollection} object.}
  \item{min.sz}{Minimum size.}
  \item{max.sz}{Maximum size.}
}

\details{
This function filters the input gene sets according to a given minimum and maximum set size.
}
\value{
A collection of gene sets that meet the given minimum and maximum set size.
}
\references{
\enc{Hänzelmann}{Hanzelmann}, S., Castelo, R. and Guinney, J.
GSVA: Gene set variation analysis for microarray and RNA-Seq data.
\emph{BMC Bioinformatics}, 14:7, 2013.
}
\author{J. Guinney}
\seealso{
  \code{\link{computeGeneSetsOverlap}}
}
\examples{

geneSets <- list(set1=as.character(1:4), set2=as.character(4:10))

filterGeneSets(geneSets, min.sz=5)

}
\keyword{Gene set}
