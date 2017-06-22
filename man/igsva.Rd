\name{igsva}
\alias{igsva}

\encoding{latin1}

\title{
Gene Set Variation Analysis
}
\description{
Starts an interactive GSVA shiny web app.
}
\usage{
igsva()
}
\details{
GSVA assesses the relative enrichment of gene sets across samples using
a non-parametric approach.  Conceptually, GSVA transforms a p-gene by n-sample
gene expression matrix into a g-geneset by n-sample pathway enrichment matrix.
This facilitates many forms of statistical analysis in the 'space' of pathways
rather than genes, providing a higher level of interpretability.

The \code{igsva()} function starts an interactive shiny web app that allows the
user to configure the arguments of the \code{\link{gsva}()} function and runs it
on the computer. Please see the manual page of the \code{\link{gsva}()} function
for a description of the arguments and their default and alternative values.

The input data may be loaded from the users workspace or by selecting a CSV file
for the expression data, and a GMT file for the gene sets data.
}
\value{
A gene-set by sample matrix of GSVA enrichment scores after
pressing the button 'Save & Close'. This result can be also downloaded as a
CSV file with the 'Download' button.
}
\references{
\enc{Hänzelmann}{Hanzelmann}, S., Castelo, R. and Guinney, J.
GSVA: Gene set variation analysis for microarray and RNA-Seq data.
\emph{BMC Bioinformatics}, 14:7, 2013.
}
\author{J. \enc{Fernández}{Fernandez} and R. Castelo}
\seealso{
  \code{\link{gsva}}
}
\examples{

\dontrun{
res <- igsva() ## this will open your browser with the GSVA shiny web app
}

}
\keyword{GSVA}
\keyword{shiny}
