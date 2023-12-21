\name{GSVA-pkg-defunct}
\alias{GSVA-pkg-defunct}
\alias{gsva,missing-method}
\alias{gsva,matrix-method}
\alias{gsva,ExpressionSet-method}
\alias{gsva,SummarizedExperiment-method}
\alias{gsva,dgCMatrix-method}
\alias{gsva,SingleCellExperiment-method}
\title{Defunct functions in package \code{GSVA}.}
\description{
The functions listed below are defunct and will be removed
in the next release.
}
\section{Details}{

Instead of \code{gsva(expr=., gset.idx.list=., method=., ...)}, use a
method-specific parameter object,
see \code{\link{plageParam}} \code{\link{zscoreParam}} \code{\link{ssgseaParam}} \code{\link{gsvaParam}},
followed by a call to the new \code{gsva()} function, see \code{\link{gsva}}.
}

\keyword{internal}
