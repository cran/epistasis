\name{plot.episim}
\alias{plot.episim}

\title{ Plot function for \code{S3} class \code{"episim"} }

\description{ Visualizes the pattern of the true graph, the adjacency matrix, precison matrix and the covariance matrix of the simulated data.}

\usage{
\method{plot}{episim}(x, layout = layout.fruchterman.reingold, ...)
}

\arguments{
	\item{x}{An object of \code{S3} class \code{"episim"}, from function \code{\link{episim}}. }
	\item{layout}{The default is \code{"layout.fruchterman.reingold"}.}
	\item{\dots}{ System reserved (No specific usage)}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com} \cr
}

\references{
P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
}

\seealso{\code{\link{episim}}}

\examples{
\dontrun{
# Generating discrete ordinal data with "genome-like" graph structure
data.sim <- episim(alpha = 0.01, beta = 0.02)
plot( data.sim )
}
}

