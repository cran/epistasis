\name{plot.episelect}
\alias{plot.episelect}

\title{ Plot function for S3 class "episelect" }

\description{
Plot the optimal graph by model selection  
}

\usage{
\method{plot}{episelect}(x, n.markers = NULL, vertex.size=NULL, vertex.label = FALSE, ...)
}

\arguments{
	\item{x}{An object with S3 class "episelect" }
	\item{n.markers}{A vector containing number of variables/markers in each group/chromosome. For example, the CviCol dataset that is provided in the package contains 5 chromosomes/ groups which the total number of markers is \eqn{p = 90}, where the first 24 markers belong into chromosome 1, the next 14 markers into chromosome 2, ..., and chromosome 5 contains 19 markers. Thus, n.mrkr = c(24,14,17,16,19). If \code{n.mrkr = NULL}, in the graph visualization all markers are represented same colour.}
	\item{vertex.size}{ The size of vertices in the graph visualization. The default value is 7.}
	\item{vertex.label}{ Assign names to the vertices. Default is FALSE.}
	\item{\dots}{System reserved (No specific usage)}
}


\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com} \cr
}

\references{
P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
}


\seealso{\code{\link{episelect}}}