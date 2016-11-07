\name{plot.epi}
\alias{plot.epi}

\title{ plot for S3 class "epi" }

\description{
Plot the graph path which is the output of two functions \code{\link{epistasis}}.
}

\usage{
\method{plot}{epi}( x, n.markers=NULL , ... )
}

\arguments{
	\item{x}{An object from "epi" class.}
	\item{n.markers}{A vector containing number of variables/markers in each group/chromosome. For example, the CviCol dataset that is provided in the package contains 5 chromosomes/ groups which the total number of markers is \eqn{p = 90}, where the first 24 markers belong into chromosome 1, the next 14 markers into chromosome 2, ..., and chromosome 5 contains 19 markers. Thus, n.mrkr = c(24,14,17,16,19). If \code{n.mrkr = NULL}, in the graph visualization all markers are represented same colour.}
	\item{\dots}{System reserved (No specific usage)}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}

\references{
P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
}


\seealso{\code{\link{epistasis}}.}

