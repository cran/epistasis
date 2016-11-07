\name{cutoffs}
\alias{cutoffs}

\title{
Cut-points 
}

\description{
Calculates cut-points of ordinal variables with respect to the Gaussian copula. 
}

\usage{
cutoffs(y)
}

\arguments{
  \item{y}{
An (\eqn{n \times p}) matrix or a \code{data.frame} corresponding to the data matrix (\eqn{n} is the sample size and \eqn{p} is the number of variables). It also could be an object of class \code{"episim"}.
}
}

\details{ The relationship between \eqn{j}th variable and \eqn{j}th latent variable is expressed through this set of cut-points.
}

\value{
\item{cutoffs}{ A \eqn{p} by \eqn{(k + 1)} matrix representing the cut-point values under the Gaussian copula, where \eqn{k} defines the number of states in the dataset.}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi <pariya.behrouzi@gmail.com>
}

\references{
P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
}

\seealso{
  \code{\link{lower.upper}} and \code{\link{epistasis-package}}.
}

\examples{
\dontrun{
D <- episim(p = 90, n = 300, k = 3)

cutoff(D$data)
}
}
