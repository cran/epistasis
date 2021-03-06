\name{lower.upper}
\alias{lower.upper}

\title{
Calculates lower band and upper band
}

\description{
Calculates lower and upper bands for each data point, using a set of cut-points which is obtained from the Gaussian copula.
}

\usage{
lower.upper(y)
}

\arguments{
	\item{y}{ An (\eqn{n \times p}) matrix or a \code{data.frame} corresponding to the data matrix (\eqn{n} is the sample size and \eqn{p} is the number of variables). It also could be an object of class \code{"episim"}.
}
}

\value{
\item{lower}{
A \eqn{n} by \eqn{p} matrix representing the lower band for each data point.
}
\item{upper}{
A \eqn{n} by \eqn{p} matrix representing the upper band for each data point.
}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi <pariya.behrouzi@gmail.com>
}

\references{
P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
}

\seealso{
  \code{\link{cutoffs}} and \code{\link{epistasis-package}}.
}

\examples{
\dontrun{
D <- episim(p = 100, n = 50, k = 3)
lower.upper(D$data)
}
}