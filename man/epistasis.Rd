\name{epistasis}
\alias{epistasis}

\title{ Detecting epistatic selection from multi-loci genotype data }

\description{
This is the main function of the \pkg{epistasis} package. Two methods are available to detect epistatic selection, including (1) approximation method, and (2) gibbs sampling within the Gaussian copula graphical model. Both methods are able to deal with missing genotypes.
}

\usage{
epistasis(y, method = "approx", rho = NULL, n.rho = NULL, rho.ratio = NULL, 
           ncores = NULL, em.iter = 10, em.tol = .001, verbose = TRUE)
}

\arguments{
	\item{y}{ An (\eqn{n \times p}) matrix or a \code{data.frame} corresponding to the data matrix (\eqn{n} is the sample size and \eqn{p} is the number of variables). It also could be an object of class \code{"episim"}. Input data can contain missing values.}
	\item{method}{ Detecting epistatic selection with two methods: "gibbs" and "approx". The default method is "approx".}
	\item{rho}{ A decreasing sequence of non-negative numbers that control the sparsity level. Leaving the input as \code{rho = NULL}, the program automatically computes a sequence of \code{rho} based on \code{n.rho} and \code{rho.ratio}. Users can also supply a decreasing sequence values to override this.}
	\item{n.rho}{ The number of regularization parameters. The default value is \code{10}.}
	\item{rho.ratio}{ Determines distance between the elements of \code{rho} sequence. A small value of \code{rho.ratio} results in a large distance between the elements of \code{rho} sequence. And a large value of \code{rho.ratio} results into a small distance between elements of \code{rho}. The default value is \code{0.3}.}
	\item{ncores}{ The number of cores to use for the calculations. Using \code{ncores = NULL} automatically detects number of available cores and runs the computations in parallel on (available cores - 1).}
	\item{em.tol}{ A criteria to stop the EM iterations. The default value is .001. }
	\item{em.iter}{ The number of EM iterations. The default value is 10.}
	\item{verbose}{Providing a detail message for tracing output. The default value is \code{TRUE}.}
}

\details{
Viability is the phenotype that is considered. This function detects the conditional dependent short- and long-range linkage disequilibrium structure of genomes and thus reveals aberrant marker-marker associations that are due to epistatic selection.
This function can be used to estimate conditional independence relationships between partially observed data that not follow Gaussianity assumption (e.g. continuous non-Gaussian, discrete, or mixed dataset).
}

\value{
An object with S3 class \code{"epi"} is returned:  
 \item{Theta}{ A list of estimated p by p precision matrices corresponding to \code{rho}.}
 \item{path}{ A list of estimated p by p adjacency matrices. This is the graph path corresponding to rho.}
 \item{Sigma}{ A list of estimated p by p covariance matrices corresponding to \code{rho}.}
 \item{ES}{ A list of estimated p by p conditional expectation corresponding to \code{rho}.}
 \item{Z}{ A list of n by p transformed data based on Gaussian copula. }
 \item{rho}{ A \code{n.rho} dimensional vector containing the penalty terms.}
 \item{loglik }{ A \code{n.rho} dimensional vector containing the maximized log-likelihood values along the graph path. }
 \item{data}{ The \eqn{n} by \eqn{p} input data matrix. }
}

\note{ This function estimates the graph path . To select an optimal graph please refer to \code{\link{episelect}}.}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainers: Pariya Behrouzi \email{pariya.behrouzi@gmail.com}\cr
}

\references{
1. P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
2. D. Witten and J. Friedman. New insights and faster computations for the graphical lasso. \emph{Journal of Computational and Graphical Statistics}, to appear, 2011. \cr
3. J. Friedman, T. Hastie and R. Tibshirani. Sparse inverse covariance estimation with the lasso, \emph{Biostatistics}, 2007.\cr
4. Guo, Jian, et al. "Graphical models for ordinal data." Journal of Computational and Graphical Statistics 24.1 (2015): 183-204. \cr
}

\seealso{\code{\link{episelect}}}

\examples{
\dontrun{
#simulate data
D <- episim(p=50, n=100, k= 3, adjacent = 3, alpha = 0.06 , beta = 0.06)
plot(D)

#epistasis path estimation using approx
out1 <-  epistasis(D$data, method="approx", n.rho= 5)
plot(out1)

#epistasis path estimation using gibbs
out2  <-  epistasis(D$data, method="gibbs", n.rho= 5, ncores= 1)
plot(out2)
}
}

