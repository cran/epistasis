\name{episelect}
\alias{episelect}

\title{ Model selection }

\description{
Estimate the optimal regularization parameter based on Extended Bayesian information criterion (ebic) at EM convergence with respect to either observed log-likelihood, or joint log-likelihood of latent and observed variables. 
}

\usage{
episelect( epi.object, ebic.gamma = 0.5, EM.select = FALSE, ncores = NULL )
}

\arguments{
	\item{epi.object}{An object with S3 class "epi"}
	\item{ebic.gamma}{The tuning parameter for ebic. The\code{ebic.gamma = 0} results in bic model selection. The default value is 0.5.}
	\item{EM.select}{Model selection criterion using ebic at the EM convergence based on either the log-likelihood of observed data (\code{EM.select = TRUE}), or the joint log-likelihood of observed and latent variables (\code{EM.select = FALSE}).}
	\item{ncores}{ The number of cores to use for the calculations. Using \code{ncores = NULL} automatically detects number of available cores and runs the computations in parallel.}
}

\details{
This function computes extended Bayesian information criteria (ebic) at EM convergence based on the observed log-likelihood, which can be written as 

\deqn{\ell_Y(\widehat{\Theta}_\lambda) = Q(\widehat{\Theta}_\lambda | \widehat{\Theta}^{(m)}) - H (\widehat{\Theta}_\lambda | \widehat{\Theta}^{(m)}),}

Where \eqn{Q} can be obtained either from \code{\link{epistasis}} and H function is 
\deqn{H(\widehat{\Theta}_\lambda | \widehat{\Theta}^{(m)}_\lambda) = E_z[\ell_{Z | Y}(\widehat{\Theta}_\lambda) | Y; \widehat{\Theta}_\lambda] = E_z[\log f(z)| Y ;\widehat{\Theta}_\lambda ] - \log p(y).}

As we are interested in graph selection, we use ebic
\deqn{ebic(\lambda) = -2 \ell(\widehat{\Theta}_\lambda) +  ( \log n + 4 \gamma \log p) df(\lambda)}

where \eqn{df} refers to the number of non-zeros offdiagonal elements of \eqn{\hat{\Theta}_\lambda}, and \eqn{\gamma \in [0, 1]}. Typical value for for ebic.gamma is 1/2, but it can also be tuned by experience.
}


\value{
An object with S3 class "episelect" is returned:
 \item{opt.path}{ The optimal graph selected from the graph path}
 \item{opt.theta}{The optimal precision matrix from the graph path}
 \item{opt.Sigma}{ The optimal covariance matrix from the graph path}
 \item{ebic.scores}{ Extended BIC scores for regularization parameter selection at the EM convergence.}
 \item{opt.index}{ The index of optimal regularization parameter. }
 \item{opt.rho}{ The selected regularization parameter. }
 
and anything else that is included in the input \code{epi.object}.
 }

\references{
1. P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
2. Ibrahim, Joseph G., Hongtu Zhu, and Niansheng Tang. "Model selection criteria for missing-data problems using the EM algorithm." Journal of the American Statistical Association (2012).
3.  D. Witten and J. Friedman. New insights and faster computations for the graphical lasso. \emph{Journal of Computational and Graphical Statistics}, to appear, 2011. \cr
4. J. Friedman, T. Hastie and R. Tibshirani. Sparse inverse covariance estimation with the lasso, \emph{Biostatistics}, 2007.\cr
}

\author{ 
  Pariya Behrouzi and Ernst C.Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com} \cr
}

\seealso{\code{\link{epistasis}}}

\examples{
\dontrun{
#generate data
D <- episim(p=20, n=10, alpha = 0.01, beta = 0.2)
#detect epistatic selection path
out <- epistasis(D$data)
#select best epistatic selection network
sel <- episelect(out)
plot(sel)
}
}