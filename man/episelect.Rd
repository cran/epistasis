\name{episelect}
\alias{episelect}

\title{ Model selection }

\description{
Estimate the optimal regularization parameter at EM convergence based on different information criteria . 
}

\usage{
episelect(epi.object, criteria = NULL, ebic.gamma = 0.5, loglik_Y = FALSE, ncores = NULL)
}

\arguments{
	\item{epi.object}{An object with S3 class "epi"}
	\item{criteria}{Model selection criteria. "ebic" and "aic" are available. BIC model selection can be calculated by fixing \code{ebic.gamma = 0}. }
	\item{ebic.gamma}{The tuning parameter for ebic. The\code{ebic.gamma = 0} results in bic model selection. The default value is 0.5.}
	\item{loglik_Y}{Model selection based on either log-likelihood of observed data (\code{loglik_Y = TRUE}), or the joint log-likelihood of observed and latent variables (\code{loglik_Y = FALSE}).}
	\item{ncores}{ The number of cores to use for the calculations. Using \code{ncores = NULL} automatically detects number of available cores and runs the computations in parallel.}
}

\details{
This function computes extended Bayesian information criteria (ebic), Bayesian information criteria, Akaike information criterion (aic) at EM convergence based on observed or joint log-likelihood.
The observed log-likelihood can be obtained through  

\deqn{\ell_Y(\widehat{\Theta}_\lambda) = Q(\widehat{\Theta}_\lambda | \widehat{\Theta}^{(m)}) - H (\widehat{\Theta}_\lambda | \widehat{\Theta}^{(m)}),}

Where \eqn{Q} can be calculated from \code{\link{epistasis}} function and H function is 
\deqn{H(\widehat{\Theta}_\lambda | \widehat{\Theta}^{(m)}_\lambda) = E_z[\ell_{Z | Y}(\widehat{\Theta}_\lambda) | Y; \widehat{\Theta}_\lambda] = E_z[\log f(z)| Y ;\widehat{\Theta}_\lambda ] - \log p(y).}

The "ebic" and "aic" model selection criteria can be obtained as follow
\deqn{ebic(\lambda) = -2 \ell(\widehat{\Theta}_\lambda) +  ( \log n + 4 \gamma \log p) df(\lambda)}

\deqn{aic(\lambda) = -2 \ell(\widehat{\Theta}_\lambda) +  2 df(\lambda)}
where \eqn{df} refers to the number of non-zeros offdiagonal elements of \eqn{\hat{\Theta}_\lambda}, and \eqn{\gamma \in [0, 1]}. Typical value for for \code{ebic.gamma} is 1/2, but it can also be tuned by experience. Fixing \code{ebic.gamma = 0 } results in bic model selection.
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
3. D. Witten and J. Friedman. New insights and faster computations for the graphical lasso. \emph{Journal of Computational and Graphical Statistics}, to appear, 2011. \cr
4. J. Friedman, T. Hastie and R. Tibshirani. Sparse inverse covariance estimation with the lasso, \emph{Biostatistics}, 2007.\cr
5. Foygel, R. and M. Drton (2010). Extended bayesian information criteria for Gaussian graphical models. In Advances in Neural Information Processing Systems, pp. 604-612.\cr
}

\author{ 
  Pariya Behrouzi and Ernst C.Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com} \cr
}

\seealso{\code{\link{epistasis}}}

\examples{
\dontrun{
#simulate data
D <- episim(p=50, n=100, k= 3, adjacent = 3, alpha = 0.06 , beta = 0.06)
plot(D)

#detect epistatic selection path
out  <-  epistasis(D$data, method="gibbs", n.rho= 5, ncores= 1)

#different graph selection methods
sel.ebic1 <- episelect(out, criteria="ebic")
plot(sel.ebic1)

sel.ebic2 <- episelect(out, criteria="ebic", loglik_Y=TRUE)
plot(sel.ebic2)

sel.aic <- episelect(out, criteria="aic")
plot(sel.aic)

sel.bic <- episelect(out, criteria="ebic", ebic.gamma = 0)
plot(sel.bic)
}
}