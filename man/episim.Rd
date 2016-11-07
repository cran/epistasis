\name{episim}
\alias{episim}

\title{ Generate discrete ordinal data  }

\description{
Generating discrete ordinal data based on underlying "genome-like" graph structure. The procedure of simulating data relies on a continues variable, which can be simulated from either multivariate normal distribution, or multivariate t-distribution with \code{d} degrees of freedom.   
}

\usage{
episim ( p = 90, n = 200, k = NULL, g = NULL, adjacent = NULL, alpha =
              NULL , beta = NULL, con.dist = "Mnorm", d = NULL, vis = FALSE)
}

\arguments{
	\item{p}{The number of variables. The default value is 90.}
	\item{n}{The number of sample size (observations). The default value is 200.}
	\item{k}{The number of states (categories). The default value is 3.}
	\item{g}{The number of groups (chromosomes) in the graph. The default value is about \eqn{p/20} if \eqn{p >= 40} and 2 if \eqn{p < 40}.}
	\item{adjacent}{ The number of adjacent variable(s) to be linked to a variable. For example, if \code{adjacent = 1} indicates a variable is linked via an edge with its adjacent variable on the left hand side, and its adjacent variable on the right hand side. The \code{adjacent = 2} defines a variable is linked via an edge with its 2 adjacent variables on its left hand side, and 2 adjacent variables on its right hand side. The default value is 1.}
	\item{alpha}{ A probability that a pair of non-adjacent variables in the same group is given an edge. The default value is 0.01.}
	\item{beta}{ A probability that variables in different groups are linked with an edge. The default value is 0.02.}
	\item{con.dist}{ The distribution of underlying continuous variable. If \code{con.dist = "Mnorm"}, a multivariate normal distribution with mean 0 is applied. If \code{con.dist = "Mt"}, the t-distribution with \code{a} degrees of freedom is applied. 
	The default distribution is \code{con.dist = "Mnorm"}.}
	\item{d}{The degrees of freedom of the continuous variable, only applicable when 
	code{con.dist = "Mt"}. The default value is 3.}
	\item{vis}{ Visualize the graph pattern and the adjacency matrix of the true graph structure. The default value is FALSE.}
}

\details{
The graph pattern is generated as below: \cr\cr
"genome-like": The \code{p} variables are evenly partitions variables into \code{g} disjoint groups; the adjacent variables within each group are linked via an edge. With a probability \code{alpha} a pair of non-adjacent variables in the same group is given an edge. Variables in different groups are linked with an edge with a probability of \code{beta}.
}

\value{
An object with S3 class "episim" is returned: \cr \cr
 \item{data}{The generated data as an \code{n} by \code{p} matrix.}
 \item{Theta}{A \code{p} by \code{p} matrix corresponding to the inverse of covariance.}
 \item{adj}{ A \code{p} by \code{p} matrix corresponding to the adjacency matrix of the true graph structure.}
 \item{Sigma}{ A \code{p} by \code{p} covariance matrix for the generated data.}
 \item{n.groups}{ The number of groups.}
 \item{groups}{A vector that indicates each variable belongs to which group.}
 \item{sparsity}{The sparsity levels of the true graph.}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi <pariya.behrouzi@gmail.com>
}

\references{
P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
}

\seealso{
\code{\link{epistasis}}, and \code{\link{epistasis-package}}
}

\examples{
\dontrun{
#genome-like graph structure
data.sim <- episim(alpha = 0.01, beta = 0.02)
#genome-like graph structure with more edges between variables in different groups
data.sim <- episim(alpha = 0.02, beta = 0.03)
}
}