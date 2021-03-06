\name{CviCol}
\alias{CviCol}
\docType{data}
\title{Arabidopsis thaliana genotype data }
\description{ 
  The genotype data of the Cvi-0 x Col-0 Recombinant Inbred Line (RIL) population.
}

\usage{data(CviCol)}
\format{
  The format is a matrix containing 90 single-nucleotide polymorphism (SNP) markers for 367 individuals.
}  

\details{
The \emph{Arabidopsis thaliana} genotype data is derived from a RIL cross between Columbia-0 (Col-0) and the Cape Verde Island (Cvi-0), where 367 individuals were genotyped for 90 genetic markers. This is a diploid population with three possible genotpe states \code{(k = 3)}, where the genotypes coded as {0, 1, 2}, where 0 and 2 represent the homozygous genotypes and 1 defines the heterozygous genotype. \cr
This data set can be used to detect epistatic selection, short- and long- range linkage disequilibrium between 90 SNP markers. 
}

\source{
It is publicly available at http://publiclines.versailles.inra.fr/page/8
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com}\cr
}

\references{
1. P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
2. Simon, M., et al. "QTL mapping in five new large RIL populations of Arabidopsis thaliana genotyped with consensus SNP markers." Genetics 178 (2008): 2253-2264.
}

\examples{
\dontrun{
data(CviCol)

# Graph path
epi <- epistasis(CviCol, method="approx", rho.ratio=0.2)
epi
plot(epi)

sel <- episelect(epi)
plot(sel)
}
}

\keyword{datasets}
