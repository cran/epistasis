\name{B73Ki11}
\alias{B73Ki11}
\docType{data}
\title{A family from Nested Association Mapping (NAM) in maize }
\description{ 
The genotype data for family B73 x ki11 Recombinant Inbred Line (RIL) from NAM population.
}

\usage{data(B73Ki11)}
\format{
  The format is a list containing two matrices. 1. the data matrix with 1106 single-nucleotide polymorphism (SNP) markers for 191 individuals. 2. Information 
  about the SNP markers regarding their location in the genome.
  
}  

\details{
The Nested Association Mapping (NAM) initiative in maize populations is designed to reveal the genetic structure of underlying complex traits in maize. As part of this study, an
inbred Ki11 maize line was crossed with the B73 reference line. This genotype data contains 1106 markers genotyped for 193 individuals. The B73 x Ki11 RIL
is a diploid population with three possible genotypes, k = 3. This data set can be used to detect epistatic selection, short- and long- range linkage disequilibrium between 1106 SNP markers. 
}

\source{
It is publicly available at http://datacommons.cyverse.org/browse/iplant/home/shared/panzea/genotypes
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com}\cr
}

\references{
1. P. Behrouzi and E. C. Wit. Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. \emph{arXiv}, 2016.\cr
2. McMullen, M. D., S. Kresovich, H. S. Villeda, P. Bradbury, H. Li, Q. Sun, S. Flint-Garcia, J. Thornsberry, C. Acharya, C. Bottoms, et al. (2009). Genetic properties of the maize nested association mapping population. Science 325 (5941), 737-740.
3. Rodgers-Melnick, E., P. J. Bradbury, R. J. Elshire, J. C. Glaubitz, C. B. Acharya, S. E. Mitchell, C. Li, Y. Li, and E. S. Buckler (2015). Recombination in diverse maize is stable, predictable, and associated with genetic load. Proceedings of the National Academy of Sciences 112 (12), 3823-3828.
}

\examples{
data(B73Ki11)
image(B73Ki11$data, ylab="Markers" , xlab= "Individuals") 
B73Ki11$info
}

\keyword{datasets}
