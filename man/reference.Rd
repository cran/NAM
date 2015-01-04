\name{reference}
\alias{reference}
\title{
Changing the Reference Genotype
}
\description{
Function changes the reference genotype. For NAM populations, this function must be used when genotypes are coded according to the reference genome intead of the standard parent.
}
\usage{
reference(gen, ref = NULL)
}
\arguments{
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}.
}
  \item{ref}{
Numeric vector of length \eqn{n} with elements coded as \code{0, 1, 2}, it represents the genotypic information of a new reference genotype. Default assumes that more frequent allele represents the reference genome.
}
}
\details{
If genotypes are coded based on the refence genome, NAM analysis are optimized by using the standard parent as reference to allele coding.
}
\value{
Returns a recoded \eqn{gen} matrix
}
\author{
Alencar Xavier
}
\examples{
data(tpod)
gen=reference(gen=gen,ref=NULL)
}
