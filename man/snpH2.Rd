\name{SNP H2}
\alias{snpH2}
\alias{plot.H2}
\title{
SNP heritability
}
\description{
Calculates the ability of markers to carry a gene.
}
\usage{
snpH2(gen)
}
\arguments{
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n} rows of observations and (\eqn{m}) columns of molecular markers. Missing values not allowed.
}
}
\value{
Numeric vector containing the heritability of each markers. 
}
\references{
Forneris, N. S., Legarra, A., Vitezica, Z. G., Tsuruta, S., Aguilar, I., Misztal, I., & Cantet, R. J. (2015). Quality Control of Genotypes Using Heritability Estimates of Gene Content at the Marker. Genetics, genetics-114.
}
\author{

Alencar Xavier

}
\examples{

data(tpod)
Heritability=snpH2(gen)
plot(Heritability,chr=chr)

}