\name{snpH2}
\alias{snpH2}
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
Forneris, N. et al (2014). Quality Control of Genotypes Using Heritability Estimates of Gene Content. In 10th World Congress on Genetics Applied to Livestock Production. Asas.
}
\author{

Alencar Xavier

}
\examples{

data(tpod)
Heritability=snpH2(gen)
plot(Heritability)

}