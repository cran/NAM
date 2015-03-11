\name{Fst}
\alias{Fst}
\title{
  Fixation Index
}
\description{
Genetic variation associated with markers distributed among subpopulations. The function generates a plot for structure diagnosis.
}
\usage{
  Fst(gen,fam)
}
\arguments{
  \item{gen}{
    Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozigous, heterozygous and reference homozygous. \code{NA} is allowed.
  }
  \item{fam}{
    Numeric vector of length (\eqn{n}) indicating which subpopulations (\eqn{i.e.} family) each observation comes from. \code{NA} is not allowed.
  }
}
\details{
Fst (Wright 1943) represents the differentiation among populations for the given locus. It is also a measure of genetic divergence between admixed populations. Neutral markers have an expected Fst 0.05, but may change according to the population.
Outlier Fst supports signatures of selection. Permutation threshold is commonly used to verify statistical significance.
}
\value{
  Wright's FST.
}
\references{
Wright, S. (1943). Isolation by distance. Genetics, 28(2), 114.
}
\author{
  Alencar Xavier and William Muir
}
\examples{
  data(tpod)
  Fst(gen=gen,fam=fam)
}
\keyword{ Fst }