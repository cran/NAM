\name{Duplicates}
\alias{cleanREP}
\title{
Manage Duplicated Genotypes
}
\description{
Identify and merge duplicate genotypes
}
\usage{
cleanREP(y,gen,fam=NULL,thr=0.95)
}
\arguments{
  \item{y}{
Numeric vector (\eqn{n}) or numeric matrix (\eqn{n} x \eqn{t}) of observations describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n} rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozygous, heterozygous and reference homozigous. \code{NA} not allowed.
}
  \item{fam}{
Numeric vector of length (\eqn{n}) indicating which subpopulations (\eqn{i.e.} family) each observation comes from. Default assumes that all observations are from the same populations.
}
  \item{thr}{
Threshold above which genotypes are considered identical. Default is 0.95, merging genotypes >95 percent identical.
}
}
\details{
The algorithm start from generating a genomic identity matrix (IBS), with pairwise percentage of identical loci among individuals. Individuals above the threshold have the phenotypes merged while keeping only one genotype.
}
\value{
List containing the inputs without replicates. Groups of replicates are replaced by a single observation with the phenotypic expected value. The algorithm keeps the genotypic information of the first individual (genotypic matrix order).
}
\author{
Alencar Xavier
}
\examples{
data(tpod)
test = cleanREP(y,gen)
}
\keyword{ NAM }
\keyword{ GWAS }
\keyword{ BLUP }
