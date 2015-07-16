\name{Report}
\alias{report}
\title{
GWAS Report
}
\description{
Run and export GWAS
}
\usage{
report(y,gen,fam=NULL,chr=NULL,exp=TRUE,reml=FALSE)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n} rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozygous, heterozygous and reference homozigous. \code{NA} is allowed.
}
  \item{fam}{
Numeric vector of length (\eqn{n}) indicating which subpopulations (\eqn{i.e.} family) each observation comes from. Default assumes that all observations are from the same populations.
}
  \item{chr}{
Numeric vector indicating the number of markers in each chromosome. The sum of \eqn{chr} must be equal to the number of columns in \eqn{gen}. Default assumes that all markers are from the same chromosome.
}
  \item{exp}{
If TRUE, it exports the GWAS report and Manhattan plot to the current directory R is set.
}
  \item{reml}{
If TRUE, it finds variance components and breeding values using REML. If false, it uses the Gibbs sampling.
}
}
\details{
This function was dedicated to provide simple output containing SNP scores and effect, exported as CSV file. Association are performed using the gwas2 function, marginal SNP effects are calculated through linear regression, joint SNP effects and variance components are computed eight by REML or Gibbs sampling.
}
\value{
It exports a CSV file with the report file and a PNG file with a Manhattan plot. In R environment, it returns a list with three objects. 1) An association study (GWAS), output of the gwas function. 2) Fitted mixed model (MixMod), output of gibbs (or reml) function, used to estimate variance components and joint effect of markers. 3) A data frame with SNP score and additive effect (Report).
}
\author{
Alencar Xavier
}
\examples{
data(tpod)
test = report(y=y,gen=gen[,1:240],fam=fam,chr=chr[1:12],exp=FALSE,reml=TRUE)
}
\keyword{ NAM }
\keyword{ GWAS }
\keyword{ BLUP }
