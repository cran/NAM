\name{GWAS}
\alias{gwas}
\alias{gwas2}
\title{
Empirical Bayes Genome Wide Association Mapping 
}
\description{
The gwas function calculates the likelihood ratio for each marker under the empirical Bayesian framework. The method also works with multiple populations. 
}
\usage{
gwas(y,gen,fam=NULL,chr=NULL,window=NULL)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozygous, heterozygous and reference homozigous. \code{NA} is allowed.
}
  \item{fam}{
Numeric vector of length (\eqn{n}) indicating which subpopulations (\eqn{i.e.} family) each observation comes from. Default assumes that all observations are from the same populations.
}
  \item{chr}{
Numeric vector indicating the number of markers in each chromosome. The sum of \eqn{chr} must be equal to the number of columns in \eqn{gen}. Default assumes that all markers are from the same chromosome.
}
  \item{window}{
Numeric. If specified, genetic distance between markers is used for moving window strategy (Wang 2015). Window must be specified in Morgans (\eqn{e.g.} 0.05 would represent 5cM). Genetic distance is calculated assuming that individuals are RILs.
}
}
\details{
Special incidence matrix is recreated to optimize the information provided by the subpopulations. Each locus is recoded as a vector with length \eqn{f} equal to number of subpopulations, or NAM families, as the interaction locus by family. For example, a locus heterozigous from an individual from subpopulation 2 is coded as [ 1, 0, 1 , ... ,\eqn{f} ], a locus homozigous for the reference allele from any subpopulation is coded as [ 2, 0, 0, ... , \eqn{f} ] and a locus homozigous for the founder allele from an individual from subpopulation 1 is coded as [ 0, 2, 0, ... ,\eqn{f} ].

The base model for genome scanning includes the fixed effect (\eqn{Xb}), the marker (\eqn{Zu}), the polygene (\eqn{g}) and the residuals (\eqn{e}). If the \eqn{window} term is specified, the model for genome scanning includes three extra terms, the left side genome ( \eqn{Zu[k-1]} ), the right side genome ( \eqn{Zu[k+1]} ) and window polygene ( \eqn{-g[k]} ). 

The polygenic term is calculated only once (Zhang et al 2010) using eigendecomposition with a GEMMA-like algorithm (Zhou ans Stephens 2012). Efficient inversion of capacitance matrix is obtained through the Woodbury matrix identities.

In order to analyze large dateset, one can avoid memory issues by using the function \eqn{gwas2} that has the same arguments, except that the argument 'window' is not implemented.
}
\value{
The function gwas returns a list containing the method deployed (\eqn{Method}), predicted parameters and statistical test (\eqn{PolyTest}), genetic map for NAM panels (\eqn{MAP}) and the marker names (\eqn{SNPs}).
}
\references{

Wang, Q. An Empirical Bayes Method for Genome-Wide Association Studies. W799/Statistical Genomics. PAG XXXII. 2015.

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide association studies. Nat. Genet. 42:355-360.

Zhou, X., & Stephens, M. (2012). Genome-wide efficient mixed-model analysis for association studies. Nature genetics, 44(7), 821-824.
}
\author{
Alencar Xavier, Tiago Pimenta, Qishan Wang and Shizhong Xu
}
\examples{
data(tpod)
test=gwas(y=y,gen=gen[,1:240],fam=fam,chr=chr[1:12])
manhattan(test,type="p",lwd=3)
}
\keyword{ NAM }
\keyword{ GWAS }
\keyword{ BLUP }
