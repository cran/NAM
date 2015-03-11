\name{REML}
\alias{reml}
\alias{MCreml}
\title{
Restricted Maximum Likelihood
}
\description{
Univariate REML estimators and variance components for a single random variable fitted by an EMMA-like algorithm.
}
\usage{
reml(y,X=NULL,Z=NULL,K=NULL)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{X}{
Formula or incidence matrix (\eqn{n} by \eqn{p}) for fixed effect. \code{NA} is not allowed.
}
  \item{Z}{
Formula or numeric matrix (\eqn{n} by \eqn{p}) that corresponds to the incidence matrix of random effect. \code{NA} is not allowed.
}
  \item{K}{
Numeric matrix (\eqn{p} by \eqn{p}). Kinship matrix for random effect with \eqn{p} parameters. \code{NA} is not allowed.
}
}
\details{
Solve mixed models with a single random effects minizing the log restricted maximum likelihood (REML) using the EMMA algorithm (Kang et al 2008). Prediction of random coefficients are performed according to VanRaden (2008). A sample strategy is also available as MCreml.
}
\value{
The function reml returns a list with variance components and heritability (VC), fixed effect coefficients and standard variations (Fixed) and estimated breeding values (EBV). 
}
\references{

Kang, H. M., Zaitlen, N. A., Wade, C. M., Kirby, A., Heckerman, D., Daly, M. J., & Eskin, E. (2008). Efficient control of population structure in model organism association mapping. Genetics, 178(3), 1709-1723.

VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of dairy science, 91(11), 4414-4423.

}
\author{
Alencar Xavier, Tiago Pimenta and Shizhong Xu
}
\examples{

# Fitting a random model
data(tpod)
FIT = reml(y=y,Z=~as.factor(fam))

# Fitting GBLUP
G = tcrossprod(gen)
G = G/mean(diag(G))
GBLUP = reml(y=y,K=G)

}
\keyword{ REML }
\keyword{ BLUP }
