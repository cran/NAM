\name{BGS}
\alias{gibbs}
\alias{plot.gibbs}
\alias{covar}
\title{
Bayesian Gibbs Sampler
}
\description{
Univariate mixed model solver through Gibbs Sampling.
}
\usage{
gibbs(y,Z=NULL,X=NULL,iK=NULL,iR=NULL,Iter=1500,Burn=500,Thin=4,DF=5,S=1,GSRU=FALSE)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{Z}{
Right-hand side formula or list of numeric matrices (\eqn{n} by \eqn{p}) with incidence matrices for random effect. \code{NA} is not allowed.
}
  \item{X}{
Right-hand side formula or incidence matrices (\eqn{n} by \eqn{p}) for fixed effect. \code{NA} is not allowed.
}
  \item{iK}{
Numeric matrix or list of numeric matrices (\eqn{p} by \eqn{p}) corresponding to the the inverse kinship matrix of each random effect with \eqn{p} parameters.
}
  \item{iR}{
Numeric matrix (\eqn{n} by \eqn{n}) corresponding to the inverse residual correlation matrix.
}
  \item{Iter}{
Integer. Number of iterations or samples to be generated.
}
  \item{Burn}{
Integer. Number of iterations or samples to be discarted.
}
  \item{Thin}{
Integer. Thinning parameter, used to save memory by storing only one every 'Thin' samples.
}
  \item{DF}{
Integer. Hyperprior degrees of freedom of variance components.
}
  \item{S}{
Integer or NULL. Hyperprior shape of variance components. If NULL, the hyperprior solution for the scale parameter is calculated as proposed by de los Campos et al. (2013).
}
  \item{GSRU}{
Logical. If TRUE, it updates the regression coefficient using Gauss-Seidel Residual Update (Legarra and Misztal 2008). Useful for p>>n, but does not work when iK or iR are provided.
}
}
\details{
Solve Gaussian mixed models in the Bayesian framework as described by Garcia-Cortes and Sorensen (1996) and Sorensen and Gianola, D. (2002) with conjugated priors. 
}
\value{
The function gibbs returns a list with variance components distribution a posteriori (Posterior.VC) and mode estimated (VC.estimate), a list with the posterior distribution of regression coefficients (Posterior.Coef) and the posterior mean (Coef.estimate), and the fitted values using the mean (Fit.mean) of posterior coefficients. 
}
\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

Legarra, A., & Misztal, I. (2008). Technical note: Computing strategies in genome-wide selection. Journal of dairy science, 91(1), 360-366.

Garcia-Cortes, L. A., and Sorensen, D. (1996). On a multivariate implementation of the Gibbs sampler. Genetics Selection Evolution, 28(1), 121-126.

Sorensen, D., & Gianola, D. (2002). Likelihood, Bayesian, and MCMC methods in quantitative genetics. Springer Science & Business Media.
}
\author{
Alencar Xavier
}
\examples{

# Fitting Mixed Model
data(tpod)
S = seq(1,350,50)
test1=gibbs(y,X=~factor(fam),Z=gen[,S],S=0.5)
plot(test1)

# Fitting GBLUP
K = tcrossprod(gen)
K = K/mean(diag(K))
iK = chol2inv(K)
test2=gibbs(y,iK=iK)
plot(test2)

# Fitting RKHS
E = dist(gen)
G = exp(-E^2/mean(E^2))
EIG = eigen(G,symmetric = TRUE)
ev = 20
U = EIG$vectors[,1:ev]
iV = diag(1/EIG$values[1:ev])
test3=gibbs(y,Z=U,iK=iV,S=NULL)
plot(test3)
}
\keyword{ BLUP }
\keyword{ BGS }