\name{Bayesian Gibbs Sampler}
\alias{gibbs}
\title{
Bayesian Gibbs Sampler
}
\description{
Univariate mixed model solver through Gibbs Sampling.
}
\usage{
gibbs(y,X=NULL,Z=NULL,iK=NULL,Iter=1500,Burn=500,Thin=4,DF=5)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{X}{
Formula or incidence matrix (\eqn{n} by \eqn{p}) for fixed effect. \code{NA} is not allowed.
}
  \item{Z}{
Formula or list of numeric matrix (\eqn{n} by \eqn{p}) with incidence matrices for random effect. \code{NA} is not allowed.
}
  \item{iK}{
Numeric matrix or list of numeric matrices (\eqn{p} by \eqn{p}) corresponding to the the inverse kinship matrix of each random effect with \eqn{p} parameters.
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
Integer. Hyper-parameter regarding degrees of freedom.
}
}
\details{
Solve Gaussian mixed models in the Bayesian framework as described by Garc?a-Cort?s and Sorensen (1996) and Sorensen and Gianola, D. (2002) with conjugated priors. The hyper-prior solution for the scale parameter of inverse Chi squared is calculated as proposed by de los Campos et al. (2013).
}
\value{
The function gibbs returns a list with variance components distribution a posteriori (Posterior.VC) and mode estimated (VC.mode), a list with the posterior distribution of regression coefficients (Posterior.Coef) and the posterior mode (Coef.mode) and mean (Coef.mean), and the fitted values using the mean (Fit.mean) and mode (Fit.mode) of posterior coefficients. 
}
\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

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
test=gibbs(y,X=gen[,S],Z=~factor(fam))
plot(test$Fit.mean,y)

# Fitting GBLUP
K = tcrossprod(gen)
K = K/mean(diag(K))
iK = chol2inv(K)
test2=gibbs(y,iK=iK)
plot(test2$Fit.mean,y)

}
\keyword{ BLUP }
\keyword{ GIBBS }