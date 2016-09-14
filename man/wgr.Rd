\name{WGR}
\alias{wgr}
\alias{ben}
\title{
Whole-genome Regression
}
\description{
Univariate model to find breeding values through regression with optional resampling techniques.
}
\usage{
wgr(y,gen,it=1500,bi=500,th=1,bag=1,rp=FALSE,iv=FALSE,pi=0,df=5,R2=0.5,
    eigK=NULL,VarK=0.95,verb=FALSE)
ben(y,gen,it=750,bi=250,th=1,bag=0.80,alpha=0.5,wpe=50,MH=FALSE,verb=TRUE)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers.
}
  \item{it}{
Integer. Number of iterations or samples to be generated.
}
  \item{bi}{
Integer. Burn-in, the number of iterations or samples to be discarted.
}
  \item{th}{
Integer. Thinning parameter, used to save memory by storing only one every 'th' samples.
}
  \item{bag}{
If different than one (= complete data), it indicates the proportion of data to be subsampled in each Markov chain. For datasets with moderate number of observations, values of bag from 0.30 to 0.60 may speed up computation without losses in predicion properties. This argument enable users to enhance MCMC through SBMC (\eqn{subssampling bootstrap Markov chain}).
}
  \item{rp}{
Logical. Use replacement for bootstrap samples when bag is different than one.
}
  \item{iv}{
Logical. Assign markers independent variance. If true, turns the default model BLUP into BayesA. For this model, the shape parameter is conjugated by a gamma with hyperpriors calculated based on the R2 rule.
}
  \item{pi}{
Value between 0 and 1. If greater than zero it activates variable selection, where markers have expected probability pi of having null effect, in other words, the lower the value of pi the more stringent the variable selection is. The model conjugates variable selection from a Beta-Binomial distribution.
}
  \item{df}{
Hyperprior degrees of freedom of variance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape as proposed by de los Campos et al. (2013).
}
  \item{eigK}{
Output of function 'eigen'. Spectral decomposition of the kernel used to compute the polygenic term.
}
  \item{VarK}{
Numeric between 0 and 1. For reduction of dimensionality. Indicates the proportion of variance explained by Eigenpairs used to fit the polygenic term.
}
  \item{verb}{
Logical. If verbose is TRUE, function displays MCMC progress bar.
}
  \item{alpha}{
Numeric between 0 and 1. Starting value of alpha parameter for the bagging elastic net model, where 0 is ridge (L2) and 1 is lasso (L1).
}
  \item{wpe}{
Weight of prediction error sum of squared for the jump function that defines alpha in the bagging elastic net model.
}
  \item{MH}{
Logical. If TRUE, the search for alpha is performed via Metropolis-Hastings. If FALSE, acceptance-rejection.
}
}
\details{
The model for the whole-genome regression is as follows:

\deqn{y = mu + Xg + u + e}

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{g} is the product of two terms (\eqn{g = bg}), \eqn{b} is the effect of an allele substitution, \eqn{d} is an indicator variable that define whether or not the marker should be included into the model, \eqn{u} is the polygenic term and \eqn{e} is the residual term.

Users can obtain four WGR methods out of this function: BRR (pi=0,iv=F), BayesA (pi=0,iv=T), BayesB (pi=0.01,iv=T) and BayesC (pi=0.01,iv=F). The full theoretical basis of each model is described by de los Campos et al. (2013).

Gibbs sampler that updates regression coefficients is adapted from GSRU algorithm (Legarra and Misztal 2008). The variable selection works through the unconditional prior algorithm proposed by Kuo and Mallick (1998). The polygenic term is solved by Bayesian algorithm of reproducing kernel Hilbert Spaces proposed by de los Campos et al. (2010).

The model for the bagging elastic net (ben) is as follows:

\deqn{y = mu + Xb + e}

Elastic net is controlled by two parameters, alpha and lambda. Lambda is analytically estimated as the ratio between residual variance and parameter variance (\eqn{\lambda = Ve/Vb}). Variance components and regression coefficients are updated as a Bayesian ridge regression. Alpha is found through random walk: A new value of alpha is proposed in each MCMC round. The acceptance of the new value depends on the algorithm (Metropolis-Hastings or Acceptance-Rejection), and the probability of a new value of alpha to be accepted is based on the Elastic-Net loss function and out-of-bag prediction error. Unfortunately, there is a computation burden associated to the fact that the regression coefficients are updated for the current and new value of alpha.

For further statistical background: \code{system(paste('open',system.file("doc","background_stat_gen.pdf",package="NAM")))}

}
\value{
The function wgr returns a list with expected value from the marker effect (\eqn{b}), probability of marker being in the model (\eqn{d}), regression coefficient (\eqn{g}), variance of each marker (\eqn{Vb}), the intercept (\eqn{mu}), the polygene (\eqn{u}) and polygenic variance (\eqn{Vk}), residual variance (\eqn{Ve}) and the fitted value (\eqn{hat}).

variance components distribution a posteriori (Posterior.VC) and mode estimated (VC.estimate), a list with the posterior distribution of regression coefficients (Posterior.Coef) and the posterior mean (Coef.estimate), and the fitted values using the mean (Fit.mean) of posterior coefficients.

}
\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

de los Campos, G., Gianola, D., Rosa, G. J., Weigel, K. A., & Crossa, J. (2010). Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods. Genetics Research, 92(04), 295-308.

Kuo, L., & Mallick, B. (1998). Variable selection for regression models. Sankhya: The Indian Journal of Statistics, Series B, 65-81.

Legarra, A., & Misztal, I. (2008). Technical note: Computing strategies in genome-wide selection. Journal of dairy science, 91(1), 360-366.

}
\author{
Alencar Xavier
}
\examples{

data(tpod)
gen = gen[,seq(1,376,5)]

# BLUP
BRR = wgr(y,gen,iv=FALSE,pi=0,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BRR$hat)

# BayesA
BA = wgr(y,gen,iv=TRUE,pi=0,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BA$hat)

# BayesB
BB = wgr(y,gen,iv=TRUE,pi=.01,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BB$hat)

# BayesC
BC = wgr(y,gen,iv=FALSE,pi=.01,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BC$hat)

# Bayes Elastic Net
BEN = ben(y,gen,it=200,bi=50)
cor(y,BEN$hat)
}