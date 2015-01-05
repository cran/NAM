\name{blup}
\alias{blup}
\title{
Best Linear Unbiased Predictor
}
\description{
Univariate BLUP estimators and variance components fitted by AI-REML.
}
\usage{
blup(y,Z=NULL,X=NULL,K=NULL,Z2=NULL,Z3=NULL,K2=NULL,K3=NULL,conv=1e-8,maxit=50,verb=TRUE)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{Z}{
Formula or numeric matrix (\eqn{n} by \eqn{p}) to incidence matrix for random effect. \code{NA} is not allowed.
}
  \item{X}{
Formula or incidence matrix (\eqn{n} by \eqn{p}) for fixed effect. \code{NA} is not allowed.
}
  \item{K}{
Numeric matrix (\eqn{p} by \eqn{p}). Kinship matrix for random effect with \eqn{p} parameters. \code{NA} is not allowed.
}
  \item{Z2}{
Formula or incidence matrix (\eqn{n} by \eqn{p}) for a second random effect. \code{NA} is not allowed.
}
  \item{Z3}{
Formula or incidence matrix (\eqn{n} by \eqn{p}) for a third random effect. \code{NA} is not allowed.
}
  \item{K2}{
Numeric matrix (\eqn{p} by \eqn{p}). Kinship matrix for the second random effect with \eqn{p} parameters. \code{NA} is not allowed.
}
  \item{K3}{
Numeric matrix (\eqn{p} by \eqn{p}). Kinship matrix for the third random effect with \eqn{p} parameters. \code{NA} is not allowed.
}
  \item{conv}{
Convergence parameter for REML.
}
  \item{maxit}{
Maximum number of iterations.
}
  \item{verb}{
Display verbose: TRUE/FALSE.
}
}
\details{
Solve mixed models with up to 3 random effects through restricted maximum likelihood (REML) using the avarage-information (AI) algorithm (Gilmour et al 1995) as decribed by Lee et al (2006). If there is only one random effect, AI-REML uses Eigen decomposition (Zhou and Stephens 2012) for improved efficiency. Coefficients are calculated using the Henderson model (animal model). Choleski and QR decomposition are used for inversions and linear solutions.
}
\value{
Returns a list with variance components (VarComp), empirical Bayes estimator (Fit), initial observations (Obs), coefficients of fixed (b) and random (u) effects.
}
\references{

Gilmour, A. R., Thompson, R., & Cullis, B. R. (1995). Average information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics, 1440-1450.

Lee, S. H., & van der Werf, J. H. (2006). An efficient variance component approach implementing an average information REML suitable for combined LD and linkage mapping with a general complex pedigree. Genetics Selection Evolution, 38(1), 1-19.

Zhou, X., & Stephens, M. (2012). Genome-wide efficient mixed-model analysis for association studies. Nature genetics, 44(7), 821-824.

}
\author{
Alencar Xavier and William Muir
}
\examples{

# Fitting a random model
data(tpod)
FIT = blup(y=y,Z=~as.factor(fam),maxit=10)

# Fitting GBLUP
G = tcrossprod(gen)
G = G/mean(diag(G))
GBLUP = blup(y=y,K=G,maxit=10)

# Fitting a Gaussian Kernel Average model
D=as.matrix(dist(gen)); D2 = D^2;
h=quantile(as.vector(D2),probs=.05)
K1 = exp(-5/h*D2)
K2 = exp(-1/h*D2)
K3 = exp(-0.2/h*D2)
KA = blup(y=y,K=K1,K2=K2,K3=K3,maxit=10)

# Fitting VanRaden's BLUP
Z=(lm(gen~1)$residuals)
VanRaden = blup(y=y,Z=Z,maxit=10)

}
\keyword{ REML }
\keyword{ BLUP }
